function [cycle_count, SRAM_Access, DRAM_Access] = ZtoXtoYtoF(Layer_param, Hardware_param, Tech_param)

%X: Direction along the width of ofmap
%Y: Direction along the height of ofmap
%Z: Direction along the channels of ifmap
%F: Direction along the channels of ofmap = Direction along the 3D filters

% This function performs the computation of the input layer using the ZXYF branch

%% Layer Specification 
filter_height = Layer_param(1);
filter_width = Layer_param(2);
ifmap_height = Layer_param(3);
ifmap_width = Layer_param(4);
Nos_of_channel = Layer_param(5);
ofmap_height = Layer_param(6);
ofmap_width = Layer_param(7);
Nos_of_filter = Layer_param(8);
stride = Layer_param(9);

%% Accelerator Specification
% bit width
bw_filter = Hardware_param(1);
bw_ifmap = Hardware_param(2);
bw_psum = Hardware_param(3);
bw_ofmap = Hardware_param(4);

% MAC array
Array_row = Hardware_param(5);
Array_column = Hardware_param(6);

% On-chip SRAM storage
SRAM_filter = Hardware_param(7); % in bit
SRAM_ifmap = Hardware_param(8); % in bit
SRAM_psum = Hardware_param(9); % in bit 

% To incorporation DRAM bandwidth induced stall in the calculation of cycle count
BW_DRAM = Hardware_param(10); %in bit/cycle
DRAM_block_size = Hardware_param(11); %in bit

%% Access energy specification 
Energy_MAC = Tech_param(1);    % Energy (in joule) for one n-bit interger/fixed point add+mul (MAC) in the respective technology node 

E_RF_to_ALU = Tech_param(2);    % Access Energy per bit in Joule
E_SRAM_to_RF = Tech_param(3);         
E_DRAM_to_SRAM = Tech_param(4);


%% Computation of the unit volume per cycle in the MAC array
unitvol_nosof_channel = min(floor(Array_row/filter_width), Nos_of_channel);   % Nos of channel in the unit volume
% If there is less Nos_of_channel than Array_row/filter_width then just process all the channels and some PE row will remain unused

unitvol_nosof_3D_filter = min(Array_column, Nos_of_filter);      % Nos of different 3D filter in the unit volume
% If there is less 3D filters than nos of array column then just process all the 3D filters and some PE column will remain unused

unit_vol_filter = (1 * filter_width * unitvol_nosof_channel) * unitvol_nosof_3D_filter; % Nos of filter elements in the unit volume
unit_vol_ifmap = (1 * filter_width * unitvol_nosof_channel);            % Nos of ifmap elements in the unit volume

used_MAC_per_column = floor(Array_row/filter_width) * filter_width;
unused_MAC_per_column = Array_row - used_MAC_per_column;  % unused MAC unit per column


%% Tree process along z-direction of ifmap

%%%%% Cycle count calculation
initial_offset_cycle = (used_MAC_per_column - 1);  %After first (used_MAC_per_column - 1) cycle, the pipeline is full
cycle_count_ideal_z = ceil(Nos_of_channel/unitvol_nosof_channel) * filter_height;   %#of cycle required to process the full z-direction with the filter-template
cycle_count_z = cycle_count_ideal_z + initial_offset_cycle;


%SRAM_filter memory requirement check;
%%%%%%% GIVE RSC LOWER BOUND ON IFMAP
% This minimum requirement is same as x-direction. While moving towards z-direction, new channels can be loaded from DRAM.
% DRAM bandwidth can result large stall if filter SRAM is too small
SRAM_req_filter_z = unit_vol_filter * filter_height * bw_filter;      %in bit, Assuming full filter height needs to fit
if (SRAM_req_filter_z <= SRAM_filter)
    disp("Full z-direction possible for filter template")
else
    % Applying a lower bound
    Min_SRAM_filter_z = (SRAM_req_filter_z)/(8 * 1024); % Minimum SRAM requirement for filter data in kB to proceed along z-direction
    disp("ALERT:Filter-SRAM is too small, increase it") %Assuming a minimum size for filter SRAM
    fprintf("Minimum SRAM memory for filter data for this layer is %f kB\n",Min_SRAM_filter_z)
    % returning inf valus for the return parameters so that this branch does not get selected and sweep fucntion continue to run for the next hardware point
    cycle_count = [inf inf inf];
    SRAM_Access = [inf inf inf];  
    DRAM_Access = [inf inf inf];
    return
end

%SRAM_psum memory requirement
Generated_psum = 1 * 1 * unitvol_nosof_3D_filter;
SRAM_req_psum_z = Generated_psum * bw_psum;   % in bit
if (SRAM_req_psum_z <= SRAM_psum)
    disp("Full z-direction possible for psum")
    Generated_psum_effective = Generated_psum;
else
    %Applying a lower bound on the psum-SRAM. Basically psum-SRAM has to hold minimum unitvol_nosof_3D_filter psum data at a time. 
    %Since z-direction is being processed first here, the ofmaps are formed and can be moved to DRAM.
    disp("ALERT: psum-SRAM is too small, increase it")
    Min_SRAM_psum_z = SRAM_req_psum_z / (8 * 1024);
    fprintf("Minimum SRAM memory for psum data is %f kB\n",Min_SRAM_psum_z)
    % returning inf valus for the return parameters so that this branch does not get selected and sweep fucntion continue to run for the next hardware point
    cycle_count = [inf inf inf];
    SRAM_Access = [inf inf inf];  
    DRAM_Access = [inf inf inf];
    return
end

%SRAM_ifmap memory requirement check
% This minimum requirement is very small. While moving towards z-direction, new channels can be loaded from DRAM.
% DRAM bandwidth can result large stall if ifmap SRAM is too small 
SRAM_req_ifmap_z = filter_width * filter_height * unitvol_nosof_channel * bw_ifmap;   % in bit
if (SRAM_req_ifmap_z <= SRAM_ifmap)
    disp("full z-direction possible for ifmap") 
else
    % Applying a lower bound
    Min_SRAM_ifmap_z = (SRAM_req_ifmap_z)/(8 * 1024); % Minimum SRAM requirement for filter data in kB to proceed along z-direction
    disp("ALERT:ifmap-SRAM is too small, increase it") %Assuming a minimum size for ifmap SRAM
    fprintf("Minimum SRAM memory for ifmap data for this layer is %f kB\n",Min_SRAM_ifmap_z)
    % returning inf valus for the return parameters so that this branch does not get selected and sweep fucntion continue to run for the next hardware point
    cycle_count = [inf inf inf];
    SRAM_Access = [inf inf inf];  
    DRAM_Access = [inf inf inf];
    return   
end


%%%%%% Access cost for filter data
gamma_DRAM_filter_z = 1; % each filter data is loaded once from DRAM
Access_DRAM_filter_z = filter_width * filter_height * Nos_of_channel * unitvol_nosof_3D_filter;  % in element, not bit 
% at each cycle "unit_vol_filter" new filter is loaded from SRAM 
Access_SRAM_filter_z = unit_vol_filter * (Nos_of_channel/unitvol_nosof_channel) * filter_height;
% the ceil is removed from cycle_count_ideal_z to reflect the fact that at the last pass, when there is less nos of channel, less nos of data will
% be loaded from SRAM and some PE will remain unused.


%%%%% Access cost for ifmap data
gamma_DRAM_ifmap_z = 1; %each ifmap data is loaded once from DRAM
Access_DRAM_ifmap_z = filter_width * filter_height * Nos_of_channel;  % in element
% at each cycle "unit_vol_ifmap" new ifmap is loaded from SRAM
Access_SRAM_ifmap_z = (Nos_of_channel/unitvol_nosof_channel) * filter_height * filter_width * unitvol_nosof_channel; % in element
% instead of cycle_count_ideal_z used it without ceil for the same reason as mentioned for the filter


% Access cost for psum data, psum does not go to DRAM, so only SRAM access
Access_SRAM_psum_z = Generated_psum_effective;  % in element, not bit
% psum write occur after finishing the z-direction only once

%This is the access energy up to a specefic stage of the decision
Weighted_Access_z = (Access_DRAM_filter_z * bw_filter + Access_DRAM_ifmap_z * bw_ifmap) * E_DRAM_to_SRAM + ...
    (Access_SRAM_filter_z * bw_filter + Access_SRAM_ifmap_z * bw_ifmap + Access_SRAM_psum_z * bw_psum) * E_SRAM_to_RF;  % Access energy (in joule)
SRAM_Access_z = (Access_SRAM_filter_z * bw_filter + Access_SRAM_ifmap_z * bw_ifmap + Access_SRAM_psum_z * bw_psum) * E_SRAM_to_RF;
DRAM_Access_z = (Access_DRAM_filter_z * bw_filter + Access_DRAM_ifmap_z * bw_ifmap) * E_DRAM_to_SRAM;

%%%%%%%%%%% Calculating DRAM bandwidth induced stall count
DRAM_filter_data_z = filter_height * filter_width * Nos_of_channel * unitvol_nosof_3D_filter * bw_filter;
DRAM_ifmap_data_z = filter_height * filter_width * Nos_of_channel * bw_ifmap;
DRAM_data_z = DRAM_filter_data_z + DRAM_ifmap_data_z; % in bit, this is the number of data loaded from DRAM to process the first ofmap element
%cycle_to_load_DRAM_data_x = ceil(DRAM_data_x/BW_DRAM)
DRAM_data_loaded_z = ceil(DRAM_data_z/DRAM_block_size) * DRAM_block_size;   % The data loaded from DRAM is an integer multiple of the DRAM block size
cycle_to_load_DRAM_data_z = ceil(DRAM_data_loaded_z/BW_DRAM);               % #of cycle to load the data from DRAM

if cycle_to_load_DRAM_data_z > cycle_count_ideal_z
    DRAM_stall_z = cycle_to_load_DRAM_data_z - cycle_count_ideal_z;
else
    DRAM_stall_z = 0;
end

% Calculating #of slack cycles from this stage
if (DRAM_stall_z == 0)
    slack_cycle_z = cycle_count_ideal_z - cycle_to_load_DRAM_data_z;  %#of extra cycles during computation when data can be loaded from DRAM
else
    slack_cycle_z = 0;   % No extra cycles for data load since data load time is higher than compute time
end


%% Tree process along z->x direction of ifmap

% filter memeory requirement check is not necessary

% ifmap memory requirement check is also not necessary, the check at depth-1 is enough

% psum memory check is not necessary.
% since z is the first direction, no effective triggers
ifmap_width_effective = ifmap_width;
ofmap_width_effective = ofmap_width;

%%%% cycle count calculation
%the same process to form the first ofmap element is repeated for all the next ofmap elements in the first row of ofmap (i.e., x-direction)
cycle_count_ideal_zx = cycle_count_ideal_z * ofmap_width_effective;   %without the initial offset cycles
cycle_count_zx = cycle_count_ideal_zx + initial_offset_cycle;


%%%%% Access cost for filter data
% Reverse formatting implementation for DRAM access count
filter_channels_to_fit_SRAM = floor(SRAM_filter/(filter_height * filter_width * bw_filter)); % #of 2D filter planes which fits into filter-SRAM
Total_filter_planes_k = unitvol_nosof_3D_filter * Nos_of_channel;    % #of total 2D filter planes in the unitvol_nosof_3D_filter filter data
if (filter_channels_to_fit_SRAM >= Total_filter_planes_k)
    disp("filter SRAM is large enough to fit RSCK volume of filter")
    gamma_DRAM_filter_zx = 1;   % filter-SRAM has enough memory, each filter data loaded from DRAM once
    Access_DRAM_filter_zx = filter_height * filter_width * Nos_of_channel * unitvol_nosof_3D_filter;
else
    disp("filter SRAM cannot fit RSCK volume of filter, will result high DRAM access for filter")
    % calculating minimum filter-SRAM size to avoid this situation and informing the user
    filter_RSCK_volume = (filter_height * filter_width * Nos_of_channel * unitvol_nosof_3D_filter * bw_filter) / (8 * 1024); %in kB
    fprintf("Minimum filter-SRAM memory for this layer to fit RSCK filter data is %f kB\n",filter_RSCK_volume)
    
    Nos_of_discarded_filter_planes = Total_filter_planes_k - filter_channels_to_fit_SRAM; % #of 2D filter planes which are discarded from SRAM
    % During the processing of first ofmap element in the x-direction, all filter are loaded once from DRAM
    First_ofmap_element_Access = filter_width * filter_height * Nos_of_channel * unitvol_nosof_3D_filter;
    % For all the next ofmap elements in the first ofmap row, the amount of discarded channels need to be loaded.
    Next_ofmap_element_Access = filter_width * filter_height * Nos_of_discarded_filter_planes;
    Access_DRAM_filter_zx = First_ofmap_element_Access + Next_ofmap_element_Access * (ofmap_width_effective - 1);
end
%the same access pattern to form the first ofmap element is repeated for all the next ofmap elements in the first row of ofmap (i.e., x-direction)
Access_SRAM_filter_zx = Access_SRAM_filter_z * ofmap_width_effective;

%%%%% Access cost for ifmap data 
ifmap_RS_template_to_fit_SRAM = floor(SRAM_ifmap / (filter_width * filter_height * bw_ifmap)); %#of ifmap RS template fits in ifmap SRAM
if (ifmap_RS_template_to_fit_SRAM >= Nos_of_channel)
    % in this case, while moving towards zx-direction, due to the sliding window, only new ifmaps can be loaded and the data not needed can be discarded
    disp("ifmap-SRAM can fit RSC volume of ifmap, hence zx direction processed by loading ifmap once from DRAM")
    gamma_DRAM_ifmap_zx = 1;  %ifmap-SRAM has enough memory, each ifmap data loaded from DRAM once
    Access_DRAM_ifmap_zx = ifmap_width_effective * filter_height * Nos_of_channel;
else
    disp("ifmap-SRAM cannot fit RSC volume of ifmap, will result high DRAM access for ifmap")
    % calculating minimum ifmap-SRAM size to avoid this situation and informing the user
    ifmap_RSC_volume = (filter_height * filter_width * Nos_of_channel * bw_ifmap) / (8 * 1024); %in kB
    fprintf("Minimum ifmap-SRAM memory for this layer to fit RSC ifmap data is %f kB\n",ifmap_RSC_volume)
    
    % Implementation with reverse order formatting
    Discarded_ifmap_RS_channel = Nos_of_channel - ifmap_RS_template_to_fit_SRAM;   % #of ifmap RS templates need to be discarded while moving towards z-direction
    First_ofmap_elem_Access_ifmap = filter_height * filter_width * Nos_of_channel;  % #of ifmap data loaded during processing the first ofmap element of first row
    % In the line below: the first term is the newly loaded ifmap data for each new ofmap element and
    % the second term is the previously discarded ifmap channels which need to be loaded again each time we process a new ofmap element
    Next_ofmap_elem_Access_ifmap = (filter_height * stride * Nos_of_channel) + (filter_height * (filter_width - stride) * Discarded_ifmap_RS_channel);
    Access_DRAM_ifmap_zx = First_ofmap_elem_Access_ifmap + Next_ofmap_elem_Access_ifmap * (ofmap_width_effective - 1); 
end
%the same access pattern to form the first ofmap element is repeated for all the next ofmap elements in the first row of ofmap (i.e., x-direction)
Access_SRAM_ifmap_zx = Access_SRAM_ifmap_z * ofmap_width_effective;


%%%%% Access cost for psum data
%the same access pattern to form the first ofmap element is repeated for all the next ofmap elements in the first row of ofmap (i.e., x-direction)
Access_SRAM_psum_zx = Access_SRAM_psum_z * ofmap_width_effective;


%This is the access energy up to this depth of the decision
Weighted_Access_zx = (Access_DRAM_filter_zx * bw_filter + Access_DRAM_ifmap_zx * bw_ifmap) * E_DRAM_to_SRAM + ...
    (Access_SRAM_filter_zx * bw_filter + Access_SRAM_ifmap_zx * bw_ifmap + Access_SRAM_psum_zx * bw_psum) * E_SRAM_to_RF;  % Access energy (in joule)
SRAM_Access_zx = (Access_SRAM_filter_zx * bw_filter + Access_SRAM_ifmap_zx * bw_ifmap + Access_SRAM_psum_zx * bw_psum) * E_SRAM_to_RF;
DRAM_Access_zx = (Access_DRAM_filter_zx * bw_filter + Access_DRAM_ifmap_zx * bw_ifmap) * E_DRAM_to_SRAM;

%%%%%%%%%%% Calculating DRAM bandwidth induced stall count
% omitting the implementation when filter-SRAM cannot fit RSCK volume of filter since in this case, this branch will incur very high stall
% and energy (have seen from the energy implementation) and will never be competitive as compared to other branches.

if (SRAM_filter >= filter_height * filter_width * Nos_of_channel * unitvol_nosof_3D_filter * bw_filter)
    % fitting RSC volume of ifmap in ifmap-SRAM is the minimum requirement
    New_loaded_filter_zx = 0;   % when moving z->x direction, the previosuly loaded filter data is reused
    New_loaded_ifmap_POE_zx = stride * filter_height * Nos_of_channel * bw_ifmap; % in bit; Additional ifmap data loaded to process each next ofmap element in x-direction
    
    DRAM_data_POE_zx = New_loaded_filter_zx + New_loaded_ifmap_POE_zx; %in bit; total new data loaded from DRAM per ofmap element (POE) in x-direction
    
    DRAM_data_loaded_POE_zx = ceil(DRAM_data_POE_zx/DRAM_block_size) * DRAM_block_size;
    cycle_to_load_DRAM_data_POE_zx = ceil(DRAM_data_loaded_POE_zx/BW_DRAM); % cycle count to load data per ofmap element in x-direction
    
    DRAM_stall_POE_zx = max(0, cycle_to_load_DRAM_data_POE_zx - cycle_count_ideal_z);
    DRAM_stall_zx = DRAM_stall_z + DRAM_stall_POE_zx * (ofmap_width_effective - 1);
    % first term: for the first element where filter dada is also loaded
    % second term: for the rest of the ofmap elements in x-direction
    % omitting the slack_cycle variable here since for the first ofmap-element we need to load more data (both filter & ifmap) than the next ofmap element (only ifmap).
    % So, if there is no stall at z, there will be no stall at x. 
    
    slack_cycle_POE_zx = max(0, cycle_count_ideal_z - cycle_to_load_DRAM_data_POE_zx);
    slack_cycle_zx = slack_cycle_z + slack_cycle_POE_zx * (ofmap_width_effective - 1);
else
    disp("filter-SRAM cannot fit RSCK volume of filter, Hence this branch will never be competitive")
end


%% Tree process along z->x->y direction of ifmap

% filter memeory requirement check is not necessary
% ifmap memory requirement check is also not necessary, the check at depth-1 is enough
% psum memory check is not necessary. ofmaps are formed and can be moved to DRAM
% no effective triggers
ifmap_height_effective = ifmap_height;
ofmap_height_effective = ofmap_height;

%%%%%% cycle count calculation
%the same process to form the first ofmap row is repeated for all the next ofmap rows
cycle_count_ideal_zxy = cycle_count_ideal_zx * ofmap_height_effective;   %without the initial offset cycles
cycle_count_zxy = cycle_count_ideal_zxy + initial_offset_cycle;

%%%%% Access cost for filter data
% Reverse formatting implementation
filter_channels_to_fit_SRAM = floor(SRAM_filter/(filter_height * filter_width * bw_filter)); % #of 2D filter planes which fits into filter-SRAM
Total_filter_planes_k = unitvol_nosof_3D_filter * Nos_of_channel;    % #of total 2D filter planes in the unitvol_nosof_3D_filter filter data
if (filter_channels_to_fit_SRAM >= Total_filter_planes_k)
    disp("filter SRAM is large enough to fit RSCK volume of filter")
    gamma_DRAM_filter_zxy = 1;   % filter-SRAM has enough memory, each filter data loaded from DRAM once
    Access_DRAM_filter_zxy = filter_height * filter_width * Nos_of_channel * unitvol_nosof_3D_filter;
else
    disp("filter SRAM cannot fit RSCK volume of filter, will result high DRAM access for filter")
    % calculating minimum filter-SRAM size to avoid this situation and informing the user
    filter_RSCK_volume = (filter_height * filter_width * Nos_of_channel * unitvol_nosof_3D_filter * bw_filter) / (8 * 1024); %in kB
    fprintf("Minimum filter-SRAM memory for this layer to fit RSCK filter data is %f kB\n",filter_RSCK_volume)
    
    Nos_of_discarded_filter_planes = Total_filter_planes_k - filter_channels_to_fit_SRAM; % #of 2D filter planes which are discarded from SRAM
    % During the processing of first ofmap element in the xy plane, all filter are loaded once from DRAM
    First_ofmap_element_Access = filter_width * filter_height * Nos_of_channel * unitvol_nosof_3D_filter;
    % For all the next ofmap elements in the xy ofmap plane, the amount of discarded channels need to be loaded.
    Next_ofmap_element_Access = filter_width * filter_height * Nos_of_discarded_filter_planes;
    Access_DRAM_filter_zxy = First_ofmap_element_Access + Next_ofmap_element_Access * ((ofmap_width_effective * ofmap_height_effective) - 1);
    % Here, (ofmap_width_effective * ofmap_height_effective) gives the nos of ofmap element in the xy plane of ofmap
end
%the same access pattern to form the first ofmap row is repeated for all the next ofmap rows 
Access_SRAM_filter_zxy = Access_SRAM_filter_zx * ofmap_height_effective;

%%%%% Access cost for ifmap data
ifmap_zR_template_to_fit_SRAM = floor(SRAM_ifmap / (Nos_of_channel * filter_height * bw_ifmap)); %#of ifmap zR template fits in ifmap SRAM
ifmap_RS_template_to_fit_SRAM = floor(SRAM_ifmap / (filter_width * filter_height * bw_ifmap)); %#of ifmap RS template fits in ifmap SRAM

if (ifmap_zR_template_to_fit_SRAM >= ifmap_width_effective) % Case-I
    disp("ifmap-SRAM can fit RWC volume of ifmap, hence zxy direction processed by loading ifmap once from DRAM")
    gamma_DRAM_ifmap_zxy = 1;  %ifmap-SRAM has enough memory, each ifmap data loaded from DRAM once
    Access_DRAM_ifmap_zxy = ifmap_width_effective * ifmap_height_effective * Nos_of_channel;
elseif ((ifmap_zR_template_to_fit_SRAM < ifmap_width_effective) && (ifmap_RS_template_to_fit_SRAM >= Nos_of_channel)) % Case-II
    disp("ifmap-SRAM can fit RSC volume of ifmap, but cannot fit RWC volume of ifmap")
    % Implementation with reverse order formatting
    Discarded_ifmap_zR_template = ifmap_width_effective - ifmap_zR_template_to_fit_SRAM;   % #of ifmap zR templates need to be discarded while moving towards x-direction
    First_ofmap_row_Access_ifmap = filter_height * ifmap_width_effective * Nos_of_channel;  % #of ifmap data loaded during processing the first ofmap row
    % In the line below: the first term is the newly loaded ifmap data for each new ofmap row and
    % the second term is the previously discarded ifmap zR templates which need to be loaded again each time we process a new ofmap row
    Next_ofmap_row_Access_ifmap = (ifmap_width_effective * stride * Nos_of_channel) + Discarded_ifmap_zR_template * (filter_height - stride) * Nos_of_channel;
    Access_DRAM_ifmap_zxy = First_ofmap_row_Access_ifmap + Next_ofmap_row_Access_ifmap * (ofmap_height_effective - 1);
else % Case-III
    disp("Case-III trigerring")
    Access_DRAM_ifmap_zxy = Access_DRAM_ifmap_zx * ofmap_height_effective;    
end
%the same access pattern to form the first ofmap row is repeated for all the next ofmap rows 
Access_SRAM_ifmap_zxy = Access_SRAM_ifmap_zx * ofmap_height_effective;

%%%%% Access cost for psum data
%the same access pattern to form the first ofmap row is repeated for all the next ofmap rows
Access_SRAM_psum_zxy = Access_SRAM_psum_zx * ofmap_height_effective;

%This is the access energy up to this depth of the decision
Weighted_Access_zxy = (Access_DRAM_filter_zxy * bw_filter + Access_DRAM_ifmap_zxy * bw_ifmap) * E_DRAM_to_SRAM + ...
    (Access_SRAM_filter_zxy * bw_filter + Access_SRAM_ifmap_zxy * bw_ifmap + Access_SRAM_psum_zxy * bw_psum) * E_SRAM_to_RF;
SRAM_Access_zxy = (Access_SRAM_filter_zxy * bw_filter + Access_SRAM_ifmap_zxy * bw_ifmap + Access_SRAM_psum_zxy * bw_psum) * E_SRAM_to_RF;
DRAM_Access_zxy = (Access_DRAM_filter_zxy * bw_filter + Access_DRAM_ifmap_zxy * bw_ifmap) * E_DRAM_to_SRAM;

%%%%%%%%%%% Calculating DRAM bandwidth induced stall count
if (SRAM_filter >= filter_height * filter_width * Nos_of_channel * unitvol_nosof_3D_filter * bw_filter)
    % no filter data need to load during the processing of y-direction, previously loaded filter data is used
    if (DRAM_stall_zx <= 0) 
        disp("CASE-1: 0 0 0")
        % This means there is no stall at z, x and each zx block can compute itself without incorporating any stall, so no stall to complete the zxy direction
        DRAM_stall_zxy = 0;
        % No need to calculate slack cycle for this case since there will be no stall while moving to F-direction as well
    else
        disp("CASE-2: 1 x x")
        % Depending on ifmap-SRAM either only additional ifmap data is loaded in each slot of next zx-y2 slot or full ifmap-z is loaded in each slot of next zx-y2.
        % For simplicity, assuming if ifmap-SRAM can fit xRz volume of ifmap then imposing only additional ifmap data to load, otherwise imposing
        % regular DRAM_ifmap_data_z ifmap to load for all next slot of zx-y2
        if (SRAM_ifmap>= ifmap_width_effective * filter_height * Nos_of_channel * bw_ifmap)
            ifmap_z_reverse = stride * filter_height * Nos_of_channel * bw_ifmap; % in bit; additional ifmap data to load to process each next slot of zx-y2
        else
            ifmap_z_reverse = DRAM_ifmap_data_z; % in bit, ifmap data to load for each next slot of zx-y2
        end
        
        cycle_to_load_data_z_reverse = ceil((ceil(ifmap_z_reverse/DRAM_block_size) * DRAM_block_size)/BW_DRAM);
        DRAM_stall_z_reverse = max(0, cycle_to_load_data_z_reverse - cycle_count_ideal_z);  %stall from ifmap data load in each slot of next zx-y2
        DRAM_stall_zx_reverse = DRAM_stall_z_reverse * ofmap_width_effective;   %stall from ifmap data load for the next zx-y2 block
              
        DRAM_stall_zxy = DRAM_stall_zx + DRAM_stall_zx_reverse * (ofmap_height_effective - 1);
        % First term: for the first zx-y1 slot
        % Next term: For each next zx-y2 slot

        % Computation of slack cycles
        % implemented the similar logic for x-z-y-F branch
        slack_cycle_z_reverse = max(0, cycle_count_ideal_z - cycle_to_load_data_z_reverse); %slack cycle from each slot of next zx-y2
        slack_cycle_zx_reverse = slack_cycle_z_reverse * ofmap_width_effective;
        slack_cycle_zxy = slack_cycle_zx + slack_cycle_zx_reverse * (ofmap_height_effective - 1);            
    end    
else
    disp("filter-SRAM cannot fit RSCK volume of filter, Hence this branch will never be competitive")
end

% DOING SLACK CYCLE COMPUTATION FOR THIS BRANCH after 3rd depth as well
% A term to decide whether to stop keeping track of slack-cycle after 3rd depth or not
if (SRAM_filter >= filter_height * filter_width * Nos_of_channel * unitvol_nosof_3D_filter * bw_filter)
    slack_decision_factor = DRAM_stall_zx/cycle_count_ideal_zxy; 
    % Since DRAM stall is not implemented when filter-SRAM cannot fit RSCK filter, DRAM_stall_zx does not exist when the above if condition is not met
end

%% Tree process along z->x->y->F direction of ifmap

% filter memeory requirement check is not necessary
% ifmap memory requirement check is also not necessary, the check at depth-1 is enough
% psum memory check is not necessary. ofmaps are formed and can be moved to DRAM
% no effective triggers

%%%%%% cycle count calculation
%the same process for the first set of 3D filters is repeated for all the next set of 3D filters
% The ceil works becasue, when not integer multiple, in the last set of 3D filters some PE columns will remain unused
cycle_count_ideal_zxyF = cycle_count_ideal_zxy * ceil(Nos_of_filter/unitvol_nosof_3D_filter);   %without the initial offset cycles
cycle_count_zxyF = cycle_count_ideal_zxyF + initial_offset_cycle;

%%%%% Access cost for filter data
% entirely new sets of 3D filter is being processed, so reverese formatting is not applicable on filter data
% For the new sets of 3D filters, same process as zxy is repeated. (for the last set, may be a little different if not integer multiple due to 
% the fact that less #of filter channels may need to be discarded, but ignoring that, wont make much of a difference I guess.)
Access_DRAM_filter_zxyF = Access_DRAM_filter_zxy * (Nos_of_filter/unitvol_nosof_3D_filter);
% The same access pattern as zxy repeats for all sets of 3D filters
Access_SRAM_filter_zxyF = Access_SRAM_filter_zxy * (Nos_of_filter/unitvol_nosof_3D_filter);


%%%%% Access cost for ifmap data
if (SRAM_ifmap >= (ifmap_width_effective * ifmap_height_effective * Nos_of_channel * bw_ifmap))
    disp("full 3D volume of ifmap fits in ifmap-SRAM")
    gamma_DRAM_ifmap_zxyF = 1;  %ifmap-SRAM has enough memory to fit the full ifmap volume, hence each ifmap data loaded from DRAM once
    Access_DRAM_ifmap_zxyF = ifmap_width_effective * ifmap_height_effective * Nos_of_channel;
else
    % Implementing reverse formatting for the ifmap data while processing towards F-direction  
    % the access pattern for the first set of 3D filters is repeated for all the next set of 3D filters
    if(ceil(Nos_of_filter/unitvol_nosof_3D_filter) == 1)  %(This part seems like not necessary since -1 in the last line inside 'else' will take care of this, THINK)
        %if all filters are covered in the first set of 3D filter pass then there is no previously loaded ifmap zR templates to subtract
        Access_DRAM_ifmap_zxyF = Access_DRAM_ifmap_zxy;
    else
        % Some 2D ifmap zR templates will already be in ifmap-SRAM after processing the first set of 3D filter
        ifmap_zR_template_to_fit_SRAM = floor(SRAM_ifmap / (Nos_of_channel * filter_height * bw_ifmap)); %#of ifmap zR template fits in ifmap SRAM
        % In the following lines: ifmap_already_in_SRAM is the amount of ifmap data which is already in the ifmap-SRAM during the processing of the previous set of 3D-filter
        % Therefore, no need to load it since performing reverese formatting       
        first_3Dset_access = Access_DRAM_ifmap_zxy;         % access cost for the first set of 3D filters
        %The following line "ifmap_already_in_SRAM" does not have the ifmap stride overlap issue because the divisor in the calculation
        %of ifmap_zR_template_to_fit_SRAM is used as multiplicator here
        ifmap_already_in_SRAM = (Nos_of_channel * filter_height * ifmap_zR_template_to_fit_SRAM);
        next_3Dset_access = (Access_DRAM_ifmap_zxy - ifmap_already_in_SRAM);    % access cost for each next set of 3D filters
        % ceil comes here because even if there is less #of 3D filter in the last pass, ifmap data need to be loaded
        Access_DRAM_ifmap_zxyF = first_3Dset_access + next_3Dset_access * ceil(Nos_of_filter/unitvol_nosof_3D_filter - 1);
    end    
end
% The same SRAM access pattern for zxy is repeated. The reason ceil comes here is that, in the last pass if there is less filter than "unitvol_nosof_3D_filter", 
% some PE column will reamin unused, but ifmap channels will need to be loaded from SRAM to RF and broadcasted to the PE columns which are being used.
Access_SRAM_ifmap_zxyF = Access_SRAM_ifmap_zxy * ceil(Nos_of_filter/unitvol_nosof_3D_filter);


%%%% Access cost for psum data
% from the zxy direction unitvol_nosof_3D_filter ofmap planes form, same process repeats for the remaining ofmap planes
Access_SRAM_psum_zxyF = Access_SRAM_psum_zxy * (Nos_of_filter/unitvol_nosof_3D_filter);
% ceil is not required here because in the last pass if there is less filter than "unitvol_nosof_3D_filter", there will be no psum out from the unused PE columns
% with z at depth-1, this is basically the ofmap volume.

%This is the access energy up to this depth of the decision
Weighted_Access_zxyF = (Access_DRAM_filter_zxyF * bw_filter + Access_DRAM_ifmap_zxyF * bw_ifmap) * E_DRAM_to_SRAM + ...
    (Access_SRAM_filter_zxyF * bw_filter + Access_SRAM_ifmap_zxyF * bw_ifmap + Access_SRAM_psum_zxyF * bw_psum) * E_SRAM_to_RF;  % Access energy (in joule)
SRAM_Access_zxyF = (Access_SRAM_filter_zxyF * bw_filter + Access_SRAM_ifmap_zxyF * bw_ifmap + Access_SRAM_psum_zxyF * bw_psum) * E_SRAM_to_RF;
DRAM_Access_zxyF = (Access_DRAM_filter_zxyF * bw_filter + Access_DRAM_ifmap_zxyF * bw_ifmap) * E_DRAM_to_SRAM; 

%%%%%%%%%%% Calculating DRAM bandwidth induced stall count
if (SRAM_filter >= filter_height * filter_width * Nos_of_channel * unitvol_nosof_3D_filter * bw_filter)
    if (DRAM_stall_zxy <= 0) 
        % This is CASE-1. This means each zxy block can do the computation without any stall. So the same zero stall repeats for all 3D filter sets
        DRAM_stall_zxyF = 0;        
    else
        % This is one of the subcases of CASE-2
        if (SRAM_filter <= filter_height * filter_width * Nos_of_channel * unitvol_nosof_3D_filter * bw_filter) || (slack_cycle_zxy <= 0)
            % There is either no slack cycle or filter-SRAM can fit exactly RSCK volume of filter
            % In this case stall from filter data load at the first element of zxy cannot be overlapped
            % the same stall from first zxy repeats for all the set of 3D filters
            DRAM_stall_zxyF = ceil(DRAM_stall_zxy * (Nos_of_filter/unitvol_nosof_3D_filter));
            %ceil not required for Nos_of_filter/unitvol_nosof_3D_filter.
            % ceil is used for the full block of xzyF

        else
            % only considering to overlap filter data at the first element since the source of this huge stall is the filter data.
            % in this branch, mainly the source of stall is the filter data, ifmap data has much less impact
            
            filter_data_z = filter_height * filter_width * Nos_of_channel * unitvol_nosof_3D_filter * bw_filter; % in bit, filter data for the next set of 3D filters
            Extra_filter_data_to_fit_SRAM = SRAM_filter - filter_data_z; %in bit, filter data exceeding RSCK volume which can fit in filter-SRAM
            
            % data loaded is the minimum between the required filter data for the next set of 3D filters and data which fits in filter-SRAM
            loaded_filter_data = min(filter_data_z,Extra_filter_data_to_fit_SRAM); 
            
            cycle_to_load_extra_filter_data = ceil((ceil(loaded_filter_data/DRAM_block_size) * DRAM_block_size)/BW_DRAM);
            
            % Number of cycles which is overlapped, hence no need to incur stall for these many cycles. This is the minimum between whatever slack cycles
            % we have from the previous zxy block and the nos of cycle required to load the data
            cycle_saved_from_stall = min(slack_cycle_zxy,cycle_to_load_extra_filter_data); 
            
            % The amount of stall from the next zxy-block is reduced by the number of cycles saved. Getting a negative value from
            % "DRAM_stall_zxy - cycle_saved_from_stall" means no stall from ifmap data as well at the first element of next zxy-block.
            % Hence the max() will take 0 stall.
            DRAM_stall_zxy_next = max(0, DRAM_stall_zxy - cycle_saved_from_stall);
            
            DRAM_stall_zxyF = DRAM_stall_zxy + ceil(DRAM_stall_zxy_next * (Nos_of_filter/unitvol_nosof_3D_filter - 1));
            % first term: for the first set of 3D filter set
            % second term: for all the next set of 3D filter set
        end
    end   
else
    disp("filter-SRAM cannot fit RSCK volume of filter, Hence this branch will never be competitive")
end

% cycle_count_ideal_zxyF;
% DRAM_stall_zxyF;
% cycle_count_zxyF_with_stall = cycle_count_zxyF + DRAM_stall_zxyF;

%% Parameters to return
% In this brach no effective can come from the psum-SRAM storage limitation since Z is the first branch

%DRAM stall to write the ofmap value back to DRAM
ofmap_data_bits = ofmap_height_effective * ofmap_width_effective * Nos_of_filter * bw_ofmap; %in bits
ofmap_data_to_write = ceil(ofmap_data_bits/DRAM_block_size) * DRAM_block_size;
cycle_to_write_ofmap = ceil(ofmap_data_to_write/BW_DRAM); % nos of stall to write the ofmap back to DRAM [WILL ADD THE OPTION FOR POOLING HERE]
Access_DRAM_ofmap_zxyF = ofmap_height_effective * ofmap_width_effective * Nos_of_filter; % in nos of element, not bit

if (SRAM_filter < filter_height * filter_width * Nos_of_channel * unitvol_nosof_3D_filter * bw_filter)
    disp ("filter-SRAM cannot fit RSCK volume of filter, Hence did not impelment DRAM_stall_zxyF and setting it to Inf")
    DRAM_stall_zxyF = Inf; % just putting an infinite value so that during total cycle comparison this branch do not get selected
end

%XY_effective = [ofmap_width_effective ofmap_height_effective ifmap_width_effective ifmap_height_effective];
cycle_count = [cycle_count_ideal_zxyF DRAM_stall_zxyF cycle_to_write_ofmap];
SRAM_Access = [Access_SRAM_filter_zxyF Access_SRAM_ifmap_zxyF Access_SRAM_psum_zxyF];  % in Nos of element, not bit
DRAM_Access = [Access_DRAM_filter_zxyF Access_DRAM_ifmap_zxyF Access_DRAM_ofmap_zxyF]; % in Nos of element, not bit

end

