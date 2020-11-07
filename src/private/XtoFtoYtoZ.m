function [XFY_effective, cycle_count, SRAM_Access, DRAM_Access] = XtoFtoYtoZ(Layer_param, Hardware_param, Tech_param)

%X: Direction along the width of ofmap
%Y: Direction along the height of ofmap
%Z: Direction along the channels of ifmap
%F: Direction along the channels of ofmap = Direction along the 3D filters

% This function performs the computation of the input layer using the XFYZ branch

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



%% Tree process along x-direction of ifmap
%SRAM_filter memory requirement check;
SRAM_req_filter_x = unit_vol_filter * filter_height * bw_filter;      %in bit, Assuming full filter height needs to fit
if (SRAM_req_filter_x <= SRAM_filter)
    disp("Full x-direction possible for filter template")
else
    Min_SRAM_filter_x = (SRAM_req_filter_x)/(8 * 1024); % Minimum SRAM requirement for filter data in kB to proceed along X-direction
    disp("ALERT:Filter-SRAM is too small, increase it") %Assuming a minimum size for filter SRAM
    fprintf("Minimum SRAM memory for filter data for this layer is %f kB\n",Min_SRAM_filter_x)
    % returning inf valus for the return parameters so that this branch does not get selected and sweep fucntion continue to run for the next hardware point
    XFY_effective = [ofmap_width ofmap_height ifmap_width ifmap_height Nos_of_filter];  % this will prevent calling the asymmetric subvolumes
    cycle_count = [inf inf inf];
    SRAM_Access = [inf inf inf];  
    DRAM_Access = [inf inf inf];
    return
end


%SRAM_psum memory requirement
Generated_psum_x = ((ifmap_width - filter_width)/stride + 1) * unitvol_nosof_3D_filter;
SRAM_req_psum_x = Generated_psum_x * bw_psum;   % in bit
if (SRAM_req_psum_x <= SRAM_psum)
    disp("Full X-direction possible by psum-SRAM")
    ifmap_width_effective = ifmap_width;
    ofmap_width_effective = ofmap_width;
    Generated_psum_effective_x = Generated_psum_x;
    flag_Xeff = 0;     % X_effective flag to indicate that psum-SRAM is not full at this level
else
    disp("Partial x-direction is being processed due to psum-SRAM size limitation")
    % calculating the ofmap width to fill the psum-SRAM
    ofmap_width_effective = floor(SRAM_psum/(unitvol_nosof_3D_filter * bw_psum));  %extent of partial ofmap width which can be processed at a time by psum-SRAM
    ifmap_width_effective = (ofmap_width_effective - 1) * stride + filter_width;
    Generated_psum_effective_x = ofmap_width_effective * unitvol_nosof_3D_filter;
    flag_Xeff = 1;      % X_effective flag to indicate that psum-SRAM is full at this level
end


%SRAM_ifmap memory requirement check
SRAM_req_ifmap_x = ifmap_width_effective * filter_height * unitvol_nosof_channel * bw_ifmap;   % in bit, SRAM requirement for ifmap to process the full x direction
if (SRAM_req_ifmap_x <= SRAM_ifmap)
    disp("X-direction obtained from psum-SRAM is possible for ifmap-SRAM")
else
    % Applying a lower bound on ifmap-SRAM requirement
    Min_SRAM_ifmap_x = (SRAM_req_ifmap_x)/(8 * 1024); % Minimum SRAM requirement for ifmap data in kB to proceed along X-direction
    disp("ALERT:ifmap-SRAM is too small, increase it")
    fprintf("Minimum SRAM memory for ifmap data for this layer is %f kB\n",Min_SRAM_ifmap_x)
    % returning inf valus for the return parameters so that this branch does not get selected and sweep fucntion continue to run for the next hardware point
    XFY_effective = [ofmap_width ofmap_height ifmap_width ifmap_height Nos_of_filter];  % this will prevent calling the asymmetric subvolumes
    cycle_count = [inf inf inf];
    SRAM_Access = [inf inf inf];  
    DRAM_Access = [inf inf inf];
    return
end


% Will keep track of "Weighted access" & "Cycle count" for each decision
initial_offset_cycle = (used_MAC_per_column - 1);  %After first (used_MAC_per_column - 1) cycle, the pipeline is full
cycle_count_x = (((ifmap_width_effective - filter_width)/stride + 1) * filter_height) + initial_offset_cycle; 
cycle_count_ideal_x = ((ifmap_width_effective - filter_width)/stride + 1) * filter_height;  %ideal cycle count to process full x-direction of ifmap omitting the first cycles to fill the pipeline
 

%%%%%% Access for filter data
gamma_DRAM_filter_x = 1;
Access_DRAM_filter_x = gamma_DRAM_filter_x * unit_vol_filter * filter_height;  % in element, not bit

gamma_SRAM_filter_x = 1;
Access_SRAM_filter_x = gamma_SRAM_filter_x * unit_vol_filter * filter_height; % each filter row is loaded from SRAM to RF only once


%%%%%% Access for ifmap data
gamma_DRAM_ifmap_x = 1;
Access_DRAM_ifmap_x = gamma_DRAM_ifmap_x * filter_height * ifmap_width_effective * unitvol_nosof_channel;  % in element

%Due to assymetric processing when filter_width = ifmap_width, the two equations inside the if/else condition are same
if (filter_width ~= 1) && (filter_height ~= 1) %for conv layers excluding 1*1 filter-sized conv layers
    if (stride == 1)
        gamma_SRAM_ifmap_x = 1;       % when stride = 1, ecah ifmap element is read from SRAM only once due to the ifmap slidding along the x-direction
        Access_SRAM_ifmap_x = gamma_SRAM_ifmap_x * filter_height * ifmap_width_effective * unitvol_nosof_channel; % in element
    end
    if (stride > 1)
       Access_SRAM_ifmap_x = cycle_count_ideal_x * filter_width * unitvol_nosof_channel; % in element
       %since in each cycle one psum out, effectively in each cycle nos of ifmap loaded is equal to "filter_width"
       %Basically in ecah cycle "unit_vol_ifmap" new ifmap is loaded from SRAM to RF
    end
end

% For 1*1 filter-sized conv layer. Basically both the equations inside the above if/else condition works for the 1*1 conv layer. 
% With R=S=1 and E = G, the two equations are same. Despite writing seperately just to be more clear about the 1*1 conv layer.
if (filter_width == 1) && (filter_height == 1) && ((stride == 1)||(stride == 2)) %for 1*1 conv with stride 2 in ResNet50
    %in ecah cycle "unit_vol_ifmap" new ifmap is loaded from SRAM to RF
    Access_SRAM_ifmap_x = cycle_count_ideal_x * filter_width * unitvol_nosof_channel; % in element
end


% Access for psum data, psum does not go to DRAM, so only SRAM access
gamma_SRAM_psum_x = 1;   
Access_SRAM_psum_x = gamma_SRAM_psum_x * (2 * Generated_psum_effective_x * (filter_height - 1) + Generated_psum_effective_x);  % in element, not bit
% First row is write only, the reamining rows of the filter template is both read-write


%This is the access energy up to a specefic stage of the decision
Weighted_Access_x = (Access_DRAM_filter_x * bw_filter + Access_DRAM_ifmap_x * bw_ifmap) * E_DRAM_to_SRAM + ...
    (Access_SRAM_filter_x * bw_filter + Access_SRAM_ifmap_x * bw_ifmap + Access_SRAM_psum_x * bw_psum) * E_SRAM_to_RF;
SRAM_Access_x = (Access_SRAM_filter_x * bw_filter + Access_SRAM_ifmap_x * bw_ifmap + Access_SRAM_psum_x * bw_psum) * E_SRAM_to_RF;

%%%%%%%%%%% Calculating DRAM bandwidth induced stall count
DRAM_filter_data_x = filter_height * filter_width * unitvol_nosof_channel * unitvol_nosof_3D_filter * bw_filter;
DRAM_ifmap_data_x = filter_height * ifmap_width_effective * unitvol_nosof_channel * bw_ifmap;
DRAM_data_x = DRAM_filter_data_x + DRAM_ifmap_data_x; % in bit, this is the number of data loaded from DRAM to process the first ofmap row
%cycle_to_load_DRAM_data_x = ceil(DRAM_data_x/BW_DRAM)

DRAM_data_loaded_x = ceil(DRAM_data_x/DRAM_block_size) * DRAM_block_size;   % The data loaded from DRAM is an integer multiple of the DRAM block size
cycle_to_load_DRAM_data_x = ceil(DRAM_data_loaded_x/BW_DRAM);               % #of cycle to load the data from DRAM

if cycle_to_load_DRAM_data_x > cycle_count_ideal_x
    DRAM_stall_x = cycle_to_load_DRAM_data_x - cycle_count_ideal_x;
else
    DRAM_stall_x = 0;
end


% Calculating #of slack cycles from this stage
if (DRAM_stall_x == 0)
    slack_cycle_x = cycle_count_ideal_x - cycle_to_load_DRAM_data_x;  %#of extra cycles during computation when data can be loaded from DRAM
else
    slack_cycle_x = 0;    % No extra cycles for data load since data load time is higher than compute time
end


%% Tree process along x->F direction of ifmap

%SRAM_psum memory requirement
% How many filter I can cover in the F-direction will depend on psum-SRAM storage since psum does not go to DRAM
Generated_psum_xF = ofmap_width_effective * Nos_of_filter;
SRAM_req_psum_xF = Generated_psum_xF * bw_psum;   % in bit
if (SRAM_req_psum_xF <= SRAM_psum)
    disp("full F-direction possible for psum-SRAM")
    Generated_psum_effective_xF = Generated_psum_xF;
    F_effective = Nos_of_filter;
    flag_Feff = 0;
else
    disp("Partial F-direction is being processed due to psum-SRAM storage requirement")
    F_effective_to_fit = floor(SRAM_psum / (ofmap_width_effective * bw_psum));
    % for the sake of regular processing pattern, F_effective is chosen to be an integer multiple of unitvol_nosof_3D_filter.
    F_effective = floor(F_effective_to_fit/unitvol_nosof_3D_filter) * unitvol_nosof_3D_filter;
    Generated_psum_effective_xF = ofmap_width_effective * F_effective;
    flag_Feff = 1;   
end


%SRAM_filter memory requirement check; This check is not necessary. 
%While moving towards F-direction, new filters can be loaded and previous filters can be discarded if necessary.
SRAM_req_filter_xF = filter_width * filter_height * unitvol_nosof_channel * F_effective * bw_filter;    
if (SRAM_req_filter_xF <= SRAM_filter)
    disp("F-effective possible for filter template")   
end

%SRAM_ifmap memory requirement check is not necessary. The same check at depth-1 is enough
%Basically, while moving towards the F-direction, same ifmap loaded for x-direction reamins stationary in SRAM

%%%% Cycle count calculation
% The same process of x will be repeated "(F_effective/unitvol_nosof_3D_filter)" times to cover up to F-effective of F-direction
% The ceil works becasue, when not integer multiple, in the last set of 3D filters some PE columns will remain unused
cycle_count_ideal_xF = ceil(F_effective/unitvol_nosof_3D_filter) * cycle_count_ideal_x;  %without the initial offset cycles
cycle_count_xF = cycle_count_ideal_xF + initial_offset_cycle;


%%%%% Access for filter data
%Loading the new 3D filters from DRAM and including the previous load, so far each filter data loaded from DRAM once
gamma_DRAM_filter_xF = 1;
Access_DRAM_filter_xF = filter_height * filter_width * unitvol_nosof_channel * F_effective;
% The access pattern of the x-direction is repeated (F_effective/unitvol_nosof_3D_filter) times
Access_SRAM_filter_xF = Access_SRAM_filter_x * (F_effective/unitvol_nosof_3D_filter);


%%%%%%% Access for ifmap data
% No new ifmap need to be loaded, ifmap data loaded from DRAM in the previous X-direction will be reused
gamma_DRAM_ifmap_xF = 0;
Access_DRAM_ifmap_xF = Access_DRAM_ifmap_x;   % in element
%The access pattern of the x-direction is repeated (F_effective/unitvol_nosof_3D_filter) times. 
%The reason ceil comes here is that, in the last pass if there is less filter than "unitvol_nosof_3D_filter", 
%some PE column will reamin unused, but ifmap channels will need to be loaded from SRAM to RF and broadcasted to the PE columns which are being used.
Access_SRAM_ifmap_xF = Access_SRAM_ifmap_x * ceil(F_effective/unitvol_nosof_3D_filter);

%%%%% Access for psum data, psum not going to DRAM
%same pattern of psum access as x-direction is repeated (F_effective/unitvol_nosof_3D_filter) times
Access_SRAM_psum_xF = Access_SRAM_psum_x * (F_effective/unitvol_nosof_3D_filter);
% ceil is not required here because in the last pass if there is less filter than "unitvol_nosof_3D_filter", there will be no psum out from the unused PE columns


%This is the access energy up to this depth of the decision
Weighted_Access_xF = (Access_DRAM_filter_xF * bw_filter + Access_DRAM_ifmap_xF * bw_ifmap) * E_DRAM_to_SRAM + ...
    (Access_SRAM_filter_xF * bw_filter + Access_SRAM_ifmap_xF * bw_ifmap + Access_SRAM_psum_xF * bw_psum) * E_SRAM_to_RF;
SRAM_Access_xF = (Access_SRAM_filter_xF * bw_filter + Access_SRAM_ifmap_xF * bw_ifmap + Access_SRAM_psum_xF * bw_psum) * E_SRAM_to_RF;
DRAM_Access_xF = (Access_DRAM_filter_xF * bw_filter + Access_DRAM_ifmap_xF * bw_ifmap) * E_DRAM_to_SRAM;


%%%%%%%%%%% Calculating DRAM bandwidth induced stall count
New_loaded_ifmap_DRAM_xF = 0;   % when moving x->F direction, the previosuly loaded ifmap data is reused
New_loaded_filter_PFS = filter_height * filter_width * unitvol_nosof_channel * unitvol_nosof_3D_filter * bw_filter; % in bit; this gives the filter data
                                                                                                                    % loaded to process each 3D filter set
DRAM_data_PFS = New_loaded_ifmap_DRAM_xF + New_loaded_filter_PFS; %in bit; total new data need to be loaded from DRAM per 3D filter set (PFS)
DRAM_data_loaded_PFS = ceil(DRAM_data_PFS/DRAM_block_size) * DRAM_block_size;
cycle_to_load_DRAM_data_PFS = ceil(DRAM_data_loaded_PFS/BW_DRAM); % cycle count to load data per 3D filter set

cycle_to_compute_PFS = cycle_count_ideal_x;                % cycle count to process the computation for each 3D filter set
if (cycle_to_load_DRAM_data_PFS > cycle_to_compute_PFS)
    DRAM_stall_PFS = cycle_to_load_DRAM_data_PFS - cycle_to_compute_PFS;  %DRAM stall per ofmap row
else
    DRAM_stall_PFS = 0;
end

DRAM_stall_F = ceil(DRAM_stall_PFS * (F_effective/unitvol_nosof_3D_filter - 1)); %excluding the first 3D filter set since it is accounted during x-direction
% ceil not used for F_effective/unitvol_nosof_3D_filter. the argument is that, in the last pass, if there is less #of unitvol_nosof_3D_filter, 
% it will take less time to load that data. ceil is used at the full block of xF to make sure the block processing stall count is integer.
% omitting the slack_cycle variable here since for r1F1 we need to load more data (both filter & ifmap) than r1F2 ... (only filter).
% So, if there is no stall at x, there will be no stall at y.

% Calculating #of slack cycles from this stage
if (DRAM_stall_F == 0)
    slack_cycle_F = (cycle_to_compute_PFS - cycle_to_load_DRAM_data_PFS) * ceil(F_effective/unitvol_nosof_3D_filter - 1); %#of extra cycles during computation when data can be loaded from DRAM
    % ceil is used here, The reason is in the last pass if there is less filter than unitvol_nosof_3D_filter, it will take less time to load
    % that set of filter data. Hence definitely in that last pass #of_slack_cycle will be higher than the standard (i.e., the first term in the equation).
    % So making the ceil will multiply with 1 for the last pass in stead of some fraction (i.e, <1) and hence closer to accurate.
else
    slack_cycle_F = 0;    % No extra cycles for data load since data load time is higher than compute time
end

DRAM_stall_xF = DRAM_stall_x + DRAM_stall_F;    % DRAM stall upto xF
slack_cycle_xF = slack_cycle_x + slack_cycle_F; % slack cycle upto xF


%% Tree process along x->F->y direction of ifmap

% filter SRAM memory requirement check is not necessary
% ifmap SRAM memory requirement check is not necessary, new ifmap rows can be loaded from DRAM

% SRAM_psum memory requirement check
% How many rows in the y-direction I can process depend on the psum-SRAM storage limit since psum does not goto DRAM
rows_tofit_psum_SRAM = floor(SRAM_psum / (ofmap_width_effective * F_effective * bw_psum)); % #of ofmap rows from each F_effective ofmap plane which fits in psum-SRAM
if (rows_tofit_psum_SRAM >= ofmap_height)
    disp("Full Y-direction possible by psum-SRAM")
    ofmap_height_effective = ofmap_height;
    ifmap_height_effective = ifmap_height;
    flag_Yeff = 0;
else
    disp("Partial Y-direction is being processed due to psum-SRAM size limitation")
    ofmap_height_effective = rows_tofit_psum_SRAM;
    ifmap_height_effective = (ofmap_height_effective - 1) * stride + filter_height;
    flag_Yeff = 1;
end

%%%%% cycle count calculation
% For each ofmap row, #of cycle required for cycle_count_ideal_xF is repeated
cycle_count_ideal_xFy = cycle_count_ideal_xF * ofmap_height_effective;    %without the initial offset cycles
cycle_count_xFy = cycle_count_ideal_xFy + initial_offset_cycle;


%%%%% Access for filter data
filter_channels_to_fit_SRAM = floor(SRAM_filter/(filter_height * filter_width * bw_filter)); % #of 2D filter planes which fits into filter-SRAM
Total_uv_filter_planes = F_effective * unitvol_nosof_channel;    % #of total 2D filter planes in the unitvol_nosof_channel of all the F_effective filters
if (filter_channels_to_fit_SRAM >= Total_uv_filter_planes)
    disp ("filter-SRAM has enough storage to process xFy direction by loading filter data once from DRAM")
    gamma_DRAM_filter_xFy = 1; % ecah filter data loaded from DRAM once
    Access_DRAM_filter_xFy = filter_width * filter_height * unitvol_nosof_channel * F_effective; % in element
else
    disp("filter-SRAM is small and requires multiple DRAM load to process xFy direction")
    Nos_of_discarded_filter_planes = Total_uv_filter_planes - filter_channels_to_fit_SRAM; % #of 2D filter planes which are discarded from SRAM
    % Implementing reverse order processing format for DRAM access
    % During the processing of first xz slice of ofmap, all filter are loaded once from DRAM
    First_xz_ofmap_Access = filter_width * filter_height * unitvol_nosof_channel * F_effective;
    % For all the next xz slice of ofmap, the amount of discarded filter planes need to be loaded.
    Next_xz_ofmap_Access = filter_width * filter_height * Nos_of_discarded_filter_planes;
    Access_DRAM_filter_xFy = First_xz_ofmap_Access + Next_xz_ofmap_Access * (ofmap_height_effective - 1);    
end
% The same SRAM access pattern for the first ofmap xz plane is repeated for all the next ofmap xz planes
Access_SRAM_filter_xFy = Access_SRAM_filter_xF * ofmap_height_effective;


%%%%%%% Access for ifmap data
% ifmap data can be loaded only once from DRAM. At one time it is enough to fit "SRAM_req_ifmap_x" ifmap data in ifmap-SRAM. While moving towards the
% Y-direction, new ifmap rows are loaded and if necessary, previous ifmap rows can be discarded.
gamma_DRAM_ifmap_xFy = 1;
Access_DRAM_ifmap_xFy = ifmap_width_effective * ifmap_height_effective * unitvol_nosof_channel;  % in element
%For ecah ofmap xz planes, same amount of SRAM ifmap access is repeated
Access_SRAM_ifmap_xFy = Access_SRAM_ifmap_xF * ofmap_height_effective;


%%%%%%% Access for psum data
%%%%% Access for psum data, psum not going to DRAM
%same pattern of psum access happen for each ofmap xz plane
Access_SRAM_psum_xFy = Access_SRAM_psum_xF * ofmap_height_effective;


%This is the access energy up to this depth of the decision
Weighted_Access_xFy = (Access_DRAM_filter_xFy * bw_filter + Access_DRAM_ifmap_xFy * bw_ifmap) * E_DRAM_to_SRAM + ...
    (Access_SRAM_filter_xFy * bw_filter + Access_SRAM_ifmap_xFy * bw_ifmap + Access_SRAM_psum_xFy * bw_psum) * E_SRAM_to_RF;  % Access energy (in joule)
SRAM_Access_xFy = (Access_SRAM_filter_xFy * bw_filter + Access_SRAM_ifmap_xFy * bw_ifmap + Access_SRAM_psum_xFy * bw_psum) * E_SRAM_to_RF;
DRAM_Access_xFy = (Access_DRAM_filter_xFy * bw_filter + Access_DRAM_ifmap_xFy * bw_ifmap) * E_DRAM_to_SRAM; 

%%%%%%%%%%% Calculating DRAM bandwidth induced stall count
if (slack_cycle_xF <= 0)  
    disp("CASE-1: 1 1")
    % CASE-1: There is no slack cycle from xF   
    % Hence, each xF block on its own, so no stall overlapped
    % However the filter data already in SRAM dont need to be loaded in the reverse order formatting,
    % Therefore, the stall from filter wont happen from those filter data
    
    % Can consider seperate stall from ifmap and filter data only at x-direction in this depth of this branch. Due to reverse order formatting 
    % at the first slot of xF-y2, only ifmap data will be loaded, no filter data need to be loaded since its already in SRAM from the processing of previous xF block
    % stall from ifmap data only at x-direction
    New_loaded_ifmap_POR = stride * ifmap_width_effective * unitvol_nosof_channel * bw_ifmap; % in bit; new ifmap data need to load to process each next ofmap row
    ifmap_data_loaded_POR = ceil(New_loaded_ifmap_POR/DRAM_block_size) * DRAM_block_size;
    cycle_to_load_ifmap_data_POR = ceil(ifmap_data_loaded_POR/BW_DRAM); % cycle count to load data per ofmap row
    DRAM_stall_ifmap_x = max(0, cycle_to_load_ifmap_data_POR - cycle_count_ideal_x);              % stall from ifmap data only at x-direction
        
    % stall from filter data only at x-direction
    New_loaded_filter_PFS = filter_height * filter_width * unitvol_nosof_channel * unitvol_nosof_3D_filter * bw_filter; % in bit; filter data loaded per 3D filter set                                                                                                              % loaded to process each 3D filter set
    filter_data_loaded_PFS = ceil(New_loaded_filter_PFS/DRAM_block_size) * DRAM_block_size;
    cycle_to_load_filter_data_PFS = ceil(filter_data_loaded_PFS/BW_DRAM); % cycle count to load filter data per 3D filter set
    DRAM_stall_filter_x = max(0, cycle_to_load_filter_data_PFS - cycle_count_ideal_x);              % stall from filter data only at x-direction

    % For simplicity using 3D filter set. For these #of 3D filter set we do not need to load the filter data
    % The implication of this is that, even though due to reverse order formatting some partial filter data from a 3D filter set is not being loaded from DRAM,
    % the systolic array is just waiting. The is the cost in performance we pay to have simplified control.
    filter_3D_set_to_fit_SRAM = floor(SRAM_filter/(filter_height * filter_width * unitvol_nosof_channel * unitvol_nosof_3D_filter * bw_filter));  
    if (SRAM_filter >= (filter_height * filter_width * unitvol_nosof_channel * F_effective * bw_filter))
        %filter-SRAM is large enough to fit (R*S*C_uv*F_eff) volume, so no filter need to be loaded again from DRAM while moving towards Y-direction
        DRAM_stall_xF_reverse = DRAM_stall_ifmap_x;  %the ifmap data is only loaded at the first slot of xF-y    
    else
        DRAM_stall_xF_reverse = DRAM_stall_ifmap_x + ceil(DRAM_stall_filter_x * (F_effective/unitvol_nosof_3D_filter - filter_3D_set_to_fit_SRAM));
        % Here the first term gives stall only from ifmap data which is loaded for the first xF-y2 slot
        % The second term gives stall from filter-data only. For all the next xF-y2 slots only filter data is loaded. However, due to reverse order formatting,
        % "filter_3D_set_to_fit_SRAM" #of 3D filter set is already in filter-SRAM and hence for those slots no need to load any filter (i.e., any data)
        % ceil not required for F_effective/unitvol_nosof_3D_filter since if there is less filter in the last pass, it will take less #of cycle to load
    end
    DRAM_stall_xFy = DRAM_stall_xF + DRAM_stall_xF_reverse * (ofmap_height_effective - 1);
    % First term: for the first xF-y1 block, no reverse order formatting applicable, so regular DRAM_stall_xF term used
    % Second term: for all the next xF-y block, reverse order formatting applicable, so DRAM_stall_xF_reverse used  
    
elseif (DRAM_stall_xF <= 0)
    disp("CASE-2: 0 0")
    %CASE-2: Each xF block can compute itself without incorporating any stall, so no stall to complete the xFy direction
    DRAM_stall_xFy = 0;
else   
    disp("CASE-3: 1 0")
    % CASE-3: There are some slack cycles from xF, as well as some stall from xF. 
    % This means the source of stall is the (ifmap+filter) data at the first xF-y1 slot.
    % The load time of only filter is not creating the stall, it can be overlapped with computation, So, no need to worry about filter-SRAM constraint. 
    % Now, due to reverse order formatting, at the first slot of xF-y2, no filter data will be loaded since the required R*S*Cuv*Fuv filter data
    % is already in filter-SRAM (min filter-SRAM requirment). so only need to check the ifmap data for the first slot of xF-y2 for the stall count.
    % And only need to check if ifmap-SRAM memory permits to load ifmap data from DRAM during those slack cycles to overlap the stall (if any).
    
    % Assumption: if the ifmap dada needed for the next xF-y2 block processing can fully fit in ifmap-SRAM, only then we will load it, 
    % otherwise wont load any ifmap data for the next xF block processing
    
    % stall from ifmap data only at x-direction
    New_loaded_ifmap_POR = stride * ifmap_width_effective * unitvol_nosof_channel * bw_ifmap; % in bit; new ifmap data need to load to process each next ofmap row
    ifmap_data_loaded_POR = ceil(New_loaded_ifmap_POR/DRAM_block_size) * DRAM_block_size;
    cycle_to_load_ifmap_data_POR = ceil(ifmap_data_loaded_POR/BW_DRAM); % cycle count to load data per ofmap row
    DRAM_stall_ifmap_x = max(0, cycle_to_load_ifmap_data_POR - cycle_count_ideal_x);              % stall from ifmap data only at x-direction
    
    if (DRAM_stall_ifmap_x <= 0)
        DRAM_stall_y = 0;      % no stall for the next xF-y2 block
        DRAM_stall_xFy = DRAM_stall_xF + DRAM_stall_y; %first term: first xF-y1 block, second term: all next xF-y2 block
    else
        if (SRAM_ifmap >= (DRAM_ifmap_data_x + New_loaded_ifmap_POR)) %There is room in SRAM_ifmap to fit the full ifmap data needed to process the next xF-y block
            if (DRAM_stall_ifmap_x <= slack_cycle_xF) 
                disp("CASE 3.1 Full stall can be overlapped")
                DRAM_stall_y = 0;
            else
                disp("CASE 3.2 Stall can be overlapped partially")
                DRAM_stall_PxF = DRAM_stall_ifmap_x - slack_cycle_xF; % stall per xF-y block
                DRAM_stall_y = DRAM_stall_PxF * (ofmap_height_effective - 1); 
                % DRAM_stall_PxF is repeated for all the next y-row, excluding the first row, it will be incorporated later by adding DRAM_stall_xF
            end
            DRAM_stall_xFy = DRAM_stall_xF + DRAM_stall_y;
        else
            % Each xy-F block on its own, so no stall overlapped and DRAM_stall_ifmap_x is repeated for all next ofmap row except the first xF-y1 block
            DRAM_stall_xFy = DRAM_stall_xF + DRAM_stall_ifmap_x * (ofmap_height_effective - 1);
        end
    end
        
end

% Assumtion: Not doing any computation of the slack cycle after xFy, For simplicity, after xFy we assume the slack variable is zero
% keeping track of slack variable becomes complicated as we move higher depth of the tree
% This assumption is based on the fact that the ratio of stall_xF/compute_cycle_xFy is very small, 
% hence it is not worth to make the implementation very complicated for this
slack_decision_factor = DRAM_stall_xF/cycle_count_ideal_xFy;

%% Tree process along x->F->y->z direction of ifmap

% filter SRAM memory requirement check is not necessary, new channels can be loaded from DRAM
% ifmap SRAM memory requirement check is not necessary, new channels can be loaded from DRAM
% psum SRAM memory requirement check is not necessary since with z-direction no new nonreducible psum is generated than the xFy direction

%%%%%%%Calculation of cycle count  
%The same process as xFy is repeated. For each ifmap plane with "unitvol_nosof_channel", same amount of cycle count is required
cycle_count_ideal_xFyz = cycle_count_ideal_xFy * ceil(Nos_of_channel/unitvol_nosof_channel); %without the initial offset cycles
cycle_count_xFyz = cycle_count_ideal_xFyz + initial_offset_cycle;

%%%%% Access for filter data
% while moving towards the xFy directions, all the computations associated with the unitvol_nosof_channel number of filter channels 
% from F_effective number of 3D filters is done. Hence for z-direction, the same process as xFy is repeated 
Access_DRAM_filter_xFyz = Access_DRAM_filter_xFy * (Nos_of_channel/unitvol_nosof_channel);
%basically here maximal filter reuse occur if R*S*|J/S|*F_effective filter volume fits in filter-SRAM, then filter data is loaded from DRAM once

%The amount of access in xFy direction is repeated (Nos_of_channel/unitvol_nosof_channel) times to cover the full z-direction
Access_SRAM_filter_xFyz = Access_SRAM_filter_xFy * (Nos_of_channel/unitvol_nosof_channel);

%%%%% Access for ifmap data
% All computations associated with the "unitvol_nosof_channel" ifmap channels are done after loading it from DRAM once. The same is true for
% all ifmap channels while moving towards z-direction
gamma_DRAM_ifmap_xFyz = 1; % so maximal utilization of ifmap data from DRAM. each ifmap data is loaded from DRAM only once
Access_DRAM_ifmap_xFyz = ifmap_width_effective * ifmap_height_effective * Nos_of_channel;  % in element
%The amount of access in xFy direction is repeated (Nos_of_channel/unitvol_nosof_channel) times to cover the full z-direction
Access_SRAM_ifmap_xFyz = Access_SRAM_ifmap_xFy * (Nos_of_channel/unitvol_nosof_channel);


%%%%% Access for psum data
%Access_SRAM_psum_xFy = For the first set of channels, the access pattern for the F_effective number of ofmap planes is same as Access_SRAM_psum_xFy
psum_xFyz_x = 2 * Generated_psum_effective_x * filter_height;    % psum SRAM access cost while processing the x-direction for next ifmap channels
psum_xFyz_xF = psum_xFyz_x * (F_effective/unitvol_nosof_3D_filter);   % psum SRAM access cost while processing the xF-direction for next ifmap channels
psum_xFyz_xFy = psum_xFyz_xF * ofmap_height_effective;  % psum SRAM access cost while processing the xFy-direction for next ifmap channels
Access_SRAM_psum_xFyz = Access_SRAM_psum_xFy + psum_xFyz_xFy * ceil((Nos_of_channel/unitvol_nosof_channel) - 1);
% ceil comes here becasue in the last pass even if there is less channel than unitvol_nosof_channel, same nos of psum read/write is required


%This is the access energy up to this depth of the decision
Weighted_Access_xFyz = (Access_DRAM_filter_xFyz * bw_filter + Access_DRAM_ifmap_xFyz * bw_ifmap) * E_DRAM_to_SRAM + ...
    (Access_SRAM_filter_xFyz * bw_filter + Access_SRAM_ifmap_xFyz * bw_ifmap + Access_SRAM_psum_xFyz * bw_psum) * E_SRAM_to_RF;  % Access energy (in joule)
SRAM_Access_xFyz = (Access_SRAM_filter_xFyz * bw_filter + Access_SRAM_ifmap_xFyz * bw_ifmap + Access_SRAM_psum_xFyz * bw_psum) * E_SRAM_to_RF;
DRAM_Access_xFyz = (Access_DRAM_filter_xFyz * bw_filter + Access_DRAM_ifmap_xFyz * bw_ifmap) * E_DRAM_to_SRAM; 


%%%%%%%%%%% Calculating DRAM bandwidth induced stall count
%The same process of xFy repeats for all the next set of channels.
DRAM_stall_xFyz = ceil(DRAM_stall_xFy * (Nos_of_channel/unitvol_nosof_channel));
% ceil not used for (Nos_of_channel/unitvol_nosof_channel), the argument is that, in the last pass, if there is less #of unitvol_nosof_channel,
% it will take less time to load that data. ceil is used for the full xFyz block


cycle_count_ideal_xFyz;
DRAM_stall_xFyz;
cycle_count_xFyz_with_stall = cycle_count_xFyz + DRAM_stall_xFyz;


%% Incorporating the calculation when some directions are processed partially
% In this brach all three effective can come from the psum-SRAM storage limitation: X-effective, F-effective, and Y-effective
% Z_effective never comes from psum-SRAM since Z-direction kills psum

%DRAM stall to write the ofmap value back to DRAM
ofmap_data_bits = ofmap_height_effective * ofmap_width_effective * F_effective * bw_ofmap; %in bits
ofmap_data_to_write = ceil(ofmap_data_bits/DRAM_block_size) * DRAM_block_size;
cycle_to_write_ofmap = ceil(ofmap_data_to_write/BW_DRAM); % nos of stall to write the ofmap back to DRAM [WILL ADD THE OPTION FOR POOLING HERE]
Access_DRAM_ofmap_xFyz = ofmap_height_effective * ofmap_width_effective * F_effective; % in nos of element, not bit


XFY_effective = [ofmap_width_effective ofmap_height_effective ifmap_width_effective ifmap_height_effective F_effective];
cycle_count = [cycle_count_ideal_xFyz DRAM_stall_xFyz cycle_to_write_ofmap];
SRAM_Access = [Access_SRAM_filter_xFyz Access_SRAM_ifmap_xFyz Access_SRAM_psum_xFyz];  % in Nos of element, not bit
DRAM_Access = [Access_DRAM_filter_xFyz Access_DRAM_ifmap_xFyz Access_DRAM_ofmap_xFyz]; % in Nos of element, not bit



end



