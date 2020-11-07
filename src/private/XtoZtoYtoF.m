function [X_effective, cycle_count, SRAM_Access, DRAM_Access] = XtoZtoYtoF(Layer_param, Hardware_param, Tech_param)

%X: Direction along the width of ofmap
%Y: Direction along the height of ofmap
%Z: Direction along the channels of ifmap
%F: Direction along the channels of ofmap = Direction along the 3D filters

% This function performs the computation of the input layer using the XZYF branch

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
    X_effective = [ofmap_width ifmap_width];  % this will prevent calling the asymmetric subvolumes
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
    X_effective = [ofmap_width ifmap_width];  % this will prevent calling the asymmetric subvolumes
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
if (filter_width ~= 1) && (filter_height ~= 1)  %for conv layers excluding 1*1 filter-sized conv layers
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


%% Tree process along x->z direction of ifmap
%%%% Cycle count calculation
cycle_count_ideal_xz = cycle_count_ideal_x * ceil(Nos_of_channel/unitvol_nosof_channel);    % required cycle count to process the full x-z direction
cycle_count_xz = cycle_count_ideal_xz + initial_offset_cycle;   % with the initial offset cycles
z_effective = Nos_of_channel;                                   % No effective triggers for Z-direction from psum-SRAM

%SRAM_filter memory requirement check; This check is not necessary at this depth
SRAM_req_filter_xz = filter_width * filter_height * z_effective * unitvol_nosof_3D_filter * bw_filter;    
if (SRAM_req_filter_xz <= SRAM_filter)
    disp("z-effective possible for filter template")   
end

%SRAM_ifmap memory requirement check; This check is not necessary at this depth
SRAM_req_ifmap_xz = ifmap_width_effective * filter_height * z_effective * bw_ifmap;   % in bit
if (SRAM_req_ifmap_xz <= SRAM_ifmap)
    disp("z-effective possible for ifmap")          
end

%SRAM_psum memory requirement
% No need to do this check since #of generated psum is same as the parent x-direction search

%%%%% Access for filter data
%Loading the new filter channels from DRAM and including the previous load
gamma_DRAM_filter_xz = 1;
Access_DRAM_filter_xz = filter_width * filter_height * z_effective * unitvol_nosof_3D_filter;   % in element
% The access pattern of the x-direction is repeated (z_effective/unitvol_nosof_channel) times
Access_SRAM_filter_xz = Access_SRAM_filter_x * (z_effective/unitvol_nosof_channel);

%%%%%%% Access for ifmap data
%Loading the new ifmap channels from DRAM and copying the previous load
gamma_DRAM_ifmap_xz = 1;
Access_DRAM_ifmap_xz = ifmap_width_effective * filter_height * z_effective;  % in element
% The access pattern of the x-direction is repeated (z_effective/unitvol_nosof_channel) times
Access_SRAM_ifmap_xz = Access_SRAM_ifmap_x * (z_effective/unitvol_nosof_channel);


%%%%% Access for psum data
Access_SRAM_psum_xz = Access_SRAM_psum_x + ceil(z_effective/unitvol_nosof_channel - 1) * (2 * Generated_psum_effective_x * filter_height);
%Access_SRAM_psum_x = for the first filter template, first row write only & next rows both read-write
%(2 * Generated_psum_effective_x * filter_height) = for the rest filter templates (except the first one) all rows have both read-write
%ceil(z_effective/unitvol_nosof_channel - 1) = The procedure repeats this number of times
% ceil comes here becasue in the last pass even if there is less channel than unitvol_nosof_channel, same nos of psum read/write is required

%This is the access energy up to this depth of the decision
Weighted_Access_xz = (Access_DRAM_filter_xz * bw_filter + Access_DRAM_ifmap_xz * bw_ifmap) * E_DRAM_to_SRAM + ...
    (Access_SRAM_filter_xz * bw_filter + Access_SRAM_ifmap_xz * bw_ifmap + Access_SRAM_psum_xz * bw_psum) * E_SRAM_to_RF;
SRAM_Access_xz = (Access_SRAM_filter_xz * bw_filter + Access_SRAM_ifmap_xz * bw_ifmap + Access_SRAM_psum_xz * bw_psum) * E_SRAM_to_RF;


%%%%%%%%%%% Calculating DRAM bandwidth induced stall count
%for x-z direction, whatever happens in x is repeated for z. Hence, the stall count and #of slack cycle is same for x & z
DRAM_stall_xz = ceil(DRAM_stall_x * (Nos_of_channel/unitvol_nosof_channel));       % DRAM stall upto xz
slack_cycle_xz = slack_cycle_x * ceil(Nos_of_channel/unitvol_nosof_channel);     % DRAM stall upto xz
%in the stall count, ceil not required for (Nos_of_channel/unitvol_nosof_channel).
%In the last pass if there is less number of channels then it will take less number of cycles to load the data
% ceil is used for the full block of xz

% ceil is used for the slack cycle count on (Nos_of_channel/unitvol_nosof_channel). 
% The reason is in the last pass if there is less filter than unitvol_nosof_channel, it will take less time to load that set of channels. 
% Hence definitely in that last pass #of_slack_cycle will be higher than the standard (i.e., the first term in the equation). 
% So making the ceil will multiply with 1 for the last pass in stead of some fraction (i.e, <1) and hence more closer to accurate.

%% Tree process along x->z->y direction of ifmap

% After Z-direction, no effective triggers from psum-SRAM
ofmap_height_effective = ofmap_height;
ifmap_height_effective = (ofmap_height_effective - 1) * stride + filter_height;


%%%% Calculation of cycle count 
%for each ofmap row, same amount of cycle is required to process a xz slice of ifmap
cycle_count_ideal_xzy = cycle_count_ideal_xz * ofmap_height_effective;
cycle_count_xzy = cycle_count_ideal_xzy + initial_offset_cycle;

%%%%% Access cost for filter data
% Reverse formatting implementation for DRAM access count
filter_channels_to_fit_SRAM = floor(SRAM_filter/(filter_height * filter_width * bw_filter)); % #of 2D filter planes which fits into filter-SRAM
Total_filter_planes_k = unitvol_nosof_3D_filter * Nos_of_channel;    % #of total 2D filter planes in the unitvol_nosof_3D_filter filter data
if (filter_channels_to_fit_SRAM >= Total_filter_planes_k)
    gamma_DRAM_filter_xzy = 1;   % filter-SRAM has enough memory, each filter data loaded from DRAM once
    Access_DRAM_filter_xzy = filter_height * filter_width * Nos_of_channel * unitvol_nosof_3D_filter;
else
    Nos_of_discarded_filter_planes = Total_filter_planes_k - filter_channels_to_fit_SRAM; % #of 2D filter planes which are discarded from SRAM
    % Implementing reverse order processing format (R1C1-->R1Cend to R2Cend--> R2C1 to R3C1-->R3Cend)
    % During the processing of first ofmap row, all filter are loaded once from DRAM
    First_ofmap_row_Access = filter_width * filter_height * Nos_of_channel * unitvol_nosof_3D_filter;
    % For all the next ofmap row, the amount of discarded channels need to be loaded.
    Next_ofmap_row_Access = filter_width * filter_height * Nos_of_discarded_filter_planes;
    Access_DRAM_filter_xzy = First_ofmap_row_Access + Next_ofmap_row_Access * (ofmap_height_effective - 1);
end

%The amount of access in xz direction is repeated ofmap_height_effective times to cover the full y-direction
Access_SRAM_filter_xzy = Access_SRAM_filter_xz * ofmap_height_effective;


%%%% Access cost for ifmap data
%Calculation on ifmap memory
ifmap_xR_template_to_fit_SRAM = floor(SRAM_ifmap / (ifmap_width_effective * filter_height * bw_ifmap)); %#of ifmap template fits in ifmap SRAM
 
if (ifmap_xR_template_to_fit_SRAM >= Nos_of_channel)
    disp("ifmap-SRAM has enough storage to process xzy direction by loading ifmap data once from DRAM")
    gamma_DRAM_ifmap_xzy = 1;  %ifmap-SRAM has enough memory, each ifmap data loaded from DRAM once
    Access_DRAM_ifmap_xzy = ifmap_width_effective * ifmap_height_effective * Nos_of_channel;
else
    disp("ifmap-SRAM does not have enough storage and multiple DRAM load required")
    % Implementation with reverse order formatting
    Discarded_ifmap_xR_channel = Nos_of_channel - ifmap_xR_template_to_fit_SRAM;   % #of ifmap templates need to be discarded while moving towards z-direction
    First_ofmap_row_Access_ifmap = filter_height * ifmap_width_effective * Nos_of_channel;  % #of ifmap data loaded during processing the first ofmap row
    % In the line below: the first term is the newly loaded ifmap data for each new ofmap row and
    % the second term is the previously discarded ifmap channels which need to be loaded again each time we process a new ofmap row
    Next_ofmap_row_Access_ifmap = (ifmap_width_effective * stride * Nos_of_channel) + (ifmap_width_effective * (filter_height - stride) * Discarded_ifmap_xR_channel);
    Access_DRAM_ifmap_xzy = First_ofmap_row_Access_ifmap + Next_ofmap_row_Access_ifmap * (ofmap_height_effective - 1);
    % Since xyz is always as good as or better than xzy I am going to prune out xzy. NO, not pruning
end   
Access_SRAM_ifmap_xzy = Access_SRAM_ifmap_xz * ofmap_height_effective;


%%%% Access cost for psum data
% from the xz direction one ofmap row forms, same process repeats for each ofmap row
Access_SRAM_psum_xzy = Access_SRAM_psum_xz * ofmap_height_effective;


%This is the access energy up to this depth of the decision
Weighted_Access_xzy = (Access_DRAM_filter_xzy * bw_filter + Access_DRAM_ifmap_xzy * bw_ifmap) * E_DRAM_to_SRAM + ...
    (Access_SRAM_filter_xzy * bw_filter + Access_SRAM_ifmap_xzy * bw_ifmap + Access_SRAM_psum_xzy * bw_psum) * E_SRAM_to_RF;
SRAM_Access_xzy = (Access_SRAM_filter_xzy * bw_filter + Access_SRAM_ifmap_xzy * bw_ifmap + Access_SRAM_psum_xzy * bw_psum) * E_SRAM_to_RF;


%%%%%%%%%%% Calculating DRAM bandwidth induced stall count
if (slack_cycle_xz <= 0)  
    disp("CASE-1: 1 1")
    % CASE-1: There is no slack cycle from xz
    % Depending on ifmap-SRAM either only additional ifmap row is loaded in each slot of next xz-y2 slot or full ifmap-x is loaded in each slot of next xz-y2.
    % For simplicity, assuming if ifmap-SRAM can fit xRz volume of ifmap then imposing only additional ifmap row to load, otherwise imposing
    % regular DRAM_ifmap_data_x ifmap to load for all next slot of xz-y2
    if (SRAM_ifmap>= ifmap_width_effective * filter_height * Nos_of_channel * bw_ifmap)
        ifmap_x_reverse = stride * ifmap_width_effective * unitvol_nosof_channel * bw_ifmap; % in bit; additional ifmap data to load to process each next slot of xz-y2
    else
        ifmap_x_reverse = DRAM_ifmap_data_x; % in bit, ifmap data to load for each next slot of xz-y2
    end
    
    % For the next xz-y2, depending on filter-SRAM, some slot will load filter data and some slot will load no filter data due to reverse order formatting.
    % For simplicity taking R*S*C_uv*F_uv volume as unit. Hence, if this amount does not fit we are imposing the full DRAM_filter_data_x load. 
    filter_3D_set_to_fit_SRAM = floor(SRAM_filter/(filter_height * filter_width * unitvol_nosof_channel * unitvol_nosof_3D_filter * bw_filter));  
    
    if (SRAM_filter >= (filter_height * filter_width * Nos_of_channel * unitvol_nosof_3D_filter * bw_filter))
        %filter-SRAM is large enough to fit R*S*C*F_uv volume of filter, so no new filter data need to be loaded to move towards Y-direction
        filter_x_reverse = 0;   % no filter data needed to load for each next slot of xz-y2
        data_x_reverse = ifmap_x_reverse + filter_x_reverse;  % total ifmap + filter data loaded for each next slot of xz-y2
        data_loaded_x_reverse = ceil(data_x_reverse/DRAM_block_size) * DRAM_block_size;
        cycle_to_load_data_x_reverse = ceil(data_loaded_x_reverse/BW_DRAM);
        DRAM_stall_x_reverse = max(0, cycle_to_load_data_x_reverse - cycle_count_ideal_x);  %DRAM stall for each next slot of xz-y2
        DRAM_stall_xz_reverse = ceil(DRAM_stall_x_reverse * (Nos_of_channel/unitvol_nosof_channel)); %DRAM stall to porcess the next ofmap row, 
        % ceil not required for (Nos_of_channel/unitvol_nosof_channel)
        
        %calculation of slack cycle
        % The slack cycles can be used only this condition is triggered, in all other cases slack cycle count does not matter
        % This slack cycle will be used at the 4th depth of the tree
        slack_cycle_x_reverse = max(0, cycle_count_ideal_x - cycle_to_load_data_x_reverse);
        slack_cycle_xz_reverse = slack_cycle_x_reverse * ceil(Nos_of_channel/unitvol_nosof_channel);
        slack_cycle_xzy = slack_cycle_xz + slack_cycle_xz_reverse * (ofmap_height_effective - 1);  % slack_cycle_xz = 0 here (see the top if condition)
               
    else
        data_x_reverse_low = ifmap_x_reverse;  %only ifmap data loaded for some slot of next xz-y2
        data_x_reverse_high = ifmap_x_reverse + DRAM_filter_data_x;  %both ifmap+filter data loaded for some slot of next xz-y2
   
        cycle_to_load_data_x_reverse_low = ceil((ceil(data_x_reverse_low/DRAM_block_size) * DRAM_block_size)/BW_DRAM);
        cycle_to_load_data_x_reverse_high = ceil((ceil(data_x_reverse_high/DRAM_block_size) * DRAM_block_size)/BW_DRAM);
        
        DRAM_stall_x_reverse_low = max(0, cycle_to_load_data_x_reverse_low - cycle_count_ideal_x);  %stall from lower data load slots
        DRAM_stall_x_reverse_high = max(0, cycle_to_load_data_x_reverse_high - cycle_count_ideal_x); %stall from higher data load slots
        
        DRAM_stall_xz_reverse = DRAM_stall_x_reverse_low * filter_3D_set_to_fit_SRAM + ...
                                ceil(DRAM_stall_x_reverse_high * (Nos_of_channel/unitvol_nosof_channel - filter_3D_set_to_fit_SRAM)); %stall to porcess the next ofmap row.
        %ceil not required for Nos_of_channel/unitvol_nosof_channel
        
        % Under this condition, slack cycle count does not matter cause filter-SRAM do not have space to hold extra filter above RSCK volume. 
        % Hence cannot use any slack cycle at 4th level of the tree
   end
    DRAM_stall_xzy = DRAM_stall_xz + DRAM_stall_xz_reverse * (ofmap_height_effective - 1);
    % First term: for the first xz-y1 slot
    % Next term: For each next xz-y2 slot (with reverse order formatting)
   
elseif (DRAM_stall_xz <= 0)
    disp("CASE-2: 0 0")
    % This means there is no stall at x, z and each xz block can compute itself without incorporating any stall, so no stall to complete the xzy direction
    DRAM_stall_xzy = 0;
    % No need to calculate slack cycle for this case since there will be no stall while moving to F-direction as well
end
    
% DOING SLACK CYCLE COMPUTATION FOR THIS BRANCH after 3rd depth as well
% A term to decide whether to stop keeping track of slack-cycle after 3rd depth or not
slack_decision_factor = DRAM_stall_xz/cycle_count_ideal_xzy;


%% Tree process along x->z->y->F direction of ifmap

% Filter SRAM memeory requirement check is not necessary since new 3D filter can be loaded from DRAM.

% ifmap SRAM memory requirement check is not necessary since the same process of xzy will be repeated for the new set of 3D filters.

% psum SRAM memory requirement check is not necessary since the full z-direction is processed for each ofmap row, the ofmaps are formed which can be 
% written back to DRAM and the same process as xzy will be repeated.


%%%%% Calculation of cycle count
% The same process of xzy will be repeated "(Nos_of_filter/unitvol_nosof_3D_filter)" times to cover the full F-direction
% The ceil works becasue, when not integer multiple, in the last set of 3D filters some PE columns will remain unused
cycle_count_ideal_xzyF = cycle_count_ideal_xzy * ceil(Nos_of_filter/unitvol_nosof_3D_filter); %without the initial offset cycles
cycle_count_xzyF = cycle_count_ideal_xzyF + initial_offset_cycle;

%%%%% Access cost for filter data
% entirely new sets of 3D filter is being processed, so reverese formatting is not applicable on filter data
% For the new sets of 3D filters, same process as xzy is repeated. (for the last set, may be a little different if not integer multiple due to 
% the fact that less #of filter channels may need to be discarded, but ignoring that, wont make much of a difference I guess.)
Access_DRAM_filter_xzyF = Access_DRAM_filter_xzy * (Nos_of_filter/unitvol_nosof_3D_filter);
% The same access pattern as xzy repeats for all sets of 3D filters
Access_SRAM_filter_xzyF = Access_SRAM_filter_xzy * (Nos_of_filter/unitvol_nosof_3D_filter);

%%%% Access cost for ifmap data
if (SRAM_ifmap >= (ifmap_width_effective * ifmap_height_effective * Nos_of_channel * bw_ifmap))
    disp("full 3D volume of ifmap fits in ifmap-SRAM")
    gamma_DRAM_ifmap_xzyF = 1;  %ifmap-SRAM has enough memory to fit the full ifmap volume, hence each ifmap data loaded from DRAM once
    Access_DRAM_ifmap_xzyF = ifmap_width_effective * ifmap_height_effective * Nos_of_channel;
else
    % Implementing reverse formatting for the ifmap data while processing towards F-direction  
    % the access pattern for the first set of 3D filters is repeated for all the next set of 3D filters   
    if(ceil(Nos_of_filter/unitvol_nosof_3D_filter) == 1)  
        %if all filters are covered in the first set of 3D filter pass then there is no previously loaded ifmap xR templates to subtract
        Access_DRAM_ifmap_xzyF = Access_DRAM_ifmap_xzy;
    else
        % Some 2D ifmap xR templates will already be in ifmap-SRAM after processing the first set of 3D filter
        ifmap_xR_template_to_fit_SRAM = floor(SRAM_ifmap / (ifmap_width_effective * filter_height * bw_ifmap)); %#of ifmap template fits in ifmap SRAM 
        % In the following lines: ifmap_already_in_SRAM is the amount of ifmap data which is already in the ifmap-SRAM during the processing of the previous set of 3D-filter
        % Therefore, no need to load it since performing reverese formatting       
        first_3Dset_access = Access_DRAM_ifmap_xzy;         % access cost for the first set of 3D filters
        %The following line "ifmap_already_in_SRAM" does not have the ifmap stride overlap issue because the divisor in the calculation
        %of ifmap_xR_template_to_fit_SRAM is used as multiplicator here
        ifmap_already_in_SRAM = (ifmap_width_effective * filter_height * ifmap_xR_template_to_fit_SRAM);
        next_3Dset_access = (Access_DRAM_ifmap_xzy - ifmap_already_in_SRAM);    % access cost for each next set of 3D filters
        % ceil comes here because even if there is less #of 3D filter in the last pass, ifmap data need to be loaded
        Access_DRAM_ifmap_xzyF = first_3Dset_access + next_3Dset_access * ceil(Nos_of_filter/unitvol_nosof_3D_filter - 1);
    end  
end
    
% The same SRAM access pattern for xzy is repeated. The reason ceil comes here is that, in the last pass if there is less filter than "unitvol_nosof_3D_filter", 
% some PE column will reamin unused, but ifmap channels will need to be loaded from SRAM to RF and broadcasted to the PE columns which are being used.
Access_SRAM_ifmap_xzyF = Access_SRAM_ifmap_xzy * ceil(Nos_of_filter/unitvol_nosof_3D_filter);

%%%% Access cost for psum data
% from the xzy direction unitvol_nosof_3D_filter ofmap planes form, same process repeats for the remaining ofmap planes
Access_SRAM_psum_xzyF = Access_SRAM_psum_xzy * (Nos_of_filter/unitvol_nosof_3D_filter);
% ceil is not required here because in the last pass if there is less filter than "unitvol_nosof_3D_filter", there will be no psum out from the unused PE columns


%This is the access energy up to this depth of the decision
Weighted_Access_xzyF = (Access_DRAM_filter_xzyF * bw_filter + Access_DRAM_ifmap_xzyF * bw_ifmap) * E_DRAM_to_SRAM + ...
    (Access_SRAM_filter_xzyF * bw_filter + Access_SRAM_ifmap_xzyF * bw_ifmap + Access_SRAM_psum_xzyF * bw_psum) * E_SRAM_to_RF;  % Access energy (in joule)
SRAM_Access_xzyF = (Access_SRAM_filter_xzyF * bw_filter + Access_SRAM_ifmap_xzyF * bw_ifmap + Access_SRAM_psum_xzyF * bw_psum) * E_SRAM_to_RF;
DRAM_Access_xzyF = (Access_DRAM_filter_xzyF * bw_filter + Access_DRAM_ifmap_xzyF * bw_ifmap) * E_DRAM_to_SRAM; 

%%%%%%%%%%% Calculating DRAM bandwidth induced stall count
if (slack_cycle_xz <= 0)  
    %"CASE-1: 1 1"
    % CASE-1: There is no slack cycle from xz
    % in this branch, mainly the source of stall is the filter data, ifmap data has much less impact
    if (SRAM_filter <= filter_height * filter_width * Nos_of_channel * unitvol_nosof_3D_filter * bw_filter) || (slack_cycle_xzy <= 0)
            % There is either no slack cycle or filter-SRAM cannot fit more than RSCK volume of filter
            % In this case stall from filter data load at the first xz-y1 block cannot be overlapped
            % the same stall from first xzy repeats for all the set of 3D filters
            DRAM_stall_xzyF = ceil(DRAM_stall_xzy * (Nos_of_filter/unitvol_nosof_3D_filter));
    else
            % only considering to overlap filter data at the first xz-y1 block since the source of major stall is the filter data. 
            % if filter data load is overlapped, stall from ifmap will be small or no stall at all. 
            filter_data_xz = filter_height * filter_width * Nos_of_channel * unitvol_nosof_3D_filter * bw_filter; % in bit, filter data for the next set of 3D filters
            Extra_filter_data_to_fit_SRAM = SRAM_filter - filter_data_xz; %in bit, filter data exceeding RSCK volume which can fit in filter-SRAM
            
            % data loaded is the minimum between the required filter data for the next set of 3D filters and data which fits in filter-SRAM
            loaded_filter_data = min(filter_data_xz, Extra_filter_data_to_fit_SRAM); 
            
            cycle_to_load_extra_filter_data = ceil((ceil(loaded_filter_data/DRAM_block_size) * DRAM_block_size)/BW_DRAM);
            
            % Number of cycles which is overlapped, hence no need to incur stall for these many cycles. This is the minimum between whatever slack cycles
            % we have from the previous xzy block and the nos of cycle required to load the data
            cycle_saved_from_stall = min(slack_cycle_xzy, cycle_to_load_extra_filter_data);
            
            % The amount of stall from the next xzy-block is reduced by the number of cycles saved. Getting a negative value from
            % "DRAM_stall_xzy - cycle_saved_from_stall" means no stall from ifmap data as well at the first xz-y1 of next xzy-block.
            % Hence the max() will take 0 stall.
            DRAM_stall_xzy_next = max(0, DRAM_stall_xzy - cycle_saved_from_stall);
            
            DRAM_stall_xzyF = DRAM_stall_xzy + ceil(DRAM_stall_xzy_next * (Nos_of_filter/unitvol_nosof_3D_filter - 1));
            % first term: for the first set of 3D filter set
            % second term: for all the next set of 3D filter set
    end
    
elseif (DRAM_stall_xz <= 0)
    % "CASE-2: 0 0")
    % This means there is no stall at x, z and each xz block can compute itself without incorporating any stall, so no stall to complete the xzy direction
    %  So the same zero stall repeats for all 3D filter sets
    DRAM_stall_xzyF = 0;
end

cycle_count_xzyF_with_stall = cycle_count_xzyF + DRAM_stall_xzyF;   

%% Incorporating the calculation when some directions are processed partially
% In this brach only one effective can come from the psum-SRAM storage limitation: X-effective
% Z_effective never comes from psum-SRAM since Z-direction kills psum
% Y and F is done after Z-direction is done, hence ofmaps are formed and psum-SRAM can be made empty by moving ofmaps to DRAM.
% Therefore Y and F effective does not occur in this branch

%DRAM stall to write the ofmap value back to DRAM
ofmap_data_bits = ofmap_height_effective * ofmap_width_effective * Nos_of_filter * bw_ofmap; %in bits
ofmap_data_to_write = ceil(ofmap_data_bits/DRAM_block_size) * DRAM_block_size;
cycle_to_write_ofmap = ceil(ofmap_data_to_write/BW_DRAM); % nos of stall to write the ofmap back to DRAM [WILL ADD THE OPTION FOR POOLING HERE]
Access_DRAM_ofmap_xzyF = ofmap_height_effective * ofmap_width_effective * Nos_of_filter; % in nos of element, not bit


X_effective = [ofmap_width_effective ifmap_width_effective];
cycle_count = [cycle_count_ideal_xzyF DRAM_stall_xzyF cycle_to_write_ofmap];
SRAM_Access = [Access_SRAM_filter_xzyF Access_SRAM_ifmap_xzyF Access_SRAM_psum_xzyF];  % in Nos of element, not bit
DRAM_Access = [Access_DRAM_filter_xzyF Access_DRAM_ifmap_xzyF Access_DRAM_ofmap_xzyF]; % in Nos of element, not bit


end
