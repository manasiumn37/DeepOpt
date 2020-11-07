function [XY_effective, cycle_count, SRAM_Access, DRAM_Access] = XtoYtoZtoF(Layer_param, Hardware_param, Tech_param)

%X: Direction along the width of ofmap
%Y: Direction along the height of ofmap
%Z: Direction along the channels of ifmap
%F: Direction along the channels of ofmap = Direction along the 3D filters

% This function performs the computation of the input layer using the XYZF branch

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

used_MAC_per_column = unitvol_nosof_channel * filter_width;
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
    XY_effective = [ofmap_width ofmap_height ifmap_width ifmap_height];  % this will prevent calling the asymmetric subvolumes
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
    XY_effective = [ofmap_width ofmap_height ifmap_width ifmap_height];  % this will prevent calling the asymmetric subvolumes
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

%% Tree process along x->y direction of ifmap
% No need to check filter memory requirement since alreday included filter height in the previous check

% SRAM-psum memory requirment check
% How many rows I can cover in the y-direction will depend on psum-SRAM storage since psum does not go to DRAM
rows_to_fit_psum_SRAM = floor(SRAM_psum / (((ifmap_width_effective - filter_width)/stride + 1) * unitvol_nosof_3D_filter * bw_psum));
if (rows_to_fit_psum_SRAM>=ofmap_height)
    disp("Full Y-direction possible by psum-SRAM")
    ofmap_height_effective = ofmap_height;
    ifmap_height_effective = (ofmap_height_effective - 1) * stride + filter_height;
    flag_Yeff = 0;     % Y_effective flag to indicate that psum-SRAM is not full at this level
else
    disp("Partial Y-direction is being processed due to psum-SRAM size limitation")
    ofmap_height_effective = rows_to_fit_psum_SRAM;
    ifmap_height_effective = (ofmap_height_effective - 1) * stride + filter_height;
    flag_Yeff = 1;     % Y_effective flag to indicate that psum-SRAM is full at this level
end


% For each ofmap row, #of cycle required is cycle_count_ideal_x
cycle_count_xy = cycle_count_ideal_x * ofmap_height_effective + initial_offset_cycle;
cycle_count_ideal_xy = cycle_count_ideal_x * ofmap_height_effective;    %without the initial offset cycles


%%%% Access for filter data
%filter data loaded from DRAM in the previous X-direction will be reused
gamma_DRAM_filter_xy = 0;
Access_DRAM_filter_xy = Access_DRAM_filter_x;  % copying the previous access result since will use it directly in the calculation of Weighted_Access
%For ecah ofmap row, same amount of SRAM filter access occur & that amount is Access_SRAM_filter_x (from the x-direction search)
Access_SRAM_filter_xy = Access_SRAM_filter_x * ofmap_height_effective;

%%%%%% Access for ifmap data
%loading the new ifmap rows except the previously loaded ones and including the previous result for direct calculation
gamma_DRAM_ifmap_xy = 1;  % This indicates no ifmap loaded from DRAM twice yet
Access_DRAM_ifmap_xy = ifmap_height_effective * ifmap_width_effective * unitvol_nosof_channel; %in element

%For ecah ofmap row, same amount of SRAM ifmap access occur & that amount is Access_SRAM_ifmap_x (from the x-direction search)
Access_SRAM_ifmap_xy = Access_SRAM_ifmap_x * ofmap_height_effective;

%%%%% Access for psum data, psum not going to DRAM
%same pattern of psum access happen for each ofmap row
Access_SRAM_psum_xy = Access_SRAM_psum_x * ofmap_height_effective;

%This is the access energy up to this depth of the decision
Weighted_Access_xy = (Access_DRAM_filter_xy * bw_filter + Access_DRAM_ifmap_xy * bw_ifmap) * E_DRAM_to_SRAM + ...
    (Access_SRAM_filter_xy * bw_filter + Access_SRAM_ifmap_xy * bw_ifmap + Access_SRAM_psum_xy * bw_psum) * E_SRAM_to_RF;
SRAM_Access_xy = (Access_SRAM_filter_xy * bw_filter + Access_SRAM_ifmap_xy * bw_ifmap + Access_SRAM_psum_xy * bw_psum) * E_SRAM_to_RF;


%%%%%%%%%%% Calculating DRAM bandwidth induced stall count
New_loaded_filter_DRAM_xy = 0;   % when moving x->y direction, the previosuly loaded filter data is reused
New_loaded_ifmap_POR = stride * ifmap_width_effective * unitvol_nosof_channel * bw_ifmap; % in bit; this gives the additional ifmap data needed to load
                                                                                          % to process each next ofmap row
DRAM_data_POR = New_loaded_filter_DRAM_xy + New_loaded_ifmap_POR; %in bit; total new data need to be loaded from DRAM per ofmap row (POR)
%cycle_to_load_DRAM_data_POR = ceil(DRAM_data_POR/BW_DRAM) % cycle count to load data per ofmap row

DRAM_data_loaded_POR = ceil(DRAM_data_POR/DRAM_block_size) * DRAM_block_size;
cycle_to_load_DRAM_data_POR = ceil(DRAM_data_loaded_POR/BW_DRAM); % cycle count to load data per ofmap row


cycle_to_compute_POR = cycle_count_ideal_x;                % cycle count to process the computation for each ofmap row
if (cycle_to_load_DRAM_data_POR > cycle_to_compute_POR)
    DRAM_stall_POR = cycle_to_load_DRAM_data_POR - cycle_to_compute_POR;  %DRAM stall per ofmap row
else
    DRAM_stall_POR = 0;
end

DRAM_stall_y = DRAM_stall_POR * (ofmap_height_effective - 1); %excluding the first row since it is accounted during x-direction
% omitting the slack_cycle variable here since for the first ofmap-row we need to load more data (both filter & ifmap) than the next ofmap rows (only ifmap).
% So, if there is no stall at x, there will be no stall at y.

% Calculating #of slack cycles from this stage
if (DRAM_stall_y == 0)
    slack_cycle_y = (cycle_to_compute_POR - cycle_to_load_DRAM_data_POR) * (ofmap_height_effective - 1); %#of extra cycles during computation when data can be loaded from DRAM
else
    slack_cycle_y = 0;    % No extra cycles for data load since data load time is higher than compute time
end

DRAM_stall_xy = DRAM_stall_x + DRAM_stall_y;    % DRAM stall upto xy
slack_cycle_xy = slack_cycle_x + slack_cycle_y; % slack cycle upto xy


%% Tree process along x->y->z direction of ifmap

%%%% Filter memory requirement check is not necessary since new filter channels can be loaded from DRAM
% Just seeing how many filter templates of "unitvol_nosof_3D_filter" fits in filter-SRAM; 
templates_tofit_filter_SRAM = floor(SRAM_filter/(unit_vol_filter * filter_height * bw_filter));

%%%% ifmap memory check is not necessary since new ifmap channels can be loaded from DRAM
% Just to see
templates_tofit_ifmap_SRAM = floor(SRAM_ifmap/(ifmap_width_effective * ifmap_height_effective * unitvol_nosof_channel * bw_ifmap));

%%%% psum memory check is not necessary since with z-direction no new nonreducible psum is generated than the xy direction

%%%% calculation of cycle count
% for each ifmap plane with "unitvol_nosof_channel", same amount of cycle count is required
cycle_count_ideal_xyz = cycle_count_ideal_xy * ceil(Nos_of_channel/unitvol_nosof_channel); %without the initial offset cycles
cycle_count_xyz = cycle_count_ideal_xyz + initial_offset_cycle;


%%%%% Access cost for filter data
gamma_DRAM_filter_xyz = 1; %each filter data is loaded once from DRAM
Access_DRAM_filter_xyz = filter_width * filter_height * Nos_of_channel * unitvol_nosof_3D_filter; % in element
%The amount of access in xy direction is repeated (Nos_of_channel/unitvol_nosof_channel) times to cover the full z-direction
Access_SRAM_filter_xyz = Access_SRAM_filter_xy * (Nos_of_channel/unitvol_nosof_channel);

%%%%% Access cost for ifmap data
gamma_DRAM_ifmap_xyz = 1; %ecah ifmap data is loaded once from DRAM
Access_DRAM_ifmap_xyz = ifmap_width_effective * ifmap_height_effective * Nos_of_channel;
%The amount of access in xy direction is repeated (Nos_of_channel/unitvol_nosof_channel) times to cover the full z-direction
Access_SRAM_ifmap_xyz = Access_SRAM_ifmap_xy * (Nos_of_channel/unitvol_nosof_channel);

%%%%% Access cost for psum data
% Access_SRAM_psum_xy = The access pattern for the first xy plane of ofmap is same as Access_SRAM_psum_xy
% psum_xyz_x, psum_xyz_xy = For the next channels, there will be a previous psum to read generated from the previous channel. Hence all row of the filter tenplate
% is both read-write when processing along the x-direction. This access pattern is repeated ((Nos_of_channel/unitvol_nosof_channel) - 1) times
psum_xyz_x = 2 * Generated_psum_effective_x * filter_height;    % psum SRAM access cost while processing the x-direction for next ifmap channels
psum_xyz_xy = psum_xyz_x * ofmap_height_effective;              % psum SRAM access cost while processing the xy-direction for next ifmap channels
Access_SRAM_psum_xyz = Access_SRAM_psum_xy + psum_xyz_xy * ceil((Nos_of_channel/unitvol_nosof_channel) - 1);
% ceil comes here becasue in the last pass even if there is less channel than unitvol_nosof_channel, same nos of psum read/write is required


%This is the access energy up to this depth of the decision
Weighted_Access_xyz = (Access_DRAM_filter_xyz * bw_filter + Access_DRAM_ifmap_xyz * bw_ifmap) * E_DRAM_to_SRAM + ...
    (Access_SRAM_filter_xyz * bw_filter + Access_SRAM_ifmap_xyz * bw_ifmap + Access_SRAM_psum_xyz * bw_psum) * E_SRAM_to_RF;
SRAM_Access_xyz = (Access_SRAM_filter_xyz * bw_filter + Access_SRAM_ifmap_xyz * bw_ifmap + Access_SRAM_psum_xyz * bw_psum) * E_SRAM_to_RF;


%%%%%%%%%%% Calculating DRAM bandwidth induced stall count
if (slack_cycle_xy <= 0)  
    disp("CASE-1: 1 1")
    % CASE-1: There is no slack cycle from xy
    % Hence, each xy block on its own, so no stall overlapped and DRAM_stall_xy is repeated for all unitvol_nosof_channel
    DRAM_stall_xyz = ceil(DRAM_stall_xy * (Nos_of_channel/unitvol_nosof_channel));
    % ceil is not required for (Nos_of_channel/unitvol_nosof_channel) since in the last pass if there is less #of channels than unitvol_nosof_channel,
    % it will take less time to load that data. ceil is used for xyz block processing
elseif (DRAM_stall_xy <= 0)
    disp("CASE-2: 0 0")
    %CASE-2: Each xy block can compute itself without incorporating any stall, so no stall to complete the xyz direction
    DRAM_stall_xyz = 0;
else 
    disp("CASE-3: 1 0")
    % CASE-3: There are some slack cycles from xy, as well as some stall from xy. 
    % This means the source of stall is the (filter+ifmap) data at the first slot of xy-z block.
    % The load time of ifmap is not creating the stall, it can be overlapped with computation,
    % So, no need to worry about ifmap-SRAM constraint   
    % Only need to check if filter-SRAM memory permits to load filter data from DRAM during those slack cycles to overlap the stall
    % Assumption: if the filter dada needed for the next xy block processing can fully fit in filter-SRAM only then we will load it, 
    % otherwise wont load any filter data for the next xy block
    filter_data_Pxy = filter_height * filter_width * unitvol_nosof_channel * unitvol_nosof_3D_filter * bw_filter; %filter data needed per xy block
    if (SRAM_filter >= (2 * filter_data_Pxy)) %There is room in SRAM_filter to fit the full filter data needed for the processing of the next xy block
        if (DRAM_stall_xy <= slack_cycle_xy)   
            disp("CASE 3.1 Full xy stall can be overlapped")
            DRAM_stall_z = 0;
        else
            disp("CASE 3.2 Stall can be overlapped partially")
            DRAM_stall_Pxy = DRAM_stall_xy - slack_cycle_xy; % stall per xy-block
            DRAM_stall_z = ceil(DRAM_stall_Pxy * (Nos_of_channel/unitvol_nosof_channel - 1)); 
            % DRAM_stall_Pxy is repeated for all the next set of channels,excluding the first set of channels, it will be incorporated later by adding DRAM_stall_xy
            % ceil is not required for Nos_of_channel/unitvol_nosof_channel since in the last pass if there is less #of channels than unitvol_nosof_channel,
            % it will take less time to load that data. ceil is used for the block processing
        end
        DRAM_stall_xyz = DRAM_stall_xy + DRAM_stall_z;
    else
        % Each xy block on its own, so no stall overlapped and DRAM_stall_xy is repeated for all unitvol_nosof_channel
        DRAM_stall_xyz = ceil(DRAM_stall_xy * (Nos_of_channel/unitvol_nosof_channel));
    end
end

% Assumtion: Not doing any computation of the slack cycle after xyz, For simplicity, after xyz we assume the slack variable is zero
% keeping track of slack variable becomes complicated as we move higher depth of the tree
% This assumption is based on the fact that the ratio of stall_xy/compute_cycle_xyz is very small, 
% hence it is not worth to make the implementation very complicated for this
slack_decision_factor = DRAM_stall_xy/cycle_count_ideal_xyz;

%% Tree process along x->y->z->F direction of ifmap

% Filter SRAM memeory requirement check is not necessary since new 3D filter can be loaded from DRAM.

% ifmap SRAM memory requirement check is not necessary since the same process of xyz will be repeated for the new set of 3D filters.

% psum SRAM memory requirement check is not necessary since after processing the full Z-direction, the ofmaps are formed which can be 
% written back to DRAM and the same process as xyz will be repeated.

%%%%% Calculation of cycle count
% The same process of xyz will be repeated "(Nos_of_filter/unitvol_nosof_3D_filter)" times to cover the full F-direction
% The ceil works becasue, when not integer multiple, in the last set of 3D filters some PE columns will remain unused
cycle_count_ideal_xyzF = cycle_count_ideal_xyz * ceil(Nos_of_filter/unitvol_nosof_3D_filter); %without the initial offset cycles
cycle_count_xyzF = cycle_count_ideal_xyzF + initial_offset_cycle;


%%%%% Access cost for filter data
gamma_DRAM_filter_xyzF = 1; % so maximal utilization of filter data from DRAM. each filter data is loaded from DRAM only once
Access_DRAM_filter_xyzF = filter_width * filter_height * Nos_of_channel * unitvol_nosof_3D_filter * (Nos_of_filter/unitvol_nosof_3D_filter); % in element

% The same access pattern as xyz repeats for all sets of 3D filters
Access_SRAM_filter_xyzF = Access_SRAM_filter_xyz * (Nos_of_filter/unitvol_nosof_3D_filter);


%%%%% Access cost for ifmap data
% Depending on the ifmap-SRAM memory, some/all ifmap channels/rows will need to be 
% discarded from SRAM and then loaded from DRAM again when moving along F-direction (i.e., new sets of 3D filters)
% ifmap xR template is the unit-amount to discard from SRAM,
% To be consistent with other branches, and for simplicity, using ifmap row as the unit to keep or discard from SRAM here
ifmap_xy_tofit_SRAM = floor(SRAM_ifmap/(ifmap_width_effective * ifmap_height_effective * bw_ifmap)); % #of ifmap channel which fits in ifmap-SRAM
ifmap_xy_reminder = rem(SRAM_ifmap,(ifmap_width_effective * ifmap_height_effective * bw_ifmap));
ifmap_xR_template_to_fit_SRAM = floor(SRAM_ifmap / (ifmap_width_effective * filter_height * bw_ifmap)); %#of ifmap xR template fits in ifmap SRAM 

if (SRAM_ifmap >= (ifmap_width_effective * ifmap_height_effective * Nos_of_channel * bw_ifmap))
    disp ("ifmap-SRAM is large enough to fit one full ifmap xyz volume")
    gamma_DRAM_ifmap_xyzF = 1; %ecah ifmap data can be loaded only once from DRAM
    Access_DRAM_ifmap_xyzF = ifmap_width_effective * ifmap_height_effective * Nos_of_channel;
else
    % During the processing of first set of 3D filters, all ifmap channels are loaded once from DRAM
    First_set_Access = ifmap_width_effective * ifmap_height_effective * Nos_of_channel;
    % For all the next set of 3D filters, the amount of discarded templates need to be loaded.The processing format for ifmap channel is 
    % first to end--> end to first--> first to end and so on. 
    if (ifmap_xy_reminder == 0) % this will trigger if ifmap-SRAM is a integer miltiple of ifmap_xy_tofit_SRAM
        %For each next set of 3D filters, "(Nos_of_channel - ifmap_xy_tofit_SRAM)" channels need to be loaded
        Next_set_Access = ifmap_width_effective * ifmap_height_effective * (Nos_of_channel - ifmap_xy_tofit_SRAM); 
    else
        nos_of_discarded_full_channel = Nos_of_channel - ifmap_xy_tofit_SRAM - 1; % nos of full xy channel which needs to be discarded
        ifmap_row_tofit_SRAM = floor(SRAM_ifmap/(ifmap_width_effective * bw_ifmap)); % #of ifmap 1D row which fits in ifmap-SRAM
        rows_in_partial_channel = ifmap_row_tofit_SRAM - (ifmap_xy_tofit_SRAM * ifmap_height_effective); % #of rows of the partial channel which fits in SRAM
        discarded_rows_in_partial_channel = ifmap_height_effective - rows_in_partial_channel;      % #of discarded rows from the partial channel
        % the following line gives the #of ifmap data which is discarded and hence loaded multiple times
        Next_set_Access = (ifmap_width_effective * ifmap_height_effective * nos_of_discarded_full_channel) + (ifmap_width_effective * discarded_rows_in_partial_channel);
    end
    % The ceil comes from the fact that, in the last pass even thought there may be less filter than "unitvol_nosof_3D_filter", the ifmap channels need to be loaded
    Access_DRAM_ifmap_xyzF = First_set_Access + Next_set_Access * ceil(Nos_of_filter/unitvol_nosof_3D_filter - 1);
    % Around 50% savings of DRAM access here by this reverse formating scheme   
end
% The same SRAM access pattern for xyz is repeated. The reason ceil comes here is that, in the last pass if there is less filter than "unitvol_nosof_3D_filter", 
% some PE column will reamin unused, but ifmap channels will need to be loaded from SRAM to RF and broadcasted to the PE columns which are being used.
Access_SRAM_ifmap_xyzF = Access_SRAM_ifmap_xyz * ceil(Nos_of_filter/unitvol_nosof_3D_filter);

%%%%% Access cost for psum data
Access_SRAM_psum_xyzF = Access_SRAM_psum_xyz * (Nos_of_filter/unitvol_nosof_3D_filter);
% ceil is not required here because in the last pass if there is less filter than "unitvol_nosof_3D_filter", there will be no psum out from the unused PE columns

%This is the access energy up to this depth of the decision
Weighted_Access_xyzF = (Access_DRAM_filter_xyzF * bw_filter + Access_DRAM_ifmap_xyzF * bw_ifmap) * E_DRAM_to_SRAM + ...
    (Access_SRAM_filter_xyzF * bw_filter + Access_SRAM_ifmap_xyzF * bw_ifmap + Access_SRAM_psum_xyzF * bw_psum) * E_SRAM_to_RF;  % Access energy (in joule)
SRAM_Access_xyzF = (Access_SRAM_filter_xyzF * bw_filter + Access_SRAM_ifmap_xyzF * bw_ifmap + Access_SRAM_psum_xyzF * bw_psum) * E_SRAM_to_RF;
DRAM_Access_xyzF = (Access_DRAM_filter_xyzF * bw_filter + Access_DRAM_ifmap_xyzF * bw_ifmap) * E_DRAM_to_SRAM; 


%%%%%%%%%%% Calculating DRAM bandwidth induced stall count
if (slack_cycle_xy <= 0)
    % CASE-1: There is no slack cycle from xy
    % In this case the ifmap data already in SRAM dont need to be loaded in the reverse order formatting,
    % Therefore, the stall from ifmap wont happen from those ifmap data
    % Assumptions: For simplicity using xy-block for ifmap row. For these #of xy block we do not need to load the ifmap data
    % The implication of this is that, even though due to reverse order formatting some partial ifmap data from a xy-block is not being loaded from DRAM,
    % the systolic array is just waiting. The is the cost in performance we pay to have simplified control.
    ifmap_xy_block_tofit_SRAM = floor(SRAM_ifmap/(ifmap_width_effective * ifmap_height_effective * unitvol_nosof_channel * bw_ifmap)); %#of ifmap xy block which fits in SRAM, 
    
    DRAM_filter_data_loaded_x = ceil(DRAM_filter_data_x/DRAM_block_size) * DRAM_block_size;  
    cycle_to_load_DRAM_filter_data_x = ceil(DRAM_filter_data_loaded_x/BW_DRAM);               % #of cycle to load the filter data from DRAM at x-direction
    DRAM_stall_filter_x = max(0, cycle_to_load_DRAM_filter_data_x - cycle_count_ideal_x);              % stall from filter data only at x-direction
    DRAM_stall_xy_reverse = DRAM_stall_filter_x + 0;  % due to reverse formatting ifmap data is not loaded, so no stall from y direction
    
    if (SRAM_ifmap >= (ifmap_width_effective * ifmap_height_effective * Nos_of_channel * bw_ifmap))
        %ifmap-SRAM is large enough to fit one full ifmap xyz volume, so no ifmap need to be loaded again from DRAM while moving towards F-direction
        DRAM_stall_xyz_reverse = ceil(DRAM_stall_xy_reverse * (Nos_of_channel/unitvol_nosof_channel));
    else
        DRAM_stall_xyz_reverse = (DRAM_stall_xy_reverse * ifmap_xy_block_tofit_SRAM) + ceil(DRAM_stall_xy * (Nos_of_channel/unitvol_nosof_channel - ifmap_xy_block_tofit_SRAM));
        % Here the first term gives lower stall which is applicable to the ifmap xy block which is already in SRAM during reverse order formatting
        % The second term gives the regular stall calculated in xy depth for the ifmap channels which are not in SRAM and has to be loaded from DRAM
    end
    
    DRAM_stall_xyzF = DRAM_stall_xyz + ceil(DRAM_stall_xyz_reverse * (Nos_of_filter/unitvol_nosof_3D_filter - 1));
    % First term: for the first set of 3D filter, no reverse order formatting applicable, so regular DRAM_stall_xyz term used
    % Second term: for the next set of 3D filters, reverse order formatting applicable, so DRAM_stall_xyz_reverse used   
else
    % CASE-2 & CASE-3
    %The same process of xyz repeats for all the next set of filters.
    % This formula works for case-2 & 3 only. Since here no stall is happening due to ifmap data,
    % even if we are doing reverse order formatting for ifmap, that do not affect the stall count
    DRAM_stall_xyzF = ceil(DRAM_stall_xyz * (Nos_of_filter/unitvol_nosof_3D_filter));
    % ceil not used for Nos_of_filter/unitvol_nosof_3D_filter, the argument is that, in the last pass, if there is less #of unitvol_nosof_3D_filter, 
    % it will take less time to load that data. ceil is used for the full xyzF block
end

DRAM_stall_xyzF;
cycle_count_ideal_xyzF;
cycle_count_xyzF_with_stall = cycle_count_xyzF + DRAM_stall_xyzF;

%% Parameters to return
% In this brach maximum two effective can come from the psum-SRAM storage limitation: X-effective and Y-effective
% Z_effective never comes from psum-SRAM since Z-direction kills psum
% F is the last level in this branch before which Z-direction is done, hence ofmaps are formed and psum-SRAM can be made empty by moving ofmaps to DRAM.
% Therefore F-effective does not occur in this branch

%DRAM stall to write the ofmap value back to DRAM
ofmap_data_bits = ofmap_height_effective * ofmap_width_effective * Nos_of_filter * bw_ofmap; %in bits
ofmap_data_to_write = ceil(ofmap_data_bits/DRAM_block_size) * DRAM_block_size;
cycle_to_write_ofmap = ceil(ofmap_data_to_write/BW_DRAM); % nos of stall to write the ofmap back to DRAM [WILL ADD THE OPTION FOR POOLING HERE]
Access_DRAM_ofmap_xyzF = ofmap_height_effective * ofmap_width_effective * Nos_of_filter; % in nos of element, not bit


XY_effective = [ofmap_width_effective ofmap_height_effective ifmap_width_effective ifmap_height_effective];
cycle_count = [cycle_count_ideal_xyzF DRAM_stall_xyzF cycle_to_write_ofmap];
SRAM_Access = [Access_SRAM_filter_xyzF Access_SRAM_ifmap_xyzF Access_SRAM_psum_xyzF];  % in Nos of element, not bit
DRAM_Access = [Access_DRAM_filter_xyzF Access_DRAM_ifmap_xyzF Access_DRAM_ofmap_xyzF]; % in Nos of element, not bit


end
