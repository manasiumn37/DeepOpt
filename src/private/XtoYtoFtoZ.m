function [XYF_effective, cycle_count, SRAM_Access, DRAM_Access] = XtoYtoFtoZ(Layer_param, Hardware_param, Tech_param)

%X: Direction along the width of ofmap
%Y: Direction along the height of ofmap
%Z: Direction along the channels of ifmap
%F: Direction along the channels of ofmap = Direction along the 3D filters

% This function performs the computation of the input layer using the XYFZ branch

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
    XYF_effective = [ofmap_width ofmap_height ifmap_width ifmap_height Nos_of_filter];  % this will prevent calling the asymmetric subvolumes
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
    XYF_effective = [ofmap_width ofmap_height ifmap_width ifmap_height Nos_of_filter];  % this will prevent calling the asymmetric subvolumes
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


%% Tree process along x->y->F direction of ifmap

% No need to check filter-SRAM memory requirement. The check in x-direction step is enough and when new set of 3D filter is being processed we can load it from DRAM

% ifmap-SRAM memory requirement is not required since the same argument from the x->y step is true

% psum-SRAM memory requirement check
% since psum does not go to DRAM, how many 3D filter I can process will depend on the psum-SRAM storage
ofmap_planes_tofit_psum_SRAM = floor(SRAM_psum/(ofmap_width_effective * ofmap_height_effective * bw_psum));
if (ofmap_planes_tofit_psum_SRAM >= Nos_of_filter)
    F_effective = Nos_of_filter; % Full ofmap volume fits into psum-SRAM
    flag_Feff = 0;   % F_effective flag to indicate that psum-SRAM is not full at this level
else
    % for the sake of regular processing pattern, F_effective is chosen to be an integer multiple of unitvol_nosof_3D_filter.
    F_effective = floor(ofmap_planes_tofit_psum_SRAM/unitvol_nosof_3D_filter) * unitvol_nosof_3D_filter;  %extent of F-direction possible to process given the psum-SRAM memory
    flag_Feff = 1;   % F_effective flag to indicate that psum-SRAM is full at this level
end


%%%%%%%Calculation of cycle count  
% for each "unitvol_nosof_3D_filter" nos of ofmap plane same amount of cycle count is required 
cycle_count_ideal_xyF = cycle_count_ideal_xy * ceil(F_effective/unitvol_nosof_3D_filter);
cycle_count_xyF = cycle_count_ideal_xyF + initial_offset_cycle;   
% The ceil works becasue, when not integer multiple, in the last set of 3D filters some PE columns will remain unused

%%%%%% Access for filter data
gamma_DRAM_filter_xyF = 1; %ecah filter data is loaded from DRAM once 
Access_DRAM_filter_xyF = filter_height * filter_width * unitvol_nosof_channel * F_effective;
%The amount of access in xy direction is repeated (F_effective/unitvol_nosof_3D_filter) times to cover the F-effective nos of 3D filter
Access_SRAM_filter_xyF = Access_SRAM_filter_xy * (F_effective/unitvol_nosof_3D_filter);


%%%%%% Access for ifmap data
if (SRAM_ifmap >= (ifmap_width_effective * ifmap_height_effective * unitvol_nosof_channel * bw_ifmap))
    disp("ifmap-SRAM has enough storage to process xyF direction by loading ifmap data once from DRAM")
    gamma_DRAM_ifmap_xyF = 1; % each ifmap data is loaded from DRAM once
    Access_DRAM_ifmap_xyF = ifmap_height_effective * ifmap_width_effective * unitvol_nosof_channel; %in element
else
    % During the processing of first set of 3D filters, all ifmap rows are loaded once from DRAM
    First_set_Access = ifmap_width_effective * ifmap_height_effective * unitvol_nosof_channel;
    % For all the next set of 3D filters, the amount of discarded rows need to be loaded. (reverse order formatting)
    % following line gives the #of ifmap 1D row from each unitvol_nosof_channel which fits in ifmap-SRAM
    ifmap_rows_to_fit_SRAM = floor(SRAM_ifmap/(ifmap_width_effective * unitvol_nosof_channel * bw_ifmap));
    discarded_ifmap_rows = ifmap_height_effective - ifmap_rows_to_fit_SRAM;  % #of discarded rows from each unitvol_nosof_channel
    Next_set_Access = ifmap_width_effective * discarded_ifmap_rows * unitvol_nosof_channel;
    % The ceil comes from the fact that, in the last pass even though there may be less filter than "unitvol_nosof_3D_filter", the ifmap rows need to be loaded
    Access_DRAM_ifmap_xyF = First_set_Access + Next_set_Access * ceil(F_effective/unitvol_nosof_3D_filter - 1);
end
% The same SRAM access pattern for xy is repeated. The reason ceil comes here is that, in the last pass if there is less filter than "unitvol_nosof_3D_filter", 
% some PE column will reamin unused, but ifmap will need to be loaded from SRAM to RF and broadcasted to the PE columns which are being used.
Access_SRAM_ifmap_xyF = Access_SRAM_ifmap_xy * ceil(F_effective/unitvol_nosof_3D_filter);


%%%%% Access cost for psum data
%the same psum access pattern for "unitvol_nosof_3D_filter" nos of xy plane of ofmap is repeated for all other pack of ofmap planes
Access_SRAM_psum_xyF = Access_SRAM_psum_xy * (F_effective/unitvol_nosof_3D_filter);
% ceil is not required here because in the last pass if there is less filter than "unitvol_nosof_3D_filter", there will be no psum out from the unused PE columns

    
%This is the access energy up to this depth of the decision
Weighted_Access_xyF = (Access_DRAM_filter_xyF * bw_filter + Access_DRAM_ifmap_xyF * bw_ifmap) * E_DRAM_to_SRAM + ...
    (Access_SRAM_filter_xyF * bw_filter + Access_SRAM_ifmap_xyF * bw_ifmap + Access_SRAM_psum_xyF * bw_psum) * E_SRAM_to_RF;  % Access energy (in joule)
SRAM_Access_xyF = (Access_SRAM_filter_xyF * bw_filter + Access_SRAM_ifmap_xyF * bw_ifmap + Access_SRAM_psum_xyF * bw_psum) * E_SRAM_to_RF;
DRAM_Access_xyF = (Access_DRAM_filter_xyF * bw_filter + Access_DRAM_ifmap_xyF * bw_ifmap) * E_DRAM_to_SRAM; 


%%%%%%%%%%% Calculating DRAM bandwidth induced stall count
if (slack_cycle_xy <= 0)  
    disp("CASE-1: 1 1")
    % CASE-1: There is no slack cycle from xy
    % Hence, each xy block on its own, so no stall overlapped
    % However the ifmap data already in SRAM dont need to be loaded in the reverse order formatting,
    % Therefore, the stall from ifmap wont happen from those ifmap data
    % For simplicity using xy-block for ifmap row. 
    % To be consistent with xyzF branch, in the case of this branch, if xy-1 fully fits in ifmpa-SRAM, no need to load ifmap data, otherwise will impose the stall.
    % The implication of this is that, even though due to reverse order formatting some ifmap data from xy-1 is not being loaded from DRAM,
    % the systolic array is just waiting. This is the cost in performance we pay to have simplified control.
    if (SRAM_ifmap >= (ifmap_width_effective * ifmap_height_effective * unitvol_nosof_channel * bw_ifmap))
    
        DRAM_filter_data_loaded_x = ceil(DRAM_filter_data_x/DRAM_block_size) * DRAM_block_size;  
        cycle_to_load_DRAM_filter_data_x = ceil(DRAM_filter_data_loaded_x/BW_DRAM);               % #of cycle to load the filter data from DRAM at x-direction
        DRAM_stall_filter_x = max(0, cycle_to_load_DRAM_filter_data_x - cycle_count_ideal_x);             % stall from filter data only at x-direction
        DRAM_stall_xy_reverse = DRAM_stall_filter_x + 0;  % due to reverse formatting ifmap data is not loaded, so no stall from y direction
    
        DRAM_stall_xyF = DRAM_stall_xy + ceil(DRAM_stall_xy_reverse * (F_effective/unitvol_nosof_3D_filter - 1));
        % First term: for the first set of 3D filter, no reverse order formatting applicable, so regular DRAM_stall_xy term used
        % Second term: for the next set of 3D filters, reverse order formatting applicable, so DRAM_stall_xyF_reverse used  
    else
        DRAM_stall_xyF = ceil(DRAM_stall_xy * (F_effective/unitvol_nosof_3D_filter));
        % ceil is not required for F_effective/unitvol_nosof_3D_filter since in the last pass if there is less #of filters than unitvol_nosof_3D_filter,
        % it will take less time to load that data. ceil is used for the xyF block processing
    end
elseif (DRAM_stall_xy <= 0)
    disp("CASE-2: 0 0")
    %CASE-2: Each xy block can compute itself without incorporating any stall, so no stall to complete the xyF direction
    DRAM_stall_xyF = 0;
else   
    disp("CASE-3: 1 0")
    % CASE-3: There are some slack cycles from xy, as well as some stall from xy. 
    % This means the source of stall is the (filter+ifmap) data at the first slot of xy-F block.
    % The load time of only ifmap is not creating the stall, it can be overlapped with computation, So, no need to worry about ifmap-SRAM constraint   
    
    % Now, due to reverse order formatting, at the first slot of xy-F2, no ifmap data will be loaded since the required R*W*Cuv ifmap data
    % is already in ifmap-SRAM (min ifmap-SRAM requirment). so only need to check the filter data for the first slot of xy-F2 for the stall count.
    % And, only need to check if filter-SRAM memory permits to load filter data from DRAM during those slack cycles to overlap the stall (if any).
 
    % Assumption: if the filter dada needed for the next xy-F2 block processing can fully fit in filter-SRAM only then we will load it, 
    % otherwise wont load any filter data for the next xy-F2 block
    
    % stall from filter data only at x-direction
    New_loaded_filter_PFS = filter_height * filter_width * unitvol_nosof_channel * unitvol_nosof_3D_filter * bw_filter; % in bit; filter data loaded per 3D filter set                                                                                                              % loaded to process each 3D filter set
    filter_data_loaded_PFS = ceil(New_loaded_filter_PFS/DRAM_block_size) * DRAM_block_size;
    cycle_to_load_filter_data_PFS = ceil(filter_data_loaded_PFS/BW_DRAM); % cycle count to load filter data per 3D filter set
    DRAM_stall_filter_x = max(0, cycle_to_load_filter_data_PFS - cycle_count_ideal_x);   % stall from filter data only at x-direction
    
    filter_data_Pxy = New_loaded_filter_PFS; %filter data needed per 3Dfilter set to process the xy-1 block
    
    if (DRAM_stall_filter_x <= 0)
        DRAM_stall_F = 0;            % no stall for the next xy-F2 block
        DRAM_stall_xyF = DRAM_stall_xy + DRAM_stall_F;  %first term: first xy-F1 block, second term: all next xy-F2 block
    else
        if (SRAM_filter >= (2 * filter_data_Pxy)) %There is room in SRAM_filter to fit the full filter data needed for the processing of the next xy-F2 block
            if (DRAM_stall_filter_x <= slack_cycle_xy)   
                disp("CASE 3.1 Full stall can be overlapped")
                DRAM_stall_F = 0;
            else
                disp("CASE 3.2 stall can be overlapped partially")
                DRAM_stall_Pxy = DRAM_stall_filter_x - slack_cycle_xy; % stall per xy-F2 block
                DRAM_stall_F = ceil(DRAM_stall_Pxy * (F_effective/unitvol_nosof_3D_filter - 1)); 
                % DRAM_stall_Pxy is repeated for all the next set of 3D filters,excluding the first set of filters, it will be incorporated later by adding DRAM_stall_xy
                % ceil is not required here since in the last pass if there is less #of filters than unitvol_nosof_3D_filter, it will take less time to load that data
            end
            DRAM_stall_xyF = DRAM_stall_xy + DRAM_stall_F;
        else
            %Each xy-F2 block on its own, so no stall overlapped and DRAM_stall_filter_x is repeated for all next unitvol_nosof_3D_filter except the first xy-F1 block
            DRAM_stall_xyF = DRAM_stall_xy + ceil(DRAM_stall_filter_x * (F_effective/unitvol_nosof_3D_filter - 1)); 
        end
        
    end
    
end

% Assumtion: Not doing any computation of the slack cycle after xyF, For simplicity, after xyF we assume the slack variable is zero
% keeping track of slack variable becomes complicated as we move higher depth of the tree
% This assumption is based on the fact that the ratio of stall_xy/compute_cycle_xyF is very small, 
% hence it is not worth to make the implementation very complicated for this
slack_decision_factor = DRAM_stall_xy/cycle_count_ideal_xyF;

%% Tree process along x->y->F->z direction of ifmap

% filter memeory requirement check is not necessary, new filter channels can be loaded from DRAM
% ifmap memory requirement chechk is no necessary, new ifmap channels can be loaded from DRAM

% psum memory requirement check is not necessary, the check in the previous xyF is enough

%%%%%%%Calculation of cycle count  
%The same process as xyF is repeated, for each ifmap plane with "unitvol_nosof_channel", same amount of cycle count is required
cycle_count_ideal_xyFz = cycle_count_ideal_xyF * ceil(Nos_of_channel/unitvol_nosof_channel); %without the initial offset cycles
cycle_count_xyFz = cycle_count_ideal_xyFz + initial_offset_cycle;

%%%%%% Access for filter data
gamma_DRAM_filter_xyFz = 1; % so maximal utilization of filter data from DRAM. each filter data is loaded from DRAM only once
Access_DRAM_filter_xyFz = filter_height * filter_width * Nos_of_channel * F_effective;
%The amount of access in xyF direction is repeated (Nos_of_channel/unitvol_nosof_channel) times to cover the full z-direction
Access_SRAM_filter_xyFz = Access_SRAM_filter_xyF * (Nos_of_channel/unitvol_nosof_channel);

%%%%% Access for ifmap data
%while moving towards the xyF directions, all the computations associated with the unitvol_nosof_channel number of ifmap channels is done
% Hence for z-direction, the same process as xyF is repeated 
Access_DRAM_ifmap_xyFz = Access_DRAM_ifmap_xyF * (Nos_of_channel/unitvol_nosof_channel);
%basically here maximal ifmap reuse occur if W*H*|J/S| ifmap volume fits in ifmap-SRAM, then ifmap data is loaded from DRAM once

%The amount of access in xyF direction is repeated (Nos_of_channel/unitvol_nosof_channel) times to cover the full z-direction
Access_SRAM_ifmap_xyFz = Access_SRAM_ifmap_xyF * (Nos_of_channel/unitvol_nosof_channel);


%%%%% Access for psum data
%Access_SRAM_psum_xyF = For the first set of channels, the access pattern for the F_effective number of ofmap planes is same as Access_SRAM_psum_xyF 
psum_xyFz_x = 2 * Generated_psum_effective_x * filter_height;    % psum SRAM access cost while processing the x-direction for next ifmap channels
psum_xyFz_xy = psum_xyFz_x * ofmap_height_effective;             % psum SRAM access cost while processing the xy-direction for next ifmap channels
psum_xyFz_xyF = psum_xyFz_xy * (F_effective/unitvol_nosof_3D_filter);  % psum SRAM access cost while processing the xyF-direction for next ifmap channels
Access_SRAM_psum_xyFz = Access_SRAM_psum_xyF + psum_xyFz_xyF * ceil((Nos_of_channel/unitvol_nosof_channel) - 1);
% ceil comes here becasue in the last pass even if there is less channel than unitvol_nosof_channel, same nos of psum read/write is required


%This is the access energy up to this depth of the decision
Weighted_Access_xyFz = (Access_DRAM_filter_xyFz * bw_filter + Access_DRAM_ifmap_xyFz * bw_ifmap) * E_DRAM_to_SRAM + ...
    (Access_SRAM_filter_xyFz * bw_filter + Access_SRAM_ifmap_xyFz * bw_ifmap + Access_SRAM_psum_xyFz * bw_psum) * E_SRAM_to_RF;  % Access energy (in joule)
SRAM_Access_xyFz = (Access_SRAM_filter_xyFz * bw_filter + Access_SRAM_ifmap_xyFz * bw_ifmap + Access_SRAM_psum_xyFz * bw_psum) * E_SRAM_to_RF;
DRAM_Access_xyFz = (Access_DRAM_filter_xyFz * bw_filter + Access_DRAM_ifmap_xyFz * bw_ifmap) * E_DRAM_to_SRAM; 

%%%%%%%%%%% Calculating DRAM bandwidth induced stall count
%The same process of xyF repeats for all the next set of channels.
DRAM_stall_xyFz = ceil(DRAM_stall_xyF * (Nos_of_channel/unitvol_nosof_channel));
% ceil not used for (Nos_of_channel/unitvol_nosof_channel), the argument is that, in the last pass, if there is less #of unitvol_nosof_channel,
% it will take less time to load that data. ceil is used for the full xyFz block

DRAM_stall_xyFz;
cycle_count_ideal_xyFz;
cycle_count_xyFz_with_stall = cycle_count_xyFz + DRAM_stall_xyFz;

%% Incorporating the calculation when some directions are processed partially
% In this brach all three effective can come from the psum-SRAM storage limitation: X-effective, Y-effective, and F-effective
% Z_effective never comes from psum-SRAM since Z-direction kills psum

%DRAM stall to write the ofmap value back to DRAM
ofmap_data_bits = ofmap_height_effective * ofmap_width_effective * F_effective * bw_ofmap; %in bits
ofmap_data_to_write = ceil(ofmap_data_bits/DRAM_block_size) * DRAM_block_size;
cycle_to_write_ofmap = ceil(ofmap_data_to_write/BW_DRAM); % nos of stall to write the ofmap back to DRAM [WILL ADD THE OPTION FOR POOLING HERE]
Access_DRAM_ofmap_xyFz = ofmap_height_effective * ofmap_width_effective * F_effective; % in nos of element, not bit


XYF_effective = [ofmap_width_effective ofmap_height_effective ifmap_width_effective ifmap_height_effective F_effective];
cycle_count = [cycle_count_ideal_xyFz DRAM_stall_xyFz cycle_to_write_ofmap];
SRAM_Access = [Access_SRAM_filter_xyFz Access_SRAM_ifmap_xyFz Access_SRAM_psum_xyFz];  % in Nos of element, not bit
DRAM_Access = [Access_DRAM_filter_xyFz Access_DRAM_ifmap_xyFz Access_DRAM_ofmap_xyFz]; % in Nos of element, not bit


end
