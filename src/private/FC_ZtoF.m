function [cycle_count, SRAM_Access, DRAM_Access] = FC_ZtoF(Layer_param, Hardware_param, Tech_param)

%X: Direction along the width of ofmap (X = 1 for a FC layer)
%Y: Direction along the height of ofmap (Y = 1 for a FC layer)
%Z: Direction along the channels of ifmap
%F: Direction along the channels of ofmap = Direction along the filters

% This function performs the computation of the input layer using the ZF branch (For FC layer without batching)

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


%% Tree process along z->F direction of computation (Implementation without batching)

%% SRAM memory requirement check, Typically SRAM size will be much larger than this minimum requirement
%SRAM_filter memory requirement check
SRAM_req_filter_zF = unit_vol_filter * bw_filter;      %in bit, only unit_vol_filter need to fit which is J*K filter
if (SRAM_req_filter_zF <= SRAM_filter)
    disp("Full zF-direction possible for filter")
else
    % Applying a lower bound
    Min_SRAM_filter_zF = (SRAM_req_filter_zF)/(8 * 1024); % Minimum SRAM requirement for filter data in kB to proceed along z-direction
    disp("ALERT:Filter-SRAM is too small, increase it") %Assuming a minimum size for filter SRAM
    fprintf("Minimum SRAM memory for filter data for this layer is %f kB\n",Min_SRAM_filter_zF)
    % returning inf valus for the return parameters so that this branch does not get selected and sweep fucntion continue to run for the next hardware point
    cycle_count = [inf inf inf];
    SRAM_Access = [inf inf inf];  
    DRAM_Access = [inf inf inf]; 
    return
end

%SRAM_psum memory requirement
Generated_psum = 1 * 1 * unitvol_nosof_3D_filter;
SRAM_req_psum_zF = Generated_psum * bw_psum;   % in bit
if (SRAM_req_psum_zF <= SRAM_psum)
    disp("Full zF-direction possible for psum")
    Generated_psum_effective = Generated_psum;
else
    %Applying a lower bound on the psum-SRAM. Basically psum-SRAM has to hold minimum unitvol_nosof_3D_filter psum data at a time. 
    %Since z-direction is being processed first here, the ofmaps are formed and can be moved to DRAM.
    disp("ALERT: psum-SRAM is too small, increase it")
    Min_SRAM_psum_zF = SRAM_req_psum_zF / (8 * 1024);
    fprintf("Minimum SRAM memory for psum data is %f kB\n",Min_SRAM_psum_zF)
    % returning inf valus for the return parameters so that this branch does not get selected and sweep fucntion continue to run for the next hardware point
    cycle_count = [inf inf inf];
    SRAM_Access = [inf inf inf];  
    DRAM_Access = [inf inf inf];
    return
end

% The RSC minimum ifmap bound is not imposed for the FC layer. For the first FC layer of VGG this can require 24.5kB of ifmap SRAM
%SRAM_ifmap memory requirement check
SRAM_req_ifmap_zF = unitvol_nosof_channel * bw_ifmap;   % in bit, only unitvol_nosof_channel ifmap need to fit which is J ifmap
if (SRAM_req_ifmap_zF <= SRAM_ifmap)
    disp("full zF-direction possible for ifmap") 
else
    % Applying a lower bound
    Min_SRAM_ifmap_zF = (SRAM_req_ifmap_zF)/(8 * 1024); % Minimum SRAM requirement for filter data in kB to proceed along z-direction
    disp("ALERT:ifmap-SRAM is too small, increase it") %Assuming a minimum size for ifmap SRAM
    fprintf("Minimum SRAM memory for ifmap data for this layer is %f kB\n",Min_SRAM_ifmap_zF)
    % returning inf valus for the return parameters so that this branch does not get selected and sweep fucntion continue to run for the next hardware point
    cycle_count = [inf inf inf];
    SRAM_Access = [inf inf inf];  
    DRAM_Access = [inf inf inf];
    return   
end


%% Data access cost (SRAM & DRAM)

%%%%% Cycle count calculation
initial_offset_cycle = (used_MAC_per_column - 1);  %After first (used_MAC_per_column - 1) cycle, the pipeline is full
cycle_count_ideal_z = ceil(Nos_of_channel/unitvol_nosof_channel);  %#of cycle required to process the full z-direction
cycle_count_ideal_zF = cycle_count_ideal_z * ceil(Nos_of_filter/unitvol_nosof_3D_filter); %#of cycle required to process the full zF-direction
cycle_count_zF = cycle_count_ideal_zF + initial_offset_cycle;


%%%%%%%%%%%%% Filter data access cost
% Without batching, each filter weight is used just for one MAC operation, Hence ecah filter is only loaded once from DRAM to SRAM and from SRAM to RF
Access_DRAM_filter_zF = Nos_of_channel * Nos_of_filter;    % in element
Access_SRAM_filter_zF = Nos_of_channel * Nos_of_filter;    % in element

%%%%%%%%%%%%%%% ifmap data access cost
%%%DRAM Access (This is the similar reverse order formatiing code from the depth-4 of Z->X->Y->F branch)
if (SRAM_ifmap >= (Nos_of_channel * bw_ifmap))
    disp("full 1D volume of ifmap fits in ifmap-SRAM")
    %ifmap-SRAM has enough memory to fit the full ifmap volume, hence each ifmap data loaded from DRAM once to process the full zF direction
    Access_DRAM_ifmap_zF = Nos_of_channel;
else
    % Implementing reverse formatting for the ifmap data while processing towards F-direction  
    % the access pattern for the first set of 3D filters is repeated for all the next set of 3D filters
    if(ceil(Nos_of_filter/unitvol_nosof_3D_filter) == 1)  %(This part seems like not necessary since -1 in the last line inside 'else' will take care of this, THINK)
        %if all filters are covered in the first set of 3D filter pass then full zF is processed by loading each ifmap data once
        Access_DRAM_ifmap_zF = Nos_of_channel;
    else 
        %Some 1D ifmap data will already be in ifmap-SRAM after processing the first set of 3D filter
        ifmap_to_fit_SRAM = floor(SRAM_ifmap / bw_ifmap); %#of ifmap data fits in ifmap SRAM; Omitting the fact that DRAM load data as a block
        % In the following lines: ifmap_already_in_SRAM is the amount of ifmap data which is already in the ifmap-SRAM during the processing of the previous set of 3D-filter
        % Therefore, no need to load it since performing reverese formatting       
        first_3Dset_access = Nos_of_channel;         % access cost for the first set of 3D filters
        ifmap_already_in_SRAM = ifmap_to_fit_SRAM;
        next_3Dset_access = (Nos_of_channel - ifmap_already_in_SRAM);    % access cost for each next set of 3D filters
        % ceil comes here because even if there is less #of 3D filter in the last pass, ifmap data need to be loaded
        Access_DRAM_ifmap_zF = first_3Dset_access + next_3Dset_access * ceil(Nos_of_filter/unitvol_nosof_3D_filter - 1);
    end
end

%%%SRAM access
% at each cycle "unit_vol_ifmap" new ifmap is loaded from SRAM to RF
Access_SRAM_ifmap_z = (Nos_of_channel/unitvol_nosof_channel) * unitvol_nosof_channel; % in element
% instead of cycle_count_ideal_z used it without ceil to reflect the fact that at the last pass, when there is less nos of channel, less nos of data will
% be loaded from SRAM and some PE will remain unused.

% The same SRAM access pattern for z is repeated. The reason ceil comes here is that, in the last pass if there is less filter than "unitvol_nosof_3D_filter", 
% some PE column will reamin unused, but ifmap channels will need to be loaded from SRAM to RF and broadcasted to the PE columns which are being used.
Access_SRAM_ifmap_zF = Access_SRAM_ifmap_z * ceil(Nos_of_filter/unitvol_nosof_3D_filter);


%%%%%%% psum data access cost
% since, z is the first processing direction each psum write occur after finishing the z-direction only once
Access_SRAM_psum_zF = Nos_of_filter;
% ofmap writeback to DRAM
Access_DRAM_ofmap_zF = Nos_of_filter;

Total_SRAM_Access_bit_zF = Access_SRAM_filter_zF * bw_filter + Access_SRAM_ifmap_zF * bw_ifmap + Access_SRAM_psum_zF * bw_psum; % in bit
Total_DRAM_Access_bit_wo_ofmap_bit_zF = Access_DRAM_filter_zF * bw_filter + Access_DRAM_ifmap_zF * bw_ifmap;   % in bit


%% DRAM induced stall count calculation

%%%%% Z-direction (depth-1)
%The code is similar to what I have at depth-1 of Z->X->Y->F branch
% to process the first element of ofmap, I need to load the full z-direction of a filter and ifmap 
DRAM_filter_data_z = Nos_of_channel * unitvol_nosof_3D_filter * bw_filter;
DRAM_ifmap_data_z = Nos_of_channel * bw_ifmap;
DRAM_data_z = DRAM_filter_data_z + DRAM_ifmap_data_z; % in bit, this is the number of data loaded from DRAM to process the first ofmap element
%cycle_to_load_DRAM_data_x = ceil(DRAM_data_x/BW_DRAM)
DRAM_data_loaded_z = ceil(DRAM_data_z/DRAM_block_size) * DRAM_block_size;   % The data loaded from DRAM is an integer multiple of the DRAM block size
cycle_to_load_DRAM_data_z = ceil(DRAM_data_loaded_z/BW_DRAM);               % #of cycle to load the data from DRAM


if cycle_to_load_DRAM_data_z > cycle_count_ideal_z
    DRAM_stall_z = cycle_to_load_DRAM_data_z - cycle_count_ideal_z;
else
    DRAM_stall_z = 0;
end
% For any feasible DRAM bandwidth, FC layer without bactching will not have any slack cycle. Hence, omitting the calculating of slack cycles


%%%%% Z->F direction (depth-2)
DRAM_filter_data_F = Nos_of_channel * unitvol_nosof_3D_filter * bw_filter; % New filter data loaded while moving towards F-direction
if (SRAM_ifmap >= (Nos_of_channel * bw_ifmap))
    %ifmap-SRAM has enough memory to fit the full ifmap volume, hence while moving towards F-direction, no need to load ifmap data,    
    DRAM_ifmap_data_F = 0;
else
    ifmap_to_fit_SRAM = floor(SRAM_ifmap / bw_ifmap);
    DRAM_ifmap_data_F = (Nos_of_channel - ifmap_to_fit_SRAM) * bw_ifmap; % newly loaded ifmap data while moving towards F-direction due to reverse order formatting
end
DRAM_data_F = DRAM_filter_data_F + DRAM_ifmap_data_F;
DRAM_data_loaded_F = ceil(DRAM_data_F/DRAM_block_size) * DRAM_block_size;   % The data loaded from DRAM is an integer multiple of the DRAM block size
cycle_to_load_DRAM_data_F = ceil(DRAM_data_loaded_F/BW_DRAM);               % #of cycle to load the data from DRAM for the next set of filters

DRAM_stall_F = max(0, cycle_to_load_DRAM_data_F - cycle_count_ideal_z);   % DARM stall for the next set of 1D filters
DRAM_stall_zF = DRAM_stall_z + ceil(DRAM_stall_F * (Nos_of_filter/unitvol_nosof_3D_filter - 1));
% first term: for the first set of 1D filter set
% second term: for all the next set of 1D filter set
% ceil not used for Nos_of_filter/unitvol_nosof_3D_filter, the argument is that, in the last pass, if there is less #of unitvol_nosof_3D_filter, 
% it will take less time to load that data. ceil is used for the full zF block


%% Parameters to return
% In this brach no effective can come from the psum-SRAM storage limitation since Z is the first branch

%DRAM stall to write the ofmap value back to DRAM
ofmap_data_bits = Nos_of_filter * bw_ofmap; %in bits
ofmap_data_to_write = ceil(ofmap_data_bits/DRAM_block_size) * DRAM_block_size;
cycle_to_write_ofmap = ceil(ofmap_data_to_write/BW_DRAM); % nos of stall to write the ofmap back to DRAM [WILL ADD THE OPTION FOR POOLING HERE]
Access_DRAM_ofmap_zF = Nos_of_filter; % in nos of element, not bit


cycle_count = [cycle_count_ideal_zF DRAM_stall_zF cycle_to_write_ofmap];
SRAM_Access = [Access_SRAM_filter_zF Access_SRAM_ifmap_zF Access_SRAM_psum_zF];  % in Nos of element, not bit
DRAM_Access = [Access_DRAM_filter_zF Access_DRAM_ifmap_zF Access_DRAM_ofmap_zF]; % in Nos of element, not bit


end