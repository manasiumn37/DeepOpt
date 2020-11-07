function [F_effective, cycle_count, SRAM_Access, DRAM_Access] = FC_FtoZ(Layer_param, Hardware_param, Tech_param)

%X: Direction along the width of ofmap (X = 1 for a FC layer)
%Y: Direction along the height of ofmap (Y = 1 for a FC layer)
%Z: Direction along the channels of ifmap
%F: Direction along the channels of ofmap = Direction along the filters

% This function performs the computation of the input layer using the FZ branch (For FC layer without batching)

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


%% Tree process along F->z direction of computation (Implementation without batching)

%% SRAM memory requirement check, Typically SRAM size will be much larger than this minimum requirement
%SRAM_filter memory requirement check
SRAM_req_filter_Fz = unit_vol_filter * bw_filter;      %in bit, only unit_vol_filter need to fit which is J*K filter
if (SRAM_req_filter_Fz <= SRAM_filter)
    disp("Full Fz-direction possible for filter")
else
    % Applying a lower bound
    Min_SRAM_filter_Fz = (SRAM_req_filter_Fz)/(8 * 1024); % Minimum SRAM requirement for filter data in kB to proceed along z-direction
    disp("ALERT:Filter-SRAM is too small, increase it") %Assuming a minimum size for filter SRAM
    fprintf("Minimum SRAM memory for filter data for this layer is %f kB\n",Min_SRAM_filter_Fz)
    % returning inf valus for the return parameters so that this branch does not get selected and sweep fucntion continue to run for the next hardware point
    F_effective = Nos_of_filter;  % this will prevent calling the asymmetric subvolumes
    cycle_count = [inf inf inf];
    SRAM_Access = [inf inf inf];  
    DRAM_Access = [inf inf inf];
    return
end

%SRAM_ifmap memory requirement check
SRAM_req_ifmap_Fz = unitvol_nosof_channel * bw_ifmap;   % in bit, only unitvol_nosof_channel ifmap need to fit which is J ifmap
if (SRAM_req_ifmap_Fz <= SRAM_ifmap)
    disp("full Fz-direction possible for ifmap") 
else
    % Applying a lower bound
    Min_SRAM_ifmap_Fz = (SRAM_req_ifmap_Fz)/(8 * 1024); % Minimum SRAM requirement for filter data in kB to proceed along z-direction
    disp("ALERT:ifmap-SRAM is too small, increase it") %Assuming a minimum size for ifmap SRAM
    fprintf("Minimum SRAM memory for ifmap data for this layer is %f kB\n",Min_SRAM_ifmap_Fz)
    % returning inf valus for the return parameters so that this branch does not get selected and sweep fucntion continue to run for the next hardware point
    F_effective = Nos_of_filter;  % this will prevent calling the asymmetric subvolumes
    cycle_count = [inf inf inf];
    SRAM_Access = [inf inf inf];  
    DRAM_Access = [inf inf inf];
    return   
end

%SRAM_psum memory requirement 
SRAM_req_psum_Fz = Nos_of_filter * bw_psum;   % in bit
if (SRAM_req_psum_Fz <= SRAM_psum)
    disp("full Fz-direction possible for psum-SRAM")
    F_effective = Nos_of_filter;
else
   disp("Partial F-direction is being processed due to psum-SRAM storage requirement")
   F_effective_to_fit = floor(SRAM_psum /bw_psum);
   % for the sake of regular processing pattern, F_effective is chosen to be an integer multiple of unitvol_nosof_3D_filter.
   F_effective = floor(F_effective_to_fit/unitvol_nosof_3D_filter) * unitvol_nosof_3D_filter;
end


%% Data access cost (SRAM & DRAM)

%%%%% Cycle count calculation
initial_offset_cycle = (used_MAC_per_column - 1);  %After first (used_MAC_per_column - 1) cycle, the pipeline is full
cycle_count_ideal_F = ceil(F_effective/unitvol_nosof_3D_filter);  %#of cycle required to process the full F-direction
cycle_count_ideal_Fz = cycle_count_ideal_F * ceil(Nos_of_channel/unitvol_nosof_channel); %#of cycle required to process the full zF-direction
cycle_count_Fz = cycle_count_ideal_Fz + initial_offset_cycle;

%%%%%%%%%%%%% Filter data access cost
% Without batching, each filter weight is used just for one MAC operation, Hence ecah filter is only loaded once from DRAM to SRAM and from SRAM to RF
Access_DRAM_filter_Fz = Nos_of_channel * F_effective;    % in element
Access_SRAM_filter_Fz = Nos_of_channel * F_effective;    % in element

%%%%%%%%%%%% ifmap data access cost
% All computations associated with the "unitvol_nosof_channel" ifmap channels are done after loading it from DRAM once. The same is true for all ifmap channels
% while moving towards z-direction. Each ifmap is also loaded from SRAM to RF only once. 
Access_DRAM_ifmap_Fz = Nos_of_channel;  % in element
Access_SRAM_ifmap_Fz = Nos_of_channel;

%%%%%%%%%%%% psum data access cost
% when going towards the F-direction, for the first set of unitvol_nosof_channel, for each filter one psum write to SRAM occur.
% For all the next set of channels, both psum read and write occur
Access_SRAM_psum_F = F_effective;
Access_SRAM_psum_Fz = F_effective + 2 * F_effective * ceil((Nos_of_channel/unitvol_nosof_channel) - 1);

% ofmap writeback to DRAM
Access_DRAM_ofmap_Fz = F_effective;

Total_SRAM_Access_bit_Fz = Access_SRAM_filter_Fz * bw_filter + Access_SRAM_ifmap_Fz * bw_ifmap + Access_SRAM_psum_Fz * bw_psum; % in bit
Total_DRAM_Access_bit_wo_ofmap_bit_Fz = Access_DRAM_filter_Fz * bw_filter + Access_DRAM_ifmap_Fz * bw_ifmap;   % in bit

%% DRAM induced stall count calculation (This is very straightforward)

%%%% F-direction
% To move towards F-direction, unitvol_nosof_channel ifmap and filter data loaded once from DRAM
DRAM_filter_data_F = unitvol_nosof_channel * F_effective * bw_filter;
DRAM_ifmap_data_F = unitvol_nosof_channel * bw_ifmap;

DRAM_data_F = DRAM_filter_data_F + DRAM_ifmap_data_F; % in bit, this is the number of data loaded from DRAM to process the first set of ifmap channels
DRAM_data_loaded_F = ceil(DRAM_data_F/DRAM_block_size) * DRAM_block_size;   % The data loaded from DRAM is an integer multiple of the DRAM block size
cycle_to_load_DRAM_data_F = ceil(DRAM_data_loaded_F/BW_DRAM);               % #of cycle to load the data from DRAM

DRAM_stall_F = max(0, cycle_to_load_DRAM_data_F - cycle_count_ideal_F);
% For any feasible DRAM bandwidth, FC layer without batching will not have any slack cycle. Hence, omitting the calculating of slack cycles

%%%% Z-direction
% DRAM_stall_F is repeated to cover the full z-direction
DRAM_stall_Fz = ceil(DRAM_stall_F * (Nos_of_channel/unitvol_nosof_channel));
% ceil not used for (Nos_of_channel/unitvol_nosof_channel), the argument is that, in the last pass, if there is less #of unitvol_nosof_channel,
% it will take less time to load that data. ceil is used for the full Fz block


%% Incorporating the calculation when some directions are processed partially
% In this brach only one effective can come from the psum-SRAM storage limitation: F-effective

%DRAM stall to write the ofmap value back to DRAM
ofmap_data_bits = F_effective * bw_ofmap; %in bits
ofmap_data_to_write = ceil(ofmap_data_bits/DRAM_block_size) * DRAM_block_size;
cycle_to_write_ofmap = ceil(ofmap_data_to_write/BW_DRAM); % nos of stall to write the ofmap back to DRAM [WILL ADD THE OPTION FOR POOLING HERE]
Access_DRAM_ofmap_Fz = F_effective; % in nos of element, not bit


F_effective;
cycle_count = [cycle_count_ideal_Fz DRAM_stall_Fz cycle_to_write_ofmap];
SRAM_Access = [Access_SRAM_filter_Fz Access_SRAM_ifmap_Fz Access_SRAM_psum_Fz];  % in Nos of element, not bit
DRAM_Access = [Access_DRAM_filter_Fz Access_DRAM_ifmap_Fz Access_DRAM_ofmap_Fz]; % in Nos of element, not bit

end
