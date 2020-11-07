function [cycle_count_layer, SRAM_Access_layer, DRAM_Access_layer] = Asymmetric_FZ(Layer_param, Hardware_param, Tech_param)

%X: Direction along the width of ofmap (X = 1 for a FC layer)
%Y: Direction along the height of ofmap (Y = 1 for a FC layer)
%Z: Direction along the channels of ifmap
%F: Direction along the channels of ofmap = Direction along the filters

% For FZ branch of a FC layer: this function is the top level wrapper to cover all the iterations to compute the full ofmap volume of a layer
% To reduce data access cost, the implementation allows processing unequal data volume across the iterations
% while iterating across different directions, the implementation also considers corner cases to emulate the actual behavior of a hardware

%% Layer Specification, 
filter_height = Layer_param(1);
filter_width = Layer_param(2);
ifmap_height = Layer_param(3);
ifmap_width = Layer_param(4);
Nos_of_channel = Layer_param(5);
ofmap_height = Layer_param(6);
ofmap_width = Layer_param(7);
Nos_of_filter = Layer_param(8);
stride = Layer_param(9);
%batch_size = Layer_param(10);

%% Calling the X->Y->Z->F branch to perform computation
[F_effective, cycle_count, SRAM_Access, DRAM_Access] = FC_FtoZ(Layer_param, Hardware_param, Tech_param);

if (F_effective == Nos_of_filter)
    disp('computation calculation done')
    % None of the effective triggers, So the computation is done
    done_flag = 1;
    
    cycle_count_ideal_Fz = cycle_count(1);
    DRAM_stall_Fz = cycle_count(2);
    cycle_to_write_ofmap = cycle_count(3);

    Access_SRAM_filter_Fz = SRAM_Access(1); % in Nos of element, not bit
    Access_SRAM_ifmap_Fz = SRAM_Access(2);
    Access_SRAM_psum_Fz = SRAM_Access(3);

    Access_DRAM_filter_Fz = DRAM_Access(1);
    Access_DRAM_ifmap_Fz = DRAM_Access(2);
    Access_DRAM_ofmap_Fz = DRAM_Access(3);
else
    %effective triggers
    done_flag = 0;
end

if (done_flag == 0)
    disp('F_effective triggers') 
    
    % Computation result for the first subvolume of ofmap
    cycle_count_ideal_Fz = cycle_count(1);
    DRAM_stall_Fz = cycle_count(2);
    cycle_to_write_ofmap = cycle_count(3);

    Access_SRAM_filter_Fz = SRAM_Access(1); % in Nos of element, not bit
    Access_SRAM_ifmap_Fz = SRAM_Access(2);
    Access_SRAM_psum_Fz = SRAM_Access(3);

    Access_DRAM_filter_Fz = DRAM_Access(1);
    Access_DRAM_ifmap_Fz = DRAM_Access(2);
    Access_DRAM_ofmap_Fz = DRAM_Access(3);
    
    % Nos of same subvolume which is repeated along the F-direction
    pass_F = floor(Nos_of_filter/F_effective);
    
    % The same subvolume repeats along the F-direction
    cycle_count_ideal_Fz = cycle_count_ideal_Fz * pass_F;
    DRAM_stall_Fz = DRAM_stall_Fz * pass_F;
    cycle_to_write_ofmap = cycle_to_write_ofmap * pass_F;
        
    Access_SRAM_filter_Fz = Access_SRAM_filter_Fz * pass_F; % in Nos of element, not bit
    Access_SRAM_ifmap_Fz = Access_SRAM_ifmap_Fz * pass_F;
    Access_SRAM_psum_Fz = Access_SRAM_psum_Fz * pass_F;

    Access_DRAM_filter_Fz = Access_DRAM_filter_Fz * pass_F;
    Access_DRAM_ifmap_Fz = Access_DRAM_ifmap_Fz * pass_F;
    Access_DRAM_ofmap_Fz = Access_DRAM_ofmap_Fz * pass_F;
    
    % The remaining part along the F-direction which is not processed yet
    new_Nos_of_filter = Nos_of_filter - F_effective * pass_F;
   
    if (new_Nos_of_filter == 0) % computation is done if there is no remaining part along the F-direction
        disp('computation calculation done')
        done_flag = 1;
    end
end

% Computation for the last subvolume
if (done_flag == 0)
    disp('Computation for the last subvolume')
    % calling the function with the updated new_Nos_of_filter from the F-effective condition
    Layer_param = [filter_height filter_width ifmap_height ifmap_width Nos_of_channel ofmap_height ofmap_width new_Nos_of_filter stride];
    [F_effective, cycle_count, SRAM_Access, DRAM_Access] = FC_FtoZ(Layer_param, Hardware_param, Tech_param);
    
    cycle_count_ideal_Fz = cycle_count_ideal_Fz + cycle_count(1);
    DRAM_stall_Fz = DRAM_stall_Fz + cycle_count(2);
    cycle_to_write_ofmap = cycle_to_write_ofmap + cycle_count(3);

    Access_SRAM_filter_Fz = Access_SRAM_filter_Fz + SRAM_Access(1); % in Nos of element, not bit
    Access_SRAM_ifmap_Fz = Access_SRAM_ifmap_Fz + SRAM_Access(2);
    Access_SRAM_psum_Fz = Access_SRAM_psum_Fz + SRAM_Access(3);

    Access_DRAM_filter_Fz = Access_DRAM_filter_Fz + DRAM_Access(1);
    Access_DRAM_ifmap_Fz = Access_DRAM_ifmap_Fz + DRAM_Access(2);
    Access_DRAM_ofmap_Fz = Access_DRAM_ofmap_Fz + DRAM_Access(3);  
end

% cycle_count_ideal_Fz       
% DRAM_stall_Fz  
% cycle_to_write_ofmap
% 
% Access_SRAM_filter_Fz
% Access_SRAM_ifmap_Fz
% Access_SRAM_psum_Fz
% 
% Access_DRAM_filter_Fz
% Access_DRAM_ifmap_Fz
% Access_DRAM_ofmap_Fz

cycle_count_layer = [cycle_count_ideal_Fz DRAM_stall_Fz cycle_to_write_ofmap];
SRAM_Access_layer = [Access_SRAM_filter_Fz Access_SRAM_ifmap_Fz Access_SRAM_psum_Fz];  % in Nos of element, not bit
DRAM_Access_layer = [Access_DRAM_filter_Fz Access_DRAM_ifmap_Fz Access_DRAM_ofmap_Fz]; % in Nos of element, not bit



end
