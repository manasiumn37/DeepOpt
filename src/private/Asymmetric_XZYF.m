function [cycle_count_layer, SRAM_Access_layer, DRAM_Access_layer] = Asymmetric_XZYF(Layer_param, Hardware_param, Tech_param)

%X: Direction along the width of ofmap
%Y: Direction along the height of ofmap
%Z: Direction along the channels of ifmap
%F: Direction along the channels of ofmap = Direction along the 3D filters

% For XZYF branch: this function is the top level wrapper to cover all the iterations to compute the full ofmap volume of a layer
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
batch_size = Layer_param(10);

%% Calling the X->Y->Z->F branch to perform computation
[X_effective, cycle_count, SRAM_Access, DRAM_Access] = XtoZtoYtoF(Layer_param, Hardware_param, Tech_param);

% Getting the value of the x effective parameters
ofmap_width_effective = X_effective(1);
ifmap_width_effective = X_effective(2);

if (ofmap_width_effective == ofmap_width)
    disp('computation calculation done')
    % None of the effective triggers, So the computation is done
    done_flag = 1;
    
    cycle_count_ideal_xzyF = cycle_count(1);
    DRAM_stall_xzyF = cycle_count(2);
    cycle_to_write_ofmap = cycle_count(3);

    Access_SRAM_filter_xzyF = SRAM_Access(1); % in Nos of element, not bit
    Access_SRAM_ifmap_xzyF = SRAM_Access(2);
    Access_SRAM_psum_xzyF = SRAM_Access(3);

    Access_DRAM_filter_xzyF = DRAM_Access(1);
    Access_DRAM_ifmap_xzyF = DRAM_Access(2);
    Access_DRAM_ofmap_xzyF = DRAM_Access(3);
else
    % effective triggers.
    done_flag = 0;
end

if (done_flag == 0) && (ofmap_width_effective < ofmap_width) % the second comparison is not required, the first one is enough
    disp('x_effective triggers') 
    
    % Computation result for the first subvolume of ifmap/ofmap
    cycle_count_ideal_xzyF = cycle_count(1);
    DRAM_stall_xzyF = cycle_count(2);
    cycle_to_write_ofmap = cycle_count(3);

    Access_SRAM_filter_xzyF = SRAM_Access(1); % in Nos of element, not bit
    Access_SRAM_ifmap_xzyF = SRAM_Access(2);
    Access_SRAM_psum_xzyF = SRAM_Access(3);

    Access_DRAM_filter_xzyF = DRAM_Access(1);
    Access_DRAM_ifmap_xzyF = DRAM_Access(2);
    Access_DRAM_ofmap_xzyF = DRAM_Access(3);
    
    % Nos of same subvolume which is repeated along the x-direction
    pass_x = floor(ofmap_width/ofmap_width_effective);
    
    % The same subvolume repeats along the x-direction
    cycle_count_ideal_xzyF = cycle_count_ideal_xzyF * pass_x;
    DRAM_stall_xzyF = DRAM_stall_xzyF * pass_x;
    cycle_to_write_ofmap = cycle_to_write_ofmap * pass_x;
        
    Access_SRAM_filter_xzyF = Access_SRAM_filter_xzyF * pass_x; % in Nos of element, not bit
    Access_SRAM_ifmap_xzyF = Access_SRAM_ifmap_xzyF * pass_x;
    Access_SRAM_psum_xzyF = Access_SRAM_psum_xzyF * pass_x;

    Access_DRAM_filter_xzyF = Access_DRAM_filter_xzyF * pass_x;
    Access_DRAM_ifmap_xzyF = Access_DRAM_ifmap_xzyF * pass_x;
    Access_DRAM_ofmap_xzyF = Access_DRAM_ofmap_xzyF * pass_x;
    
    % The remaining part along the x-direction which is not processed yet
    new_ofmap_width = ofmap_width - ofmap_width_effective * pass_x;
    new_ifmap_width = (new_ofmap_width - 1) * stride + filter_width;
   
    if (new_ofmap_width == 0) % computation is done if there is no remaining part along the x-direction
        disp('computation calculation done')
        done_flag = 1;
    end
end

% Computation for the last subvolume
if (done_flag == 0)
    disp('Computation for the last subvolume')
    % calling the function with the updated ifmap/ofmap width from the x-effective condition
    Layer_param = [filter_height filter_width ifmap_height new_ifmap_width Nos_of_channel ofmap_height new_ofmap_width Nos_of_filter stride];
    [XY_effective, cycle_count, SRAM_Access, DRAM_Access] = XtoZtoYtoF(Layer_param, Hardware_param, Tech_param);
    
    cycle_count_ideal_xzyF = cycle_count_ideal_xzyF + cycle_count(1);
    DRAM_stall_xzyF = DRAM_stall_xzyF + cycle_count(2);
    cycle_to_write_ofmap = cycle_to_write_ofmap + cycle_count(3);

    Access_SRAM_filter_xzyF = Access_SRAM_filter_xzyF + SRAM_Access(1); % in Nos of element, not bit
    Access_SRAM_ifmap_xzyF = Access_SRAM_ifmap_xzyF + SRAM_Access(2);
    Access_SRAM_psum_xzyF = Access_SRAM_psum_xzyF + SRAM_Access(3);

    Access_DRAM_filter_xzyF = Access_DRAM_filter_xzyF + DRAM_Access(1);
    Access_DRAM_ifmap_xzyF = Access_DRAM_ifmap_xzyF + DRAM_Access(2);
    Access_DRAM_ofmap_xzyF = Access_DRAM_ofmap_xzyF + DRAM_Access(3);  
end

% cycle_count_ideal_xzyF       
% DRAM_stall_xzyF  
% cycle_to_write_ofmap
% 
% Access_SRAM_filter_xzyF
% Access_SRAM_ifmap_xzyF
% Access_SRAM_psum_xzyF
% 
% Access_DRAM_filter_xzyF
% Access_DRAM_ifmap_xzyF
% Access_DRAM_ofmap_xzyF


cycle_count_layer = [cycle_count_ideal_xzyF DRAM_stall_xzyF cycle_to_write_ofmap];
SRAM_Access_layer = [Access_SRAM_filter_xzyF Access_SRAM_ifmap_xzyF Access_SRAM_psum_xzyF];  % in Nos of element, not bit
DRAM_Access_layer = [Access_DRAM_filter_xzyF Access_DRAM_ifmap_xzyF Access_DRAM_ofmap_xzyF]; % in Nos of element, not bit


end
