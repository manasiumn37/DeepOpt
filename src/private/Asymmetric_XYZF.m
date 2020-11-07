function [cycle_count_layer, SRAM_Access_layer, DRAM_Access_layer] = Asymmetric_XYZF(Layer_param, Hardware_param, Tech_param)

%X: Direction along the width of ofmap
%Y: Direction along the height of ofmap
%Z: Direction along the channels of ifmap
%F: Direction along the channels of ofmap = Direction along the 3D filters

% For XYZF branch: this function is the top level wrapper to cover all the iterations to compute the full ofmap volume of a layer
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
[XY_effective, cycle_count, SRAM_Access, DRAM_Access] = XtoYtoZtoF(Layer_param, Hardware_param, Tech_param);

% Getting the value of the x,y effective parameters
ofmap_width_effective = XY_effective(1);
ofmap_height_effective = XY_effective(2);
ifmap_width_effective = XY_effective(3);
ifmap_height_effective = XY_effective(4);


if (ofmap_width_effective == ofmap_width) && (ofmap_height_effective == ofmap_height)
    disp('computation calculation done')
    % None of the effective triggers, So the computation is done
    done_flag = 1;
    
    cycle_count_ideal_xyzF = cycle_count(1);
    DRAM_stall_xyzF = cycle_count(2);
    cycle_to_write_ofmap = cycle_count(3);

    Access_SRAM_filter_xyzF = SRAM_Access(1); % in Nos of element, not bit
    Access_SRAM_ifmap_xyzF = SRAM_Access(2);
    Access_SRAM_psum_xyzF = SRAM_Access(3);

    Access_DRAM_filter_xyzF = DRAM_Access(1);
    Access_DRAM_ifmap_xyzF = DRAM_Access(2);
    Access_DRAM_ofmap_xyzF = DRAM_Access(3);
else
    % Either both or one of the effective triggers. Lets start with assigining zero to all the computation variables
    done_flag = 0;
    
    cycle_count_ideal_xyzF = 0;
    DRAM_stall_xyzF = 0;
    cycle_to_write_ofmap = 0;

    Access_SRAM_filter_xyzF = 0; 
    Access_SRAM_ifmap_xyzF = 0;
    Access_SRAM_psum_xyzF = 0;

    Access_DRAM_filter_xyzF = 0;
    Access_DRAM_ifmap_xyzF = 0;
    Access_DRAM_ofmap_xyzF = 0;
end

new_ofmap_width = ofmap_width;
new_ofmap_height = ofmap_height;
new_ifmap_width = (new_ofmap_width - 1) * stride + filter_width;
new_ifmap_height = (new_ofmap_height - 1) * stride + filter_height;

if (ofmap_width_effective < ofmap_width) % X_effective triggers, so Y_effective = 1
    disp('x_effective triggers, so y_effective = 1') 
    
    % Computation result for the first subvolume of ifmap/ofmap
    cycle_count_ideal_xyzF = cycle_count(1);
    DRAM_stall_xyzF = cycle_count(2);
    cycle_to_write_ofmap = cycle_count(3);

    Access_SRAM_filter_xyzF = SRAM_Access(1); % in Nos of element, not bit
    Access_SRAM_ifmap_xyzF = SRAM_Access(2);
    Access_SRAM_psum_xyzF = SRAM_Access(3);

    Access_DRAM_filter_xyzF = DRAM_Access(1);
    Access_DRAM_ifmap_xyzF = DRAM_Access(2);
    Access_DRAM_ofmap_xyzF = DRAM_Access(3);
    
    % Nos of same subvolume which is repeated along the x-direction
    pass_x = floor(ofmap_width/ofmap_width_effective);
    
    % Since y-effective = 1, the same subvolume also repeats along the full y-direction
    cycle_count_ideal_xyzF = cycle_count_ideal_xyzF * ofmap_height * pass_x;
    DRAM_stall_xyzF = DRAM_stall_xyzF * ofmap_height * pass_x;
    cycle_to_write_ofmap = cycle_to_write_ofmap * ofmap_height * pass_x;
        
    Access_SRAM_filter_xyzF = Access_SRAM_filter_xyzF * ofmap_height * pass_x; % in Nos of element, not bit
    Access_SRAM_ifmap_xyzF = Access_SRAM_ifmap_xyzF * ofmap_height * pass_x;
    Access_SRAM_psum_xyzF = Access_SRAM_psum_xyzF * ofmap_height * pass_x;

    Access_DRAM_filter_xyzF = Access_DRAM_filter_xyzF * ofmap_height * pass_x;
    Access_DRAM_ifmap_xyzF = Access_DRAM_ifmap_xyzF * ofmap_height * pass_x;
    Access_DRAM_ofmap_xyzF = Access_DRAM_ofmap_xyzF * ofmap_height * pass_x;
    
    % The remaining part along the x-direction which is not processed yet
    new_ofmap_width = ofmap_width - ofmap_width_effective * pass_x;
    new_ifmap_width = (new_ofmap_width - 1) * stride + filter_width;
   
    if (new_ofmap_width == 0) % computation is done if there is no remaining part along the x-direction
        disp('computation calculation done')
        done_flag = 1;
    end
end

% This condition is triggered when the computation is not done yet and hence full y-direction is not covered yet
if (done_flag == 0) && (ofmap_height_effective < ofmap_height)
    disp('y_effective triggers')
    % calling the function with the updated ifmap/ofmap width if the x-effective condition updateds them or with the original ifmap/ofmap width if not updated
    Layer_param = [filter_height filter_width ifmap_height new_ifmap_width Nos_of_channel ofmap_height new_ofmap_width Nos_of_filter stride];
    [XY_effective, cycle_count, SRAM_Access, DRAM_Access] = XtoYtoZtoF(Layer_param, Hardware_param, Tech_param);
    
    ofmap_width_effective = XY_effective(1);
    ofmap_height_effective = XY_effective(2);
    ifmap_width_effective = XY_effective(3);
    ifmap_height_effective = XY_effective(4);
    
    % Nos of same subvolume which is repeated along the y-direction
    pass_y = floor(ofmap_height/ofmap_height_effective);
    
    % Computation result of the subvolume
    cycle_count_ideal_xyzF_sub = cycle_count(1);
    DRAM_stall_xyzF_sub = cycle_count(2);
    cycle_to_write_ofmap_sub = cycle_count(3);
    
    Access_SRAM_filter_xyzF_sub = SRAM_Access(1); % in Nos of element, not bit
    Access_SRAM_ifmap_xyzF_sub = SRAM_Access(2);
    Access_SRAM_psum_xyzF_sub = SRAM_Access(3);

    Access_DRAM_filter_xyzF_sub = DRAM_Access(1);
    Access_DRAM_ifmap_xyzF_sub = DRAM_Access(2);
    Access_DRAM_ofmap_xyzF_sub = DRAM_Access(3);
    
    % The same subvolume repeats along the y-direction
    % Adding the new computation results with the previously computed results
    % Note that if the x-effective condition is not triggered then the first terms are set to zero beforehand
    cycle_count_ideal_xyzF = cycle_count_ideal_xyzF + cycle_count_ideal_xyzF_sub * pass_y;
    DRAM_stall_xyzF = DRAM_stall_xyzF + DRAM_stall_xyzF_sub * pass_y;
    cycle_to_write_ofmap = cycle_to_write_ofmap + cycle_to_write_ofmap_sub * pass_y;

    Access_SRAM_filter_xyzF = Access_SRAM_filter_xyzF + Access_SRAM_filter_xyzF_sub * pass_y; % in Nos of element, not bit
    Access_SRAM_ifmap_xyzF = Access_SRAM_ifmap_xyzF + Access_SRAM_ifmap_xyzF_sub * pass_y;
    Access_SRAM_psum_xyzF = Access_SRAM_psum_xyzF + Access_SRAM_psum_xyzF_sub * pass_y;

    Access_DRAM_filter_xyzF = Access_DRAM_filter_xyzF + Access_DRAM_filter_xyzF_sub * pass_y;
    Access_DRAM_ifmap_xyzF = Access_DRAM_ifmap_xyzF + Access_DRAM_ifmap_xyzF_sub * pass_y;
    Access_DRAM_ofmap_xyzF = Access_DRAM_ofmap_xyzF + Access_DRAM_ofmap_xyzF_sub * pass_y;

    % The remaining part along the y-direction which is not processed yet
    new_ofmap_height = ofmap_height - ofmap_height_effective * pass_y;
    new_ifmap_height = (new_ofmap_height - 1) * stride + filter_height;
    
    if (new_ofmap_height == 0) % computation is done if there is no remaining part along the y-direction
        disp('computation calculation done')
        done_flag = 1;
    end  
end

% Computation for the last subvolume
if (done_flag == 0)
    disp('Computation for the last subvolume')
    % calling the function with the updated ifmap/ofmap width if the x-effective condition updateds them or with the original ifmap/ofmap width if not updated
    % and with the updated height from the y-effective condition
    Layer_param = [filter_height filter_width new_ifmap_height new_ifmap_width Nos_of_channel new_ofmap_height new_ofmap_width Nos_of_filter stride];
    [XY_effective, cycle_count, SRAM_Access, DRAM_Access] = XtoYtoZtoF(Layer_param, Hardware_param, Tech_param);
    
    cycle_count_ideal_xyzF = cycle_count_ideal_xyzF + cycle_count(1);
    DRAM_stall_xyzF = DRAM_stall_xyzF + cycle_count(2);
    cycle_to_write_ofmap = cycle_to_write_ofmap + cycle_count(3);

    Access_SRAM_filter_xyzF = Access_SRAM_filter_xyzF + SRAM_Access(1); % in Nos of element, not bit
    Access_SRAM_ifmap_xyzF = Access_SRAM_ifmap_xyzF + SRAM_Access(2);
    Access_SRAM_psum_xyzF = Access_SRAM_psum_xyzF + SRAM_Access(3);

    Access_DRAM_filter_xyzF = Access_DRAM_filter_xyzF + DRAM_Access(1);
    Access_DRAM_ifmap_xyzF = Access_DRAM_ifmap_xyzF + DRAM_Access(2);
    Access_DRAM_ofmap_xyzF = Access_DRAM_ofmap_xyzF + DRAM_Access(3);
    
end
        
% cycle_count_ideal_xyzF       
% DRAM_stall_xyzF  
% cycle_to_write_ofmap
% 
% Access_SRAM_filter_xyzF
% Access_SRAM_ifmap_xyzF
% Access_SRAM_psum_xyzF
% 
% Access_DRAM_filter_xyzF
% Access_DRAM_ifmap_xyzF
% Access_DRAM_ofmap_xyzF

cycle_count_layer = [cycle_count_ideal_xyzF DRAM_stall_xyzF cycle_to_write_ofmap];
SRAM_Access_layer = [Access_SRAM_filter_xyzF Access_SRAM_ifmap_xyzF Access_SRAM_psum_xyzF];  % in Nos of element, not bit
DRAM_Access_layer = [Access_DRAM_filter_xyzF Access_DRAM_ifmap_xyzF Access_DRAM_ofmap_xyzF]; % in Nos of element, not bit

        
end       
        
        
        
        
 


