function [cycle_count_layer, SRAM_Access_layer, DRAM_Access_layer] = Asymmetric_XYFZ(Layer_param, Hardware_param, Tech_param)

%X: Direction along the width of ofmap
%Y: Direction along the height of ofmap
%Z: Direction along the channels of ifmap
%F: Direction along the channels of ofmap = Direction along the 3D filters

% For XYFZ branch: this function is the top level wrapper to cover all the iterations to compute the full ofmap volume of a layer
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

%% Calling the X->Y->F->Z branch to perform computation
[XYF_effective, cycle_count, SRAM_Access, DRAM_Access] = XtoYtoFtoZ(Layer_param, Hardware_param, Tech_param);

% Getting the value of the x,y,F effective parameters
ofmap_width_effective = XYF_effective(1);
ofmap_height_effective = XYF_effective(2);
ifmap_width_effective = XYF_effective(3);
ifmap_height_effective = XYF_effective(4);
F_effective = XYF_effective(5);
% When F_effective triggers, F_effective >= Nos of column in the systolic array

if (ofmap_width_effective == ofmap_width) && (ofmap_height_effective == ofmap_height) && (F_effective == Nos_of_filter)
    disp('computation calculation done')
    % None of the effective triggers, So the computation is done
    done_flag = 1;
    
    cycle_count_ideal_xyFz = cycle_count(1);
    DRAM_stall_xyFz = cycle_count(2);
    cycle_to_write_ofmap = cycle_count(3);

    Access_SRAM_filter_xyFz = SRAM_Access(1); % in Nos of element, not bit
    Access_SRAM_ifmap_xyFz = SRAM_Access(2);
    Access_SRAM_psum_xyFz = SRAM_Access(3);

    Access_DRAM_filter_xyFz = DRAM_Access(1);
    Access_DRAM_ifmap_xyFz = DRAM_Access(2);
    Access_DRAM_ofmap_xyFz = DRAM_Access(3);
else
    % One or more effective had triggered. Lets start with assigining zero to all the computation variables
    done_flag = 0;
    
    cycle_count_ideal_xyFz = 0;
    DRAM_stall_xyFz = 0;
    cycle_to_write_ofmap = 0;

    Access_SRAM_filter_xyFz = 0; 
    Access_SRAM_ifmap_xyFz = 0;
    Access_SRAM_psum_xyFz = 0;

    Access_DRAM_filter_xyFz = 0;
    Access_DRAM_ifmap_xyFz = 0;
    Access_DRAM_ofmap_xyFz = 0;
end

new_ofmap_width = ofmap_width;
new_ofmap_height = ofmap_height;
new_ifmap_width = (new_ofmap_width - 1) * stride + filter_width;
new_ifmap_height = (new_ofmap_height - 1) * stride + filter_height;
new_Nos_of_filter = Nos_of_filter;

if (ofmap_width_effective < ofmap_width) % X_effective triggers, so Y_effective = 1 and F_effective = unitvol_nosof_3D_filter
    disp('x_effective triggers, so y_effective = 1 and F_effective = unitvol_nosof_3D_filter') 
    
    % Computation result for the first subvolume of ifmap/ofmap
    cycle_count_ideal_xyFz = cycle_count(1);
    DRAM_stall_xyFz = cycle_count(2);
    cycle_to_write_ofmap = cycle_count(3);

    Access_SRAM_filter_xyFz = SRAM_Access(1); % in Nos of element, not bit
    Access_SRAM_ifmap_xyFz = SRAM_Access(2);
    Access_SRAM_psum_xyFz = SRAM_Access(3);

    Access_DRAM_filter_xyFz = DRAM_Access(1);
    Access_DRAM_ifmap_xyFz = DRAM_Access(2);
    Access_DRAM_ofmap_xyFz = DRAM_Access(3);
    
    % Nos of same subvolume which is repeated along the x-direction
    pass_x = floor(ofmap_width/ofmap_width_effective);
    
    % Since y-effective = 1, the same subvolume also repeats along the full y-direction
    cycle_count_ideal_xyFz = cycle_count_ideal_xyFz * ofmap_height * pass_x;
    DRAM_stall_xyFz = DRAM_stall_xyFz * ofmap_height * pass_x;
    cycle_to_write_ofmap = cycle_to_write_ofmap * ofmap_height * pass_x;
        
    Access_SRAM_filter_xyFz = Access_SRAM_filter_xyFz * ofmap_height * pass_x; % in Nos of element, not bit
    Access_SRAM_ifmap_xyFz = Access_SRAM_ifmap_xyFz * ofmap_height * pass_x;
    Access_SRAM_psum_xyFz = Access_SRAM_psum_xyFz * ofmap_height * pass_x;

    Access_DRAM_filter_xyFz = Access_DRAM_filter_xyFz * ofmap_height * pass_x;
    Access_DRAM_ifmap_xyFz = Access_DRAM_ifmap_xyFz * ofmap_height * pass_x;
    Access_DRAM_ofmap_xyFz = Access_DRAM_ofmap_xyFz * ofmap_height * pass_x;
    
    % The remaining part along the x-direction which is not processed yet
    new_ofmap_width = ofmap_width - ofmap_width_effective * pass_x;
    new_ifmap_width = (new_ofmap_width - 1) * stride + filter_width;
    new_Nos_of_filter = F_effective;
   
    if (new_ofmap_width == 0) % computation is done if there is no remaining part along the x-direction
        disp('computation calculation done for the full xy direction')
        %done_flag = 1;
    elseif (new_ofmap_width > 0) && (ofmap_height_effective == ofmap_height)
        % This condition triggers only for the FC layer with batching when Y = 1. In this situation the
        % following original Y-effective condition will not trigger since Y = 1. This code will compute the remaining X-directional subvolume
        disp('Computation for the last X directional subvolume')
        % calling the function with the updated ifmap/ofmap width from the x-effective condition
        % Also using the original F_effective to call the function
        Layer_param = [filter_height filter_width ifmap_height new_ifmap_width Nos_of_channel ofmap_height new_ofmap_width F_effective stride];
        [XYF_effective, cycle_count, SRAM_Access, DRAM_Access] = XtoYtoFtoZ(Layer_param, Hardware_param, Tech_param);
        
        cycle_count_ideal_xyFz = cycle_count_ideal_xyFz + cycle_count(1);
        DRAM_stall_xyFz = DRAM_stall_xyFz + cycle_count(2);
        cycle_to_write_ofmap = cycle_to_write_ofmap + cycle_count(3);

        Access_SRAM_filter_xyFz = Access_SRAM_filter_xyFz + SRAM_Access(1); % in Nos of element, not bit
        Access_SRAM_ifmap_xyFz = Access_SRAM_ifmap_xyFz + SRAM_Access(2);
        Access_SRAM_psum_xyFz = Access_SRAM_psum_xyFz + SRAM_Access(3);

        Access_DRAM_filter_xyFz = Access_DRAM_filter_xyFz + DRAM_Access(1);
        Access_DRAM_ifmap_xyFz = Access_DRAM_ifmap_xyFz + DRAM_Access(2);
        Access_DRAM_ofmap_xyFz = Access_DRAM_ofmap_xyFz + DRAM_Access(3);
        
    end
end

% This condition is triggered when the computation in xy is not done yet and hence full y-direction is not covered yet
if (new_ofmap_width > 0) && (ofmap_height_effective < ofmap_height)
    disp('y_effective triggers so F_effective = unitvol_nosof_3D_filter')
    % calling the function with the updated ifmap/ofmap width if the x-effective condition updateds them or with the original ifmap/ofmap width if not updated
    % Also calling the function with F_effective if the x-effective condition triggers so that F-direction of the ofmap volume is fixed ...
    % for regularity in the 3D volume. Otherwise calling it with the original Nos_of_filter
    Layer_param = [filter_height filter_width ifmap_height new_ifmap_width Nos_of_channel ofmap_height new_ofmap_width new_Nos_of_filter stride];
    [XYF_effective, cycle_count, SRAM_Access, DRAM_Access] = XtoYtoFtoZ(Layer_param, Hardware_param, Tech_param);
    
    ofmap_width_effective = XYF_effective(1);
    ofmap_height_effective = XYF_effective(2);
    ifmap_width_effective = XYF_effective(3);
    ifmap_height_effective = XYF_effective(4);
    %F_effective = XYF_effective(5); do not update F_effective to keep F_effective fixed as obtained before since we do not know whether with
    %the updated ofmap_width from x-condition, F-effective will remain same or not
    
    % Nos of same subvolume which is repeated along the y-direction
    pass_y = floor(ofmap_height/ofmap_height_effective);
    
    % Computation result of the subvolume
    cycle_count_ideal_xyFz_sub = cycle_count(1);
    DRAM_stall_xyFz_sub = cycle_count(2);
    cycle_to_write_ofmap_sub = cycle_count(3);
    
    Access_SRAM_filter_xyFz_sub = SRAM_Access(1); % in Nos of element, not bit
    Access_SRAM_ifmap_xyFz_sub = SRAM_Access(2);
    Access_SRAM_psum_xyFz_sub = SRAM_Access(3);

    Access_DRAM_filter_xyFz_sub = DRAM_Access(1);
    Access_DRAM_ifmap_xyFz_sub = DRAM_Access(2);
    Access_DRAM_ofmap_xyFz_sub = DRAM_Access(3);
    
    % The same subvolume repeats along the y-direction
    % Adding the new computation results with the previously computed results
    % Note that if the x-effective condition is not triggered then the first terms are set to zero beforehand
    cycle_count_ideal_xyFz = cycle_count_ideal_xyFz + cycle_count_ideal_xyFz_sub * pass_y;
    DRAM_stall_xyFz = DRAM_stall_xyFz + DRAM_stall_xyFz_sub * pass_y;
    cycle_to_write_ofmap = cycle_to_write_ofmap + cycle_to_write_ofmap_sub * pass_y;

    Access_SRAM_filter_xyFz = Access_SRAM_filter_xyFz + Access_SRAM_filter_xyFz_sub * pass_y; % in Nos of element, not bit
    Access_SRAM_ifmap_xyFz = Access_SRAM_ifmap_xyFz + Access_SRAM_ifmap_xyFz_sub * pass_y;
    Access_SRAM_psum_xyFz = Access_SRAM_psum_xyFz + Access_SRAM_psum_xyFz_sub * pass_y;

    Access_DRAM_filter_xyFz = Access_DRAM_filter_xyFz + Access_DRAM_filter_xyFz_sub * pass_y;
    Access_DRAM_ifmap_xyFz = Access_DRAM_ifmap_xyFz + Access_DRAM_ifmap_xyFz_sub * pass_y;
    Access_DRAM_ofmap_xyFz = Access_DRAM_ofmap_xyFz + Access_DRAM_ofmap_xyFz_sub * pass_y;

    % The remaining part along the y-direction which is not processed yet
    new_ofmap_height = ofmap_height - ofmap_height_effective * pass_y;
    new_ifmap_height = (new_ofmap_height - 1) * stride + filter_height;
    
    if (new_ofmap_height == 0) % computation is done for the full xy direction if there is no remaining part along the y-direction
        disp('computation calculation done for the full xy direction')
        %done_flag = 1;
    else
        % Computation for the last subvolume of the xy-direction.
        % This situation of last subvolume of xy comes only when y-effective triggers. Hence can write inside the y-effective condition
        disp('Computation for the last y directional subvolume')
        % calling the function with the updated ifmap/ofmap width if the x-effective condition updateds them or with the original ifmap/ofmap width if not updated
        % and with the updated height from the y-effective condition.
        % Also using the original F_effective to call the function
        Layer_param = [filter_height filter_width new_ifmap_height new_ifmap_width Nos_of_channel new_ofmap_height new_ofmap_width F_effective stride];
        [XYF_effective, cycle_count, SRAM_Access, DRAM_Access] = XtoYtoFtoZ(Layer_param, Hardware_param, Tech_param);
    
        cycle_count_ideal_xyFz = cycle_count_ideal_xyFz + cycle_count(1);
        DRAM_stall_xyFz = DRAM_stall_xyFz + cycle_count(2);
        cycle_to_write_ofmap = cycle_to_write_ofmap + cycle_count(3);

        Access_SRAM_filter_xyFz = Access_SRAM_filter_xyFz + SRAM_Access(1); % in Nos of element, not bit
        Access_SRAM_ifmap_xyFz = Access_SRAM_ifmap_xyFz + SRAM_Access(2);
        Access_SRAM_psum_xyFz = Access_SRAM_psum_xyFz + SRAM_Access(3);

        Access_DRAM_filter_xyFz = Access_DRAM_filter_xyFz + DRAM_Access(1);
        Access_DRAM_ifmap_xyFz = Access_DRAM_ifmap_xyFz + DRAM_Access(2);
        Access_DRAM_ofmap_xyFz = Access_DRAM_ofmap_xyFz + DRAM_Access(3);
    end
end

% Now covering the full F_direction
% if y-effective triggers, F-effective will trigger too and F_effective = unitvol_nosof_3D_filter
% Here, done_flag will never be zero inside x & y effective condition, besides F-effective is not updated inside x & y effective condition
if (F_effective < Nos_of_filter)
    disp('F-effective triggers')
    pass_F = floor(Nos_of_filter/F_effective);
    
    if (cycle_count_ideal_xyFz == 0) % This means neither x or y effective triggered, only F-effective triggered, Hence the variable values are zero
        % Since made the variable values zero, calling the function again with the original layer parameters to obtain the computation result for full xy dimension
        Layer_param = [filter_height filter_width ifmap_height ifmap_width Nos_of_channel ofmap_height ofmap_width Nos_of_filter stride];
        [XYF_effective, cycle_count, SRAM_Access, DRAM_Access] = XtoYtoFtoZ(Layer_param, Hardware_param, Tech_param);
        
        cycle_count_ideal_xyFz = cycle_count(1);
        DRAM_stall_xyFz = cycle_count(2);
        cycle_to_write_ofmap = cycle_count(3);

        Access_SRAM_filter_xyFz = SRAM_Access(1); % in Nos of element, not bit
        Access_SRAM_ifmap_xyFz = SRAM_Access(2);
        Access_SRAM_psum_xyFz = SRAM_Access(3);

        Access_DRAM_filter_xyFz = DRAM_Access(1);
        Access_DRAM_ifmap_xyFz = DRAM_Access(2);
        Access_DRAM_ofmap_xyFz = DRAM_Access(3);    
    end
    
    % The full xy direction is done and the same computation repeats along the F-direction. This full xy computation results comes either from the above if
    % condition when only F-effective triggers or from the above y-effective condition when x or y or both triggers
    cycle_count_ideal_xyFz = cycle_count_ideal_xyFz * pass_F;
    DRAM_stall_xyFz = DRAM_stall_xyFz * pass_F;
    cycle_to_write_ofmap = cycle_to_write_ofmap * pass_F;
        
    Access_SRAM_filter_xyFz = Access_SRAM_filter_xyFz * pass_F; % in Nos of element, not bit
    Access_SRAM_ifmap_xyFz = Access_SRAM_ifmap_xyFz * pass_F;
    Access_SRAM_psum_xyFz = Access_SRAM_psum_xyFz * pass_F;

    Access_DRAM_filter_xyFz = Access_DRAM_filter_xyFz * pass_F;
    Access_DRAM_ifmap_xyFz = Access_DRAM_ifmap_xyFz * pass_F;
    Access_DRAM_ofmap_xyFz = Access_DRAM_ofmap_xyFz * pass_F;
   
    % The remaining part along the F-direction which is not processed yet
    new_Nos_of_filter = Nos_of_filter - F_effective * pass_F;
    
    if (new_Nos_of_filter == 0) % computation is done for the full xyF direction if there is no remaining part along the F-direction
        disp('computation calculation done for the full xyF direction')
        done_flag = 1;
    else
        % Computation for the remaining part in the F-direction
        % This situation of last F directional subvolume comes only when F-effective triggers. Hence can write inside the F-effective condition
        disp('Computation for the last F directional subvolume')
        % F-effective will not trigger any more. Now we have to go through the same process where only x and y is candidate effective parameters
        % Calling the function with updated layer parameter to do that
        Layer_param = [filter_height filter_width ifmap_height ifmap_width Nos_of_channel ofmap_height ofmap_width new_Nos_of_filter stride];
        [cycle_count_sub, SRAM_Access_sub, DRAM_Access_sub] = Asym_XYeffect_XYFZ(Layer_param, Hardware_param, Tech_param);   
    
        cycle_count_ideal_xyFz = cycle_count_ideal_xyFz + cycle_count_sub(1);
        DRAM_stall_xyFz = DRAM_stall_xyFz + cycle_count_sub(2);
        cycle_to_write_ofmap = cycle_to_write_ofmap + cycle_count_sub(3);

        Access_SRAM_filter_xyFz = Access_SRAM_filter_xyFz + SRAM_Access_sub(1); % in Nos of element, not bit
        Access_SRAM_ifmap_xyFz = Access_SRAM_ifmap_xyFz + SRAM_Access_sub(2);
        Access_SRAM_psum_xyFz = Access_SRAM_psum_xyFz + SRAM_Access_sub(3);

        Access_DRAM_filter_xyFz = Access_DRAM_filter_xyFz + DRAM_Access_sub(1);
        Access_DRAM_ifmap_xyFz = Access_DRAM_ifmap_xyFz + DRAM_Access_sub(2);
        Access_DRAM_ofmap_xyFz = Access_DRAM_ofmap_xyFz + DRAM_Access_sub(3);
    
    end 
end



% cycle_count_ideal_xyFz       
% DRAM_stall_xyFz  
% cycle_to_write_ofmap
% 
% Access_SRAM_filter_xyFz
% Access_SRAM_ifmap_xyFz
% Access_SRAM_psum_xyFz
% 
% Access_DRAM_filter_xyFz
% Access_DRAM_ifmap_xyFz
% Access_DRAM_ofmap_xyFz

cycle_count_layer = [cycle_count_ideal_xyFz DRAM_stall_xyFz cycle_to_write_ofmap];
SRAM_Access_layer = [Access_SRAM_filter_xyFz Access_SRAM_ifmap_xyFz Access_SRAM_psum_xyFz];  % in Nos of element, not bit
DRAM_Access_layer = [Access_DRAM_filter_xyFz Access_DRAM_ifmap_xyFz Access_DRAM_ofmap_xyFz]; % in Nos of element, not bit


end
