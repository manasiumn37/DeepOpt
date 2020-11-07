function [Combined_Table, Table_LOS] = LOS_Computation (Hardware_param, Tech_param, Pmat_flag, Network_flag)

%This file generates the Performance Matric (Pmat) of a network while using a fixed branch (FS) vs. layer-specific optimal branches (LOS) to compute the full network.

%% Obtaining the the Network Parameters
[Layer, filter_size_Net, Nos_of_Filter_Net, Ifmap_size_Net, Nos_of_Channel_Net, Stride_Net, Ofmap_size_Net, Batch_size_Net, fc_flag]...
                                                                                                                = Network_Parameters (Network_flag);
                                                                                                            
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

%% Technology Parameters
Energy_Add_perbit = Tech_param(1);      % per bit
Energy_Mul_element = Tech_param(2);     % for each bw_filter/bw_ifmap bit element
E_RF_to_ALU = Tech_param(3);            % per bit
E_fsSRAM_to_RF = Tech_param(4);         % per bit
E_isSRAM_to_RF = Tech_param(5);         % per bit
E_psSRAM_to_RF = Tech_param(6);         % per bit
E_DRAM_to_SRAM = Tech_param(7);         % per bit


%% Performing the computation with layer-specific best branch
for i = 1:1:length(Layer)
    % Specification of a single layer
    filter_height = filter_size_Net(i);
    filter_width = filter_size_Net(i);
    ifmap_height = Ifmap_size_Net(i);
    ifmap_width = Ifmap_size_Net(i);
    Nos_of_channel = Nos_of_Channel_Net(i);
    ofmap_height = Ofmap_size_Net(i);
    ofmap_width = Ofmap_size_Net(i);
    Nos_of_filter = Nos_of_Filter_Net(i);
    stride = Stride_Net(i);
    batch_size = Batch_size_Net(i);    

    % layer parameters together to pass
    Layer_param = [filter_height filter_width ifmap_height ifmap_width Nos_of_channel ofmap_height ofmap_width Nos_of_filter stride batch_size];
    
    if (fc_flag(i) == 0)  % Conv Layers
        disp('Calling the function to compute using X->Y->Z->F branch') %SRAM/DRAM access are in Nos of element, not bit
        [cycle_count_xyzF, SRAM_Access_xyzF, DRAM_Access_xyzF] = Asymmetric_XYZF(Layer_param, Hardware_param, Tech_param);
        Total_cycle_xyzF = sum(cycle_count_xyzF);
        Total_weighted_access_xyzF = SRAM_Access_xyzF(1) * bw_filter * E_fsSRAM_to_RF + SRAM_Access_xyzF(2) * bw_ifmap * E_isSRAM_to_RF...
                                   + SRAM_Access_xyzF(3) * bw_psum * E_psSRAM_to_RF...
                                   + (DRAM_Access_xyzF(1) * bw_filter + DRAM_Access_xyzF(2) * bw_ifmap + DRAM_Access_xyzF(3) * bw_ofmap) * E_DRAM_to_SRAM; %in Joule
        
        disp('Calling the function to compute using X->Z->Y->F branch')
        [cycle_count_xzyF, SRAM_Access_xzyF, DRAM_Access_xzyF] = Asymmetric_XZYF(Layer_param, Hardware_param, Tech_param);
        Total_cycle_xzyF = sum(cycle_count_xzyF);
        Total_weighted_access_xzyF = SRAM_Access_xzyF(1) * bw_filter * E_fsSRAM_to_RF + SRAM_Access_xzyF(2) * bw_ifmap * E_isSRAM_to_RF...
                                   + SRAM_Access_xzyF(3) * bw_psum * E_psSRAM_to_RF...
                                   + (DRAM_Access_xzyF(1) * bw_filter + DRAM_Access_xzyF(2) * bw_ifmap + DRAM_Access_xzyF(3) * bw_ofmap) * E_DRAM_to_SRAM; %in Joule
        
        disp('Calling the function to compute using X->Y->F->Z branch')
        [cycle_count_xyFz, SRAM_Access_xyFz, DRAM_Access_xyFz] = Asymmetric_XYFZ(Layer_param, Hardware_param, Tech_param);
        Total_cycle_xyFz = sum(cycle_count_xyFz);
        Total_weighted_access_xyFz = SRAM_Access_xyFz(1) * bw_filter *  E_fsSRAM_to_RF + SRAM_Access_xyFz(2) * bw_ifmap *  E_isSRAM_to_RF...
                                   + SRAM_Access_xyFz(3) * bw_psum *  E_psSRAM_to_RF...
                                   + (DRAM_Access_xyFz(1) * bw_filter + DRAM_Access_xyFz(2) * bw_ifmap + DRAM_Access_xyFz(3) * bw_ofmap) * E_DRAM_to_SRAM; %in Joule
        
        disp('Calling the function to compute using X->F->Y->Z branch')
        [cycle_count_xFyz, SRAM_Access_xFyz, DRAM_Access_xFyz] = Asymmetric_XFYZ(Layer_param, Hardware_param, Tech_param);
        Total_cycle_xFyz = sum(cycle_count_xFyz);
        Total_weighted_access_xFyz = SRAM_Access_xFyz(1) * bw_filter * E_fsSRAM_to_RF + SRAM_Access_xFyz(2) * bw_ifmap * E_isSRAM_to_RF...
                                   + SRAM_Access_xFyz(3) * bw_psum * E_psSRAM_to_RF...
                                   + (DRAM_Access_xFyz(1) * bw_filter + DRAM_Access_xFyz(2) * bw_ifmap + DRAM_Access_xFyz(3) * bw_ofmap) * E_DRAM_to_SRAM; %in Joule

        disp('Calling the function to compute using Z->X->Y->F branch')
        [cycle_count_zxyF, SRAM_Access_zxyF, DRAM_Access_zxyF] = ZtoXtoYtoF(Layer_param, Hardware_param, Tech_param);
        Total_cycle_zxyF = sum(cycle_count_zxyF);
        Total_weighted_access_zxyF = SRAM_Access_zxyF(1) * bw_filter * E_fsSRAM_to_RF + SRAM_Access_zxyF(2) * bw_ifmap * E_isSRAM_to_RF...
                                   + SRAM_Access_zxyF(3) * bw_psum * E_psSRAM_to_RF...
                                   + (DRAM_Access_zxyF(1) * bw_filter + DRAM_Access_zxyF(2) * bw_ifmap + DRAM_Access_zxyF(3) * bw_ofmap) * E_DRAM_to_SRAM; %in Joule
                                 
        %%%%%% Determining best branch in terms of the performance metric (i.e.,E, D, ED, E^2*D, and E*D^2)
        branch_name = ["xyzF", "xzyF", "xyFz", "xFyz", "zxyF"];
        Cycle_all_branch = [Total_cycle_xyzF Total_cycle_xzyF Total_cycle_xyFz Total_cycle_xFyz Total_cycle_zxyF];
        WeightedAccess_all_branch = [Total_weighted_access_xyzF Total_weighted_access_xzyF Total_weighted_access_xyFz...
                                                                        Total_weighted_access_xFyz Total_weighted_access_zxyF];
        
        % MAC computation energy & Register file access energy for a layer, this is fixed for all branches
        Nos_of_MAC_layer = filter_height * filter_width * Nos_of_channel * ofmap_height * ofmap_width * Nos_of_filter;
        MAC_Energy_layer = Nos_of_MAC_layer * (Energy_Add_perbit * bw_psum + Energy_Mul_element); % in Joule
        % for each MAC, one filter element read, one ifmap element read, one psum element read, and one psum element write occur from/to the RF level
        RF_Energy_layer = (E_RF_to_ALU * bw_filter + E_RF_to_ALU * bw_ifmap + 2 * E_RF_to_ALU * bw_psum) * Nos_of_MAC_layer ; % in Joule, 
        
        % Total energy including MAC and RF energy with the SRAM and DRAM access energy
        WeightedAccess_all_branch = WeightedAccess_all_branch + (MAC_Energy_layer + RF_Energy_layer); % in Joule
        
        if (Network_flag == 1) % Applicable only for AlexNet due to the grouping in C2, C4, C5
            if (i == 2) || (i == 4) || (i == 5)   % multiplying the result of C2, C4, C5 with 2
                Cycle_all_branch = [Total_cycle_xyzF Total_cycle_xzyF Total_cycle_xyFz Total_cycle_xFyz Total_cycle_zxyF].*2;
                WeightedAccess_all_branch = WeightedAccess_all_branch .* 2;
            end
        end
        
        %Chaning performance metric based on Pmat_flag
        if (Pmat_flag == 1)
            Pmat_all_branch = WeightedAccess_all_branch.*Cycle_all_branch;                                  % when the performance metric is Energy * Delay
        elseif (Pmat_flag == 2)
            Pmat_all_branch = (WeightedAccess_all_branch.*WeightedAccess_all_branch).*Cycle_all_branch;     % when the performance metric is Energy^2 * Delay
        elseif (Pmat_flag == 3)
            Pmat_all_branch = WeightedAccess_all_branch.*Cycle_all_branch.*Cycle_all_branch;                % when the performance metric is Energy * Delay^2
        elseif (Pmat_flag == 4)
            Pmat_all_branch = WeightedAccess_all_branch;                                                    % when the performance metric is Energy
        elseif (Pmat_flag == 5)
            Pmat_all_branch = Cycle_all_branch;                                                             % when the performance metric is Delay
        end
        
        [Sorted_Pmat_all_branch(i,:), sort_index_Pm]= sort(Pmat_all_branch);
        Sorted_branch_Pm(i,:) = branch_name(sort_index_Pm);
                
        % Ordered result according to branch_name (i.e., in "xyzF", "xzyF", "xyFz", "xFyz", "zxyF" order), not sorted.
        % This is to compute the fixed branch performance
        Ordered_Pmat_all_branch(i,:) = Pmat_all_branch;

                                     
    elseif (fc_flag(i) == 1) % FC Layers
        if (batch_size == 1) % without batching
            disp('FC layer without batching')
            
            disp('Calling the function to compute using Z->F branch')
            [cycle_count_zF, SRAM_Access_zF, DRAM_Access_zF] = FC_ZtoF(Layer_param, Hardware_param, Tech_param);
            Total_cycle_zF = sum(cycle_count_zF);
            Total_weighted_access_zF = SRAM_Access_zF(1) * bw_filter * E_fsSRAM_to_RF + SRAM_Access_zF(2) * bw_ifmap * E_isSRAM_to_RF...
                                     + SRAM_Access_zF(3) * bw_psum * E_psSRAM_to_RF...
                                     + (DRAM_Access_zF(1) * bw_filter + DRAM_Access_zF(2) * bw_ifmap + DRAM_Access_zF(3) * bw_ofmap) * E_DRAM_to_SRAM; %in Joule

            disp('Calling the function to compute using F->Z branch')
            [cycle_count_Fz, SRAM_Access_Fz, DRAM_Access_Fz] = Asymmetric_FZ(Layer_param, Hardware_param, Tech_param);
            Total_cycle_Fz = sum(cycle_count_Fz);
            Total_weighted_access_Fz = SRAM_Access_Fz(1) * bw_filter * E_fsSRAM_to_RF + SRAM_Access_Fz(2) * bw_ifmap * E_isSRAM_to_RF...
                                     + SRAM_Access_Fz(3) * bw_psum * E_psSRAM_to_RF...
                                     + (DRAM_Access_Fz(1) * bw_filter + DRAM_Access_Fz(2) * bw_ifmap + DRAM_Access_Fz(3) * bw_ofmap) * E_DRAM_to_SRAM; %in Joule
                                   
            %%%%% Determining best branch in terms of the performance metric (i.e.,E, D, ED, E^2*D, and E*D^2)
            %(putting inf & N/A for the last 3 locations to match the metric from conv layer)                    
            branch_name = ["zF", "Fz", "N/A", "N/A", "N/A"];
            Cycle_all_branch = [Total_cycle_zF Total_cycle_Fz inf inf inf];
            WeightedAccess_all_branch = [Total_weighted_access_zF Total_weighted_access_Fz inf inf inf];
            
            % MAC computation energy & Register file access energy for a layer, this is fixed for all branches
            Nos_of_MAC_layer = filter_height * filter_width * Nos_of_channel * ofmap_height * ofmap_width * Nos_of_filter;
            MAC_Energy_layer = Nos_of_MAC_layer * (Energy_Add_perbit * bw_psum + Energy_Mul_element); % in Joule
            % for each MAC, one filter element read, one ifmap element read, one psum element read, and one psum element write occur from/to the RF level
            RF_Energy_layer = (E_RF_to_ALU * bw_filter + E_RF_to_ALU * bw_ifmap + 2 * E_RF_to_ALU * bw_psum) * Nos_of_MAC_layer ; % in Joule, 

            % Total energy including MAC and RF energy with the SRAM and DRAM access energy
            WeightedAccess_all_branch = WeightedAccess_all_branch + (MAC_Energy_layer + RF_Energy_layer); % in Joule
            
            %Chaning performance metric based on Pmat_flag
            if (Pmat_flag == 1)
                Pmat_all_branch = WeightedAccess_all_branch.*Cycle_all_branch;                                  % when the performance metric is Energy * Delay
            elseif (Pmat_flag == 2)
                Pmat_all_branch = (WeightedAccess_all_branch.*WeightedAccess_all_branch).*Cycle_all_branch;     % when the performance metric is Energy^2 * Delay
            elseif (Pmat_flag == 3)
                Pmat_all_branch = WeightedAccess_all_branch.*Cycle_all_branch.*Cycle_all_branch;                % when the performance metric is Energy * Delay^2
            elseif (Pmat_flag == 4)
                Pmat_all_branch = WeightedAccess_all_branch;                                                    % when the performance metric is Energy
            elseif (Pmat_flag == 5)
                Pmat_all_branch = Cycle_all_branch;                                                             % when the performance metric is Delay
            end
        
            [Sorted_Pmat_all_branch(i,:), sort_index_Pm]= sort(Pmat_all_branch);
            Sorted_branch_Pm(i,:) = branch_name(sort_index_Pm);
                        
            % Ordered result according to branch_name (i.e., in "xyzF", "xzyF", "xyFz", "xFyz", "zxyF" order), not sorted.
            % This is to compute the fixed branch performance
            % For FC layer putting the result of zF = xyzF, xzyF, zxyF and Fz = xyFz, xFyz to preserve the branch order
            Ordered_cycle_branches = [Total_cycle_zF Total_cycle_zF Total_cycle_Fz Total_cycle_Fz Total_cycle_zF];
            Ordered_WAccess_branches = [Total_weighted_access_zF Total_weighted_access_zF Total_weighted_access_Fz ...
                                                                        Total_weighted_access_Fz Total_weighted_access_zF];
            
            % Total ordered energy including MAC and RF energy with the SRAM and DRAM access energy                                                       
            Ordered_WAccess_branches = Ordered_WAccess_branches + (MAC_Energy_layer + RF_Energy_layer); % in Joule
            
            %Chaning performance metric based on Pmat_flag  
            if (Pmat_flag == 1)
                Ordered_Pmat_all_branch(i,:) = Ordered_WAccess_branches.*Ordered_cycle_branches;                             %performance metric is Energy * Delay
            elseif (Pmat_flag == 2)
                Ordered_Pmat_all_branch(i,:) = Ordered_WAccess_branches.*Ordered_WAccess_branches.*Ordered_cycle_branches;   %performance metric is Energy^2 * Delay
            elseif (Pmat_flag == 3)
                Ordered_Pmat_all_branch(i,:) = Ordered_WAccess_branches.*Ordered_cycle_branches.*Ordered_cycle_branches;     %performance metric is Energy * Delay^2
            elseif (Pmat_flag == 4)
                Ordered_Pmat_all_branch(i,:) = Ordered_WAccess_branches;                                                     %performance metric is Energy
            elseif (Pmat_flag == 5)
                Ordered_Pmat_all_branch(i,:) = Ordered_cycle_branches;                                                       %performance metric is Delay
            end
            
                                                                               
        elseif (batch_size > 1) % with batching
            disp('FC layer with batching')
            % FC layer with batching is equvalent to 1*1 conv layer with X = N, Y = 1
            ifmap_width = batch_size;
            ofmap_width = batch_size;
            Layer_param = [filter_height filter_width ifmap_height ifmap_width Nos_of_channel ofmap_height ofmap_width Nos_of_filter stride batch_size];
            
            disp('Calling the function to compute using N->Z->F branch which is equivalent to calling X->Y->Z->F branch')
            % Using X->Z->Y->F to do the computation of N->Z->F gives the same access energy. However the DRAM stall count is slightly heigher. It happens because
            % in this branch, the ifmap reverse order formatting is omitted during stall count for simplicity. Hence, better to use X->Y->Z->F.
            [cycle_count_NzF, SRAM_Access_NzF, DRAM_Access_NzF] = Asymmetric_XYZF(Layer_param, Hardware_param, Tech_param);
            cycle_count_NzF_1N = ceil(cycle_count_NzF./batch_size); % per image result
            SRAM_Access_NzF_1N = SRAM_Access_NzF./batch_size;
            DRAM_Access_NzF_1N = DRAM_Access_NzF./batch_size;
            
            Total_cycle_NzF_1N = sum(cycle_count_NzF_1N);
            Total_weighted_access_NzF_1N = SRAM_Access_NzF_1N(1) * bw_filter * E_fsSRAM_to_RF + SRAM_Access_NzF_1N(2) * bw_ifmap * E_isSRAM_to_RF...
                                         + SRAM_Access_NzF_1N(3) * bw_psum * E_psSRAM_to_RF...
                                         + (DRAM_Access_NzF_1N(1) * bw_filter + DRAM_Access_NzF_1N(2) * bw_ifmap + DRAM_Access_NzF_1N(3) * bw_ofmap) *  E_DRAM_to_SRAM; %in Joule

            disp('Calling the function to compute using N->F->Z branch which is equivalent to calling X->Y->F->Z branch')
            % Can also use X->F->Y->Z
            [cycle_count_NFz, SRAM_Access_NFz, DRAM_Access_NFz] = Asymmetric_XYFZ(Layer_param, Hardware_param, Tech_param);
            cycle_count_NFz_1N = ceil(cycle_count_NFz./batch_size); % per image result
            SRAM_Access_NFz_1N = SRAM_Access_NFz./batch_size;
            DRAM_Access_NFz_1N = DRAM_Access_NFz./batch_size;    
            
            Total_cycle_NFz_1N = sum(cycle_count_NFz_1N);
            Total_weighted_access_NFz_1N = SRAM_Access_NFz_1N(1) * bw_filter * E_fsSRAM_to_RF + SRAM_Access_NFz_1N(2) * bw_ifmap * E_isSRAM_to_RF...
                                         + SRAM_Access_NFz_1N(3) * bw_psum * E_psSRAM_to_RF...
                                         + (DRAM_Access_NFz_1N(1) * bw_filter + DRAM_Access_NFz_1N(2) * bw_ifmap + DRAM_Access_NFz_1N(3) * bw_ofmap) *  E_DRAM_to_SRAM; %in Joule

            %Determining best branch in terms of energy & performance (putting inf & N/A for the last 3 locations to match the metric from conv layer)                      
            branch_name = ["NzF-1N", "NFz-1N", "N/A", "N/A", "N/A"];
            Cycle_all_branch = [Total_cycle_NzF_1N Total_cycle_NFz_1N inf inf inf];
            WeightedAccess_all_branch = [Total_weighted_access_NzF_1N Total_weighted_access_NFz_1N inf inf inf];
            
            % MAC computation energy & Register file access energy for a layer, this is fixed for all branches
            Nos_of_MAC_layer = filter_height * filter_width * Nos_of_channel * ofmap_height * ofmap_width * Nos_of_filter;
            MAC_Energy_layer = Nos_of_MAC_layer * (Energy_Add_perbit * bw_psum + Energy_Mul_element); % in Joule
            % for each MAC, one filter element read, one ifmap element read, one psum element read, and one psum element write occur from/to the RF level
            RF_Energy_layer = (E_RF_to_ALU * bw_filter + E_RF_to_ALU * bw_ifmap + 2 * E_RF_to_ALU * bw_psum) * Nos_of_MAC_layer ; % in Joule, 

            % Total energy including MAC and RF energy with the SRAM and DRAM access energy
            WeightedAccess_all_branch = WeightedAccess_all_branch + (MAC_Energy_layer + RF_Energy_layer); % in Joule
            
            %Chaning performance metric based on Pmat_flag
            if (Pmat_flag == 1)
                Pmat_all_branch = WeightedAccess_all_branch.*Cycle_all_branch;                                  % when the performance metric is Energy * Delay
            elseif (Pmat_flag == 2)
                Pmat_all_branch = (WeightedAccess_all_branch.*WeightedAccess_all_branch).*Cycle_all_branch;     % when the performance metric is Energy^2 * Delay
            elseif (Pmat_flag == 3)
                Pmat_all_branch = WeightedAccess_all_branch.*Cycle_all_branch.*Cycle_all_branch;                % when the performance metric is Energy * Delay^2
            elseif (Pmat_flag == 4)
                Pmat_all_branch = WeightedAccess_all_branch;                                                    % when the performance metric is Energy
            elseif (Pmat_flag == 5)
                Pmat_all_branch = Cycle_all_branch;                                                             % when the performance metric is Delay
            end
        
            [Sorted_Pmat_all_branch(i,:), sort_index_Pm]= sort(Pmat_all_branch);
            Sorted_branch_Pm(i,:) = branch_name(sort_index_Pm);
                         
            % Ordered result according to branch_name (i.e., in "xyzF", "xzyF", "xyFz", "xFyz", "zxyF" order), not sorted.
            % This is to compute the fixed branch performance
            % For FC layer putting the result of NzF-1N = xyzF, xzyF, zxyF and NFz-1N = xyFz, xFyz to preserve the branch order
            Ordered_cycle_branches = [Total_cycle_NzF_1N Total_cycle_NzF_1N Total_cycle_NFz_1N Total_cycle_NFz_1N Total_cycle_NzF_1N];
            Ordered_WAccess_branches = [Total_weighted_access_NzF_1N Total_weighted_access_NzF_1N Total_weighted_access_NFz_1N ...
                                                                        Total_weighted_access_NFz_1N Total_weighted_access_NzF_1N];
                                                                    
            % Total ordered energy including MAC and RF energy with the SRAM and DRAM access energy                                                       
            Ordered_WAccess_branches = Ordered_WAccess_branches + (MAC_Energy_layer + RF_Energy_layer); % in Joule
                                                                    
            %Chaning performance metric based on Pmat_flag
            if (Pmat_flag == 1)
                Ordered_Pmat_all_branch(i,:) = Ordered_WAccess_branches.*Ordered_cycle_branches;                             %performance metric is Energy * Delay
            elseif (Pmat_flag == 2)
                Ordered_Pmat_all_branch(i,:) = Ordered_WAccess_branches.*Ordered_WAccess_branches.*Ordered_cycle_branches;   %performance metric is Energy^2 * Delay
            elseif (Pmat_flag == 3)
                Ordered_Pmat_all_branch(i,:) = Ordered_WAccess_branches.*Ordered_cycle_branches.*Ordered_cycle_branches;     %performance metric is Energy * Delay^2
            elseif (Pmat_flag == 4)
                Ordered_Pmat_all_branch(i,:) = Ordered_WAccess_branches;                                                     %performance metric is Energy
            elseif (Pmat_flag == 5)
                Ordered_Pmat_all_branch(i,:) = Ordered_cycle_branches;                                                       %performance metric is Delay
            end
  
            
        end        
    end    
end

Sorted_Pmat_all_branch;
Sorted_branch_Pm;
Ordered_Pmat_all_branch;

%% Result using layer-specific best branching
Best_Pm_branch_tbl = "Hybrid"; 
Best_Pmat_all_layer = Sorted_Pmat_all_branch(:,1);
Best_Pm_tbl = sum(Best_Pmat_all_layer);

% Layer specific optimal branches
Best_branch_all_layer = Sorted_branch_Pm(:,1);
Layer_No = [Layer; "Total"];
Best_branch = [Best_branch_all_layer; "Hybrid"]; % The last one is for the "Total" row
Best_Pmat = [Best_Pmat_all_layer; Best_Pm_tbl];

Table_LOS = table(Layer_No, Best_branch, Best_Pmat);  % Table showing the optimal branch to compute each layer under LOS


%% Result using fixed xyzF branching
Pmat_all_layer_xyzF = Ordered_Pmat_all_branch(:,1);
branch_xyzF_tbl = "xyzF";   
Pmat_xyzF_tbl = sum(Pmat_all_layer_xyzF);

% penalty as compared to the best case
Percent_penalty_xyzF_tbl = ((Pmat_xyzF_tbl - Best_Pm_tbl)/Best_Pm_tbl)*100;

% variable names to create the table
bit_width_filter = bw_filter;
bw_ifmap_ofmap = bw_ifmap;
bit_width_psum = bw_psum;

PE_Row = Array_row;
PE_Column = Array_column;

SRAM_filter_KB = SRAM_filter/(1024 * 8);  % size in KB
SRAM_ifmap_KB = SRAM_ifmap/(1024 * 8);    % size in KB
SRAM_psum_KB = SRAM_psum/(1024 * 8);      % size in KB

Bandwidth_DRAM = BW_DRAM;                 % in bit/cycle

Fixed_branch = branch_xyzF_tbl;
Fixed_Pmat = Pmat_xyzF_tbl;
Best_branch = Best_Pm_branch_tbl;
Best_Pmat = Best_Pm_tbl;
Percent_penalty = Percent_penalty_xyzF_tbl;

Table_xyzF = table(bit_width_filter, bw_ifmap_ofmap, bit_width_psum, PE_Row, PE_Column, SRAM_filter_KB, SRAM_ifmap_KB, SRAM_psum_KB, Bandwidth_DRAM,...
                    Fixed_branch, Fixed_Pmat, Best_branch, Best_Pmat, Percent_penalty);


%% Result using fixed xzyF branching
Pmat_all_layer_xzyF = Ordered_Pmat_all_branch(:,2);
branch_xzyF_tbl = "xzyF";   
Pmat_xzyF_tbl = sum(Pmat_all_layer_xzyF);

% penalty as compared to the best case
Percent_penalty_xzyF_tbl = ((Pmat_xzyF_tbl - Best_Pm_tbl)/Best_Pm_tbl)*100;

% variable names to create the table
bit_width_filter = bw_filter;
bw_ifmap_ofmap = bw_ifmap;
bit_width_psum = bw_psum;

PE_Row = Array_row;
PE_Column = Array_column;

SRAM_filter_KB = SRAM_filter/(1024 * 8);  % size in KB
SRAM_ifmap_KB = SRAM_ifmap/(1024 * 8);    % size in KB
SRAM_psum_KB = SRAM_psum/(1024 * 8);      % size in KB

Bandwidth_DRAM = BW_DRAM;                 % in bit/cycle

Fixed_branch = branch_xzyF_tbl;
Fixed_Pmat = Pmat_xzyF_tbl;
Best_branch = Best_Pm_branch_tbl;
Best_Pmat = Best_Pm_tbl;
Percent_penalty = Percent_penalty_xzyF_tbl;

Table_xzyF = table(bit_width_filter, bw_ifmap_ofmap, bit_width_psum, PE_Row, PE_Column, SRAM_filter_KB, SRAM_ifmap_KB, SRAM_psum_KB, Bandwidth_DRAM,...
                    Fixed_branch, Fixed_Pmat, Best_branch, Best_Pmat, Percent_penalty);
                

%% Result using fixed xyFz branching
Pmat_all_layer_xyFz = Ordered_Pmat_all_branch(:,3);
branch_xyFz_tbl = "xyFz";
Pmat_xyFz_tbl = sum(Pmat_all_layer_xyFz);

% penalty as compared to the best case
Percent_penalty_xyFz_tbl = ((Pmat_xyFz_tbl - Best_Pm_tbl)/Best_Pm_tbl)*100;

% variable names to create the table
bit_width_filter = bw_filter;
bw_ifmap_ofmap = bw_ifmap;
bit_width_psum = bw_psum;

PE_Row = Array_row;
PE_Column = Array_column;

SRAM_filter_KB = SRAM_filter/(1024 * 8);  % size in KB
SRAM_ifmap_KB = SRAM_ifmap/(1024 * 8);    % size in KB
SRAM_psum_KB = SRAM_psum/(1024 * 8);      % size in KB

Bandwidth_DRAM = BW_DRAM;                 % in bit/cycle

Fixed_branch = branch_xyFz_tbl;
Fixed_Pmat = Pmat_xyFz_tbl;
Best_branch = Best_Pm_branch_tbl;
Best_Pmat = Best_Pm_tbl;
Percent_penalty = Percent_penalty_xyFz_tbl;

Table_xyFz = table(bit_width_filter, bw_ifmap_ofmap, bit_width_psum, PE_Row, PE_Column, SRAM_filter_KB, SRAM_ifmap_KB, SRAM_psum_KB, Bandwidth_DRAM,...
                    Fixed_branch, Fixed_Pmat, Best_branch, Best_Pmat, Percent_penalty);


%% Result using fixed xFyz branching
Pmat_all_layer_xFyz = Ordered_Pmat_all_branch(:,4);
branch_xFyz_tbl = "xFyz";   
Pmat_xFyz_tbl = sum(Pmat_all_layer_xFyz);

% penalty as compared to the best case
Percent_penalty_xFyz_tbl = ((Pmat_xFyz_tbl - Best_Pm_tbl)/Best_Pm_tbl)*100;

% variable names to create the table
bit_width_filter = bw_filter;
bw_ifmap_ofmap = bw_ifmap;
bit_width_psum = bw_psum;

PE_Row = Array_row;
PE_Column = Array_column;

SRAM_filter_KB = SRAM_filter/(1024 * 8);  % size in KB
SRAM_ifmap_KB = SRAM_ifmap/(1024 * 8);    % size in KB
SRAM_psum_KB = SRAM_psum/(1024 * 8);      % size in KB

Bandwidth_DRAM = BW_DRAM;                 % in bit/cycle

Fixed_branch = branch_xFyz_tbl;
Fixed_Pmat = Pmat_xFyz_tbl;
Best_branch = Best_Pm_branch_tbl;
Best_Pmat = Best_Pm_tbl;
Percent_penalty = Percent_penalty_xFyz_tbl;

Table_xFyz = table(bit_width_filter, bw_ifmap_ofmap, bit_width_psum, PE_Row, PE_Column, SRAM_filter_KB, SRAM_ifmap_KB, SRAM_psum_KB, Bandwidth_DRAM,...
                    Fixed_branch, Fixed_Pmat, Best_branch, Best_Pmat, Percent_penalty);


%% Result using fixed zxyF branching
Pmat_all_layer_zxyF = Ordered_Pmat_all_branch(:,5);
branch_zxyF_tbl = "zxyF";   
Pmat_zxyF_tbl = sum(Pmat_all_layer_zxyF);

% penalty as compared to the best case
Percent_penalty_zxyF_tbl = ((Pmat_zxyF_tbl - Best_Pm_tbl)/Best_Pm_tbl)*100;

% variable names to create the table
bit_width_filter = bw_filter;
bw_ifmap_ofmap = bw_ifmap;
bit_width_psum = bw_psum;

PE_Row = Array_row;
PE_Column = Array_column;

SRAM_filter_KB = SRAM_filter/(1024 * 8);  % size in KB
SRAM_ifmap_KB = SRAM_ifmap/(1024 * 8);    % size in KB
SRAM_psum_KB = SRAM_psum/(1024 * 8);      % size in KB

Bandwidth_DRAM = BW_DRAM;                 % in bit/cycle

Fixed_branch = branch_zxyF_tbl;
Fixed_Pmat = Pmat_zxyF_tbl;
Best_branch = Best_Pm_branch_tbl;
Best_Pmat = Best_Pm_tbl;
Percent_penalty = Percent_penalty_zxyF_tbl;

Table_zxyF = table(bit_width_filter, bw_ifmap_ofmap, bit_width_psum, PE_Row, PE_Column, SRAM_filter_KB, SRAM_ifmap_KB, SRAM_psum_KB, Bandwidth_DRAM,...
                    Fixed_branch, Fixed_Pmat, Best_branch, Best_Pmat, Percent_penalty);

Combined_Table = [Table_xyzF; Table_xzyF; Table_xyFz; Table_xFyz; Table_zxyF];

end


