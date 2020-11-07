clc;
clear all;
close all;

%% Selecting CNN topology
% 1: AlexNet
% 2: VGG-16
% 3: SqueezeNet-v1.1
% 4: GoogleNet-v1
% 5: ResNet-50
Network_flag = 1;  % Flag to choose the network topology; can take values 1, 2, 3, 4, 5


%% Accelerator Specification
% bit width
bw_filter = 8;
bw_ifmap = bw_filter; %using same bitwidth for filter & ifmap since the current parameter for mul-energy does not support different widths for the two operands
bw_psum = 32;
bw_ofmap = bw_ifmap;

% MAC array
Array_row = 32;
Array_column = 64;

% On-chip SRAM storage
% Use any size between 4 to 512 kB with a step of 1kB. Currently, CACTI data is not available in the repository outside this range. 
% User have to extract data from CACTI if they wish to use a SRAM size outside the above mentioned range 
SRAM_filter = 16 * 1024 * 8; % in bit 
SRAM_ifmap = 32 * 1024 * 8; % in bit 
SRAM_psum = 8 * 1024 * 8; % in bit 

% To incorporation DRAM bandwidth induced stall in the calculation of cycle count
BW_DRAM = 128; %in bit/cycle
DRAM_block_size = 128; %in bit


%% Access energy specification 
s = (65/45)*(1/0.9)^2;                      %Technology scaling parameter; 65nm Technology node, Vdd = 1V; 45nm Technology node, Vdd = 0.9V

% to obtain 16-bit add and mul, the direct 8-bit add and mul data from energy computing paper are linearly and quadratically sacled, respectively
Energy_MAC16 = (0.03 * 2 + 0.2 * 4)*1e-12 * s;   %Energy (in joule) for one 16-bit interger/fixed point add+mul (MAC) in 65nm node

E_RF_to_ALU = Energy_MAC16/16;               %RF to ALU access Energy per bit in Joule (Using data from Eyeriss paper)       
E_DRAM_to_SRAM = (200*Energy_MAC16)/16;      %DRAM to SRAM access Energy per bit in Joule (Using data from Eyeriss paper)

% Base case
Energy_Add32 = 0.1 * 1e-12 * s; %Energy (in joule) for one 32-bit interger/fixed point addition in 65nm node      (direct data from energy computing paper used)
Energy_Mul8 = 0.2 * 1e-12 * s;  %Energy (in joule) for one 8-bit interger/fixed point multiplication in 65nm node (direct data from energy computing paper used)

%Energy parameters are quadratically scaled for multiplication and linearly scaled for addition and memory access 
Energy_Add_perbit = Energy_Add32/32; %Energy (in joule) for one-bit interger/fixed point addition in 65nm node

%%% from CACTI to Eyeriss projection for the data access energy from SRAM to RF
load SRAMSize           % SRAM size in kB
load E_EySRAM_perbit    % SRAM to RF access energy per bit in Joule, 65 nm

if (SRAM_filter/(8*1024) < SRAMSize(1)) || (SRAM_ifmap/(8*1024) < SRAMSize(1)) || (SRAM_psum/(8*1024) < SRAMSize(1)) ||...
   (SRAM_filter/(8*1024) > SRAMSize(end)) || (SRAM_ifmap/(8*1024) > SRAMSize(end)) || (SRAM_psum/(8*1024) > SRAMSize(end))   
    disp("At least one of the SRAM size is outside the range of available CACTI data")
    return
end

if (floor(SRAM_filter/(8*1024)) ~= (SRAM_filter/(8*1024))) || (floor(SRAM_ifmap/(8*1024)) ~= (SRAM_ifmap/(8*1024))) ||...
        (floor(SRAM_psum/(8*1024)) ~= (SRAM_psum/(8*1024)))
    disp("CACTI data is not available for Fractional kB of SRAM size")  
    return
end

%% Performance Metric (PM)
% 1: Performance metric is Energy * Delay (Delay denotes the per image inference time of a network)
% 2: Performance metric is Energy^2 * Delay
% 3: Performance metric is Energy * Delay^2
% 4: Performance metric is Energy
% 5: Performance metric is Delay
Pmat_flag = 1;     % flag to determine the performance metric; can take values 1, 2, 3, 4, 5

%RF sizes should be less than their respective SRAM
Array_RF_ifmap = Array_row * bw_ifmap; % in bit, Total reg file storage for ifmap in a column,same ifmap vector get copied to all column
Array_RF_filter = Array_row * Array_column * bw_filter; % in bit
Array_RF_psum = Array_column * bw_psum;  % in bit

% will not run the simulation if the respective SRAM sizes are not larger than the register file storage in the PE array
% There is a minimum requirement of filter and ifmap SRAM inside the code of each branch as well.
% For psum-SRAM this if condition will help since for x-> branches, currently there is no constraint placed on psum-SRAM.
if (SRAM_filter >= Array_RF_filter) && (SRAM_ifmap >= Array_RF_ifmap) && (SRAM_psum >= Array_RF_psum)
    % Hardware parameters together to pass
    Hardware_param = [bw_filter bw_ifmap bw_psum bw_ofmap Array_row Array_column SRAM_filter SRAM_ifmap SRAM_psum BW_DRAM DRAM_block_size];
    
    %Preparing the remaining Tech Param                                                                                
    %Obtaining SRAM energy from CACTI projection data
    fsindex = find(SRAM_filter == SRAMSize * 1024 * 8); % index of corresponding filter-SRAM
    isindex = find(SRAM_ifmap == SRAMSize * 1024 * 8);
    psindex = find(SRAM_psum == SRAMSize * 1024 * 8);

    E_fsSRAM_to_RF = E_EySRAM_perbit(fsindex);  % access energy per bit in joule for filter-SRAM
    E_isSRAM_to_RF = E_EySRAM_perbit(isindex);  % access energy per bit in joule for ifmap-SRAM
    E_psSRAM_to_RF = E_EySRAM_perbit(psindex);  % access energy per bit in joule for psum-SRAM

    %Calculating Energy for Multiplication using quadratic scaling
    scaling_factor = bw_filter/8;   % 8-bit MUL is my base case
    Energy_Mul_element = Energy_Mul8 * (scaling_factor^2); %Energy(in joule) for one bw_filter bit interger/fixed point mul in 65nm                                                                       

    % Technology parameters together to pass
    Tech_param = [Energy_Add_perbit, Energy_Mul_element, E_RF_to_ALU, E_fsSRAM_to_RF, E_isSRAM_to_RF, E_psSRAM_to_RF, E_DRAM_to_SRAM];                                                                                                                     

    % Calling the function to perform the computation
    [Combined_Table, Table_LOS] = LOS_Computation (Hardware_param, Tech_param, Pmat_flag, Network_flag);
else
    disp("At least one of the SRAMs is too small, increase the size")
end

% Outputs:
Table_LOS           % Table showing the optimal "Best_branch" to compute each layer and its associated PM value

Combined_Table      % Table showing the PM value for each of the fixed branches (Fixed_Pmat), 
                    % PM value for layer-specific optimal scheduling (Best_Pmat), 
                    % and Percent_penalty of using each of the fixed branches as compared to layer-specific optimal scheduling
                    

if (Network_flag == 1)
    %Performance metric to write based on Pmat_flag
    if (Pmat_flag == 1)
        writetable(Combined_Table,'Energy*Delay_AlexNet.csv')
        writetable(Table_LOS,'Energy*Delay_AlexNet_LOS.csv')
    elseif (Pmat_flag == 2)
        writetable(Combined_Table,'Energy^2*Delay_AlexNet.csv')
        writetable(Table_LOS,'Energy^2*Delay_AlexNet_LOS.csv')
    elseif (Pmat_flag == 3)
        writetable(Combined_Table,'Energy*Delay^2_AlexNet.csv')
        writetable(Table_LOS,'Energy*Delay^2_AlexNet_LOS.csv')
    elseif (Pmat_flag == 4)
        writetable(Combined_Table,'Energy_AlexNet.csv')
        writetable(Table_LOS,'Energy_AlexNet_LOS.csv')
    elseif (Pmat_flag == 5)
        writetable(Combined_Table,'Delay_AlexNet.csv')
        writetable(Table_LOS,'Delay_AlexNet_LOS.csv')
    end
    
elseif (Network_flag == 2)
    %Performance metric to write based on Pmat_flag
    if (Pmat_flag == 1)
        writetable(Combined_Table,'Energy*Delay_VGG.csv')
        writetable(Table_LOS,'Energy*Delay_VGG_LOS.csv')
    elseif (Pmat_flag == 2)
        writetable(Combined_Table,'Energy^2*Delay_VGG.csv')
        writetable(Table_LOS,'Energy^2*Delay_VGG_LOS.csv')
    elseif (Pmat_flag == 3)
        writetable(Combined_Table,'Energy*Delay^2_VGG.csv')
        writetable(Table_LOS,'Energy*Delay^2_VGG_LOS.csv')
    elseif (Pmat_flag == 4)
        writetable(Combined_Table,'Energy_VGG.csv')
        writetable(Table_LOS,'Energy_VGG_LOS.csv')
    elseif (Pmat_flag == 5)
        writetable(Combined_Table,'Delay_VGG.csv')
        writetable(Table_LOS,'Delay_VGG_LOS.csv')
    end
    
elseif (Network_flag == 3)
    %Performance metric to write based on Pmat_flag
    if (Pmat_flag == 1)
        writetable(Combined_Table,'Energy*Delay_SqNet.csv')
        writetable(Table_LOS,'Energy*Delay_SqNet_LOS.csv')
    elseif (Pmat_flag == 2)
        writetable(Combined_Table,'Energy^2*Delay_SqNet.csv')
        writetable(Table_LOS,'Energy^2*Delay_SqNet_LOS.csv')
    elseif (Pmat_flag == 3)
        writetable(Combined_Table,'Energy*Delay^2_SqNet.csv')
        writetable(Table_LOS,'Energy*Delay^2_SqNet_LOS.csv')
    elseif (Pmat_flag == 4)
        writetable(Combined_Table,'Energy_SqNet.csv')
        writetable(Table_LOS,'Energy_SqNet_LOS.csv')
    elseif (Pmat_flag == 5)
        writetable(Combined_Table,'Delay_SqNet.csv')
        writetable(Table_LOS,'Delay_SqNet_LOS.csv')
    end
    
elseif (Network_flag == 4)
    %Performance metric to write based on Pmat_flag
    if (Pmat_flag == 1)
        writetable(Combined_Table,'Energy*Delay_GNet.csv')
        writetable(Table_LOS,'Energy*Delay_GNet_LOS.csv')
    elseif (Pmat_flag == 2)
        writetable(Combined_Table,'Energy^2*Delay_GNet.csv')
        writetable(Table_LOS,'Energy^2*Delay_GNet_LOS.csv')
    elseif (Pmat_flag == 3)
        writetable(Combined_Table,'Energy*Delay^2_GNet.csv')
        writetable(Table_LOS,'Energy*Delay^2_GNet_LOS.csv')
    elseif (Pmat_flag == 4)
        writetable(Combined_Table,'Energy_GNet.csv')
        writetable(Table_LOS,'Energy_GNet_LOS.csv')
    elseif (Pmat_flag == 5)
        writetable(Combined_Table,'Delay_GNet.csv')
        writetable(Table_LOS,'Delay_GNet_LOS.csv')
    end
    
elseif (Network_flag == 5)
    %Performance metric to write based on Pmat_flag
    if (Pmat_flag == 1)
        writetable(Combined_Table,'Energy*Delay_ResNet.csv')
        writetable(Table_LOS,'Energy*Delay_ResNet_LOS.csv')
    elseif (Pmat_flag == 2)
        writetable(Combined_Table,'Energy^2*Delay_ResNet.csv')
        writetable(Table_LOS,'Energy^2*Delay_ResNet_LOS.csv')
    elseif (Pmat_flag == 3)
        writetable(Combined_Table,'Energy*Delay^2_ResNet.csv')
        writetable(Table_LOS,'Energy*Delay^2_ResNet_LOS.csv')
    elseif (Pmat_flag == 4)
        writetable(Combined_Table,'Energy_ResNet.csv')
        writetable(Table_LOS,'Energy_ResNet_LOS.csv')
    elseif (Pmat_flag == 5)
        writetable(Combined_Table,'Delay_ResNet.csv')
        writetable(Table_LOS,'Delay_ResNet_LOS.csv')
    end
    
end
    
    
