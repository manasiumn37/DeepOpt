clc;
clear all;
close all;

%%% This file performs an exhaustive search to find the optimal allocation of hardware resources based on the user-specified area budget and performance metric
%%% It usually takes several hours to complete the search and to return the optimal hardware specification

%% Setecting CNN topology
% 1: AlexNet
% 2: VGG-16
% 3: SqueezeNet-v1.1
% 4: GoogleNet-v1
% 5: ResNet-50
Network_flag = 4;  % Flag to choose the network topology; can take values 1, 2, 3, 4, 5

%% Initial Data & Parameters
Area_SRAM_per_byte = 32.545670202636720; %SRAM area per byte in um^2; from TSMC 65 nm memory compiler

% Obtaining area of a PE (ecah PE contains 8-bit multiplier and 32 bit adder)
%from post place-and-route implementation of a systolic array with each PE having 2kB of SRAM; PE logics constitute 10% of the area while memory constitues 90% area
Area_SRAM_2kb = Area_SRAM_per_byte * 2 * 1024;    % Area of a 2kB SRAM in um^2
Area_per_PE = 0.1 * (Area_SRAM_2kb / 0.9);        % in um^2

% Imposing a minimum boundary on each SRAM
min_fSRAM_size = 4; %  in kB
min_iSRAM_size = 4; %  in kB
min_pSRAM_size = 4; %  in kB

min_fSRAM_area = min_fSRAM_size * 1024 * Area_SRAM_per_byte; %  in um^2
min_iSRAM_area = min_iSRAM_size * 1024 * Area_SRAM_per_byte; %  in um^2
min_pSRAM_area = min_pSRAM_size * 1024 * Area_SRAM_per_byte; %  in um^2


min_row = 8; % minimum #of row in the PE array; for AlexNet use min_row = 12  (minimum requirement of row for AlexNet is 11)
min_col = 8; % minimum #of column in the PE array

Area_budget =  4000 * 4000; %in um^2

%Maximum boundary based on the area budget
Max_fSRAM_Area = Area_budget - min_iSRAM_area - min_pSRAM_area - min_row * min_col * Area_per_PE;
Max_iSRAM_Area = Area_budget - min_fSRAM_area - min_pSRAM_area - min_row * min_col * Area_per_PE;
Max_pSRAM_Area = Area_budget - min_iSRAM_area - min_fSRAM_area - min_row * min_col * Area_per_PE;

Max_fSRAM_size = (Max_fSRAM_Area/Area_SRAM_per_byte)/1024;   % in kB
Max_iSRAM_size = (Max_iSRAM_Area/Area_SRAM_per_byte)/1024;   % in kB
Max_pSRAM_size = (Max_pSRAM_Area/Area_SRAM_per_byte)/1024;   % in kB

Max_nos_of_PE = floor((Area_budget - min_iSRAM_area - min_fSRAM_area - min_pSRAM_area)/Area_per_PE);
Max_row = floor(Max_nos_of_PE/min_col);
Max_col = floor(Max_nos_of_PE/min_row);

stepj = 2; %step size for #of PE row
stepk = 2; %step size for #of PE column

% row_vect_size = length(min_row:stepj:Max_row); % to estimate simulation time
% col_vect_size = length(min_col:stepk:Max_col);

if (Max_fSRAM_Area < min_fSRAM_area) || (Max_iSRAM_Area < min_iSRAM_area) || (Max_pSRAM_Area < min_pSRAM_area) || (Max_row < min_row) || (Max_col < min_col)
    disp ("Area budget is too small, increase it")
    return
end

%Creating the array of SRAM sizes which are power of 2 and within the min & max SRAM sizes
% filter SRAM
minf_pow = ceil(log(min_fSRAM_size)/log(2));
maxf_pow = floor(log(Max_fSRAM_size)/log(2));
fSRAM_array = 2.^(minf_pow : maxf_pow);       % in kB
% ifmap SRAM
mini_pow = ceil(log(min_iSRAM_size)/log(2));
maxi_pow = floor(log(Max_iSRAM_size)/log(2));
iSRAM_array = 2.^(mini_pow : maxi_pow);       % in kB
% psum SRAM
minp_pow = ceil(log(min_pSRAM_size)/log(2));
maxp_pow = floor(log(Max_pSRAM_size)/log(2));
pSRAM_array = 2.^(minp_pow : maxp_pow);       % in kB


%% Some Accelerator Specification
% bit width
% The technology parameters used for this search are for 8-bit mul and 32-bit add. The area per PE is also for these bit widths
% Hence if the following bit-widths are changed the technology parameters and area per PE will also have to be adjusted accordingly.
bw_filter = 8;
bw_ifmap = 8;
bw_psum = 32;
bw_ofmap = 8;

% To incorporation DRAM bandwidth induced stall in the calculation of cycle count
BW_DRAM = 128; %in bit/cycle
DRAM_block_size = 128; %in bit

%% Some Access energy specification 
s = (65/45)*(1/0.9)^2;                      %Technology scaling parameter; 65nm Technology node, Vdd = 1V; 45nm Technology node, Vdd = 0.9V

% to obtain 16-bit add and mul, the direct 8-bit add and mul data from energy computing paper are linearly and quadratically sacled, respectively
Energy_MAC16 = (0.03 * 2 + 0.2 * 4)*1e-12 * s;   %Energy (in joule) for one 16-bit interger/fixed point add+mul (MAC) in 65nm node

Energy_Add32 = 0.1 * 1e-12 * s; %Energy (in joule) for one 32-bit interger/fixed point addition in 65nm node      (direct data from energy computing paper used)
Energy_Mul8 = 0.2 * 1e-12 * s;  %Energy (in joule) for one 8-bit interger/fixed point multiplication in 65nm node (direct data from energy computing paper used)

E_RF_to_ALU = Energy_MAC16/16;               %RF to ALU access Energy per bit in Joule (Using data from Eyeriss paper)      
E_DRAM_to_SRAM = (200*Energy_MAC16)/16;      %DRAM to SRAM access Energy per bit in Joule (Using data from Eyeriss paper)

%%% from CACTI to Eyeriss projection for the data access energy from SRAM to RF
load SRAMSize           % SRAM size in kB
load E_EySRAM_perbit    % SRAM access energy per bit in Joule, 65 nm

if (fSRAM_array(1) < SRAMSize(1)) || (iSRAM_array(1) < SRAMSize(1)) || (pSRAM_array(1) < SRAMSize(1)) ||...
   (fSRAM_array(end) > SRAMSize(end)) || (iSRAM_array(end) > SRAMSize(end)) || (pSRAM_array(end) > SRAMSize(end))   
    disp("either min or max SRAM size is outside the available CACTI data")
    return
end

%% Performance Metric (PM)
% 1: Performance metric is Energy * Delay (Delay denotes the per image inference time of a network)
% 2: Performance metric is Energy^2 * Delay
% 3: Performance metric is Energy * Delay^2
% 4: Performance metric is Energy
% 5: Performance metric is Delay
Pmat_flag = 1;     % flag to determine the performance metric

%% Performing the hardware sweep
Best_Pmat = inf;
Best_Area = inf;

sweep_point = 1;
for j = min_row : stepj : Max_row
    for k = min_col : stepk : Max_col
        for fs = 1:1:length(fSRAM_array)
            for is = 1:1:length(iSRAM_array)
                for ps = 1:1:length(pSRAM_array)
                    
                    max_jk = max(j,k);
                    min_jk = min(j,k);
                    if (min_jk >= 0.5 * max_jk) % This if statement is to have a spec where row (col) is within 50% of col(row)
                    
                        Chip_Area = (j*k*Area_per_PE) + (fSRAM_array(fs) + iSRAM_array(is) + pSRAM_array(ps)) * 1024 * Area_SRAM_per_byte;  % in um^2 
                        if Chip_Area <= Area_budget  
                            % MAC array
                            Array_row = j;
                            Array_column = k;

                            % On-chip SRAM storage
                            SRAM_filter = fSRAM_array(fs) * 1024 * 8; % in bit 
                            SRAM_ifmap = iSRAM_array(is) * 1024 * 8; % in bit 
                            SRAM_psum = pSRAM_array(ps) * 1024 * 8; % in bit 

                            %These sizes should be less than their respective SRAM
                            Array_RF_ifmap = Array_row * bw_ifmap; % in bit, Total reg file storage for ifmap in a column,same ifmap vector get copied to all column
                            Array_RF_filter = Array_row * Array_column * bw_filter; % in bit
                            Array_RF_psum = Array_column * bw_psum;  % in bit

                            % will not run the combinations if the respective SRAM sizes are not larger than the register file storage in the PE array
                            % There is a minimum requirement of filter and ifmap SRAM inside the code of each branch as well.
                            % For psum-SRAM this if condition will help since for x-> branches, there is no constraint on psum-SRAM.
                            if (SRAM_filter >= Array_RF_filter) && (SRAM_ifmap >= Array_RF_ifmap) && (SRAM_psum >= Array_RF_psum)
                                Spec = [j k fSRAM_array(fs) iSRAM_array(is) pSRAM_array(ps)];
                                % Hardware parameters together to pass
                                Hardware_param = [bw_filter bw_ifmap bw_psum bw_ofmap Array_row Array_column SRAM_filter SRAM_ifmap SRAM_psum BW_DRAM DRAM_block_size];

                                %Obtaining SRAM energy from CACTI projection data
                                fsindex = find(SRAM_filter == SRAMSize * 1024 * 8); % index of corresponding filter-SRAM
                                isindex = find(SRAM_ifmap == SRAMSize * 1024 * 8);
                                psindex = find(SRAM_psum == SRAMSize * 1024 * 8);

                                E_fsSRAM_to_RF = E_EySRAM_perbit(fsindex);  % access energy per bit in joule for filter-SRAM
                                E_isSRAM_to_RF = E_EySRAM_perbit(isindex);  % access energy per bit in joule for ifmap-SRAM
                                E_psSRAM_to_RF = E_EySRAM_perbit(psindex);  % access energy per bit in joule for psum-SRAM

                                % Technology parameters together to pass
                                Tech_param = [Energy_Add32, Energy_Mul8, E_RF_to_ALU, E_fsSRAM_to_RF, E_isSRAM_to_RF, E_psSRAM_to_RF, E_DRAM_to_SRAM];

                                % Calling the function to perform the computation
                                Pmat_L2Lbest = Opt_Hardware_Allocation (Hardware_param, Tech_param, Pmat_flag, Network_flag); %performance metric based on layer to layer best branching
                                
                                if (Network_flag == 1)
                                    All_Res_Alex(sweep_point,:) = [Spec Chip_Area Pmat_L2Lbest];
                                elseif (Network_flag == 2)
                                    All_Res_VGG(sweep_point,:) = [Spec Chip_Area Pmat_L2Lbest];
                                elseif (Network_flag == 3)
                                    All_Res_SqNet(sweep_point,:) = [Spec Chip_Area Pmat_L2Lbest];
                                elseif (Network_flag == 4)
                                    All_Res_GNet(sweep_point,:) = [Spec Chip_Area Pmat_L2Lbest];
                                elseif (Network_flag == 5)
                                    All_Res_ResNet(sweep_point,:) = [Spec Chip_Area Pmat_L2Lbest];
                                end
   
                                sweep_point = sweep_point + 1;

                                if (Best_Pmat > Pmat_L2Lbest)
                                    Best_Pmat = Pmat_L2Lbest;
                                    Best_Area = Chip_Area;     % in um^2
                                    Best_Spec = Spec;                                
                                end

                                if (Best_Pmat == Pmat_L2Lbest)
                                    if (Best_Area > Chip_Area)
                                        Best_Pmat = Pmat_L2Lbest;
                                        Best_Area = Chip_Area;     % in um^2
                                        Best_Spec = Spec;    
                                    end
                                end


                            end 
                        end
                    end 
                    
                end
            end
        end
    end
end

%Outputs:
Best_Pmat  % Value of performance metric at the Optimal hardeare Spec.
Best_Area  % Area at the Optimal hardeare Spec.
Best_Spec  % Optimal hardware Spec.

%% Unit of the output Performance Metric
% Energy*Delay in Joule * cycle
% Energy^2*Delay in Joule^2 * cycle
% Energy*Delay^2 in Joule * cycle^2
% Energy in Joule
% Delay in cycle

disp ("computation done for all the loop")
%Storing performance metric based on Pmat_flag
if (Pmat_flag == 1)
    Best_row = Best_Spec(1);
    Best_col = Best_Spec(2);
    Best_fSRAM_kB = Best_Spec(3);
    Best_iSRAM_kB = Best_Spec(4);
    Best_pSRAM_kB = Best_Spec(5);
    Best_EnergyDelay = Best_Pmat;
    Best_Area_um2 = Best_Area;
    Best_Res_Network = table(Best_row, Best_col, Best_fSRAM_kB, Best_iSRAM_kB, Best_pSRAM_kB, Best_EnergyDelay, Best_Area_um2)
    if (Network_flag == 1)
        writetable(Best_Res_Network,'Best_Energy*Delay_Alex.csv');
    elseif (Network_flag == 2)
        writetable(Best_Res_Network,'Best_Energy*Delay_VGG.csv');
    elseif (Network_flag == 3)
        writetable(Best_Res_Network,'Best_Energy*Delay_SqNet.csv');
    elseif (Network_flag == 4)
        writetable(Best_Res_Network,'Best_Energy*Delay_GNet.csv');
    elseif (Network_flag == 5)
        writetable(Best_Res_Network,'Best_Energy*Delay_ResNet.csv');
    end
    
elseif (Pmat_flag == 2)
    Best_row = Best_Spec(1);
    Best_col = Best_Spec(2);
    Best_fSRAM_kB = Best_Spec(3);
    Best_iSRAM_kB = Best_Spec(4);
    Best_pSRAM_kB = Best_Spec(5);
    Best_Energy2Delay = Best_Pmat;
    Best_Area_um2 = Best_Area;
    Best_Res_Network = table(Best_row, Best_col, Best_fSRAM_kB, Best_iSRAM_kB, Best_pSRAM_kB, Best_Energy2Delay, Best_Area_um2)
    if (Network_flag == 1)
        writetable(Best_Res_Network,'Best_Energy^2*Delay_Alex.csv');
    elseif (Network_flag == 2)
        writetable(Best_Res_Network,'Best_Energy^2*Delay_VGG.csv');
    elseif (Network_flag == 3)
        writetable(Best_Res_Network,'Best_Energy^2*Delay_SqNet.csv');
    elseif (Network_flag == 4)
        writetable(Best_Res_Network,'Best_Energy^2*Delay_GNet.csv');
    elseif (Network_flag == 5)
        writetable(Best_Res_Network,'Best_Energy^2*Delay_ResNet.csv');
    end
   
elseif (Pmat_flag == 3)
    Best_row = Best_Spec(1);
    Best_col = Best_Spec(2);
    Best_fSRAM_kB = Best_Spec(3);
    Best_iSRAM_kB = Best_Spec(4);
    Best_pSRAM_kB = Best_Spec(5);
    Best_EnergyDelay2 = Best_Pmat;
    Best_Area_um2 = Best_Area;
    Best_Res_Network = table(Best_row, Best_col, Best_fSRAM_kB, Best_iSRAM_kB, Best_pSRAM_kB, Best_EnergyDelay2, Best_Area_um2)
    if (Network_flag == 1)
        writetable(Best_Res_Network,'Best_Energy*Delay^2_Alex.csv'); 
    elseif (Network_flag == 2)
        writetable(Best_Res_Network,'Best_Energy*Delay^2_VGG.csv'); 
    elseif (Network_flag == 3)
        writetable(Best_Res_Network,'Best_Energy*Delay^2_SqNet.csv'); 
    elseif (Network_flag == 4)
        writetable(Best_Res_Network,'Best_Energy*Delay^2_GNet.csv'); 
    elseif (Network_flag == 5)
        writetable(Best_Res_Network,'Best_Energy*Delay^2_ResNet.csv'); 
    end
   
elseif (Pmat_flag == 4)
    Best_row = Best_Spec(1);
    Best_col = Best_Spec(2);
    Best_fSRAM_kB = Best_Spec(3);
    Best_iSRAM_kB = Best_Spec(4);
    Best_pSRAM_kB = Best_Spec(5);
    Best_Energy_J = Best_Pmat;
    Best_Area_um2 = Best_Area;
    Best_Res_Network = table(Best_row, Best_col, Best_fSRAM_kB, Best_iSRAM_kB, Best_pSRAM_kB, Best_Energy_J, Best_Area_um2)   
    if (Network_flag == 1)
        writetable(Best_Res_Network,'Best_Energy_Alex.csv'); 
    elseif (Network_flag == 2)
        writetable(Best_Res_Network,'Best_Energy_VGG.csv'); 
    elseif (Network_flag == 3)
        writetable(Best_Res_Network,'Best_Energy_SqNet.csv'); 
    elseif (Network_flag == 4)
        writetable(Best_Res_Network,'Best_Energy_GNet.csv'); 
    elseif (Network_flag == 5)
        writetable(Best_Res_Network,'Best_Energy_ResNet.csv'); 
    end
       
elseif (Pmat_flag == 5)
    Best_row = Best_Spec(1);
    Best_col = Best_Spec(2);
    Best_fSRAM_kB = Best_Spec(3);
    Best_iSRAM_kB = Best_Spec(4);
    Best_pSRAM_kB = Best_Spec(5);
    Best_Delay = Best_Pmat;
    Best_Area_um2 = Best_Area;
    Best_Res_Network = table(Best_row, Best_col, Best_fSRAM_kB, Best_iSRAM_kB, Best_pSRAM_kB, Best_Delay, Best_Area_um2)
    if (Network_flag == 1)
        writetable(Best_Res_Network,'Best_Delay_Alex.csv');  
    elseif (Network_flag == 2)
        writetable(Best_Res_Network,'Best_Delay_VGG.csv');  
    elseif (Network_flag == 3)
        writetable(Best_Res_Network,'Best_Delay_SqNet.csv');  
    elseif (Network_flag == 4)
        writetable(Best_Res_Network,'Best_Delay_GNet.csv');  
    elseif (Network_flag == 5)
        writetable(Best_Res_Network,'Best_Delay_ResNet.csv');  
    end
    
end


% Saving result for all hardware point in the search
if (Network_flag == 1)
    save All_Res_Alex All_Res_Alex       % saving all res for a single performance metric
elseif (Network_flag == 2)
    save All_Res_VGG All_Res_VGG         % saving all res for a single performance metric
elseif (Network_flag == 3)
    save All_Res_SqNet All_Res_SqNet     % saving all res for a single performance metric
elseif (Network_flag == 4)
    save All_Res_GNet All_Res_GNet       % saving all res for a single performance metric
elseif (Network_flag == 5)
    save All_Res_ResNet All_Res_ResNet   % saving all for a single performance metric
end



