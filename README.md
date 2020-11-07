# DeepOpt: 
-A Tool for layer-specific and hardware-specific optimized scheduling of CNN Workload for ASIC-based systolic accelerator.

-A Tool to optimally allocate hardware resources under user-specified area budget and performance metric (i.e., Energy, Delay, Energy×Delay, Energy^2×Delay, and Energy×Delay^2).

If you use any part of this project please cite:

S. D. Manasi and S. S. Sapatnekar, “DeepOpt: Optimized Scheduling of CNN Workloads for ASIC-based Systolic Deep Learning Accelerators”, Asia and South Pacific Design Automation Conference (ASP-DAC), 2021 [Accepted].

The directory named "src" and "data" contain all the required files to use the tool.

## Guideline to use DeepOpt

**Files in the "data" directory:**

(1) .txt files: 

-These files contain the shape parameters of all the layers for SqueezeNet-v1.1, GoogleNet-v1, and ResNet-50.

-There are two files for each network: one with the description of the format of each row and another without any description.

-Used as input by the Main.m files in the "src" directory.

(2) .mat files:

-These files contain a range of SRAM sizes and associated data access energy extracted from CACTI

-Used as input by the Main.m files in the "src" directory.

**Files in the "src" directory:**

(1) LOS_Main.m: The main file which computes the performance metric (PM) of an input network using layer-specific optimal scheduling (LOS)

-Input: CNN topology, Hardware Specification, Technology parameters 

-Input: Performance Metric, one out of five PMs -- Energy, Delay, Energy×Delay, Energy^2×Delay, and Energy×Delay^2

-Output: Name of the optimal scheduling branch for each layer

-Output: Value of the PM using LOS and its improvement over the five fixed scheduling (FS) schemes

(2) LOS_Computation.m: function to compute the PM of a network (used by LOS_Main.m)

(3) Network_Parameters.m: function to obtain the parameter of all the networks (AlexNet, VGG-16, SqueezeNet-v1.1, GoogleNet-v1, and ResNet-50).

(4) OHA_Main.m: The Main file to determine the optimal allocation of hardware resources based on a user-specified area budget and performance metric.

(5) Opt_Hardware_Allocation.m: function to compute PM of a network using LOS (used by OHA_Main.m)

(6) "private" directory: functions of the implementations of all the scheduling schemes.

**Note:**
In order to run the main.m files, place all the files of the "data" directory inside the "private" folder of the "src" directory
