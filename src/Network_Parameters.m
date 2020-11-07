function [Layer, filter_size_Net, Nos_of_Filter_Net, Ifmap_size_Net, Nos_of_Channel_Net, Stride_Net, Ofmap_size_Net, Batch_size_Net, fc_flag]...
                                                                                                                        = Network_Parameters (Network_flag)
%Function for the Network Parameters

%For Network_flag:
% 1: AlexNet
% 2: VGG-16
% 3: SqueezeNet-v1.1
% 4: GoogleNet-v1
% 5: ResNet-50

if (Network_flag == 1)
    %% AlexNet Layer Input Parameters (The layer input parameters are directly provided here instead of a text file since the number of layers are small)
    Layer = [1 2 3 4 5 6 7 8]';                                      % conv & FC layer index number
    filter_size_Net = [11 5 3 3 3 1 1 1]';                           % Height/Width of filter
    Nos_of_Filter_Net = [96 128 384 192 128 4096 4096 1000]';        % #of total 3D filters
    Ifmap_size_Net = [227 31 15 15 15 1 1 1]';                       % Height/Width of padded ifmap
    Nos_of_Channel_Net = [3 48 256 192 192 9216 4096 4096]';         % #of total channels in the ifmap/filter
    Stride_Net = [4 1 1 1 1 1 1 1]';                                 % Convolution stride
    Ofmap_size_Net = [55 27 13 13 13 1 1 1]';                        % Height/Width of ofmap
    Batch_size_Net = [1 1 1 1 1 16 16 16]';                          % Batch size, The option to use a batch size > 1 is added only for FC layer
    %Batch_size_Net = [1 1 1 1 1 1 1 1]';

    fc_flag = [0 0 0 0 0 1 1 1]';    % the flag is 0 for conv layers and 1 for FC layers

    % For the FC layers, ifmap/ofmap & filter sizes are 1;
    % C2, C4, C5 needs to be multiplied by 2. In AlexNet, to reduce #of computation output channels are grouped into 2. This is equivalent to running them seperately. 
    % The "Nos_of_Filter" params of C2, C4, C5 are adjusted accordingly.
 
elseif (Network_flag == 2)
    %% VGG-16 Layer Input Parameters (The layer input parameters are directly provided here instead of a text file since the number of layers are small)
    Layer = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16]';                                           % layer index number
    filter_size_Net = [3 3 3 3 3 3 3 3 3 3 3 3 3 1 1 1]';                                        % Height/Width of filter
    Nos_of_Filter_Net = [64 64 128 128 256 256 256 512 512 512 512 512 512 4096 4096 1000]';     % #of total 3D filters
    Ifmap_size_Net = [226 226 114 114 58 58 58 30 30 30 16 16 16 1 1 1]';                        % Height/Width of padded ifmap
    Nos_of_Channel_Net = [3 64 64 128 128 256 256 256 512 512 512 512 512 25088 4096 4096]';     % #of total channels in the ifmap/filter
    Stride_Net = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]';                                             % Convolution stride
    Ofmap_size_Net = [224 224 112 112 56 56 56 28 28 28 14 14 14 1 1 1]';                        % Height/Width of ofmap
    Batch_size_Net = [1 1 1 1 1 1 1 1 1 1 1 1 1 16 16 16]';                                 % Batch size,The option to use a batch size > 1 is added only for FC layer
    %Batch_size_Net = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]';                                        

    fc_flag = [0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1]';

    % For the FC layers, ifmap/ofmap & filter sizes are 1;
    
elseif (Network_flag == 3)
    %% SqueezeNet-v11 Layer Input Parameters
    load 'SqueezeNet_v11_Input_File.txt';  %.txt file which contains the layer input parameters for SqueezeNet-v11

    Layer = SqueezeNet_v11_Input_File(:,1);                  % layer index number
    Ifmap_size_unpad = SqueezeNet_v11_Input_File(:,2);       % Height/Width of ifmap without padding
    Nos_of_Channel_Net = SqueezeNet_v11_Input_File(:,3);     % #of total channels in the ifmap/filter
    filter_size_Net = SqueezeNet_v11_Input_File(:,4);        % Height/Width of filter
    Nos_of_Filter_Net = SqueezeNet_v11_Input_File(:,5);      % #of total 3D filters
    Ofmap_size_Net = SqueezeNet_v11_Input_File(:,6);         % Height/Width of ofmap
    Stride_Net = SqueezeNet_v11_Input_File(:,7);             % Convolution stride
    Pad = SqueezeNet_v11_Input_File(:,8);                    % Amount of padding in the ifmap
    Ifmap_size_Net = Ifmap_size_unpad + (2.*Pad);            % Height/Width of padded ifmap
    Batch_size_Net = ones(length(Layer),1);                  % Batch size, all 1 since there is no FC layer in SqueezeNet-v1.1
                                                             % The option to use a batch size > 1 is added only for FC layer

    fc_flag = zeros(length(Layer),1);                        % No FC layer in SqueezeNet

elseif (Network_flag == 4)
    %% GoogleNet-v1 Layer Input Parameters
    load 'GoogleNet_v1_Input_File.txt';      %.txt file which contains the layer input parameters for GoogleNet-v1

    Layer = GoogleNet_v1_Input_File(:,1);                  % layer index number
    Ifmap_size_unpad = GoogleNet_v1_Input_File(:,2);       % Height/Width of ifmap without padding
    Nos_of_Channel_Net = GoogleNet_v1_Input_File(:,3);     % #of total channels in the ifmap/filter
    filter_size_Net = GoogleNet_v1_Input_File(:,4);        % Height/Width of filter
    Nos_of_Filter_Net = GoogleNet_v1_Input_File(:,5);      % #of total 3D filters
    Ofmap_size_Net = GoogleNet_v1_Input_File(:,6);         % Height/Width of ofmap
    Stride_Net = GoogleNet_v1_Input_File(:,7);             % Convolution stride
    Pad = GoogleNet_v1_Input_File(:,8);                    % Amount of padding in the ifmap
    Ifmap_size_Net = Ifmap_size_unpad + (2.*Pad);          % Height/Width of padded ifmap
    Batch_size_Net = ones(length(Layer),1);                % Batch size 
    Batch_size_Net(end) = 16;                              % Only the last layer is FC. Hence batching > 1 is used only for the last layer
                                                           % The option to use a batch size > 1 is added only for FC layer

    fc_flag = zeros(length(Layer),1);
    fc_flag(end) = 1;                                      % Only the last layer is FC
    
    
    
elseif (Network_flag == 5)
    %% ResNet-50 Layer Input Parameters
    load 'ResNet50_Input_File.txt';  %.txt file which contains the layer input parameters for ResNet-50

    Layer = ResNet50_Input_File(:,1);                  % layer index number
    Ifmap_size_unpad = ResNet50_Input_File(:,2);       % Height/Width of ifmap without padding
    Nos_of_Channel_Net = ResNet50_Input_File(:,3);     % #of total channels in the ifmap/filter
    filter_size_Net = ResNet50_Input_File(:,4);        % Height/Width of filter
    Nos_of_Filter_Net = ResNet50_Input_File(:,5);      % #of total 3D filters
    Ofmap_size_Net = ResNet50_Input_File(:,6);         % Height/Width of ofmap
    Stride_Net = ResNet50_Input_File(:,7);             % Convolution stride
    Pad = ResNet50_Input_File(:,8);                    % Amount of padding in the ifmap
    Ifmap_size_Net = Ifmap_size_unpad + (2.*Pad);      % Height/Width of padded ifmap
    Batch_size_Net = ones(length(Layer),1);            % Batch size
    Batch_size_Net(end) = 16;                          % Only the last layer is FC. Hence batching > 1 is used only for the last layer
                                                       % The option to use a batch size > 1 is added only for FC layer

    fc_flag = zeros(length(Layer),1);  
    fc_flag(end) = 1;                                  % Only the last layer is FC
    
    
else
    disp("Invalid Network flag")
    return
end


end
