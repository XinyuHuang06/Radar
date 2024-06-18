clc; clear; close all;
% % Intialization
OutFolderPath = 'Output';
% Load data.
InitialParameter = loadjson('./Data/initial_Parameter.json');
N = round(InitialParameter.signal.T*InitialParameter.signal.fs);
% add function definition
addpath("function/");
Max_ItersNum = InitialParameter.Max_ItersNum;
Threshold = InitialParameter.Threshold;
% A progress bar for code runs.
forwaitbar = forWaitbar(Max_ItersNum);
% Generate the parameter packets
ParameterPackets = GenerateParameter(N, InitialParameter.signal.M, InitialParameter.Flag.Sparse, InitialParameter.Threshold);
% % Generate xr, br, cr.
x = exp(1j*2*pi*rand(N,1));
xr = [real(x);imag(x)]; 
b = exp(1j*2*pi*rand(N,1));
br = [real(b);imag(b)]; 
c = generator_LFM(InitialParameter.signal.fs,InitialParameter.signal.fc,InitialParameter.signal.B,InitialParameter.signal.T);
cr = [real(c);imag(c)];
% % Intialize the parameter of SDR.
% % SDR Iterations
for i_m = 1:Max_ItersNum
    % % ADMM update
    
    forwaitbar.show_bar; % progress bar
end
clear *_k
% % 
PlotAndExport(DataSetPackets.packets, DataRecordPack, InitialParameter, OutFolderPath);
rmpath('function/');
profile viewer