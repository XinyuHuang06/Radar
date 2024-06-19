% Example:
% :param :
% :return :
% detailed description:
%------------------------------------------------------------------------------
% Created by: Xinyu Huang.
% On: 15/05/2024.
% Copyright (C) 2024 Xinyu Huang (learning_huang@163.com).
% All Rights Reserved.
% Unauthorized copying of this file, via any medium is strictly prohibited.
% Proprietary and confidential.
%------------------------------------------------------------------------------
% Jsonlab
clc;clear; close all;
tic
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
% x = exp(1j*2*pi*rand(N,1));
% xr = [real(x);imag(x)]; 
% b = exp(1j*2*pi*rand(N,1));
% br = [real(b);imag(b)]; 
c = generator_LFM(InitialParameter.signal.fs,InitialParameter.signal.fc,InitialParameter.signal.B,InitialParameter.signal.T);
cr = [real(c);imag(c)];
xr = cr;
br = cr;
% % Intialize the parameter of ADMM.
% The coef of Lagrange terms.
lambda_0 = 0*ones(2*N,1); 
lambda_1 = 0*ones(N,1); 
% The coef of penalty terms.
rho_0 = 1e-2;
rho_1 = 1e-2*ones(N,1);
% r: The barrier function parameter h: The 
r = 0.01;
h = -0.01;
% The similarty parameter.
vartheta = 0.1;
% % Data Record.
DataRecordPack = DataRecord(Max_ItersNum);
% The initial parameter Caculation
DataSetPackets = DataSet(xr, br, cr, lambda_0, lambda_1, rho_0, rho_1, r, h, vartheta, N);
DataRecordPack.UpdateDataSet(CaculateTargetFun(DataSetPackets.packets, ParameterPackets), 1);
% % 
figure(1)
% % ADMM Iterations
for i_m = 1:Max_ItersNum
    % % ADMM update
    xr = Update_xr(DataSetPackets.packets, ParameterPackets);% % Step 1 , Solving the x_r. 
    DataSetPackets.update(xr,'xr');
    br = Update_br(DataSetPackets.packets, ParameterPackets);% % Step 2 , Solving the b_r
    DataSetPackets.update(br,'br');
    h = Update_h(DataSetPackets.packets, ParameterPackets);% % Step 3 , Solving the h
    DataSetPackets.update(h,'h');
    [rho_0, rho_1] = Update_rho(DataSetPackets, ParameterPackets);% % Step 4 , Solving the rho_0 and rho_1
    [lambda_0, lambda_1] = Update_lambda(DataSetPackets.packets, ParameterPackets);% % Step 5 , Solving the lambda_0 and lambda_1
    DataSetPackets.update(lambda_0,'lambda_0');DataSetPackets.update(lambda_1,'lambda_1');
    % % Other
    DataRecordPack.UpdateDataSet(CaculateTargetFun(DataSetPackets.packets, ParameterPackets), i_m + 1);
    forwaitbar.show_bar; % progress bar
end
clear *_k
% % 
PlotAndExport(DataSetPackets.packets, DataRecordPack, InitialParameter, OutFolderPath);
rmpath('function/');
toc