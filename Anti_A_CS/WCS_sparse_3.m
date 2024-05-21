% Example:
% :param :
% :return :
% detailed description:
%------------------------------------------------------------------------------
% Created by: Xinyu Huang.
% On: 16/04/2024.
% Copyright (C) 2024 Xinyu Huang (learning_huang@163.com).
% All Rights Reserved.
% Unauthorized copying of this file, via any medium is strictly prohibited.
% Proprietary and confidential.
%------------------------------------------------------------------------------
profile on
clc;clear;close all;
% add function definition
addpath("function/");
% % Parameter Setting
Flag_PlotExport = 1;
Flag_Sparse = 1;
Flag_PrintLog = 0;
Max_ItersNum = 20; % 最大迭代次数
Threshold = 20;
% % Data IntializationS
N = 128;
M = 8;
fs = 2e6;
T = N/fs;
fc = 1e5;
temp_B = 4e5;
c = generator_LFM(fs,fc,temp_B,T);
cr = [real(c);imag(c)];
Data = Generate_data(N, M, Flag_Sparse, Threshold);
x = exp(1j*2*pi*rand(N,1));
xr = [real(x);imag(x)]; 
% br = xr;
% xr = cr;
b = exp(1j*2*pi*rand(N,1));
br = [real(b);imag(b)]; 
xr = xr/sqrt(N);
br = br/sqrt(N);
cr = cr/sqrt(N);
% % ADMM参数初始化
lambda_0 = 1*ones(2*N,1);
lambda_1 = 1*ones(N,1); 
rho_0 = 1;
rho_1 = 1*ones(N,1); % 罚函数项系数
r = 0.001; % 障碍函数系数
h = -0.01; % 相似度约束辅助变量
xi_inc = 3;
xi_dec = 2;
vartheta = (0.1/sqrt(N)); % 参考信号相关度阈值
% % Data Record
% % The initial parameter Caculation
xr_init = xr;
TarS = struct('sum',zeros(Max_ItersNum+1,1),'Tar',zeros(Max_ItersNum+1,1),'Lar_1',zeros(Max_ItersNum+1,1),'pen_1',zeros(Max_ItersNum+1,1),...
                    'Lar_2',zeros(Max_ItersNum+1,1),'pen_2',zeros(Max_ItersNum+1,1));
TarS_2 = struct('sum',zeros(Max_ItersNum+1,1),'Tar',zeros(Max_ItersNum+1,1),'Lar_1',zeros(Max_ItersNum+1,1),'pen_1',zeros(Max_ItersNum+1,1),...
    'Lar_2',zeros(Max_ItersNum+1,1),'pen_2',zeros(Max_ItersNum+1,1));
[TarS.sum(1),TarS.Tar(1), TarS.Lar_1(1), TarS.pen_1(1), TarS.Lar_2(1), TarS.pen_2(1)] = ...
    Debug_fun(xr, br, cr, rho_0, rho_1, lambda_0, lambda_1, h, vartheta, Data);
forwaitbar = forWaitbar(Max_ItersNum);
parameter_record = struct('rho_0',cell(Max_ItersNum+1,1));
% parameter_record.rho_0{1} = rho_0;
% % ADMM
% fig1 = figure;
FNr = Data.FNr;
for i_m = 1:Max_ItersNum
    % clf(fig1);
    % hold on;plot(xr);plot(br);plot(cr);hold off;legend('xr','br','cr');SetDrawStyle(fig1);
    xr_k = xr;
    br_k = br;
    h_k = h;
    lambda_0_k = lambda_0;
    lambda_1_k = lambda_1;
    rho_0_k = rho_0;
    rho_1_k = rho_1;
    % % Step 1 , Solving the x_r. 
    [A1,BT1,~] = Caculate_matrix_xr(xr, br, cr, rho_0, rho_1, lambda_0, lambda_1, h, vartheta, Data);
    [xr] = Update_xr(A1,BT1,xr);
    % out1 = xr_k'*A1*xr_k + BT1*xr_k;
    % out2 = xr'*A1*xr + BT1*xr;
    % % Step 2 , Solving the b_r
    [br] = Update_br(xr, br, cr, rho_0, rho_1, lambda_0, lambda_1, h, vartheta, Data);
    % Step 3 , Solving the h
    h = Update_h(xr, br, cr, rho_1, lambda_1, vartheta, r, Data);
    % Step 4 , Solving the u, v_n
    [lambda_0, lambda_1] = Update_lambda(xr, br, cr, rho_0, rho_1, lambda_0, lambda_1, h, vartheta, Data);
    % Step 5 , Update rho_0 rho_1
    rho_0 = rho_0*1.1;
    rho_1 = rho_1*1.1;
    % [rho_0,rho_1,delta_0,delta_1] = Update_rho(rho_0_k,rho_1_k,xr_k,br_k,h_k,xr,br,cr,vartheta,h_k,chi_matrix,Flag_Sparse,xi_inc,xi_dec);
    % 参数记录
    forwaitbar.show_bar;
    % [A1,BT1,Tar_seq] = Caculate_matrix_xr(N, FNr,xr, br, cr, omega_alpha_1, rho_0, rho_1, lambda_0, lambda_1, h, vartheta, Taf_1, Taf_2, chi_matrix, Flag_Sparse);
    % TarS_2.sum(i_m) = Tar_seq(1);TarS_2.Tar(i_m) = Tar_seq(2);TarS_2.Lar_1(i_m) = Tar_seq(3);TarS_2.pen_1(i_m) = Tar_seq(4);TarS_2.Lar_2(i_m) = Tar_seq(5);TarS_2.pen_2(i_m) = Tar_seq(6);
    % [TarS.sum(i_m+1),TarS.Tar(i_m+1), TarS.Lar_1(i_m+1), TarS.pen_1(i_m+1), TarS.Lar_2(i_m+1), TarS.pen_2(i_m+1)] = ...
    %     Debug_fun(N, FNr, xr, br, cr, omega_alpha_1, rho_0, rho_1, lambda_0, lambda_1, h, vartheta, Taf_1, Taf_2, chi_matrix, Flag_Sparse);
    % parameter_record.rho_0{i_m+1} = rho_0;
end
clear *_k
% s_xr_init = RadarSignal(xr_init,'xr_init');s_xr_post = RadarSignal(xr,'xr_post');s_cr = RadarSignal(cr,'cr');
% print_log(TarS,s_xr_init,s_xr_post,s_cr,Flag_PrintLog);
% Plot_cruve(xr,br,cr,fs,N,M,TarS,Flag_PlotExport);
% % Output_data('./output_files/test.xlsx', TarS);
% sum(abs(cr - xr) < 0.1/sqrt(N))
rmpath('function/');
profile viewer