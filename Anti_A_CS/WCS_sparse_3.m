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
clc;clear;close all;
% add function definition
addpath("function/");
% % Parameter Setting
flag_PlotAndExport = 1;
flag_Sparse = 1;
flag_PrintLog = 0;
flag_Tars = 1;
num_Mator = 1;
maxnum = 50; % 最大迭代次数
threshhold = 20;
% % Data IntializationS
N = 128;
M = 8;
fs = 2e6;
T = N/fs;
fc = 1e5;
temp_B = 4e5;
c = generator_LFM(fs,fc,temp_B,T);
cr = [real(c);imag(c)];
[chi_matrix, Taf_1, Taf_2, FNr, omega_alpha_1] = Generate_data(N, M, flag_Sparse, threshhold);
x = exp(1j*2*pi*rand(N,1));
xr = [real(x);imag(x)]; 
% br = xr;
% xr = cr;
b = exp(1j*2*pi*rand(N,1))/sqrt(N);
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
xi_inc = 1.5;
xi_dec = 1.2;
vartheta = (0.2/sqrt(N)); % 参考信号相关度阈值
% % Data Record
% % The initial parameter Caculation
xr_init = xr;
TarS = struct('sum',zeros(maxnum+1,1),'Tar',zeros(maxnum+1,1),'Lar_1',zeros(maxnum+1,1),'pen_1',zeros(maxnum+1,1),...
                    'Lar_2',zeros(maxnum+1,1),'pen_2',zeros(maxnum+1,1));
TarS_2 = struct('sum',zeros(maxnum+1,1),'Tar',zeros(maxnum+1,1),'Lar_1',zeros(maxnum+1,1),'pen_1',zeros(maxnum+1,1),...
    'Lar_2',zeros(maxnum+1,1),'pen_2',zeros(maxnum+1,1));
[TarS.sum(1),TarS.Tar(1), TarS.Lar_1(1), TarS.pen_1(1), TarS.Lar_2(1), TarS.pen_2(1)] = ...
    Debug_fun(N, FNr, xr, br, cr, omega_alpha_1, rho_0, rho_1, lambda_0, lambda_1, h, vartheta, Taf_1, Taf_2, chi_matrix, flag_Sparse);
forwaitbar = forWaitbar(maxnum);
parameter_record = struct('rho_0',cell(maxnum+1,1));
% parameter_record.rho_0{1} = rho_0;
% % ADMM
% fig1 = figure;
for i_m = 1:maxnum
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
    [A1,BT1,~] = Caculate_matrix_xr(N, FNr,xr, br, cr, omega_alpha_1, rho_0, rho_1, lambda_0, lambda_1, h, vartheta, Taf_1, Taf_2, chi_matrix, flag_Sparse);
    [xr] = Update_xr(A1,BT1,xr);
    % % Step 2 , Solving the b_r
    [br] = Update_br(N, FNr, xr, br, cr, omega_alpha_1, rho_0, rho_1, lambda_0, lambda_1, h, vartheta, Taf_1, Taf_2, chi_matrix, flag_Sparse);
    % Step 3 , Solving the h
    h = Update_h(N, xr, br, cr, rho_1, lambda_1, vartheta, r, chi_matrix, flag_Sparse);
    % Step 4 , Solving the u, v_n
    [lambda_0, lambda_1] = Update_lambda(N, xr, br, cr, rho_0, rho_1, lambda_0, lambda_1, h, vartheta, chi_matrix, flag_Sparse);
    % Step 5 , Update rho_0 rho_1
    [rho_0,rho_1,delta_0,delta_1] = Update_rho(rho_0_k,rho_1_k,xr_k,br_k,h_k,xr,br,cr,vartheta,h,chi_matrix,flag_Sparse,xi_inc,xi_dec);
    % 参数记录
    forwaitbar.show_bar;
    % [A1,BT1,Tar_seq] = Caculate_matrix_xr(N, FNr,xr, br, cr, omega_alpha_1, rho_0, rho_1, lambda_0, lambda_1, h, vartheta, Taf_1, Taf_2, chi_matrix, flag_Sparse);
    % TarS_2.sum(i_m) = Tar_seq(1);TarS_2.Tar(i_m) = Tar_seq(2);TarS_2.Lar_1(i_m) = Tar_seq(3);TarS_2.pen_1(i_m) = Tar_seq(4);TarS_2.Lar_2(i_m) = Tar_seq(5);TarS_2.pen_2(i_m) = Tar_seq(6);
    [TarS.sum(i_m+1),TarS.Tar(i_m+1), TarS.Lar_1(i_m+1), TarS.pen_1(i_m+1), TarS.Lar_2(i_m+1), TarS.pen_2(i_m+1)] = ...
        Debug_fun(N, FNr, xr, br, cr, omega_alpha_1, rho_0, rho_1, lambda_0, lambda_1, h, vartheta, Taf_1, Taf_2, chi_matrix, flag_Sparse);
    % parameter_record.rho_0{i_m+1} = rho_0;
end
clear *_k
s_xr_init = RadarSignal(xr_init,'xr_init');s_xr_post = RadarSignal(xr,'xr_post');s_cr = RadarSignal(cr,'cr');
print_log(TarS,s_xr_init,s_xr_post,s_cr,flag_PrintLog);
Plot_cruve(xr,br,cr,fs,N,M,TarS,flag_PlotAndExport);
Output_data('./output_files/test.xlsx', TarS);
rmpath('function/');