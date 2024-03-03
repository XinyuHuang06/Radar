% 本程序包含对所有Commons内signal_evaluate与signal_evaluate函数的测试
%% 信号生成类函数测试--signal_generator

% gnerator_LFM--LFM信号
T = 10.24*1e-6; %
B = 40*1e6; % 40 MHz
fs = 100*1e6; % 100 MHz
fc = 10*1e6; % 10MHz
[signal_LFM, t] = generator_LFM(fs,fc,B,T,'type',0); % 信号生成
% NLFM_generator.m--NLFM信号

% OFDM_generator.m--OFDM信号

% % PSL_generator.m--PSK信号
fs = 4*1e6;
fc = 1*1e5;
Rb = 1*1e5;
code_seq = round(rand(1,64));
[signal_BPSK,t] = generator_PSK(fs,fc,Rb,code_seq,2,'noiseF',1,'SNR',10);

%% 信号评估类函数测试--signal_evaluate

% sidelobe--信号旁瓣、PSL、ISL、主瓣宽度(3dB)分析

% AF_Analysis--模糊函数分析

% CS_Analysis--二阶循环谱分析

% CWD分析

% FS_Analysis--频谱分析

% 

%% 其它函数

% SetDrawStyle.m--设置绘图风格
% Q.m