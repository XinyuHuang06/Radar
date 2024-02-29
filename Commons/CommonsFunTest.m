% 本程序包含对所有Commons内signal_evaluate与signal_evaluate函数的测试
%% 信号生成类函数测试--signal_evaluate
% gnerator_LFM--LFM信号
T = 5*1e-7;
B = 8*1e6;
fs = 20*1e6;
fc = 0;
[signal_LFM, t] = generator_LFM(fs,fc,B,T,'type',1,'ratio',0.5,'num_pulse',3); % 信号生成
% NLFM_generator--NLFM信号

% OFDM_generator--OFDM信号

% % PSL_generator--PSK信号
fs = 4*1e6;
fc = 1*1e5;
Rb = 1*1e5;
code_seq = round(rand(1,64));
[signal_BPSK,t] = generator_PSK(fs,fc,Rb,code_seq,2,'noiseF',1,'SNR',10);

%% 信号评估类函数测试--signal-evaluate
% sidelobe--信号旁瓣、PSL、ISL、主瓣宽度(3dB)分析

% AF_Analysis--模糊函数分析

% 低截获性能评估
% CS_Analysis--二阶循环谱分析

% CWD分析

% FS_Analysis--频谱分析

% 

%% 其它函数
% SetDrawStyle--设置