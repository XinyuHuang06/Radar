% 本程序包含对所有Commons内signal_evaluate与signal_evaluate函数的测试
%% 信号生成类函数测试--signal_generator

% gnerator_LFM--LFM信号
T = 5.12*1e-6; %
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

% Evaluate_eps.m--测试信号Epilison值
seq_SNR = 0:5:50;
Evaluate_eps(signal_LFM,128,seq_SNR, 1000);

% sidelobe--信号旁瓣、PSL、ISL、主瓣宽度(3dB)分析

% Analysis_CS_DFSM--使用DFSM方法进行二阶循环谱分析

% CS_Analysis--二阶循环谱分析
% Analysis_CS_FAM(fs,signal_LFM,256);
% [signal_PRPC, t] = generator_PRPC(fs, 1, 512);
% Analysis_CS_FAM(fs,signal_PRPC,256);
% CWD分析

% Analysis_Sidelobe--旁瓣分析

% Evaluate_eps--epsilon分析
seq_SNR = 0:5:50;
DFTm = dftmtx(length(signal_LFM));
Evaluate_eps(real(signal_LFM), 16, seq_SNR);

%% 其它函数
% SetDrawStyle()--设置绘图风格
% HEX2RGB()--<'XXXXXX'>HEX格式转RGB<[12,34,56]>格式
% RGB2HEX()--RGB<[12,34,56]>格式转<'XXXXXX'>HEX格式
% Q()--
% Qinv()--
