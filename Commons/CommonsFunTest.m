% 本程序包含对所有Commons内signal_evaluate与signal_evaluate函数的测试
clear
%% 信号生成类函数测试--signal_evaluate
% gnerator_LFM--LFM信号
% 
T = 5*1e-6;
B = 8*1e6;
fs = 20*1e6;
fc = 0;
[signal_LFM, t] = generator_LFM(fs,fc,B,T); % 信号生成

% OFDM_generator--OFDM信号
%
% % generator_PSK--PSK信号
% %
% fs = 4*1e6;
% fc = 1*1e5;
% Rb = 1*1e5;
% code_seq = round(rand(1,64));
% [signal_BPSK, t] = generator_PSK(fs,fc,Rb,code_seq,2,'noiseF',1,'SNR',10);

%% 信号评估类函数测试--signal-evaluate
% Analysis_AF--模糊函数分析

% Analysis_CS_DFSM--使用DFSM方法进行二阶循环谱分析

% Analysis_CS_FAM--使用FAM方法进行二阶循环谱分析

% Analysis_Fs--频谱分析

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