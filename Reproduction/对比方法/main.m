% Example: 
% :param : 
% :return : 
% detailed description: 
%------------------------------------------------------------------------------
% Created by: Xinyu Huang.
% On: 19/09/2024.
% Copyright (C) 2024 Xinyu Huang (learning_huang@163.com).
% All Rights Reserved.
% Unauthorized copying of this file, via any medium is strictly prohibited.
% Proprietary and confidential.
%------------------------------------------------------------------------------


% Pre
addpath("../WCFS_SDR/function")
addpath("../LPI Radar Waveform Design With Desired Cyclic Spectrum and Pulse Compression Properties/function")
addpath("../LPI Radar Waveform Design With Desired Cyclic Spectrum and Pulse Compression Properties")
%% addpath为临时添加路径函数，在每次程序释放时将被清除
fs = 128*1e6;
B = 20*1e6;
fc = 0;
T  = 2*1e-6;
N  = ceil(T*fs);
M  = 4;
t  = (0:1/fs:(N-1)/fs)';
f_ordinate = zeros(N,N);
threshold = 30;
alfa_ordinate = zeros(N,N);
for i_1 = 1:N
    i_2 = 1:N;
    f_ordinate(i_1,i_2) = (i_1 + i_2 -1)/(2*N) - 1/2;
    alfa_ordinate(i_1,i_2) = (i_1 - i_2 + N - 1)/N - 1;
end
% % 绘制典型信号的循环谱分布
%% 0. 随机信号
Num = 10; % 样本数
CS_RS = zeros(N,N);
for i = 1:Num
    signal_random = exp(1j*2*pi*rand(N,1));
    CS_RS = CS_RS + CF_diag(signal_random, N, M)/Num;
end
CS_RS = abs(CS_RS);
normalize_value = max(max(CS_RS));
CS_RS = CS_RS/max(max(CS_RS));
target = diag(abs(fliplr(CS_RS)));
target(N/2) = (CS_RS(N/2,N/2) + CS_RS(N/2+1,N/2+1))/2 ;
plot(alfa_ordinate(:,N/2),target);
xlabel('\alpha');
ylabel('Normailized');
title('The profile(f=0)');
SetDrawStyle;
%% 1. LFM
signal_LFM = generator_LFM(fs,fc,B,T);
CS_LFM = CF_diag(signal_LFM, N, M);
CS_LFM = abs(CS_LFM);
%% 2. BPSK
% Nc = 4;
% Rb = fs/Nc;
% phase_seq = randi([0 1], N/Nc, 1);
% [signal_BPSK, ~] = generator_PSK(fs, fc, Rb, phase_seq, 2);
% CF_BPSK = CF_diag(signal_BPSK, N, M);
% CF_BPSK = abs(CF_BPSK);
% %% 3. QPSK
% Rb = fs/4;
% Nc = 4;
% phase_seq = randi([0 3], N/Nc, 1);
% [signal_QPSK, ~] = generator_PSK(fs, fc, Rb, phase_seq, 4);
% CF_QPSK = CF_diag(signal_QPSK, N, M);
% CF_QPSK = abs(CF_QPSK);
% ADMM
% % 生成数据文件
bool_generatemat = true;
if(bool_generatemat)
    
else

end
epsilonset = [0.1,0.2,0.3];
SDR_set = cell(3,1);
ADMM_set = cell(3,1);
DataPath = ["output/parameter_01.mat" "output/parameter_02.mat" "output/parameter_03.mat"];
OutPath = ["output/ADMM01" "output/ADMM02" "output/ADMM03"];
parfor i = 1:3
    epsilon = epsilonset(i);
    SDR_set{i} = WCFS_SDR(N,threshold,signal_LFM,epsilon);
    WCFS_ADMM(DataPath(i), OutPath(i));
end


% SDR

