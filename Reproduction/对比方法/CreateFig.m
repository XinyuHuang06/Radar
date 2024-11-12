% % 对比图 1
% 读取波形数据
clc;clear;close all;
fs = 128*1e6;
B = 40*1e6;
fc = 0;
T  = 2*1e-6;
signal_LFM = generator_LFM(fs,fc,B,T);
addpath("../WCFS_SDR/function")
addpath("../LPI Radar Waveform Design With Desired Cyclic Spectrum and Pulse Compression Properties/function")
addpath("../LPI Radar Waveform Design With Desired Cyclic Spectrum and Pulse Compression Properties")
load("./out_set.mat");
SDR_set = outset;
N = 256;
M = 4;
% 
ADMM_Set = cell(3,1);
for i=1:3
	str = ['./output/ADMM0',num2str(i),'/Record.mat'];
	temp_set = load(str);
	xr = temp_set.DataRecordPack.ParaRecord.xr{end};
	xr = xr(1:N) + 1j*xr(N+1:end);
	ADMM_Set{i} = xr;
end
% % 绘制对比图
% Calculate the maxvalue
CS_SDR_set = cell(3,1);
CS_ADMM_set = cell(3,1);
parfor i = 1:3
	temp_CS1 = CF_diag(ADMM_Set{i}, N, M);
    temp_CS1 = abs(temp_CS1);
    CS_ADMM_set{i} = temp_CS1;
	temp_CS2 = CF_diag(SDR_set{i}, N, M);
    temp_CS2 = abs(temp_CS2);
	CS_SDR_set{i} = temp_CS2;
end
%%
CS_Stand = CF_diag(signal_LFM, N, M);

for i=1:3
    figure
    % plot(real(ADMM_Set{i}));hold on;
    plot(real(SDR_set{i}));hold on;
    plot(real(signal_LFM));
end
%%
for i=1:3
    figure
    temp_CS = CS_SDR_set{i};
    target = diag(abs(fliplr(temp_CS)));
    plot(target);hold on;
    temp_CS = CS_ADMM_set{i};
    target = diag(abs(fliplr(temp_CS)));
    plot(target);hold on;
    temp_CS = CS_Stand;
    target = diag(abs(fliplr(temp_CS)));
    plot(target);hold on;
    legend('SDR','ADMM',"LFM")
end

%% Test the moduls
figure
for i = 1:3
    ADMM_xr = ADMM_Set{i};
    SDR_xr = SDR_set{i};
    plot(abs(ADMM_xr));hold on;
    plot(abs(SDR_xr));hold on;
end
legend('ADMM01','SDR01','ADMM02','SDR02','ADMM03','SDR03');