% Example: 
% :param : 
% :return : 
% detailed description: 
%------------------------------------------------------------------------------
% Created by: Xinyu Huang.
% On: 13/09/2024.
% Copyright (C) 2024 Xinyu Huang (learning_huang@163.com).
% All Rights Reserved.
% Unauthorized copying of this file, via any medium is strictly prohibited.
% Proprietary and confidential.
%------------------------------------------------------------------------------
 
% MTD Detection
% The Parameter of target

% 参数设置
numTargets = 5; % 目标数量
rangeMax = 1000; % 最大距离（米）
velocityMax = 50; % 最大速度（米/秒）
numPulses = 64; % 脉冲数量
fs = 1e6; % 采样频率（Hz）
fc = 10e9; % 载波频率（Hz）
bw = 1e6; % 带宽（Hz）
pulseWidth = 1e-6; % 脉冲宽度（秒）

% 生成LFM信号
t = 0:1/fs:pulseWidth-1/fs;
lfmPulse = chirp(t, fc-bw/2, pulseWidth, fc+bw/2);

% 随机初始化目标位置和速度
targetRanges = rangeMax * rand(1, numTargets);
targetVelocities = velocityMax * (2*rand(1, numTargets) - 1);

% 模拟回波信号
echoSignal = zeros(numPulses, length(lfmPulse));
for i = 1:numTargets
delay = 2 * targetRanges(i) / 3e8; % 距离引起的延迟
dopplerShift = 2 * targetVelocities(i) * fc / 3e8; % 多普勒频移
delayedPulse = circshift(lfmPulse, round(delay * fs));
for j = 1:numPulses
echoSignal(j, :) = echoSignal(j, :) + delayedPulse .* exp(1j * 2 * pi * dopplerShift * (j-1) * pulseWidth);
end
end

% 添加噪声
noise = 0.1 * (randn(size(echoSignal)) + 1j * randn(size(echoSignal)));
receivedSignal = echoSignal + noise;

% 进行匹配滤波
matchedFilter = conj(fliplr(lfmPulse));
rangeDopplerMap = zeros(numPulses, length(lfmPulse));
for i = 1:numPulses
pulseSegment = receivedSignal(i, :);
rangeDopplerMap(i, :) = abs(ifft(fft(pulseSegment) .* fft(matchedFilter, length(pulseSegment))));
end

% 绘图显示距离-速度二维图像
figure;
imagesc(20*log10(abs(rangeDopplerMap)));
xlabel('Range Bins');
ylabel('Doppler Bins');
title('Range-Doppler Map');
colorbar;