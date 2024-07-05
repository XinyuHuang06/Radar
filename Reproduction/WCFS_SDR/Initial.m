% Example: 
% :param : 
% :return : 
% detailed description: 
%------------------------------------------------------------------------------
% Created by: Xinyu Huang.
% On: 19/06/2024.
% Copyright (C) 2024 Xinyu Huang (learning_huang@163.com).
% All Rights Reserved.
% Unauthorized copying of this file, via any medium is strictly prohibited.
% Proprietary and confidential.
%------------------------------------------------------------------------------
% 初始化
JsonPath = "data/initial_Parameter.json";
addpath("function/");
initial_parameter = loadjson(JsonPath);
T = initial_parameter.signal.T; % 信号时宽
CV = initial_parameter.signal.CV; % 码元传输速率
fs = initial_parameter.signal.fs; % 采样速率
fc = initial_parameter.signal.fc; % 载波速率
threshold = initial_parameter.Threshold;
NPC = T*CV; % 码元数目
N = T*fs; % 样本数
M = fs/CV; % 码元所占样本数
% 初始化码元序列
x_PC = 2*pi*rand(NPC,1); % 随机初始化
index = kron((1:NPC)',ones(M,1));
t = (0:1/fs:(N-1)/fs)';
x = exp(1j*2*pi*fc*t).*exp(1j*x_PC(index)); % 初始化序列
S = GenerateMatrixS(N, 4, threshold);

Max_iter = 100;

for i = 1:Max_iter
    % Coordinate-Descent
    for iter_CD = 1:N
        x_n = x(iter_CD);
        x_Non_n = x; x_Non_n(iter_CD) = 0;
        e_n = zeros(N,1); e_n(iter_CD) = 1;
        x_n_k2 = exp(1j*(pi/2-angle(x_Non_n'*S*e_n)));
        x(iter_CD) = x_n_k2;
    end
end
% 直接忽略掉秩一约束 然后使用CVX对其进行求解，其次，再使用相应



% 