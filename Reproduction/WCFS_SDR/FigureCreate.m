% Example: 
% :param : 
% :return : 
% detailed description: 
%------------------------------------------------------------------------------
% Created by: Xinyu Huang.
% On: 08/07/2024.
% Copyright (C) 2024 Xinyu Huang (learning_huang@163.com).
% All Rights Reserved.
% Unauthorized copying of this file, via any medium is strictly prohibited.
% Proprietary and confidential.
%------------------------------------------------------------------------------
addpath("function/")
fs = 128*1e6;
B = 20*1e6;
fc = 0;
T  = 2*1e-6;
N  = ceil(T*fs);
M  = 4;
t  = (0:1/fs:(N-1)/fs)';
f_ordinate = zeros(N,N);
alfa_ordinate = zeros(N,N);
for i_1 = 1:N
    i_2 = 1:N;
    f_ordinate(i_1,i_2) = (i_1 + i_2 -1)/(2*N) - 1/2;
    alfa_ordinate(i_1,i_2) = (i_1 - i_2 + N - 1)/N - 1;
end
% % 绘制典型信号的循环谱分布
%% 0. 随机信号
Num = 100; % 样本数
CS_RS = zeros(N,N);
parfor i = 1:Num
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
ylabel('归一化幅度');
title('循环谱剖面(f=0)');
SetDrawStyle;
%% 1. LFM
signal_LFM = generator_LFM(fs,fc,B,T);
CS_LFM = CF_diag(signal_LFM, N, M);
CS_LFM = abs(CS_LFM);
%% 2. BPSK
Nc = 4;
Rb = fs/Nc;
phase_seq = randi([0 1], N/Nc, 1);
[signal_BPSK, ~] = generator_PSK(fs, fc, Rb, phase_seq, 2);
CF_BPSK = CF_diag(signal_BPSK, N, M);
CF_BPSK = abs(CF_BPSK);

%% 3. QPSK
Rb = fs/4;
Nc = 4;
phase_seq = randi([0 3], N/Nc, 1);
[signal_QPSK, ~] = generator_PSK(fs, fc, Rb, phase_seq, 4);
CF_QPSK = CF_diag(signal_QPSK, N, M);
CF_QPSK = abs(CF_QPSK);






%% 绘制使用SDR方法抗循环谱分析的结果
N = 256;
M = 4;
threshold = 30;
fs = 128*1e6;
fc = 0e6;
B = 40e6;
T = 2e-6;

% % 2. 相似度约束矩阵
c = generator_LFM(fs,fc,B,T); % 参考信号
epsilon = 0.3; % 相似性约束
delta = (1-epsilon/2)^2;
delta = N^2*delta;
R0 = c*c';

% % 4. 循环谱约束矩阵
S = GenerateMatrixS(N, 4, threshold);
S = (S + S')/2;
x = exp(1j*2*pi*rand(N,1));
x0 = x;
x = c;
X0 = x*x';

%% SDP问题求解
fig1 = figure;
temp_pro = max(target);
target = target/temp_pro;
plot(alfa_ordinate(:,N/2),target); hold on;
CS = CF_diag(c, N, M);
CS = abs(CS);
CS = CS/max(max(CS));
profile1 = diag(abs(fliplr(CS)));
plot(alfa_ordinate(:,N/2),profile1); hold on;
epsilonset = [0.1, 0.2, 0.3]; % 相似性约束
outset = cell(length(epsilonset),1);

for i1 = 1:3
    epsilon = epsilonset(i1);   
    delta = (1-epsilon/2)^2;
    delta = N^2*delta;

    cvx_begin sdp quiet
        cvx_solver SDPT3
            variable X(N,N) complex hermitian
            minimize(real(trace(S*X)));
            subject to
                X >= 0;
                diag(X) == 1;
                (trace(R0*X)) >= delta;
    cvx_end
    % % % 向量求解
    [eigenvector,D] = eig(X);
    eigenvalue = diag(D);
    eigenvalue_threshold = max(eigenvalue)*1e-4;
    rank_X = sum(eigenvalue > eigenvalue_threshold);
    [value,num] = max(eigenvalue);
    x1 = sqrt(value)*eigenvector(:,num);
    outset{i1} = x1;
    CS = CF_diag(x1, N, M);
    CS = abs(CS);
    CS = CS/max(max(CS));
    profile1 = diag(abs(fliplr(CS)));
    plot(alfa_ordinate(:,N/2),profile1); hold on;
end
hold off;
legend(["Noise","cr(LFM)","xr(\delta=0.1)","xr(\delta=0.2)","xr(\delta=0.3)"]);
xlabel('\alpha');ylabel('Amplify');
SetDrawStyle;
exportgraphics(fig1,"./output/fig2.pdf");
% % Gaussian Randomzation


%% Plot Figure 1
figure
subplot(221)
mesh(f_ordinate, alfa_ordinate, CS_RS);
xlabel('f', 'Interpreter', 'tex');
ylabel('\alpha', 'Interpreter', 'tex');
view(0,90);
title('噪声信号循环谱特性');

subplot(222)
mesh(f_ordinate, alfa_ordinate, CS_LFM);
xlabel('f', 'Interpreter', 'tex');
ylabel('\alpha', 'Interpreter', 'tex');
view(0,90);
title('LFM信号循环谱特性');

subplot(223)
mesh(f_ordinate, alfa_ordinate, CF_BPSK);
xlabel('f', 'Interpreter', 'tex');
ylabel('\alpha', 'Interpreter', 'tex');
view(0,90);
title('BPSK信号循环谱特性');

subplot(224)
mesh(f_ordinate, alfa_ordinate, CF_QPSK);
xlabel('f', 'Interpreter', 'tex');
ylabel('\alpha', 'Interpreter', 'tex');
view(0,90);
title('QPSK信号循环谱特性');



%% 旁瓣特性
Analysis_Sidelobe(c, c, 'bool_draw', 1);hold on;
for i1 = 1:length(epsilonset)
    x1 = outset{i1};
    Analysis_Sidelobe(x1, x1, 'bool_draw', 1);hold on;
end
legend(["参考信号(LFM)","xr(\delta=0.1)","xr(\delta=0.2)","xr(\delta=0.3)"]);
SetDrawStyle;
ylim([-60 0]);


