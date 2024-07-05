clc; clear; close all;
addpath("function/");
N = 128;
M = 4;
threshold = 10;
fs = 128*1e6;
fc = 0e6;
B = 40e6;
T = 1e-6;
noise = exp(1j*2*pi*fc);

num_Mator = 10000;
iter_waitbar = forWaitbar(num_Mator);
test_Out = 0;
parfor i = 1:num_Mator
    x =  exp(1j*2*pi*rand(N,1));
    % x = randn(N,1) + 1j*randn(N,1);
    CS = CF_diag(real(x), N, M);
    test_Out = test_Out + CS/num_Mator;
end
test_Out = abs(test_Out);
test_a = diag(abs(fliplr(test_Out)));plot(test_a);
% % 1. 检测概率约束矩阵 假设噪声为复高斯分布，均值为0，方差为sigma
% sigma2 = 0.5;
% M = diag(ones(N,1)*(sigma2 + 1j*sigma2));
% R = inv(M);

% % 2. 相似度约束矩阵
c = generator_LFM(fs,fc,B,T); % 参考信号
epsilon = 0.3; % 相似性约束
delta = (1-epsilon/2)^2;
delta = N^2*delta;
R0 = c*c';

% % 3. 阻带约束矩阵
% 阻带约束定义
% s = [0 0.0617;...
%     0.0988 0.2469;...
%     0.2593,0.2840;...
%     0.3086, 0.3827;...
%     0.4074, 0.4938;...
%     0.5185, 0.5556;...
%     0.9383, 1.0000];
% s = s-1/2;
% ws = [0.1 0.2 0.3 0.4 0.5 0.6 0.7];
% ws = ws/sum(ws);
s = [0.4 0.8];
ws = 1;
k = size(s, 1);
Es = 0.2; % 阻带能量约束定义
Rss = zeros(N);
for ik = 1:k
    Rs = zeros(N);
    for m = 1:N
        for l = 1:N
            if m==l
                Rs(m,l) = s(k,2) - s(k,1);
            else
                Rs(m,l) = (exp(1j*2*pi*s(k,2)*(m-l))-exp(1j*2*pi*s(k,1)*(m-l)))/(1j*2*pi*(m-l));
            end
        end
    end
    Rss = Rss + ws(ik)*Rs;
end
Rss = Rss/N;
% % 4. 循环谱约束矩阵
S = GenerateMatrixS(N, 4, threshold);
S = (S + S')/2;
x = exp(1j*2*pi*rand(N,1));
x0 = x;
x = c;

figure
epsilonset = [0.1, 0.2, 0.3]; % 相似性约束
out = CF_diag(c, N, M);
plot(diag(abs(fliplr(out)))); hold on;
Temp = cell(0,1);
for i1 = 1:3
    epsilon = epsilonset(i1);
    delta = (1-epsilon/2)^2;
    delta = N^2*delta;

    max_iter = 1000;
    atest = zeros(max_iter,4);
    X = x*x';
    % real(trace(S*x*x'))
    trace(Rss*X)
    cvx_begin sdp
        cvx_solver SDPT3
            variable X(N,N) complex hermitian
            minimize(real(trace(S*X)));
            subject to
                diag(X) == 1;
                trace(R0*X) >= delta;
                trace(Rss*X) <= Es;
                X >= 0;
            Temp{end+1} = X;
    cvx_end
    trace(Rss*X)
    [V,D] = eig(X);
    [value,num] = max(diag(D));
    x = sqrt(value)*V(:,num);
    out = CF_diag(x, N, M);
    plot(diag(abs(fliplr(out)))); hold on;
end
hold off;
legend(["c","xr(\delta=0.1)","xr(\delta=0.2)","xr(\delta=0.3)"]);
SetDrawStyle