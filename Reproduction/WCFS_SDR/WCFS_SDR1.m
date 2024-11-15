clc; clear; close all;
addpath("function/");
N = 256;
M = 4;
threshold = 20;
fs = 128*1e6;
fc = 0e6;
B = 40e6;
T = 2e-6;
noise = exp(1j*2*pi*fc);

% num_Mator = 100;
% iter_waitbar = forWaitbar(num_Mator);
% test_Out = 0;
% parfor i = 1:num_Mator
%     x =  exp(1j*2*pi*rand(N,1));
%     % x = randn(N,1) + 1j*randn(N,1);
%     CS = CF_diag(real(x), N, M);
%     test_Out = test_Out + CS/num_Mator;
% end
% test_Out = abs(test_Out);
% test_a = diag(abs(fliplr(test_Out)));plot(test_a);


% % 1. 检测概率约束矩阵 假设噪声为复高斯分布，均值为0，方差为sigma
% sigma2 = 0.5;
% M = diag(ones(N,1)*(sigma2 + 1j*sigma2));
% R = inv(M);

% % 2. 相似度约束矩阵
c = generator_LFM(fs,fc,B,T); % 参考信号
epsilon = 0.2; % 相似性约束
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
s = [0.35 0.385;...
    0.55 0.575];
% s = [-1 1];
ws = [0.5;0.5];
k = size(s, 1);
Es = 1e-10; % 阻带能量约束定义
Rss = zeros(N);
Rss2 = zeros(N);
for ik = 1:k
    Rs = zeros(N);
    for m = 1:N
        for l = 1:N
            if m==l
                Rs(m,l) = s(ik,2) - s(ik,1);
            else
                Rs(m,l) = (exp(1j*2*pi*s(ik,2)*(m-l))-exp(1j*2*pi*s(ik,1)*(m-l)))/(1j*2*pi*(m-l));
            end
        end
    end
    if ik == 2
        Rss2 = Rs;
    elseif ik == 1
        Rss = Rs;
    end
end
Rss = Rss/N;
Rss2 = Rss2/N;
% % 4. 循环谱约束矩阵
S = GenerateMatrixS(N, 4, threshold);
S = (S + S')/2;
x = exp(1j*2*pi*rand(N,1));
x0 = x;
x = c;
X0 = x*x';

cvx_begin sdp
    cvx_solver SDPT3
        variable X(N,N) complex hermitian semidefinite
        minimize(real(trace(S*X)));
        subject to
            diag(X) == 1;
            (trace(R0*X)) >= delta;
            (trace(Rss*X)) <= Es;
            trace(Rss2*X) <= 1e-10;
            % X >= 0;
        % Temp{end+1} = X;
cvx_end
[eigenvector,D] = eig(X);
eigenvalue = diag(D);
eigenvalue_threshold = max(eigenvalue)*1e-4;
rank_X = sum(eigenvalue > eigenvalue_threshold);
[value,num] = max(eigenvalue);
x1 = sqrt(value)*eigenvector(:,num);
%%
% NumL = 1000; % 生成Num_L组随机向量
% randCN = zeros(NumL, size(X,1));
% randCN_out = zeros(NumL, 1);
% L = chol(X, 'lower');
% for i = 1:NumL
%     z = (randn(size(X, 1), 1) + 1i * randn(size(X, 1), 1)) / sqrt(2);
%     temp_randCN = L * z;
%     randCN(i,:) = temp_randCN;
%     randCN_out(i) = temp_randCN'*S*temp_randCN;
% end
% [~, index] = min(randCN_out);
% x1 = randCN(index, :);


Analysis_ES(x1,fs,N*2,'bool_draw', 1);title('Spectrum')
real(trace(Rss*X))
real(trace(Rss2*X))
% epsilonset = [0.1, 0.2, 0.3]; % 相似性约束
% Temp = cell(0,1);
% for i1 = 1:1
%     epsilon = epsilonset(i1);
%     delta = (1-epsilon/2)^2;
%     delta = N^2*delta;
% 
%     max_iter = 1000;
%     atest = zeros(max_iter,4);
%     X = x*x';
%     % real(trace(S*x*x'))
%     trace(Rss*X)
%     cvx_begin sdp
%         cvx_solver SDPT3
%             variable X(N,N) complex hermitian
%             minimize(real(trace(S*X)));
%             subject to
%                 diag(X) == 1;
%                 % trace(R0*X) >= delta;
%                 % trace(Rss*X) <= Es;
%                 trace(X) == N;
%                 X >= 0;
%             Temp{end+1} = X;
%     cvx_end
% 
%     trace(Rss*X)
%     [eigenvector,D] = eig(X);
%     eigenvalue = diag(D);
%     eigenvalue_threshold = max(eigenvalue)*1e-4;
%     rank_X = sum(eigenvalue > eigenvalue_threshold);
%     if rank_X == 1
%         [value,num] = max(eigenvalue);
%         x = sqrt(value)*eigenvector(:,num);
% 
%     elseif rank_X == 2
%         % % Slove the new problem
%         x = zeros(N,1);
%         for i_X = 1:rank_X
%             x = x + 0;
% 
%         end
%     elseif rank_X >= 3
%         % % Slove the new problem
% 
%     end
% 
% 
%     [Xs,r,rdel,rdelb]=fcndcmps12(R0,Rss,X);
% 
% 
%     out = CF_diag(x, N, M);
%     plot(diag(abs(fliplr(out)))); hold on;
% end
% figure
% out = CF_diag(c, N, M);
% plot(diag(abs(fliplr(out)))); hold on;
% hold off;
% legend(["c","xr(\delta=0.1)","xr(\delta=0.2)","xr(\delta=0.3)"]);
% SetDrawStyle
% figure
% Analysis_ES(x,fs,N,'bool_draw', 1);title('Spectrum')