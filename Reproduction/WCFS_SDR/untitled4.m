clc; clear; close all;
addpath("function/");
N = 128;
M = 4;
threshold = 10;
fs = 64*1e6;
fc = 0e6;
B = 30e6;
T = 2e-6;
Es = 1; % 阻带能量约束定义
c = generator_LFM(fs,fc,B,T);
epsilon = 0.9; % 相似性约束
delta = N^2*(1-epsilon/2)^2;
% 初始化矩阵S
S = GenerateMatrixS(N, M, threshold);
S = (S + S')/2;
% S = S/max(abs(S),[],'all');
% 初始化矩阵R0
R0 = c*c';
% 初始化矩阵Rs
x = exp(1j*2*pi*rand(N,1));
x0 = x;
x = c;
% FNr
FN = dftmtx(N); 
Fshift = zeros(N);
Fshift(1:N/2,N/2+1:end) = eye(N/2);
Fshift(N/2+1:end,1:N/2) = eye(N/2);
FN = Fshift*FN;

p = 0:N-M;
q = 0:1:N-M;
m = -M/2:M/2-1;
fftx = FN*c;
out = zeros(N-M+1,(M-N)/2+1);
for ip = p
    for iq = q
        q1 = ip + iq + m;
        q2 = ip - iq + m;
        t1 = (0 <=q1 & q1 <= N-1) & (0 <= q2 & q2 <= N-1);
        q1 = q1 + 1;
        q2 = q2 + 1;
        out(ip+1,iq+1) = sum((fftx(q2(t1))).*fftx(q1(t1)));
    end
end
out2 = fliplr(out);
mesh(abs([out2,out]));
out3 = abs([out2,out]);
test = out3(:,N+1);
test = test/max(test);
plot(test)
