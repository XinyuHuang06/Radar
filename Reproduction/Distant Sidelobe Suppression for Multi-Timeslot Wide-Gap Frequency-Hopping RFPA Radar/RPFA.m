clc; clear; close all;
% 导入函数定义
addpath(".\function");

% 参数初始化
fC = 3*1e9; % 中心频率 3 GHz
M = 1024; % 脉冲数
TW = 1*1e-6; % 脉冲宽度 1 us
TR = 10*1e-6; % 平均PRI间隔 10 us
TJ = 6*1e-6; % 最大PRI捷变范围 6 us
PRI_PointsNum = 192; % PRI 捷变点数
% PRI捷变分布 均匀分布
B = 64*1e6; % 载频捷变范围 64 MHz
Fre_PointsNum = 2048; % 载频捷变点数
% 载频捷变分布 高斯分布
N = 3; % 相邻脉冲数
G = 10*1e6; % 最小频率间隔 10 MHz
fS = 128*1e6; % 采样频率 128 MHz
phi0 = 0; % 初始相位
TS = 1/(fS); % 采样间隔 且有 TS <= 1/(2*B)
% 探测相关参数
Delta_r = 10; % 距离探测精度
Delta_v = 10; % 速度探测精度
max_l = 100; % 最大距离探测单元
max_k = 100; % 最大速度探测单元
fpass = 1*1e6; % 低通滤波器通带贷款
% 其它参数
c = 3*1e8; % 光速

fnSeq = generate_fnSeq(fC, B, G, N, M); % 生成载波频率序列
tnSeq = generate_tnSeq(M,TJ); % 生成波形起始时间偏移量

% 目标信息
Target_num = 4;
TargetMatrix = zeros(3,Target_num); % 三个维度指：所处距离单元，所处速度单元，散射系数
TargetMatrix(:,1) = [10;10;0.01]; % 第一个目标，位于第(10,10)个单元上

% 绘制目标回波
n = 100; % 当前第n个信号的回波
m = 100; % 信号编号 m
p = 1:ceil(TW/TS)*2; % 回波采样点编号
value = smlk(m,TargetMatrix(1,1),TargetMatrix(2,1),n,p,TS,Delta_r,Delta_v,c,phi0,TW,tnSeq,fnSeq);
value2 = value.*exp(1j*(2*pi*10*1e6*(1:length(value))*TS)+phi0);
plot(real(value));hold on; % 因为信号
plot(real(value2));hold off;
figure
Analysis_ES((value),fS,length(value));
% 频移后回波

hold on;
n = 100; % 当前第n个信号的回波
m = 30; % 信号编号 m

value = smlk(m,TargetMatrix(1,1),TargetMatrix(2,1),n,p,TS,Delta_r,Delta_v,c,phi0,TW,tnSeq,fnSeq);
Analysis_ES((value2),fS,length(value));
hold off 

legend(["s_{m,l,k}(n,p)  m=n","s_{m,l,k}(n,p)  m\neq n"]);
% plot(real(value))
hp = h();
% 绘制滤波器的时域响应与频域响应
figure
subplot(121);
plot(hp);
subplot(122);
[aaa,aaaw] = freqz(hp, 1, length(hp));
plot(abs(aaa))
hat_smlk = conv(value,hp);
hat_snlk = conv(value2,hp);
% hat_smlk = hat_smlk(1:end-length(hp));
% hat_snlk = hat_snlk(1:end-length(hp));
% hat_smlk = filter(hp, 1, value);
% hat_snlk = filter(hp, 1, value2);
figure 
Analysis_ES(hat_smlk, fS,length(value));hold on;
Analysis_ES(hat_snlk, fS,length(value));hold off;
legend(["s_{m,l,k}(n,p)  m=n","s_{m,l,k}(n,p)  m\neq n"]);

% 相位编码信号
CR = 128*1e6; % 码元传输速率 128 MHz
TC = 1/CR; % 单码元传输间隔
MW = TW/TC; % 单脉冲码元传输个数



% 移除函数定义
rmpath(".\function");