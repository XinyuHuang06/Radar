% 生成导引矢量
clc;clear;close all;
N = 128; % 相参处理脉冲数目
M = 8; % 频率间隔数目
rng(1);
dn = randi(M,[1,N])-1; % 随机生成频率序列
fc = 10*1e9; % 载波频率
Tr = 1*1e-6; % 脉冲重复周期
c  = 3*1e8; % 光速
delta_f = 8*1e6; % 子频带带宽
P = M; % 高精度距离分辨单元数
Q = N; % 高精度速度分辨单元数

delta_R = c/(2*M*delta_f); % 高精度距离分辨力
delta_v = c/(2*fc*N*Tr); % 高精度速度分辨力
target_coef_dis = zeros(P*Q,1); % 高分辨单元对应距离相位因子
target_coef_veo = zeros(P*Q,1); % 高分辨单元对应速度相位因子
for i = 1:Q
    for j = 1:P
        % 子单元距离与速度
        R = (j-1)*delta_R+0.5*delta_R; 
        v = (i-1)*delta_v+0.5*delta_v;
        target_coef_dis( (i-1)*P+(j-1)+1 ) = -4*pi*delta_f*R/c;
        target_coef_veo( (i-1)*P+(j-1)+1 ) = -4*pi*fc*v*Tr/c;
    end
end

target = zeros(P*Q, 1); % 高精度分辨单元反射系数矩阵
target_index_seq = [10,105,180,341,456,542,786,890,963];
for i = 1:length(target_index_seq)
    target(target_index_seq(i)) = 1;
end
target_coef_RCS = zeros(P*Q,1);
for i = 1:Q
    for j = 1:P
        target_coef_RCS( (i-1)*P+(j-1)+1 ) = target( (i-1)*P+(j-1)+1 );
    end
end

ones_seq = 1 + dn.*delta_f/fc;
phi_matrix = zeros(N,P*Q); % 观测矩阵生成

for l = 1:P*Q
    phi_matrix(:,l) = exp(1j*dn*target_coef_dis(l) + 1j*target_coef_veo(l).*ones_seq.*(0:1:N-1));
end

% 由目标信息生成对应回波
theory_y = phi_matrix*target;
% 求解x
% % 利用匹配滤波求解x
x_1 = phi_matrix'*theory_y;
% % 对角加载算法
lambda = 1;
x_2 = (phi_matrix'*phi_matrix + lambda*eye(P*Q))\phi_matrix'*theory_y;
% % % l2范数最小化
x_3 = phi_matrix'*inv(phi_matrix*phi_matrix')*theory_y;
% % % l1范数最小化--基追踪算法
% 使用 cvx 求解基追踪问题
cvx_begin
    variable x_4(P*Q)
    minimize(norm(x_4, 1))
    subject to
    theory_y == phi_matrix*x_4;
cvx_end
target_index_seq = target_index_seq + 0.5;
figure
subplot(221);
plot(abs(my_normalize_complex(x_1)));hold on;
plot(target_index_seq, ones(size(target_index_seq)), 'o');hold off;
xlabel('距离-速度分辨单元');ylabel('目标反射系数');
legend(['匹配滤波',"真值"]);
subplot(222);
plot(abs(my_normalize_complex(x_2)));hold on;
plot(target_index_seq, ones(size(target_index_seq)), 'o');hold off;
xlabel('距离-速度分辨单元');ylabel('目标反射系数');
legend(['对角加载',"真值"]);
subplot(223);
plot(abs(my_normalize_complex(x_3)));hold on;
plot(target_index_seq, ones(size(target_index_seq)), 'o');hold off;
xlabel('距离-速度分辨单元');ylabel('目标反射系数');
legend(['l_2范数最小化',"真值"]);
subplot(224);
plot(abs(my_normalize_complex(x_4)));hold on;
plot(target_index_seq, ones(size(target_index_seq)), 'o');hold off;
xlabel('距离-速度分辨单元');ylabel('目标反射系数');
legend(['l_1范数最小化',"真值"]);
SetDrawStyle();
exportgraphics(gcf,'output.jpg');
% % 有噪情况下各算法恢复情况

% % 转化为二维图像，显示无噪声情况下恢复效果

% 由目标信息生成对应回波
SNR = 5; % 信噪比
noise_coef = 1/10^(SNR/10); % 考虑恒定功率发射波形
theory_y = phi_matrix*target + noise_coef*(1j*randn(N,1)+randn(N,1));
% 求解x
% % 利用匹配滤波求解x
x_1 = phi_matrix'*theory_y;
% % 对角加载算法
lambda = 1;
x_2 = (phi_matrix'*phi_matrix + lambda*eye(P*Q))\phi_matrix'*theory_y;
% % % l2范数最小化
x_3 = phi_matrix'*inv(phi_matrix*phi_matrix')*theory_y;
% % % l1范数最小化--lasso方法
% 使用 cvx 求解基lasso问题
cvx_begin
    variable x_4(P*Q)
    minimize(norm(theory_y - phi_matrix*x_4, 2))
    subject to
    norm(x_4, 1) < length(target_index_seq)+2;
cvx_end

figure
subplot(221);
plot(abs(my_normalize_complex(x_1)));hold on;
plot(target_index_seq, ones(size(target_index_seq)), 'o');hold off;
xlabel('距离-速度分辨单元');ylabel('目标反射系数');
legend(['匹配滤波',"真值"]);
subplot(222);
plot(abs(my_normalize_complex(x_2)));hold on;
plot(target_index_seq, ones(size(target_index_seq)), 'o');hold off;
xlabel('距离-速度分辨单元');ylabel('目标反射系数');
legend(['对角加载',"真值"]);
subplot(223);
plot(abs(my_normalize_complex(x_3)));hold on;
plot(target_index_seq, ones(size(target_index_seq)), 'o');hold off;
xlabel('距离-速度分辨单元');ylabel('目标反射系数');
legend(['l_2范数最小化',"真值"]);
subplot(224);
plot(abs(my_normalize_complex(x_4)));hold on;
plot(target_index_seq, ones(size(target_index_seq)), 'o');hold off;
xlabel('距离-速度分辨单元');ylabel('目标反射系数');
legend(['l_1范数最小化',"真值"]);
SetDrawStyle();
exportgraphics(gcf,'output_noise.jpg');
% % 考虑模型失配情况下各个算法恢复情况

function y = my_normalize_complex(x)
    y = x/max(abs(x));
end

function out = inv_seq(seq, M, N)
    out = 0;
end

