% Example:
% :param :
% :return :
% detailed description: 并行化优化低截获雷达波形程序--（低维度优化/内存充裕时运行）
%------------------------------------------------------------------------------
% Created by: Xinyu Huang.
% On: 26/03/2024.
% Copyright (C) 2024 Xinyu Huang (learning_huang@163.com).
% All Rights Reserved.
% Unauthorized copying of this file, via any medium is strictly prohibited.
% Proprietary and confidential.
%------------------------------------------------------------------------------

% % 并行功能初始化
% NumCores =  feature('numcores'); % 获取当前CPU内核数
% delete(gcp('nocreate')); % 关闭已有的并行化资源池
% parpool(NumCores); % 开启并行化计算资源池

% % 数据初始化
N = 64; % 信号长度
M = 16; % 窗长
x = exp(1j*2*pi*rand(N,1))/sqrt(N); % 目标信号，随机初始化
xr = [real(x);imag(x)]; 
br = xr; % 辅助变量
FN = dftmtx(N); % 离散傅里叶变换矩阵
FNr = complex2real(FN);
% % cr 参考信号生成
T = N*1e-6;
fs = 1e6;
fc = 1e5;
B = 2e5;
c = generator_LFM(fs,fc,B,T)/sqrt(N);
cr = [real(c);imag(c)];
vartheta = 0.1; % 参考信号相关度阈值
h = -0.1; % 相似度约束辅助变量
chi = zeros(2*N,2*N,N); % 相似性约束矩阵
for i_temp = 1:N
    matrix_temp = zeros(N,N);
    matrix_temp(i_temp,i_temp) = 1;
    chi(:,:,i_temp) = complex2real(matrix_temp);
end
omega_alpha = zeros((N-M)/2,1); % 循环谱加权向量
omega_alpha(1:40) = 0.01;
omega_alpha(41:end) = 10;
omega_alpha_1 = ones((N-M)*(N-M),1);
% 循环谱计算矩阵 
alpha = 1:N-M;
cf = 1:N-M;
Taf_1 = zeros(2*N,2*N,(N-M)*(N-M));
Taf_2 = zeros(2*N,2*N,(N-M)*(N-M));
for i_temp = 1:(N-M)*(N-M)
    temp_alpha = fix(i_temp/(N-M))+1;
    temp_cf = mod(i_temp,N-M)+1;
    temp_taf = zeros(N,N);
    temp_taf(temp_cf:temp_cf+M-1, temp_alpha:temp_alpha+M-1) = diag(ones(M,1));
    temp_taf_1 = [real(temp_taf),-imag(temp_taf);imag(temp_taf),real(temp_taf)];
    temp_taf_2 = [imag(temp_taf),-real(temp_taf);real(temp_taf),imag(temp_taf)];
    Taf_1(:,:,i_temp) = temp_taf_1;
    Taf_2(:,:,i_temp) = temp_taf_2;
end
% % ADMM参数初始化
lambda_0 = 0.01*ones(2*N,1); % 拉格朗日乘子系数1
lambda_1 = 0.01*ones(N,1); % 拉格朗日乘子系数2
rho_0 = 0.01; % 罚函数项系数1
rho_1 = 0.01*ones(N,1); % 罚函数项系数2
maxnum = 200; % 最大迭代次数
Tar_out = zeros(maxnum,1);
t = 8; % 障碍函数参数

% % ADMM 求解
for i_maxnum = 1:maxnum
    
    % % Step 1 , Solving the x_r. 
    A1_temp1 = pagemtimes(pagemtimes(FNr,"transpose",Taf_1,"transpose"),FNr*br); % 子项1.1
    A1_temp1 = pagemtimes(A1_temp1,"none",A1_temp1,"transpose");
    A1_11 = sum(bsxfun(@times, A1_temp1, reshape(omega_alpha_1, 1, 1, [])),3); % 加权求和
    A1_temp2 = pagemtimes(pagemtimes(FNr,"transpose",Taf_2,"transpose"),FNr*br); % 子项1.2
    A1_temp2 = pagemtimes(A1_temp2,"none",A1_temp2,"transpose");
    A1_12 = sum(bsxfun(@times, A1_temp2, reshape(omega_alpha_1, 1, 1, [])),3); % 加权求和
    A1_temp3 = pagemtimes(chi,"transpose",br-cr,"none");
    A1_temp3 = pagemtimes(A1_temp3,"none",A1_temp3,"transpose");
    A1_2 = sum(bsxfun(@times, A1_temp3, reshape(rho_1, 1, 1, [])),3);
    A1_3 = rho_0/2*eye(2*N,2*N);
    A1 = A1_11 + A1_12 + A1_2 + A1_3;
    % Caculate the matrix B_1^T
    BT1_1 = lambda_0';
    BT1_2 = -rho_0*br';
    BT1_3_temp1 = pagemtimes(br-cr,"transpose",chi,"transpose");
    BT1_3 = sum(bsxfun(@times, BT1_3_temp1, reshape(lambda_1, 1, 1, [])), 3);
    BT1_4_temp1 = pagemtimes(chi,"none",br-cr,"none");
    BT1_4_temp1 = pagemtimes(cr,"transpose",pagemtimes(BT1_4_temp1,"none",BT1_4_temp1,"transpose"),"none");
    BT1_4 = -sum(bsxfun(@times, BT1_4_temp1, reshape(rho_1, 1, 1, [])), 3);
    BT1_5_temp1 = pagemtimes(br-cr,"transpose",chi,"none");
    BT1_5 = -(vartheta+h)*sum(bsxfun(@times, BT1_5_temp1, reshape(rho_1, 1, 1, [])), 3);
    BT1 = BT1_1 + BT1_2 + BT1_3 + BT1_4 + BT1_5;
    Tar_out(i_maxnum) = xr'*A1*xr+BT1*xr;
    if i_maxnum == 1
        fprintf("迭代开始\n初始目标函数值:%d \n",xr'*A1*xr+BT1*xr);
    end
    % 如何求解此时的lambda_xr 与 xr 
    options = optimoptions(@fmincon,'Algorithm','interior-point',...
    'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
    'HessianFcn',@(x,lambda)quadhess(x,lambda,Q,H));
    x0 = xr; % 初始点
    [xr,fval,exitflag,output] = fmincon(@(x) x'*A1*x + BT1*x, xr,[],[],[],[],[],[],...
        @(x) x'*x,options);
    % temp_xr = xr;
    % lambda_xr = 1; % ???
    % rho_xr = 10;
    % temp_xr = -pinv(2*A1+2*lambda_xr*eye(size(A1)))*BT1';
    % xr = temp_xr;
    % xr = xr/norm(xr);
    % % Step 2 , Solving the b_r
    A2_temp1 = pagemtimes(pagemtimes(FNr,"transpose",Taf_1,"transpose"),FNr*xr); % 子项1.1
    A2_temp1 = pagemtimes(A2_temp1,"none",A2_temp1,"transpose");
    A2_11 = sum(bsxfun(@times, A2_temp1, reshape(omega_alpha_1, 1, 1, [])),3); % 加权求和
    A2_temp2 = pagemtimes(pagemtimes(FNr,"transpose",Taf_2,"transpose"),FNr*xr); % 子项1.2
    A2_temp2 = pagemtimes(A2_temp2,"none",A2_temp2,"transpose");
    A2_12 = sum(bsxfun(@times, A2_temp2, reshape(omega_alpha_1, 1, 1, [])),3); % 加权求和
    A2_1 = A2_11 + A2_12;
    A2_temp3 = pagemtimes(chi,"transpose",xr-cr,"none");
    A2_temp3 = pagemtimes(A2_temp3,"none",A2_temp3,"transpose");
    A2_2 = sum(bsxfun(@times, A2_temp3, reshape(rho_1, 1, 1, [])),3);
    A2_3 = rho_0/2*eye(2*N,2*N);
    A2 = A2_1 + A2_2 + A2_3;
    BT2_1 = -lambda_0';
    BT2_2 = rho_0*xr';
    BT2_3_temp1 = pagemtimes(xr-cr,"transpose",chi,"none");
    BT2_3 = sum(bsxfun(@times,BT2_3_temp1,reshape(lambda_1,1,1,[])),3);
    BT2_4_temp1 = pagemtimes(chi,"none",xr-cr,"none");
    BT2_4_temp1 = pagemtimes(cr,"transpose",pagemtimes(BT2_4_temp1,"none",BT2_4_temp1,"transpose"),"none");
    BT2_4 = -sum(bsxfun(@times, BT1_4_temp1, reshape(rho_1, 1, 1, [])), 3);
    BT2_5_temp1 = pagemtimes(xr-cr,"transpose",chi,"none");
    BT2_5 = -(vartheta+h)*sum(bsxfun(@times, BT2_5_temp1, reshape(rho_1, 1, 1, [])), 3);
    BT2 = BT2_1 + BT2_2 + BT2_3 + BT2_4 + BT2_5;
    br = -inv(2*A2)*BT2';
    br = br/norm(br);
    % Step 3 , Solving the h
    A3 = sum(rho_1)/2;
    B3_1 = -sum(lambda_1);
    B3_temp_1 = pagemtimes(xr-cr,"transpose",chi,"none");
    B3_temp_1 = pagemtimes(B3_temp_1,"none",br-cr,"none");
    B3_2 = -sum(bsxfun(@times,B3_temp_1,reshape(rho_1,1,1,[])));
    B3_3 = vartheta*sum(rho_1);
    B3 = B3_1 + B3_2 + B3_3;
    h = (-B3+sqrt(B3^2+8*A3/t))/(4*A3);
    % Step 4 , Solving the u, v_n
    lambda_0 = lambda_0 + rho_0*(xr-br);
    for i_n = 1:N
        lambda_1(i_n) = lambda_1(i_n) + rho_1(i_n)*((xr-cr)'*chi(:,:,i_n)*(br-cr)-h-vartheta);
    end

    if i_maxnum == maxnum
        fprintf("迭代结束\n最终目标函数值:%d \n",xr'*A1*xr+BT1*xr);
    end
    
end
figure;
plot(Tar_out);
exportgraphics(gcf, './output_files/out_tar_value_iters.pdf','ContentType', 'vector');

figure;
Analysis_Sidelobe(xr(1:N),xr(1:N),'bool_draw',1);hold on;
Analysis_Sidelobe(cr(1:N),cr(1:N),'bool_draw',1);hold off;
legend('xr','cr');
exportgraphics(gcf, './output_files/out_sidelobe.pdf','ContentType', 'vector')

figure;
plot(xr(1:N));hold on;
plot(cr(1:N));hold off;
legend('xr','cr');
exportgraphics(gcf, './output_files/xr_and_cr.pdf','ContentType', 'vector');

figure 
Analysis_CS_DFSM(xr(1:N),fs,fs/N,M);
exportgraphics(gcf, './output_files/xr_CS.pdf','ContentType', 'vector');

figure 
Analysis_CS_DFSM(cr(1:N),fs,fs/N,M);
exportgraphics(gcf, './output_files/cr_CS.pdf','ContentType', 'vector');