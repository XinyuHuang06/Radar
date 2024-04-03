% Example:
% :param :
% :return :
% detailed description: 并行化优化低截获雷达波形程序
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
clc
clear

% % 数据初始化
N = 128; % 信号长度
M = 16; % 窗长
% 目标信号，随机初始化
x = exp(1j*2*pi*rand(N,1))/sqrt(N);  
xr = [real(x);imag(x)]; 
% 辅助变量
br = xr; 
% 离散傅里叶变换矩阵
FN = dftmtx(N); 
FNr = complex2real(FN);
% cr 参考信号生成
T = N*1e-6;
fs = 1e6;
fc = 1e5;
B = 2e5;
c = generator_LFM(fs,fc,B,T)/sqrt(N);
cr = [real(c);imag(c)];
vartheta = 0.1; % 参考信号相关度阈值
h = -0.1; % 相似度约束辅助变量
% 相似性约束矩阵
chi_sparse = cell(N,1);
for i_temp = 1:N
    matrix_temp = zeros(N,N);
    matrix_temp(i_temp,i_temp) = 1;
    chi_sparse{i_temp} = sparse(complex2real(matrix_temp));
end
% 循环谱加权向量
omega_alpha = zeros((N-M)/2,1);
omega_alpha(1:40) = 0.01;
omega_alpha(41:end) = 10;

omega_alpha_1 = ones((N-M)*(N-M),1);
omega_alpha_1 = sparse(omega_alpha_1);

% % ADMM参数初始化
% 拉格朗日乘子系数
lambda_0 = 0.1*ones(2*N,1);
lambda_1 = 0.1*ones(N,1); 
rho_0 = 0.1;
rho_1 = 0.1*ones(N,1); % 罚函数项系数
% 最大迭代次数
maxnum = 100; 
t = 8;


% 循环谱计算矩阵 
alpha = 1:N-M;
cf = 1:N-M;
Taf_1_sparse = cell((N-M)*(N-M),1);
Taf_2_sparse = cell((N-M)*(N-M),1);

for i_temp = 1:(N-M)*(N-M)
    temp_alpha = fix(i_temp/(N-M))+1;
    temp_cf = mod(i_temp,N-M)+1;
    temp_taf = zeros(N,N);
    temp_taf(temp_cf:temp_cf+M-1, temp_alpha:temp_alpha+M-1) = diag(ones(M,1));
    temp_taf_1 = [real(temp_taf),-imag(temp_taf);imag(temp_taf),real(temp_taf)];
    temp_taf_2 = [imag(temp_taf),-real(temp_taf);real(temp_taf),imag(temp_taf)];
    Taf_1_sparse{i_temp} = sparse(temp_taf_1);
    Taf_2_sparse{i_temp} = sparse(temp_taf_2);
end

clear temp_* *_temp
Tar_out = zeros(maxnum,1);
% % ADMM 求解
for i_maxnum = 1:maxnum
    % % Step 1 , Solving the x_r. 
    % 
    % A1_temp1 = pagemtimes(pagemtimes(FNr,"transpose",Taf_1,"transpose"),FNr*br); % 子项1.1
    % A1_temp1 = pagemtimes(A1_temp1,"none",A1_temp1,"transpose");
    % A1_11 = sum(bsxfun(@times, A1_temp1, reshape(omega_alpha_1, 1, 1, [])),3); % 加权求和
    % A1_temp2 = pagemtimes(pagemtimes(FNr,"transpose",Taf_2,"transpose"),FNr*br); % 子项1.2
    % A1_temp2 = pagemtimes(A1_temp2,"none",A1_temp2,"transpose");
    % A1_12 = sum(bsxfun(@times, A1_temp2, reshape(omega_alpha_1, 1, 1, [])),3); % 加权求和
    % A1_temp3 = pagemtimes(chi,"transpose",br-cr,"none");
    % A1_temp3 = pagemtimes(A1_temp3,"none",A1_temp3,"transpose");
    % A1_2 = sum(bsxfun(@times, A1_temp3, reshape(rho_1, 1, 1, [])),3);
    A1_1 = 0;
    A1_2 = 0;
    for i_temp = 1:(N-M)^2
        A1_temp1 = FNr'*(Taf_1_sparse{i_temp})'*(FNr*br);
        A1_temp2 = FNr'*(Taf_2_sparse{i_temp})'*(FNr*br);
        A1_1 = A1_1 + omega_alpha_1(i_temp)*(A1_temp1*A1_temp1'+A1_temp2*A1_temp2');
    end
    for i_temp = 1:N
        A1_temp3 = (chi_sparse{i_temp})'*(br-cr);
        A1_2 = A1_2 + rho_1(i_temp)*(A1_temp3*A1_temp3');
    end
    A1_3 = rho_0/2*eye(2*N,2*N);
    A1 = A1_1 + A1_2 + A1_3;
    % Caculate the matrix B_1^T
    BT1_1 = lambda_0';
    BT1_2 = -rho_0*br';

    % BT1_3_temp1 = pagemtimes(br-cr,"transpose",chi,"transpose");
    % BT1_3 = sum(bsxfun(@times, BT1_3_temp1, reshape(lambda_1, 1, 1, [])), 3);
    % 
    % BT1_4_temp1 = pagemtimes(chi,"none",br-cr,"none");
    % BT1_4_temp1 = pagemtimes(cr,"transpose",pagemtimes(BT1_4_temp1,"none",BT1_4_temp1,"transpose"),"none");
    % BT1_4 = -sum(bsxfun(@times, BT1_4_temp1, reshape(rho_1, 1, 1, [])), 3);
    % 
    % BT1_5_temp1 = pagemtimes(br-cr,"transpose",chi,"none");
    % BT1_5 = -(vartheta+h)*sum(bsxfun(@times, BT1_5_temp1, reshape(rho_1, 1, 1, [])), 3);
    BT1_3 = 0;
    BT1_4 = 0;
    BT1_5 = 0;
    for i_temp = 1:N
        BT1_3 = BT1_3 + lambda_1(i_temp)*(br-cr)'*(chi_sparse{i_temp}');
        BT1_4_temp1 = (br-cr)'*(chi_sparse{i_temp})';
        BT1_4 = BT1_4 + rho_1(i_temp)*cr'*(BT1_4_temp1*BT1_4_temp1');
        BT1_5 = -(vartheta+h)*rho_1(i_temp)*(br-cr)'*(chi_sparse{i_temp});
    end
    BT1 = BT1_1 + BT1_2 + BT1_3 + BT1_4 + BT1_5;
    % 如何求解此时的lambda_xr 与 xr 
    Tar_out(i_maxnum) = xr'*A1*xr+BT1*xr;
    if i_maxnum == 1
        fprintf("迭代开始\n初始目标函数值:%d \n",xr'*A1*xr+BT1*xr);
    end
    lambda_xr = 0; % ???
    % 
    xr = -inv(2*A1+2*lambda_xr*eye(size(A1)))*BT1';
    xr = xr/norm(xr);
    
    % % Step 2 , Solving the b_r
    % A2_temp1 = pagemtimes(pagemtimes(FNr,"transpose",Taf_1,"transpose"),FNr*xr); % 子项1.1
    % A2_temp1 = pagemtimes(A2_temp1,"none",A2_temp1,"transpose");
    % A2_11 = sum(bsxfun(@times, A2_temp1, reshape(omega_alpha_1, 1, 1, [])),3); % 加权求和
    % A2_temp2 = pagemtimes(pagemtimes(FNr,"transpose",Taf_2,"transpose"),FNr*xr); % 子项1.2
    % A2_temp2 = pagemtimes(A2_temp2,"none",A2_temp2,"transpose");
    % A2_12 = sum(bsxfun(@times, A2_temp2, reshape(omega_alpha_1, 1, 1, [])),3); % 加权求和
    % A2_1 = A2_11 + A2_12;
    % A2_temp3 = pagemtimes(chi,"transpose",xr-cr,"none");
    % A2_temp3 = pagemtimes(A2_temp3,"none",A2_temp3,"transpose");
    % A2_2 = sum(bsxfun(@times, A2_temp3, reshape(rho_1, 1, 1, [])),3);
    % A2_3 = rho_0/2*eye(2*N,2*N);
    A2_1 = 0;
    A2_2 = 0;
    for i_temp = 1:(N-M)^2
        A2_temp1 = FNr'*(Taf_1_sparse{i_temp})'*(FNr*xr);
        A2_temp2 = FNr'*(Taf_2_sparse{i_temp})'*(FNr*xr);
        A2_1 = A2_1 + omega_alpha_1(i_temp)*(A2_temp1*A2_temp1'+A2_temp2*A2_temp2');
    end
    for i_temp = 1:N
        A2_temp3 = (chi_sparse{i_temp})'*(xr-cr);
        A2_2 = A2_2 + rho_1(i_temp)*(A2_temp3*A2_temp3');
    end
    A2_3 = rho_0/2*eye(2*N,2*N);
    A2 = A2_1 + A2_2 + A2_3;
    
    BT2_1 = -lambda_0';
    BT2_2 = rho_0*xr';
    % 
    % BT2_3_temp1 = pagemtimes(xr-cr,"transpose",chi,"none");
    % BT2_3 = sum(bsxfun(@times,BT2_3_temp1,reshape(lambda_1,1,1,[])),3);
    % 
    % BT2_4_temp1 = pagemtimes(chi,"none",xr-cr,"none");
    % BT2_4_temp1 = pagemtimes(cr,"transpose",pagemtimes(BT2_4_temp1,"none",BT2_4_temp1,"transpose"),"none");
    % BT2_4 = -sum(bsxfun(@times, BT1_4_temp1, reshape(rho_1, 1, 1, [])), 3);
    % 
    % BT2_5_temp1 = pagemtimes(xr-cr,"transpose",chi,"none");
    % BT2_5 = -(vartheta+h)*sum(bsxfun(@times, BT2_5_temp1, reshape(rho_1, 1, 1, [])), 3);
    BT2_3 = 0;
    BT2_4 = 0;
    BT2_5 = 0;
    for i_temp = 1:N
        BT2_3 = BT2_3 + lambda_1(i_temp)*(xr-cr)'*chi_sparse{i_temp};
        BT2_4_temp1 = chi_sparse{i_temp}*(xr-cr);
        BT2_4 = -rho_1(i_temp)*cr'*(BT2_4_temp1*BT2_4_temp1');
        BT2_5 = -(vartheta+h)*rho_1(i_temp)*(xr-cr)'*(chi_sparse{i_temp});
    end

    BT2 = BT2_1 + BT2_2 + BT2_3 + BT2_4 + BT2_5;

    br = -inv(2*A2)*BT2';
    br = br/norm(br);
    % Step 3 , Solving the h
    A3 = sum(rho_1)/2;
    B3_1 = -sum(lambda_1);

    % B3_temp_1 = pagemtimes(xr-cr,"transpose",chi,"none");
    % B3_temp_1 = pagemtimes(B3_temp_1,"none",br-cr,"none");
    % B3_2 = -sum(bsxfun(@times,B3_temp_1,reshape(rho_1,1,1,[])));
    B3_2 = 0;
    for i_temp = 1:N
        B3_2 = -rho_1(i_temp)*(xr-cr)'*chi_sparse{i_temp}*(br-cr);
    end
    B3_3 = vartheta*sum(rho_1);
    B3 = B3_1 + B3_2 + B3_3;
    h = (-B3+sqrt(B3^2+8*A3/t))/(4*A3);
    % Step 4 , Solving the u, v_n
    lambda_0 = lambda_0 + rho_0*(xr-br);
    % for i_n = 1:N
    %     lambda_1(i_n) = lambda_1(i_n) + rho_1(i_n)*((xr-cr)'*chi(:,:,i_n)*(br-cr)-h-vartheta);
    % end
    for i_n = 1:N
        lambda_1(i_n) = lambda_1(i_n) + rho_1(i_n)*((xr-cr)'*chi_sparse{i_n}*(br-cr)-h-vartheta);
    end
    if i_maxnum == 1
        fprintf("%f",xr'*A1*xr + BT1*xr);
    elseif i_maxnum == maxnum
        fprintf("%f",xr'*A1*xr + BT1*xr);
        
    end
end
