% Example:
% :param :
% :return :
% detailed description: 并行化优化低截获雷达波形程序--稀疏矩阵版本（高维度与内存不足情况下运行）
%------------------------------------------------------------------------------
% V 0.0.1
% Created by: Xinyu Huang.
% On: 26/03/2024.
% V 0.0.2
% Modified By: Xinyu Huang
% On: 10/04/2024
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
close all
% % 数据初始化

N = 128; % 信号长度
M = 4; % 窗长
x = exp(1j*2*pi*rand(N,1))/sqrt(N);  % 目标信号，随机初始化
xr = [real(x);imag(x)]; 
br = xr; % 辅助变量
FN = dftmtx(N); % 离散傅里叶变换矩阵
Fshift = zeros(N);
Fshift(1:N/2,N/2+1:end) = eye(N/2);
Fshift(N/2+1:end,1:N/2) = eye(N/2);
FN = Fshift*FN;
FNr = complex2real(FN);
% % cr 参考信号生成
fs = 5e6;
T = N/fs;
fc = 1e5;
temp_B = 2e5;
c = generator_LFM(fs,fc,temp_B,T)/sqrt(N);
cr = [real(c);imag(c)];
vartheta = 0.1/sqrt(N); % 参考信号相关度阈值

% 相似性约束矩阵
chi_sparse = cell(N,1);
for i_temp = 1:N
    matrix_temp = zeros(N,N);
    matrix_temp(i_temp,i_temp) = 1;
    chi_sparse{i_temp} = sparse(complex2real(matrix_temp));
end

% 循环谱计算矩阵 20240410 Modified
% 考虑到循环谱计算矩阵是对称的，在这里我们可以仅仅计算接近一半的数值，由N^2 -> (N+1)*N/2，减小接近一半的计算量与存储空间
omega_alpha_1 = ones((N+1)*N/2,1); % 循环谱加权向量
temp_m = -M/2+1:M/2;
Taf_1_sparse = cell((N+1)*N/2,1);
Taf_2_sparse = cell((N+1)*N/2,1);
i_temp = 0;

for temp_1 = 1:N
    for temp_2 = temp_1:N
        i_temp = i_temp + 1;
        temp_B = max(1-temp_1, 1-temp_2);
        temp_A = min (N-temp_1, N-temp_2);
        temp_n = temp_m((temp_m<=temp_A) & (temp_m>=temp_B));
        if isempty(temp_n)
            temp_taf = zeros(N,N);
        else
            temp_p = temp_1 + temp_n;
            temp_q = temp_2 + temp_n;
            temp_taf = zeros(N,N);
            temp_taf(temp_p(1):temp_p(end),temp_q(1):temp_q(end)) = eye(length(temp_q));
        end
        temp_taf_1 = [real(temp_taf),-imag(temp_taf);imag(temp_taf),real(temp_taf)];
        temp_taf_2 = [imag(temp_taf),-real(temp_taf);real(temp_taf),imag(temp_taf)];
        Taf_1_sparse{i_temp} = sparse(temp_taf_1);
        Taf_2_sparse{i_temp} = sparse(temp_taf_2);
        if abs(temp_1-temp_2) < 20 
            omega_alpha_1(i_temp) = 0.01;
        else
            omega_alpha_1(i_temp) = 10;
        end
    end
end

% % The initial parameter Caculation
s_xr_init = struct('signal',[],'ISL',[],'ISLR',[],'PSL',[],'PSLR',[],'PAPR',[]);
s_xr_init.PAPR = Evaluate_PAPR(xr(1:N));
[~ ,s_xr_init.PSL, s_xr_init.ISL, s_xr_init.PSLR, s_xr_init.ISLR] = Analysis_Sidelobe(xr(1:N),xr(1:N));

% % ADMM参数初始化
% 拉格朗日乘子系数
h = -0.1; % 相似度约束辅助变量
lambda_0 = 1*ones(2*N,1);
lambda_1 = 1*ones(N,1); 
rho_0 = 1;
rho_1 = 10*ones(N,1); % 罚函数项系数
maxnum = 50; % 最大迭代次数
t = 8; % 障碍函数系数


Tar_out = zeros(maxnum,1);
% % ADMM 求解
forwaitbar = forWaitbar(maxnum);
for i_maxnum = 1:maxnum
    % % Step 1 , Solving the x_r. 
    A1_1 = 0;
    A1_2 = 0;
    for i_temp = 1:(N+1)*N/2
            A1_1_1 = FNr'*(Taf_1_sparse{i_temp})'*(FNr*br);
            A1_1_2 = FNr'*(Taf_2_sparse{i_temp})'*(FNr*br);
            A1_1 = A1_1 + omega_alpha_1(i_temp)*(A1_1_1*A1_1_1'+A1_1_2*A1_1_2')/((N+1)*N*10/2);
    end
    for i_temp = 1:N
        A1_temp_3= (chi_sparse{i_temp})'*(br-cr);
        A1_2 = A1_2 + rho_1(i_temp)*(A1_temp_3*A1_temp_3');
    end
    A1_3 = rho_0/2*eye(2*N,2*N);
    A1 = A1_1 + A1_2 + A1_3;
    
    BT1_1 = lambda_0'; % Caculate the matrix B_1^T
    BT1_2 = -rho_0*br';
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
    Tar_out(i_maxnum) = xr'*A1_1*xr;

    temp_H = cell(1);temp_k= cell(1);temp_d= cell(1);
    temp_H{1} = eye(2*N)*2;
    temp_k{1} = zeros(2*N,1);
    temp_d{1} = -1;
    options = optimoptions(@fmincon,'Algorithm','interior-point',...
        'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
        'HessianFcn',@(x,lambda)quadhess(x,lambda,A1/2,temp_H),'Display','off');
    % x0 = xr; % 初始点
    [xr,~,~,~] = fmincon(@(x) quadobj(x,A1/2,BT1',0), xr,[],[],[],[],[],[],...
        @(x) quadconstr(x,temp_H,temp_k,temp_d),options);
    % if sum(abs(br - xr)> vartheta) > 1
    %     fprintf([num2str(sum(abs(br - xr)> vartheta)),'\n']);
    % end
    % % Step 2 , Solving the b_r
    A2_1 = 0;
    A2_2 = 0;
    for i_temp = 1:(N+1)*N/2
        A2_temp1 = FNr'*(Taf_1_sparse{i_temp})'*(FNr*xr);
        A2_temp2 = FNr'*(Taf_2_sparse{i_temp})'*(FNr*xr);
        A2_1 = A2_1 + omega_alpha_1(i_temp)*(A2_temp1*A2_temp1'+A2_temp2*A2_temp2'/((N+1)*N*10/2));
    end
    for i_temp = 1:N
        A2_temp3 = (chi_sparse{i_temp})'*(xr-cr);
        A2_2 = A2_2 + rho_1(i_temp)*(A2_temp3*A2_temp3');
    end
    A2_3 = rho_0/2*eye(2*N,2*N);
    A2 = A2_1 + A2_2 + A2_3;
    
    BT2_1 = -lambda_0';
    BT2_2 = rho_0*xr';
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
    
    

    % br = br / norm(br);
    temp_H = cell(1);temp_k= cell(1);temp_d= cell(1);
    temp_H{1} = eye(2*N)*2;
    temp_k{1} = zeros(2*N,1);
    temp_d{1} = -1;
    options = optimoptions(@fmincon,'Algorithm','interior-point',...
        'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
        'HessianFcn',@(x,lambda)quadhess(x,lambda,A1/2,temp_H),'Display','off');
    % x0 = xr; % 初始点
    [br,fval,exitflag,output] = fmincon(@(x) quadobj(x,A2/2,BT2',0), br,[],[],[],[],[],[],...
        @(x) quadconstr(x,temp_H,temp_k,temp_d),options);   
    % br = -inv(2*A2)*BT2';
    % Test_1 = br'*A2*br + BT2*br;
    % Test_2 = br'*A2*br + BT2*br;
    % br = br/norm(br);
    % Test_3 = br'*A2*br + BT2*br;
    % disp(Test_1);disp(Test_2);disp(Test_3);
    % fprintf('1:%.2d.  2:%.2d   . 3:%.2d \n',Test_1,Test_2,Test_3);
    % Step 3 , Solving the h
    A3 = sum(rho_1)/2;
    B3_1 = -sum(lambda_1);
    B3_2 = 0;
    for i_temp = 1:N
        B3_2 = -rho_1(i_temp)*(xr-cr)'*chi_sparse{i_temp}*(br-cr);
    end
    B3_3 = vartheta*sum(rho_1);
    B3 = B3_1 + B3_2 + B3_3;
    h = (-B3-sqrt(B3^2+8*A3/t))/(4*A3);
    % Step 4 , Solving the u, v_n
    lambda_0 = lambda_0 + rho_0*(xr-br);
    for i_n = 1:N
        lambda_1(i_n) = lambda_1(i_n) + rho_1(i_n)*((xr-cr)'*chi_sparse{i_n}*(br-cr)-h-vartheta);
    end
    forwaitbar.show_bar;
end
clear temp_* *_temp A1_* BT1_* A2_* BT2_* B3_*


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
Analysis_CS_DFSM(fs,xr(1:N),fs/N,M,'bool_draw',1);
exportgraphics(gcf, './output_files/xr_CS.pdf','ContentType', 'vector');
Analysis_CS_DFSM(fs,cr(1:N),fs/N,M,'bool_draw',1);
exportgraphics(gcf, './output_files/cr_CS.pdf','ContentType', 'vector');


% % The post-parameter after optimization
s_xr_post = struct('signal',[],'ISL',[],'ISLR',[],'PSL',[],'PSLR',[],'PAPR',[]);
s_xr_post.PAPR = Evaluate_PAPR(xr(1:N));
[~ ,s_xr_post.PSL, s_xr_post.ISL, s_xr_post.PSLR, s_xr_post.ISLR] = Analysis_Sidelobe(xr(1:N),xr(1:N));
% % The parameter of signal c Caculation
s_cr = struct('signal',[],'ISL',[],'ISLR',[],'PSL',[],'PSLR',[],'PAPR',[]);
s_cr.PAPR = Evaluate_PAPR(cr(1:N));
[~ ,s_cr.PSL, s_cr.ISL, s_cr.PSLR, s_cr.ISLR] = Analysis_Sidelobe(cr(1:N),cr(1:N));


% The Qudratic Target Function
function [y,grady] = quadobj(x,Q,f,c)
    y = 1/2*x'*Q*x + f'*x + c;
    if nargout > 1
        grady = Q*x + f;
    end
end

% The Qudratic Constraint Function
function [y,yeq,grady,gradyeq] = quadconstr(x,H,k,d)
    jj = length(H);  % jj is the number of equality constraints
    yeq = zeros(1,jj);
    for i = 1:jj
        yeq(i) = 1/2*x'*H{i}*x + k{i}'*x + d{i};
    end
    y = [];
    if nargout > 2
        gradyeq = zeros(length(x),jj);
        for i = 1:jj
            gradyeq(:,i) = H{i}*x + k{i};
        end
    end
    grady = [];
end

% 定义Hessian函数
function hess = quadhess(x,lambda,Q,H)
    hess = Q;
    jj = length(H);  % jj is the number of equality constraints
    for i = 1:jj
        hess = hess + lambda.eqnonlin(i)*H{i};
    end
end