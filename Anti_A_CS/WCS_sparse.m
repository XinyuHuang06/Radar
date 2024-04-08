% Example:
% :param :
% :return :
% detailed description: 并行化优化低截获雷达波形程序--稀疏矩阵版本（高维度与内存不足情况下运行）
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
N = 8; % 信号长度
M = 2; % 窗长
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
omega_alpha_1 = ones(N^2,1); % 循环谱加权向量
m = -M/2+1:M/2;
Taf_1_sparse = cell(N^2,1);
Taf_2_sparse = cell(N^2,1);
for iter_1 = 1:N
    for iter_2 = 1:N
        i_temp = (iter_1-1)*N+iter_2;
        B = max(1-iter_1, 1-iter_2);
        A = min (N-iter_1, N-iter_2);
        n = m((m<=A) & (m>=B));
        if isempty(n)
            temp_taf = zeros(N,N);
        else
            taf_1 = zeros(N,1);
            taf_2 = zeros(N,1);
            temp_p = iter_1 + n;
            temp_q = iter_2 + n;
            taf_1(temp_p) = 1;
            taf_2(temp_q) = 1;
            taf_1 = diag(taf_1);
            taf_2 = diag(taf_2);
            temp_taf = taf_1*taf_2;
        end
        temp_taf_1 = [real(temp_taf),-imag(temp_taf);imag(temp_taf),real(temp_taf)];
        temp_taf_2 = [imag(temp_taf),-real(temp_taf);real(temp_taf),imag(temp_taf)];
        Taf_1_sparse{i_temp} = uint8(temp_taf_1);
        Taf_2_sparse{i_temp} = uint8(temp_taf_2);
        if abs(iter_1-iter_2) < 10 
            omega_alpha_1(i_temp) = 0.01;
        else
            omega_alpha_1(i_temp) = 10;
        end
    end
end


% 初始化一个空的 cell 数组来存储结果
result = {};

% 创建一个副本用于搜索
searchArray = Taf_1_sparse;

% 使用一个 for 循环来遍历 Taf_1_sparse 中的每个元素
for i = 1:length(Taf_1_sparse)
    % 使用 cellfun 和 isequal 函数来找出与当前元素相同的所有元素的下标
    idx = find(cellfun(@(x) isequal(x, Taf_1_sparse{i}), searchArray));
    
    % 如果找到了多于一个的相同元素，那么就将这些下标添加到结果中
    if length(idx) > 1
        result{end+1} = idx;
        % 将已经找到的重复元素从搜索范围中移除
        searchArray(idx(2:end)) = {[]};
    end
end

% 打印结果
disp(result);




clear temp_* *_temp
Tar_out = zeros(maxnum,1);
% % ADMM 求解
for i_maxnum = 1:maxnum
    % % Step 1 , Solving the x_r. 
    A1_1 = 0;
    A1_2 = 0;
    for i_temp = 1:N^2
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

    % lambda_xr = 0; % ???
    % while true
    %     xr = -inv(2*A1+2*lambda_xr*eye(size(A1)))*BT1';
    %     if abs(xr'*xr-1) < 1
    %         break;
    %     else
    %         continue;
    %     end
    % end
    % % xr = xr/norm(xr);
    temp_H = cell(1);temp_k= cell(1);temp_d= cell(1);
    temp_H{1} = eye(2*N)*2;
    temp_k{1} = zeros(2*N,1);
    temp_d{1} = -1;
    options = optimoptions(@fmincon,'Algorithm','interior-point',...
        'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
        'HessianFcn',@(x,lambda)quadhess(x,lambda,A1/2,temp_H));
    x0 = xr; % 初始点

    [xr,fval,exitflag,output] = fmincon(@(x) quadobj(x,A1/2,BT1',0), xr,[],[],[],[],[],[],...
        @(x) quadconstr(x,temp_H,temp_k,temp_d),options);
    % % Step 2 , Solving the b_r
    A2_1 = 0;
    A2_2 = 0;
    for i_temp = 1:N^2
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
    % br = br/norm(br);
    % Step 3 , Solving the h
    A3 = sum(rho_1)/2;
    B3_1 = -sum(lambda_1);
    B3_2 = 0;
    for i_temp = 1:N
        B3_2 = -rho_1(i_temp)*(xr-cr)'*chi_sparse{i_temp}*(br-cr);
    end
    B3_3 = vartheta*sum(rho_1);
    B3 = B3_1 + B3_2 + B3_3;
    h = (-B3+sqrt(B3^2+8*A3/t))/(4*A3);
    % Step 4 , Solving the u, v_n
    lambda_0 = lambda_0 + rho_0*(xr-br);
    for i_n = 1:N
        lambda_1(i_n) = lambda_1(i_n) + rho_1(i_n)*((xr-cr)'*chi_sparse{i_n}*(br-cr)-h-vartheta);
    end
    if i_maxnum == maxnum
        fprintf("%f",xr'*A1*xr + BT1*xr);
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

Analysis_CS_DFSM(fs,xr(1:N),fs/N,M,'bool_draw',1);
exportgraphics(gcf, './output_files/xr_CS.pdf','ContentType', 'vector');

Analysis_CS_DFSM(fs,cr(1:N),fs/N,M,'bool_draw',1);
exportgraphics(gcf, './output_files/cr_CS.pdf','ContentType', 'vector');

% 定义目标函数
function [y,grady] = quadobj(x,Q,f,c)
    y = 1/2*x'*Q*x + f'*x + c;
    if nargout > 1
        grady = Q*x + f;
    end
end

% 定义约束函数
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