% Example:
% :param :
% :return :
% detailed description:
%------------------------------------------------------------------------------
% Created by: Xinyu Huang.
% On: 20/05/2024.
% Copyright (C) 2024 Xinyu Huang (learning_huang@163.com).
% All Rights Reserved.
% Unauthorized copying of this file, via any medium is strictly prohibited.
% Proprietary and confidential.
%------------------------------------------------------------------------------
clc; clear; close all;
N = 256;
M = 4;
max_iternums = 100;
% 
FN = dftmtx(N); 
Fshift = zeros(N);
Fshift(1:N/2,N/2+1:end) = eye(N/2);
Fshift(N/2+1:end,1:N/2) = eye(N/2);
FN = Fshift*FN;
taf = zeros(N,N);
temp_m = -M/2+1:M/2;


x = exp(1j*2*pi*rand(N,1));
for temp_1 = 1:N
    for temp_2 = 1:N
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
            if abs(temp_1 - temp_2) > -10 && abs(temp_1 - temp_2) < 10
                taf = taf + 0.01*temp_taf;
            else
                taf = taf + temp_taf;
            end
        end
    end
end 
S = FN'*taf*FN;
if ~issymmetric(taf)
    fprintf('为非对称矩阵！\n'); % 若为非对称矩阵终止脚本。
    return;
end
% if ~ishermitian(S)
%     fprintf('为非共轭对称矩阵！\n'); % 若为非共轭堆成矩阵终止脚本。
%     return;
% end
% % 
Analysis_CS_DFSM(1e6, x, 1e6/100, M, 'bool_draw', 1)
target_seq = zeros(max_iternums,1);
for i_m = 1:max_iternums
    target_seq(i_m) = real(x'*S*x);
    for k = 1:N

        % Coordinate-Descent Iterations
        x_nonk = x;
        x_nonk(k) = 0;
        e_k = zeros(N,1);
        e_k(k) = 1;
        x_k = exp(-1j*angle(x_nonk'*S*e_k));

        x(k) = x_k;
        % a_k = real(S(k,k));
        % S_k_n = S(k,:); S_k_n(k) = 0;
        % S_m_k = S(:,k); S_m_k(k) = 0;
        % temp_x = x; 
        % b_k = real(S_k_n*temp_x) + real(temp_x'*S_m_k);
        % c_k = imag(S_k_n*temp_x) + imag(temp_x'*S_m_k);
        % S_m_n = S; S_m_n(k,:) = 0; S_m_n(:,k) = 0;
        % d_k = real(temp_x'*S_m_n*temp_x);
        % if a_k > 0
        %     % a_k > 0 
        %     PointSet = roots([4*a_k, -2*a_k, 0, -2*b_k, 2*a_k-4*c_k, 2*b_k]);
        %     PointSet = PointSet(real(PointSet)==PointSet);
        %     temp_value = zeros(size(PointSet));
        %     for i_point = 1:length(PointSet)
        %         temp_value(i_point) = f_beta_k(a_k,b_k,c_k,d_k,PointSet(i_point));
        %     end
        %     [~,max_index] = min(temp_value);
        %     beta_k = PointSet(max_index); % 找到最优的beta_K值
        % elseif a_k < 0
        %     beta_k = -100000;
        % end
    
        % beta_k1 = tan(angle(x(k))/2);
        % beta_k2 = beta_k;
 

        % temp_x = x;
        % x(k) = exp(1j*2*atan(beta_k));
        % disp([i_m,k]);
        % fprintf(['fvalue--Pre:%d\n','After:%d\n'],...
        %     f_beta_k(a_k,b_k,c_k,d_k,beta_k1),f_beta_k(a_k,b_k,c_k,d_k,beta_k2));
        % fprintf(['Tarvalue--Pre:%d\n','After:%d\n'],...
        %     real(temp_x'*S*temp_x),real(x'*S*x));
        % % % Dinkebach 迭代
        % lambda = 10;
        % delta = 1e-6;
        % while(true)
        %     beta_k = argmin_f_beta_k(lambda);
        %     if f_beta_k(a_k,b_k,c_k,d_k,beta_k) <= delta

        %     end

        % end

    end
end

Analysis_CS_DFSM(1e6, x, 1e6/100, M, 'bool_draw', 1)

function value = f_beta_k(a_k,b_k,c_k,d_k,beta_k)
    value = g_beta_k(a_k,b_k,c_k,d_k,beta_k)/h_beta_k(beta_k);
end

function value = g_beta_k(a_k,b_k,c_k,d_k,beta_k)
    value = a_k*beta_k^4 + (2*a_k - c_k + d_k)*beta_k^2 +...
            2*b_k*beta_k + (a_k + c_k +d_k);
end

function value = h_beta_k(beta_k)
    value = (1+beta_k^2);
end