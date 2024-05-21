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

N = 128;
M =4;
max_iternums = 100;
% 
FN = dftmtx(N); 
Fshift = zeros(N);
Fshift(1:N/2,N/2+1:end) = eye(N/2);
Fshift(N/2+1:end,1:N/2) = eye(N/2);
FN = Fshift*FN;
taf = zeros(N,N);
temp_m = -M/2+1:M/2;
for temp_1 = 1:N
    for temp_2 = temp_1:N
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
            taf = taf + temp_taf;
        end
    end
end   

% % 
S = FN'*taf*FN;

for i_m = 1:max_iternums
    for k = 1:N
        % Coordinate-Descent Iterations
        a_k = real(S(k,k));
        S_k_n = S(k,:); S_k_n(k) = 0;
        S_m_k = S(:,k); S_m_k(k) = 0;
        b_k = real(S_k_n*x) + real(x'*S_m_k);
        c_k = imag(S_k_n*x) + imag(x'*S_m_k);
        S_m_n = S; S_m_n(k,:) = 0; S_m_n(:,k) = 0;
        d_k = real(x'*S_m_n*x);
        % % Dinkebach 迭代
        lambda = 10;
        delta = 1e-6;
        while(true)
            beta_k = argmin_f_beta_k(lambda);
            if f_beta_k(a_k,b_k,c_k,d_k,beta_k) <= delta

            end

        end

    end
end



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