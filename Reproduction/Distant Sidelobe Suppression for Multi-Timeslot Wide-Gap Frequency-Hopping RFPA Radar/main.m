% Example:
% :param :
% :return :
% detailed description:
%------------------------------------------------------------------------------
% Created by: Xinyu Huang.
% On: 21/05/2024.
% Copyright (C) 2024 Xinyu Huang (learning_huang@163.com).
% All Rights Reserved.
% Unauthorized copying of this file, via any medium is strictly prohibited.
% Proprietary and confidential.
%------------------------------------------------------------------------------
% % Initialization parameter
% 雷达参数
addpath("./function/")
global parapack;
parapack = parameter;

% % 2
% n = 6; % 第n个传输脉冲, 0 <= n <= M-1
% l_0_p = floor( c/(2*Delta_r)*(p*T_S - T_W - T_M) );
% l_1_p = ceil( c/(2*Delta_r)*(p*T_S - T_R + T_J + T_M) );
% l_2_p = floor( c/(2*Delta_r)*(p*T_S + T_R - T_J - T_W -T_M));
% L_W   = ceil( c/(2*Delta_r)*(T_W + 2*T_M) ) + 1;

N_1 = min(N, M-1-n);
N_2 = min(N, n);


P = 1024; % LPF滤波器点数
h_p = fir1(P-1, 0.5); % 低通滤波，截止频率fs/4

% CS解析数据


rmpath("./function/")



function out = s_m_t(m,t)
    t_m = m*parapack.T_R + parapack.T_m_Seq(m);
    if (t - parapack.t_m)/parapack.T_W >= 0 && (t - t_m)/parapack.T_W <= 1 
        out = exp(1j*(2*pi*parapack.f_m_Seq(m)*t + parapack.phi_0));
    else
        out = zeros(size(t));
    end
end

function out = s_m_l_k(n,p,m,l,k)
    t_n = n*parapack.T_R + parapack.T_m_Seq(n);
    t_n_p = t_n + p*parapack.T_S;
    out = s_m_t(m, t_n_p-2*(l*parapack.Delta_r-k*parapack.Delta_v*t_n_p)/parapack.c);
end

