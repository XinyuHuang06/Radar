function [signal_NLFM, t_grid] = generator_NLFM(varargin)
% Example: 
% :param fs:
% :param fc:
% :param 
% :return signal_NLFM:
% :return t:
% 依据相位逗留原理设计NLFM信号
%------------------------------------------------------------------------------
% Created by: Xinyu Huang.
% On: 20/10/2023.
% Copyright (C) 2023 Xinyu Huang (learning_huang@163.com).
% All Rights Reserved.
% Unauthorized copying of this file, via any medium is strictly prohibited.
% Proprietary and confidential.
%------------------------------------------------------------------------------
    in_par = inputParser;
    addOptional(in_par, 'fs', 0);
    addOptional(in_par, 'fc', 0);
    addOptional(in_par, 'B', 0);
    addOptional(in_par, 'N', 0);
    addOptional(in_par, 'fun', 0);
    parse(in_par, varargin{:});
    fs = in_par.Results.fs;
    fc = in_par.Results.fc;
    B = in_par.Results.B;
    N = in_par.Results.N;
    fun = in_par.Results.fun;

    T0 = N/fs;
    t_grid = (0:1/fs:(N-1)/fs)';
    fun = fun(:)/max(fun); % 转化为列向量并归一化
    fun = fun*B;
    m_fun = tril(ones(N))*fun; % 积分
    m_fun = fun;
    carrier = exp(1j*2*pi*fc*t_grid);
    signal_NLFM = carrier.*exp(1j*2*pi*m_fun.*t_grid);

%传统基于逗留相位原理的窗函数反求法-NLFM信号
clear all;clc ;close all 

fs=100e6;%100Mhz
fc=10e6;%10Mhz
T=10.24e-6;
B=40e6;
N=T*fs;
t=-T/2:1/fs:T/2-1/fs;
f=(-N/2:N/2-1)*B/N;
%% LFM信号
K=B/T;
s0=exp(j*2*pi*fc+j*pi*K*t.^2);


%% NLFM信号
Tf=T*f/B+(0.426*T/pi).*sin(2*pi*f/B);%海明窗
% Tf=T*f/B+(T/(2*pi))*sin(2*pi*f/B);%汉宁窗
% Tf=T*f/B+((2*T)/(3*pi))*sin(2*pi*f/B)+(T/(12*pi))*sin(4*pi*f/B);%余弦四次方窗
% figure(1) 
% subplot(211) 
% plot(f,Tf);grid on 
% title('群时延函数Tf/f') 
% xlabel('frequency/hz')
% ylabel('time/s')
f_t=interp1(Tf,f,t,'spline'); %三次样条插值法
% % [p,~,mu]=polyfit(Tf,f,20);  %多项式拟合
% % f_t=polyval(p,t,[],mu);
% subplot(212) 
% plot(t,f_t);
% grid on ;
% title('时频函数f/Tf') 
% xlabel('time/s')
% ylabel('frequency/hz')
% te=t(1); 
sumi=0; 
%数值积分
for i=1:N 
%     dt=t(i)-te; 
%     te=t(i); 
    f_t(i)=f_t(i)*(1/fs)+sumi; 
    sumi=f_t(i); 
    phi(i)=2*pi*f_t(i);%相位
end 
% phi=2*pi*f_t;%相位


s=exp(1j*phi);%S形NLFM信号
 Z = cmplxambiguity('true', 1024, 10,s,s);
    
end