clc; clear; close all;
% 导入函数定义
addpath("./function");
% 参数初始化
flag = 1;
if flag
    % 进行参数初始化
    Parameter = Initialization_Parameter();
    save('./data/Parameter.mat','Parameter');
else
    % 加载参数文件
    load(".\data\Parameter.mat");
end
% 生成信号
m = 100;
n = 100;
p = (1:600)';
num_tar = 1;
value = smlk(m, n, p, num_tar, Parameter);
value2 = value.*exp(1j*2*pi*10*1e6*(p-1)*Parameter.TS);
figure
plot(real(value)); hold on;
plot(real(value2)); hold off;
figure
Analysis_ES((value), Parameter.fS, round(Parameter.TW*Parameter.fS),'normalized', 0);hold on;
Analysis_ES((value2), Parameter.fS, round(Parameter.TW*Parameter.fS),'normalized', 0);hold off;
legend(["m=n","m\neqn"])
figure
lf_value = conv(Parameter.lf_h,value);
Analysis_ES(lf_value, Parameter.fS, round(Parameter.TW*Parameter.fS),'normalized', 0);hold on;
lf_value2 = conv(Parameter.lf_h,value2);
Analysis_ES(lf_value2, Parameter.fS, round(Parameter.TW*Parameter.fS),'normalized', 0);hold off;
% 编码更新
% 对于相位编码信号，不太可能存在那样的频谱...
N1 = min(Parameter.N, Parameter.M-1-n);
N2 = min(Parameter.N, n);
Rm = Rm(N1, N2, m, Parameter);
Update_pc(m, Rm, 1e-6, Parameter);
value = smlk(m, n, p, num_tar, Parameter);
figure
Analysis_ES(value, Parameter.fS, round(Parameter.TW*Parameter.fS),'normalized', 0);
title('更新后')
% 移除函数定义
rmpath(".\function");