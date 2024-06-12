clc; clear; close all;
% 导入函数定义
addpath(".\function");
% 参数初始化
flag = 1;
if flag
    % 进行参数初始化
    Parameter = Initialization_Parameter();
    save('./data/Parameter.mat','Parameter');
else
    % 加载参数文件
    load("data\Parameter.mat");
end
% 生成信号
m = 100;
n = 100;
p = 1:128;
num_tar = 1;
value = smlk(m, n, p, num_tar, Parameter);
plot(real(value))
% 移除函数定义
rmpath(".\function");