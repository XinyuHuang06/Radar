% Example: Parameter = Initialization_Parameter()
% :param :
% :return :
% detailed description: 通用参数结构体设置，调用函数时可以使用此参数包传递通用参数，参数设置在此函数内进行
%------------------------------------------------------------------------------
% Created by: Xinyu Huang.
% On: 12/06/2024.
% Copyright (C) 2024 Xinyu Huang (learning_huang@163.com).
% All Rights Reserved.
% Unauthorized copying of this file, via any medium is strictly prohibited.
% Proprietary and confidential.
%------------------------------------------------------------------------------
function varargout = Initialization_Parameter()
    % 参数初始化
    fC = 3*1e9; % 中心频率 3 GHz
    M = 1024; % 脉冲数
    TW = 1*1e-6; % 脉冲宽度 1 us
    TR = 10*1e-6; % 平均PRI间隔 10 us
    TJ = 6*1e-6; % 最大PRI捷变范围 6 us
    PRI_PointsNum = 192; % PRI 捷变点数
    % PRI捷变分布 均匀分布
    B = 64*1e6; % 载频捷变范围 64 MHz
    Fre_PointsNum = 2048; % 载频捷变点数
    % 载频捷变分布 高斯分布
    N = 3; % 相邻脉冲数
    G = 10*1e6; % 最小频率间隔 10 MHz
    fS = 128*1e6; % 采样频率 128 MHz
    phi0 = 0; % 初始相位
    TS = 1/(fS); % 采样间隔 且有 TS <= 1/(2*B)
    % 探测相关参数
    Delta_r = 10; % 距离探测精度
    Delta_v = 10; % 速度探测精度
    max_l = 100; % 最大距离探测单元
    max_k = 100; % 最大速度探测单元
    % 其它参数
    c = 3*1e8; % 光速
    % 目标信息
    Target_num = 4;
    TargetMatrix = zeros(3,Target_num); % 三个维度指：所处距离单元，所处速度单元，散射系数
    TargetMatrix(:,1) = [10;10;0.01]; % 第一个目标，位于第(10,10)个单元上
    % 生成低通滤波器
    fpass = 1*1e6;  % 通带频率 1 MHz
    fstop = 2.5*1e6; % 阻带频率 2.5 MHz
    Rp = 1; % 通带波纹 1 dB
    As = 30; % 阻带衰减 30 dB
    lf_h = h(fpass,fstop,fS,Rp,As);
    % 生成波形随机参数
    fnSeq = generate_fnSeq(fC, G, B, N, M);
    tnSeq = generate_tnSeq(M,TJ,TR);
    % 相位编码相关参数
    flag_PC = 0; % 是否使用相位编码 bool type
    CR = 128*1e6; % 码元传输速率 128 MHz
    TC = 1/CR; % 单码元传输间隔
    MW = round(TW/TC); % 单脉冲码元传输个数
    PhaseCode = exp(1j*2*pi*rand(MW,1)); % 初始相位编码序列
    % 初始化结构体
    Parameter.fC = fC;
    Parameter.M = M;
    Parameter.TW = TW;
    Parameter.TR = TR;
    Parameter.TJ = TJ;
    Parameter.PRI_PointsNum = PRI_PointsNum;
    Parameter.B = B;
    Parameter.Fre_PointsNum = Fre_PointsNum;
    Parameter.N = N;
    Parameter.G = G;
    Parameter.phi0 = phi0;
    Parameter.TS = TS;
    Parameter.Delta_r = Delta_r;
    Parameter.Delta_v = Delta_v;
    Parameter.max_l = max_l;
    Parameter.max_k = max_k;
    Parameter.c = c;
    Parameter.TargetMatrix = TargetMatrix;
    Parameter.fpass = fpass;
    Parameter.fstop = fstop;
    Parameter.Rp = Rp;
    Parameter.As = As;
    Parameter.lf_h = lf_h;
    Parameter.fnSeq = fnSeq;
    Parameter.tnSeq = tnSeq;
    % 相位编码信号
    Parameter.flag_PC = flag_PC;
    PhaseCodeParameter.CR = CR;
    PhaseCodeParameter.TC = TC;
    PhaseCodeParameter.MW = MW;
    PhaseCodeParameter.PhaseCode = PhaseCode;
    Parameter.PhaseCodeParameter = PhaseCodeParameter;
    varargout{1} = Parameter;
end