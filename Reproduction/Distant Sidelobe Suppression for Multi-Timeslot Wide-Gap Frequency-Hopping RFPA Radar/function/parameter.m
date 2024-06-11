function parapack = parameter()
    parapack.f_s = 128*1e6; % 采样频率: 128 MHz
    parapack.T_W = 1e-6; % 脉冲宽度: 1 us
    parapack.M = 1024; % 传输脉冲总数: 1024
    parapack.T_R = 10e-6; % 平均PRI间隔: 10 us
    parapack.T_J = 6e-6; % 最大PRI捷变值: 6 us
    parapack.f_c = 3e9; % 基础载波频率: 3 GHz
    parapack.B = 0.5*f_c; % 载波捷变范围: 64 MHz
    parapack.G = 10e6; % 载频捷变时，每个脉冲间的最小频率间隔 10MHz
    parapack.ChipRate = 128e6; % Chip Rate 使用相位编码调制的码率
    parapack.T_C = 1/ChipRate; 
    parapack.M_W = parapack.T_W/parapack.T_C; % 每个脉冲对应的码片数
    parapack.L = 1000; % 最大距离探测单元数
    parapack.K = 100; % 最大速度探测单元数
    parapack.Delta_r = 10; % 距离探测精度
    parapack.Delta_v = 10; % 速度探测精度
    parapack.c = 3*1e8; % 光速
    parapack.phi_0 = 0; % 初始相位
    parapack.T_S = 1/(2*parapack.B); % 采样间隔
    parapack.Num_CP = 10; % 相干处理脉冲数目
    parapack.T_M = 2*parapack.Num_CP*parapack.T_R*(parapack.K*parapack.Delta_v)/parapack.c; % CPI内最大距离徙动对应的时延
    parapack.NumPRI = 192;
    parapack.NumFre = 2048;
    % 目标参数
    parapack.TargetCoef = zeros(parapack.L,parapack.M); % 目标散射系数矩阵
    parapack.TargetV = zeros(parapack.L,parapack.M); % 目标径向速度矩阵
    parapack.TargetLocation = zeros(parapack.L,parapack.M); % 目标径向距离矩阵
    parapack.N = floor( (2*parapack.L*parapack.Delta_r/parapack.c + parapack.T_J + parapack.T_W + 2*parapack.T_M)/parapack.T_R ); % 保持宽间隔需求所需的时隙数目
    % % 发射波形初始化
    parapack.T_m_Seq = rand(parapack.M,1)*parapack.T_J;            % 发射脉冲起始时间偏移量[0-T_J], 均匀分布
    parapack.f_m_Seq = parapack.f_c ; % 发射脉冲载波脉冲偏移量f_c + [-B/2,B/2], 高斯分布，在这里，载波偏移量需要满足载波频率间隔约束
    parapack.SamplePoints = ceil(parapack.T_W*parapack.f_s); % 时宽采样点数
end