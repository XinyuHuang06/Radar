addpath("function/");
% Load data.
InitialParameter = loadjson('./Data/initial_Parameter.json');
N = round(InitialParameter.signal.T*InitialParameter.signal.fs);
Max_ItersNum = InitialParameter.Max_ItersNum;
Threshold = InitialParameter.Threshold;
% Generate the parameter packets
ParameterPackets = GenerateParameter(N, InitialParameter.signal.M, InitialParameter.Flag.Sparse, InitialParameter.Threshold);
% % Generate xr, br, cr.
% x = exp(1j*2*pi*rand(N,1));
% xr = [real(x);imag(x)]; pwd
% b = exp(1j*2*pi*rand(N,1));
% br = [real(b);imag(b)]; 
c = generator_LFM(InitialParameter.signal.fs,InitialParameter.signal.fc,InitialParameter.signal.B,InitialParameter.signal.T);
cr = [real(c);imag(c)];
xr = cr;
br = cr;
% % Intialize the parameter of ADMM.
% The coef of Lagrange terms.
lambda_0 = 0*ones(2*N,1); 
lambda_1 = 0*ones(N,1); 
% The coef of penalty terms.
rho_0 = 1;
rho_1 = 1*ones(N,1);
% r: The barrier function parameter h: The 
r = 0.01;
h = -0.01;
% The similarty parameter.
vartheta = 0.3;
% % Data Record.
DataRecordPack = DataRecord(Max_ItersNum);
% The initial parameter Caculation
DataSetPackets = DataSet(xr, br, cr, lambda_0, lambda_1, rho_0, rho_1, r, h, vartheta, N);
DataRecordPack.UpdateTarRecord(CaculateTargetFun(DataSetPackets.packets, ParameterPackets), 1);
% % Save the mat file
save("data/parameter_03.mat", 'ParameterPackets', 'DataRecordPack', 'DataSetPackets', 'InitialParameter');
rmpath('function/');