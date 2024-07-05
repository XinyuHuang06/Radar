% Load data.
JsonPath = './data/initial_Parameter.json';
InitialParameter = loadjson(JsonPath);
Max_ItersNum = InitialParameter.Max_ItersNum;
N = round(InitialParameter.signal.T*InitialParameter.signal.fs);
% Generate the parameter packets
ParameterPackets = GenerateParameter(N, InitialParameter.signal.M, InitialParameter.Flag.Sparse, InitialParameter.Threshold);
% % Generate xr, br, cr.
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
r = 0.001;
h = -0.001;
% The similarty parameter.



% % Data Record.
DataRecordPack = DataRecord(Max_ItersNum);
% The initial parameter Caculation
% % Save the mat file

vartheta = 1;
DataSetPackets = DataSet(xr, br, cr, lambda_0, lambda_1, rho_0, rho_1, r, h, vartheta, N);
DataRecordPack.UpdateTarRecord(1, CaculateTargetFun(DataSetPackets.packets, ParameterPackets));
DataRecordPack.InitialConstantParameter(JsonPath);
save("data/parameter_01.mat", 'ParameterPackets', 'DataRecordPack', 'DataSetPackets');
vartheta = 1;
DataSetPackets = DataSet(xr, br, cr, lambda_0, lambda_1, rho_0, rho_1, r, h, vartheta, N);
DataRecordPack.UpdateTarRecord(1, CaculateTargetFun(DataSetPackets.packets, ParameterPackets));
DataRecordPack.InitialConstantParameter(JsonPath);
save("data/parameter_02.mat", 'ParameterPackets', 'DataRecordPack', 'DataSetPackets');
DataRecordPack.UpdateTarRecord(1, CaculateTargetFun(DataSetPackets.packets, ParameterPackets));
DataRecordPack.InitialConstantParameter(JsonPath);
vartheta = 1;
DataSetPackets = DataSet(xr, br, cr, lambda_0, lambda_1, rho_0, rho_1, r, h, vartheta, N);
save("data/parameter_03.mat", 'ParameterPackets', 'DataRecordPack', 'DataSetPackets');