% 整合后的复现代码
% 进行数据初始化

% % 01. 数据初始化
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
vartheta = 0.1;
DataSetPackets = DataSet(xr, br, cr, lambda_0, lambda_1, rho_0, rho_1, r, h, vartheta, N);
DataRecordPack.UpdateTarRecord(1, CaculateTargetFun(DataSetPackets.packets, ParameterPackets));
DataRecordPack.InitialConstantParameter(JsonPath);
save("data/parameter_01.mat", 'ParameterPackets', 'DataRecordPack', 'DataSetPackets');
vartheta = 0.2;
DataSetPackets = DataSet(xr, br, cr, lambda_0, lambda_1, rho_0, rho_1, r, h, vartheta, N);
DataRecordPack.UpdateTarRecord(1, CaculateTargetFun(DataSetPackets.packets, ParameterPackets));
DataRecordPack.InitialConstantParameter(JsonPath);
save("data/parameter_02.mat", 'ParameterPackets', 'DataRecordPack', 'DataSetPackets');
DataRecordPack.UpdateTarRecord(1, CaculateTargetFun(DataSetPackets.packets, ParameterPackets));
DataRecordPack.InitialConstantParameter(JsonPath);
vartheta = 0.3;
DataSetPackets = DataSet(xr, br, cr, lambda_0, lambda_1, rho_0, rho_1, r, h, vartheta, N);
save("data/parameter_03.mat", 'ParameterPackets', 'DataRecordPack', 'DataSetPackets');



% % 00. 函数及其它函数定义
function ParameterPackets = GenerateParameter(N, M, Flag_Sparse, threshhold)
    % FNr
    FN = dftmtx(N); 
    Fshift = zeros(N);
    Fshift(1:N/2,N/2+1:end) = eye(N/2);
    Fshift(N/2+1:end,1:N/2) = eye(N/2);
    FN = Fshift*FN;
    FNr = complex2real(FN);
    % omega_alpha_1
    omega_alpha_1 = ones((N+1)*N/2,1);
    if Flag_Sparse
        chi_matrix = cell(N,1);
        for i_temp = 1:N
            matrix_temp = zeros(N,N);
            matrix_temp(i_temp,i_temp) = 1;
            chi_matrix{i_temp} = sparse(complex2real(matrix_temp));
        end
        % 循环谱计算矩阵 20240410 Modified   
        temp_m = -M/2+1:M/2;
        Taf_1 = cell((N+1)*N/2,1);
        Taf_2 = cell((N+1)*N/2,1);
        Taf = cell((N+1)*N/2,1);
        i_temp = 0;
        for temp_1 = 1:N
            for temp_2 = temp_1:N
                i_temp = i_temp + 1;
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
                end
                temp_taf_1 = [real(temp_taf),-imag(temp_taf);imag(temp_taf),real(temp_taf)];
                temp_taf_2 = [imag(temp_taf),-real(temp_taf);real(temp_taf),imag(temp_taf)];
                Taf_1{i_temp} = sparse(temp_taf_1);
                Taf_2{i_temp} = sparse(temp_taf_2);
                Taf{i_temp} = sparse(temp_taf);
                if abs(N + 1 - temp_1-temp_2) < threshhold
                    omega_alpha_1(i_temp) = 0.01;
                else
                    omega_alpha_1(i_temp) = 10;
                end
            end
        end    
    else
        chi_matrix = zeros(2*N,2*N,N); % 相似性约束矩阵
        for i_temp = 1:N
            matrix_temp = zeros(N,N);
            matrix_temp(i_temp,i_temp) = 1;
            chi_matrix(:,:,i_temp) = complex2real(matrix_temp);
        end
        Taf_1 = zeros(2*N,2*N,(N+1)*N/2);
        Taf_2 = zeros(2*N,2*N,(N+1)*N/2);
        Taf = zeros(N,N,(N+1)*N/2);
        temp_m = -M/2+1:M/2;
        i_temp = 0;
        for temp_1 = 1:N
            for temp_2 = temp_1:N
                i_temp = i_temp + 1;
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
                end
                temp_taf_1 = [real(temp_taf),-imag(temp_taf);imag(temp_taf),real(temp_taf)];
                temp_taf_2 = [imag(temp_taf),-real(temp_taf);real(temp_taf),imag(temp_taf)];
                Taf_1(:,:,i_temp) = sparse(temp_taf_1);
                Taf_2(:,:,i_temp) = sparse(temp_taf_2);
                Taf(:,:,i_temp) = temp_taf;
                if abs(temp_1-temp_2) < threshhold
                    omega_alpha_1(i_temp) = 0.01;
                else
                    omega_alpha_1(i_temp) = 10;
                end
            end
        end  
    end
    % Intialize Data Struct
    ParameterPackets.CHIMatrix = chi_matrix;
    ParameterPackets.Taf1 = Taf_1;
    ParameterPackets.Taf2 = Taf_2;
    ParameterPackets.FNr  = FNr;
    ParameterPackets.OmegaAlpha = omega_alpha_1;
    ParameterPackets.Taf  = Taf;
    ParameterPackets.FlagSparse = Flag_Sparse;
    ParameterPackets.N = N;
    ParameterPackets.FN = FN;
end

% % 