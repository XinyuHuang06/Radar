% Example:
% :param :
% :return :
% detailed description:
%------------------------------------------------------------------------------
% Created by: Xinyu Huang.
% On: 15/05/2024.
% Copyright (C) 2024 Xinyu Huang (learning_huang@163.com).
% All Rights Reserved.
% Unauthorized copying of this file, via any medium is strictly prohibited.
% Proprietary and confidential.
%------------------------------------------------------------------------------
% Jsonlab
function WCFS_ADMM(DataPath, OutPath)
    % % Intialization
    load(DataPath);
    Max_ItersNum = DataRecordPack.ConstantPara.Max_ItersNum;
    % % ADMM Iterations
    for i_m = 1:Max_ItersNum
        if i_m == 1
            DataRecordPack.UpdateTarRecord(1, CaculateTargetFun(DataSetPackets.packets, ParameterPackets));
            DataRecordPack.UpdateParaRecord(1, DataSetPackets.packets.xr, 'xr', DataSetPackets.packets.lambda_0, 'lambda_0', DataSetPackets.packets.lambda_1, 'lambda_1', DataSetPackets.packets.h, 'h', DataSetPackets.packets.rho_0, 'rho_0', DataSetPackets.packets.rho_1, 'rho_1');
        end
        % % ADMM update
        xr = Update_xr(DataSetPackets.packets, ParameterPackets);% % Step 1 , Solving the x_r. 
        DataSetPackets.update(xr,'xr');
        br = Update_br(DataSetPackets.packets, ParameterPackets);% % Step 2 , Solving the b_r
        DataSetPackets.update(br,'br');
        h = Update_h(DataSetPackets.packets, ParameterPackets);% % Step 3 , Solving the h
        DataSetPackets.update(h,'h');
        [lambda_0, lambda_1] = Update_lambda(DataSetPackets.packets, ParameterPackets);% % Step 4 , Solving the lambda_0 and lambda_1
        DataSetPackets.update(lambda_0,'lambda_0',lambda_1,'lambda_1');
        % [rho_0, rho_1] = Update_rho(DataSetPackets, ParameterPackets);% % Step 4 , Solving the rho_0 and rho_1
        % DataSetPackets.update(rho_0, 'rho_0');
        % DataSetPackets.update(rho_1, 'rho_1');
        % % Other
        DataRecordPack.UpdateTarRecord(i_m+1, CaculateTargetFun(DataSetPackets.packets, ParameterPackets));
        DataRecordPack.UpdateParaRecord(i_m+1, xr, 'xr');
        DataRecordPack.UpdateParaRecord(i_m+1, lambda_0, 'lambda_0', lambda_1, 'lambda_1');
        DataRecordPack.UpdateParaRecord(i_m+1, h, 'h');
        % DataRecordPack.UpdateParaRecord(i_m+1, rho_0, 'rho_0', rho_1, 'rho_1');
        % % Stop 
        if DataRecordPack.JudgeStop(i_m+1, 1e-6)
            break;
        end
    end
    if exist(OutPath,"dir") == 0
        mkdir(OutPath);
    end
    % If the parfor are used, don't directly use the save function.
    parsave(strcat(OutPath,"/","Record.mat"), DataRecordPack);
    PlotAndExport(DataSetPackets.packets, DataRecordPack, OutPath);
end