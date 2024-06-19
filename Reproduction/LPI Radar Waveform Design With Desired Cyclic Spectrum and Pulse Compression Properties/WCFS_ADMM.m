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
    % % ADMM Iterations
    for i_m = 1:Max_ItersNum
        % % ADMM update
        xr = Update_xr(DataSetPackets.packets, ParameterPackets);% % Step 1 , Solving the x_r. 
        DataSetPackets.update(xr,'xr');
        br = Update_br(DataSetPackets.packets, ParameterPackets);% % Step 2 , Solving the b_r
        DataSetPackets.update(br,'br');
        h = Update_h(DataSetPackets.packets, ParameterPackets);% % Step 3 , Solving the h
        DataSetPackets.update(h,'h');
        [lambda_0, lambda_1] = Update_lambda(DataSetPackets.packets, ParameterPackets);% % Step 4 , Solving the lambda_0 and lambda_1
        DataSetPackets.update(lambda_0,'lambda_0',lambda_1,'lambda_1');

        [rho_0, rho_1] = Update_rho(DataSetPackets, ParameterPackets);% % Step 4 , Solving the rho_0 and rho_1
        DataSetPackets.update(rho_0, 'rho_0', rho_1, 'rho_1');
        % % Other
        DataRecordPack.UpdateTarRecord(CaculateTargetFun(DataSetPackets.packets, ParameterPackets), i_m + 1);
        DataRecordPack.UpdateParaRecord(xr, 'xr', rho_0, 'rho_0', rho_1, 'rho_1', lambda_0, 'lambda_0', lambda_1, 'lambda_1');
    end
    PlotAndExport(DataSetPackets.packets, DataRecordPack, InitialParameter, OutPath);
    save(strcat(OutPath,"/","Record.mat"), "DataRecordPack");
end