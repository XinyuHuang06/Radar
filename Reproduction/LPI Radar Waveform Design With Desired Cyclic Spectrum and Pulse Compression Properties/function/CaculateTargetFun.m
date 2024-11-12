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

function [TarRecord] = CaculateTargetFun(DataSet, Data)
    % Unpacket Data
    FNr = Data.FNr;
    flag_Sparse = Data.FlagSparse;
    chi_matrix = Data.CHIMatrix;
    Taf_1 = Data.Taf1;
    Taf_2 = Data.Taf2;
    omega_alpha_1 = Data.OmegaAlpha;
    % Unpacket DataSet
    xr = DataSet.xr;
    br = DataSet.br;
    cr = DataSet.cr;
    lambda_0 = DataSet.lambda_0;
    lambda_1 = DataSet.lambda_1;
    rho_0 = DataSet.rho_0;
    rho_1 = DataSet.rho_1;
    h = DataSet.h;
    vartheta = DataSet.vartheta;


    TargetFun = Caculate_TarFunCValue(FNr, xr, Taf_1, Taf_2, flag_Sparse, omega_alpha_1);
    TargetLagrange_1 = lambda_0'*(xr-br);
    TargetPenalty_1 = rho_0/2*(norm(xr-br))^2;
    if flag_Sparse
        Tar_4_temp = cellfun(@(lambda,chi) lambda*((xr-cr)'*chi*(br-cr)-vartheta-h), num2cell(lambda_1), chi_matrix,"UniformOutput",false);
        Tar_5_temp = cellfun(@(rho,chi) 1/2*rho*(norm((xr-cr)'*chi*(br-cr)-vartheta-h)^2), num2cell(rho_1),chi_matrix,"UniformOutput",false);
        TargetLagrange_2 = sum(cell2mat(Tar_4_temp));
        TargetPenalty_2 = sum(cell2mat(Tar_5_temp));
    else
        Tar_45_temp = pagemtimes((pagemtimes(xr-cr,"transpose",chi_matrix,"none")),"none",br-cr,"none")-vartheta-h;
        Tar_45_temp = reshape(Tar_45_temp,N,1);
        TargetLagrange_2 = sum(Tar_45_temp);
        TargetPenalty_2 = 0.5*sum(Tar_45_temp.^2);
    end
    SUM = TargetFun + TargetLagrange_1 + TargetLagrange_2 + TargetPenalty_1 + TargetPenalty_2;
    TarRecord.SUM = SUM;
    TarRecord.TargetFun = TargetFun;
    TarRecord.TargetLagrange_1=TargetLagrange_1;
    TarRecord.TargetLagrange_2=TargetLagrange_2;
    TarRecord.TargetPenalty_1=TargetPenalty_1;
    TarRecord.TargetPenalty_2=TargetPenalty_2;
end

function TarFunCValue = Caculate_TarFunCValue(FNr, xr, Taf_1, Taf_2, flag_Sparse, omega_alpha_1)
    TarFunCValue = 0;
    if flag_Sparse
        for iter_i = 1:length(omega_alpha_1)
            TarFunCValue = TarFunCValue + omega_alpha_1(iter_i)*(((FNr*xr)'*Taf_1{iter_i}*FNr*xr)^2 + ((FNr*xr)'*Taf_2{iter_i}*FNr*xr)^2);
        end
    else
        TarFunCValue_temp = pagemtimes(FNr,"none",xr,"none");
        temp_1 = pagemtimes(pagemtimes(TarFunCValue_temp,"transpose",Taf_1,"none"),"none",TarFunCValue_temp,"none");
        temp_1 = temp_1.^2;
        temp_2 = pagemtimes(pagemtimes(TarFunCValue_temp,"transpose",Taf_2,"none"),"none",TarFunCValue_temp,"none");
        temp_2 = temp_2.^2;
        TarFunCValue = TarFunCValue + sum( bsxfun( @times,temp_1,reshape(omega_alpha_1, 1, 1, []) ), 3);
        TarFunCValue = TarFunCValue + sum( bsxfun( @times,temp_2,reshape(omega_alpha_1, 1, 1, []) ), 3);
    end
end