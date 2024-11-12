function [lambda_0, lambda_1] = Update_lambda(DataSet, Data)
    % Unpacket Data
    flag_Sparse = Data.FlagSparse;
    chi_matrix = Data.CHIMatrix;
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
    N = DataSet.N;

    lambda_0 = lambda_0 + rho_0*(xr-br);
    lambda_1 = zeros(size(lambda_1));
    if flag_Sparse
        for i_n = 1:N
            lambda_1(i_n) = lambda_1(i_n) + rho_1(i_n)*((xr-cr)'*chi_matrix{i_n}*(br-cr)-h-vartheta);
        end
    else
        lambda_1 = reshape(pagemtimes(xr-cr,"transpose",pagemtimes(chi_matrix,"none",br-cr,"none"),"none")-h-vartheta, [], 1).*rho_1 + lambda_1;
    end
end