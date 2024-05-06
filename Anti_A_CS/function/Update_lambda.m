function [lambda_0, lambda_1] = Update_lambda(N, xr, br, cr, rho_0, rho_1, lambda_0_in, lambda_1_in, h, vartheta, chi_matrix, flag_Sparse)
    lambda_0 = lambda_0_in + rho_0*(xr-br);
    lambda_1 = zeros(size(lambda_1_in));
    if flag_Sparse
        for i_n = 1:N
            lambda_1(i_n) = lambda_1_in(i_n) + rho_1(i_n)*((xr-cr)'*chi_matrix{i_n}*(br-cr)-h-vartheta);
        end
    else
        lambda_1 = reshape(pagemtimes(xr-cr,"transpose",pagemtimes(chi_matrix,"none",br-cr,"none"),"none")-h-vartheta, [], 1).*rho_1 + lambda_1_in;
    end
end