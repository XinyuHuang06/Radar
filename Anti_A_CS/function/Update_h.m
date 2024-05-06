function h = Update_h(N, xr, br, cr, rho_1, lambda_1, vartheta, r, chi_matrix, flag_Sparse)
    A3 = sum(rho_1)/2;
    BT3_1 = -sum(lambda_1);
    if flag_Sparse
        BT3_2 = 0;
        for i_temp = 1:N
            BT3_2 = -rho_1(i_temp)*(xr-cr)'*chi_matrix{i_temp}*(br-cr);
        end
    else
        BT3_2_temp = pagemtimes(xr-cr,"transpose",pagemtimes(chi_matrix,"none",br-cr,"none"),"none");
        BT3_2 = sum(bsxfun(@times, BT3_2_temp, reshape(-rho_1, 1, 1, [])), 3);
    end
    BT3_3 = vartheta*sum(rho_1);
    BT3 = BT3_1 + BT3_2 + BT3_3;
    h = (-BT3-sqrt(BT3^2+8*A3*r))/(4*A3);
end