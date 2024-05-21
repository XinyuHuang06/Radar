function [rho_0,rho_1,delta_0,delta_1] = Update_rho(rho_0_k,rho_1_k,xr_k,br_k,h_k,xr,br,cr,vartheta,h,chi_matrix,flag_Sparse,xi_inc,xi_dec)
    delta_0 = pen1_der_rho0(xr,br) - pen1_der_rho0(xr_k,br_k);
    if delta_0 < 0
        rho_0 = rho_0_k * xi_inc;
    elseif delta_0 > 0
        rho_0 = rho_0_k / xi_dec;
    else
        rho_0 = rho_0_k;
    end
    pen_2 = pen2_der_rho1(xr,br,cr,vartheta,h,chi_matrix,flag_Sparse);
    pen_2_k = pen2_der_rho1(xr_k,br_k,cr,vartheta,h_k,chi_matrix,flag_Sparse);
    % % 更新策略1
    delta_1 = pen_2 - pen_2_k; 
    temp_multi = ones(size(rho_1_k));
    temp_multi(delta_1 < 0) = 1 * xi_inc;
    temp_multi(delta_1 > 0) = 1 / xi_dec;
    rho_1 = rho_1_k.*temp_multi;

    % % 更新策略2
    % delta_1 = sum(pen_2-pen_2_k);
    % if delta_1 < 0
    %     rho_1 = rho_1_k * xi_inc;
    % elseif delta_1 > 0
    %     rho_1 = rho_1_k / xi_dec;
    % else
    %     rho_1 = rho_1_k;
    % end
    % fid = fopen("output_files/rho_1.log","a+");
    % fprintf(fid, '%.2d\n',delta_1);
    % fclose(fid);
end

function [pen_1] = pen1_der_rho0(xr,br)
    pen_1 = 1/2*(norm(xr-br))^2;
end

function [pen_2] = pen2_der_rho1(xr,br,cr,vartheta,h,chi_matrix,flag_Sparse)
    N = length(xr)/2;
    if flag_Sparse
        pen_2 = cellfun(@(chi) 1/2*(((xr-cr)'*chi*(br-cr)-vartheta-h).^2),chi_matrix,"UniformOutput",false);
        pen_2 = cell2mat(pen_2);
    else
        pen_2 = pagemtimes((pagemtimes(xr-cr,"transpose",chi_matrix,"none")),"none",br-cr,"none")-vartheta-h;
        pen_2 = reshape(pen_2.^2,N,1);
    end
end