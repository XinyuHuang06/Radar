function [Tar_sum,Tar_1, Tar_2, Tar_3, Tar_4, Tar_5] = Debug_fun(N, FNr, xr, br, cr, omega_alpha_1, rho_0, rho_1, lambda_0, lambda_1, h, vartheta, Taf_1, Taf_2, chi_matrix, flag_Sparse)
    Tar_1 = Caculate_TarFunCValue(FNr, xr, Taf_1, Taf_2, flag_Sparse, omega_alpha_1);
    Tar_2 = lambda_0'*(xr-br);
    % Tar_2 = norm(xr-br);
    Tar_3 = rho_0/2*(norm(xr-br))^2;
    if flag_Sparse
        Tar_4_temp = cellfun(@(lambda,chi) lambda*((xr-cr)'*chi*(br-cr)-vartheta-h), num2cell(lambda_1), chi_matrix,"UniformOutput",false);
        Tar_5_temp = cellfun(@(rho,chi) 1/2*rho*(norm((xr-cr)'*chi*(br-cr)-vartheta-h)^2), num2cell(rho_1),chi_matrix,"UniformOutput",false);
        % Tar_4_temp = cellfun(@(chi) ((xr-cr)'*chi*(br-cr)-vartheta-h), chi_matrix,"UniformOutput",false);
        % Tar_5_temp = cellfun(@(chi) 1/2*(norm((xr-cr)'*chi*(br-cr)-vartheta-h)^2),chi_matrix,"UniformOutput",false);
        Tar_4 = sum(cell2mat(Tar_4_temp));
        Tar_5 = sum(cell2mat(Tar_5_temp));
    else
        Tar_45_temp = pagemtimes((pagemtimes(xr-cr,"transpose",chi_matrix,"none")),"none",br-cr,"none")-vartheta-h;
        Tar_45_temp = reshape(Tar_45_temp,N,1);
        % Tar_4 = lambda_1'*Tar_45_temp;
        % Tar_5 = (rho_1/2)'*(Tar_45_temp.^2);
        Tar_4 = sum(Tar_45_temp);
        Tar_5 = 0.5*sum(Tar_45_temp.^2);
    end
    Tar_sum = Tar_1 + Tar_2 + Tar_3 + Tar_4 + Tar_5;
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