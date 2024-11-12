function [A1,BT1,Tar_seq] = Caculate_matrix_xr(xr_in, br, cr, rho_0, rho_1, lambda_0, lambda_1, h, vartheta, Data)
    omega_alpha_1 = Data.omega_alpha;
    Taf_1 = Data.Taf_1;
    Taf_2 = Data.Taf_2;
    chi_matrix = Data.chi_matrix;
    FNr = Data.FNr;
    flag_Sparse = Data.sparse;
    N = Data.N;
    if flag_Sparse
        % % A1
        XA1_1_1 = cellfun( @(x) FNr'*(x)'*(FNr*br), Taf_1, "UniformOutput",false);
        XA1_1_2 = cellfun( @(x) FNr'*(x)'*(FNr*br), Taf_2, "UniformOutput",false);
        A1_1_temp = (cellfun( @(x1,x2,omega_alpha) omega_alpha*(x1*x1'+x2*x2'), XA1_1_1, XA1_1_2, num2cell(omega_alpha_1), "UniformOutput",false))';
        A1_1_temp = reshape([A1_1_temp{:}],2*N, 2*N, (N+1)*N/2);
        A1_1 = sum(A1_1_temp, 3);
        XA1_2 = cellfun( @(chi,rho) rho/2*chi*(br-cr)*(br-cr)'*chi', chi_matrix, num2cell(rho_1), "UniformOutput",false);
        XA1_2 = XA1_2';
        A1_2_temp = reshape( [XA1_2{:}], 2*N, 2*N, N);
        A1_2 = sum(A1_2_temp, 3);
        A1_3 = rho_0/2*eye(2*N,2*N);
        A1 = A1_1 + A1_2 + A1_3;
        % % B1
        BT1_1 = lambda_0'; % Caculate the matrix B_1^T
        BT1_2 = -rho_0*br';
        XBT1_3 = cellfun(@(lambda,chi) lambda*(br-cr)'*chi', num2cell(lambda_1), chi_matrix, "UniformOutput", false);
        BT1_3 = (sum(reshape([XBT1_3{:}]',2*N,N), 2)');
        XBT1_4 = cellfun( @(rho,chi) -rho*cr'*(chi*(br-cr)*(br-cr)'*chi'), num2cell(rho_1), chi_matrix, "UniformOutput", false);
        BT1_4 = (sum(reshape([XBT1_4{:}]',2*N,N), 2)');
        XBT1_5 = cellfun(@(rho,chi) -(vartheta+h)*rho*(br-cr)'*chi', num2cell(rho_1), chi_matrix,"UniformOutput",false);
        BT1_5 =  (sum(reshape([XBT1_5{:}]',2*N,N), 2)');
        
        % BT1 = BT1_1 + BT1_2 + BT1_3 + BT1_4 + BT1_5;
        BT1 = BT1_2 + BT1_3 + BT1_4 + BT1_5;
        % % Target Value Test
        % C1
        C1_1 = -lambda_0'*br;
        C1_2 = rho_0/2*(br'*br);
        C1_3 = sum(cellfun(@(lambda,chi) lambda*(-cr'*chi*(br-cr)-vartheta-h),num2cell(lambda_1),chi_matrix));
        C1_4 = sum(cellfun(@(rho,chi) rho/2*((vartheta+h)^2+cr'*chi*(br-cr)*(br-cr)'*chi'*cr),num2cell(rho_1),chi_matrix));
        C1 = C1_1 + C1_2 + C1_3 + C1_4;
    else
        A1_temp1 = pagemtimes(pagemtimes(FNr,"transpose",Taf_1,"transpose"),FNr*br);
        A1_temp1 = pagemtimes(A1_temp1,"none",A1_temp1,"transpose");
        A1_11 = sum(bsxfun(@times, A1_temp1, reshape(omega_alpha_1, 1, 1, [])),3); 
        A1_temp2 = pagemtimes(pagemtimes(FNr,"transpose",Taf_2,"transpose"),FNr*br);
        A1_temp2 = pagemtimes(A1_temp2,"none",A1_temp2,"transpose");
        A1_12 = sum(bsxfun(@times, A1_temp2, reshape(omega_alpha_1, 1, 1, [])),3);
        A1_1 = A1_11 + A1_12;
        A1_temp3 = pagemtimes(chi_matrix,"transpose",br-cr,"none");
        A1_temp3 = pagemtimes(A1_temp3,"none",A1_temp3,"transpose");
        A1_2 = sum(bsxfun(@times, A1_temp3, reshape(rho_1/2, 1, 1, [])),3);
        A1_3 = rho_0/2*eye(2*N,2*N);
        A1 = A1_1 + A1_2 + A1_3;
        BT1_1 = lambda_0'; % Caculate the matrix B_1^T
        BT1_2 = -rho_0*br';
        BT1_3_temp1 = pagemtimes(br-cr,"transpose",chi_matrix,"transpose");
        BT1_3 = sum(bsxfun(@times, BT1_3_temp1, reshape(lambda_1, 1, 1, [])), 3);
        BT1_4_temp1 = pagemtimes(chi_matrix,"none",br-cr,"none");
        BT1_4_temp1 = pagemtimes(cr,"transpose",pagemtimes(BT1_4_temp1,"none",BT1_4_temp1,"transpose"),"none");
        BT1_4 = sum(bsxfun(@times, BT1_4_temp1, reshape(-rho_1, 1, 1, [])), 3);
        BT1_5_temp1 = pagemtimes(br-cr,"transpose",chi_matrix,"transpose");
        BT1_5 = -(vartheta+h)*sum(bsxfun(@times, BT1_5_temp1, reshape(rho_1, 1, 1, [])), 3);
        % BT1 = BT1_1 + BT1_2 + BT1_3 + BT1_4 + BT1_5;
        BT1 = BT1_2 + BT1_3 + BT1_4 + BT1_5;
    end
    Tar_seq = zeros(1,6);
    Tar_seq(1) = xr_in'*A1*xr_in+BT1*xr_in+C1;
    Tar_seq(2) = xr_in'*A1_1*xr_in;
    Tar_seq(3) = BT1_1*xr_in + C1_1;
    Tar_seq(4) = xr_in'*A1_3*xr_in + BT1_2*xr_in + C1_2;
    Tar_seq(5) = C1_3 + BT1_3*xr_in;
    Tar_seq(6) = C1_4 + xr_in'*A1_2*xr_in + (BT1_4+BT1_5)*xr_in;
end