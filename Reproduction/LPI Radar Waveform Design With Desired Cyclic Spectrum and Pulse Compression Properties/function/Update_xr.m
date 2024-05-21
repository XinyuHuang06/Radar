function xr_out = Update_xr(DataSet, Data)
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
    N = DataSet.N;
    % Update xr
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
        
        BT1 = BT1_1 + BT1_2 + BT1_3 + BT1_4 + BT1_5;
        % BT1 = BT1_2 + BT1_3 + BT1_4 + BT1_5;
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
        BT1 = BT1_1 + BT1_2 + BT1_3 + BT1_4 + BT1_5;
        % BT1 = BT1_2 + BT1_3 + BT1_4 + BT1_5;
    end
    % Constant Modules Constraints 1/2*xT*H*x + k*x + d = 0 => xT*(2*I)*x + 0*x -1 = 0
    N = length(xr)/2;
    temp_H = cell(1);temp_k= cell(1);temp_d= cell(1);
    temp_H{1} = eye(2*N)*2;temp_k{1} = zeros(2*N,1);temp_d{1} = -N;
    options = optimoptions(@fmincon,'Algorithm','interior-point',...
        'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
        'HessianFcn',@(x,lambda)quadhess(lambda,A1*2,temp_H),'Display','off');
    [xr_out,~,~,~] = fmincon(@(x) quadobj(x,A1*2,BT1',0), xr,[],[],[],[],[],[],...
        @(x) quadconstr(x,temp_H,temp_k,temp_d),options);
    % Caculate the new xr
    % xr_out = (A1 + 2*temp_lambda)/B;
    % xr'*A1*xr + BT1*xr
    % xr_out'*A1*xr_out + BT1*xr_out
end


function [y,grady] = quadobj(x,Q,f,c)
    y = 1/2*x'*Q*x + f'*x + c;
    if nargout > 1
        grady = Q*x + f;
    end
end

function [y,yeq,grady,gradyeq] = quadconstr(x,H,k,d)
    jj = length(H);  % jj is the number of equality constraints
    yeq = zeros(1,jj);
    for i = 1:jj
        yeq(i) = 1/2*x'*H{i}*x + k{i}'*x + d{i};
    end
    y = [];
    if nargout > 2
        gradyeq = zeros(length(x),jj);
        for i = 1:jj
            gradyeq(:,i) = H{i}*x + k{i};
        end
    end
    grady = [];
end

function hess = quadhess(lambda,Q,H)
    hess = Q;
    jj = length(H);  % jj is the number of equality constraints
    for i = 1:jj
        hess = hess + lambda.eqnonlin(i)*H{i};
    end
end