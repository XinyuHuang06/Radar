function [br] = Update_br(N, FNr, xr, br_in, cr, omega_alpha_1, rho_0, rho_1, lambda_0, lambda_1, h, vartheta, Taf_1, Taf_2, chi_matrix, flag_Sparse)
    if flag_Sparse
        XA2_1_1 = cellfun( @(x) FNr'*(x)'*(FNr*xr), Taf_1, "UniformOutput",false);
        XA2_1_2 = cellfun( @(x) FNr'*(x)'*(FNr*xr), Taf_2, "UniformOutput",false);
        A2_1_temp = (cellfun( @(x1,x2,omega_alpha) omega_alpha*(x1*x1'+x2*x2'), XA2_1_1, XA2_1_2, num2cell(omega_alpha_1), "UniformOutput",false))';
        A2_1_temp = reshape([A2_1_temp{:}],2*N, 2*N, (N+1)*N/2);
        A2_1 = sum(A2_1_temp, 3);
        XA2_2 = cellfun( @(chi,rho) rho/2*chi*(xr-cr)*(xr-cr)'*chi', chi_matrix, num2cell(rho_1), "UniformOutput",false); % 2024/04/23 修正1/2系数
        A2_2_temp = reshape( [XA2_2{:}], 2*N, 2*N, N);
        A2_2 = sum(A2_2_temp, 3);
        A2_3 = rho_0/2*eye(2*N,2*N);
        A2 = A2_1 + A2_2 + A2_3;
        BT2_1 = -lambda_0';
        BT2_2 = -rho_0*xr';
        XBT2_3 = cellfun( @(lambda,chi) lambda*(xr-cr)'*chi, num2cell(lambda_1), chi_matrix, "UniformOutput", false);
        BT2_3 = (sum(reshape([XBT2_3{:}]',2*N,N), 2)');
        XBT2_4 = cellfun(@(rho,chi) -rho*cr'*((chi*(xr-cr))*(xr-cr)'*chi'), num2cell(rho_1), chi_matrix, "UniformOutput", false);
        BT2_4 = (sum(reshape([XBT2_4{:}]',2*N,N), 2)');
        XBT2_5 = cellfun(@(rho,chi) -(vartheta+h)*rho*(xr-cr)'*chi, num2cell(rho_1), chi_matrix,"UniformOutput",false);
        BT2_5 = (sum(reshape([XBT2_5{:}]',2*N,N), 2)');
        % BT2 = BT2_1 + BT2_2 + BT2_3 + BT2_4 + BT2_5;
        BT2 = BT2_2 + BT2_3 + BT2_4 + BT2_5;
    else
        A2_temp1 = pagemtimes(pagemtimes(FNr,"transpose",Taf_1,"transpose"),FNr*xr); % 子项1.1
        A2_temp1 = pagemtimes(A2_temp1,"none",A2_temp1,"transpose");
        A2_11 = sum(bsxfun(@times, A2_temp1, reshape(omega_alpha_1, 1, 1, [])),3); % 加权求和
        A2_temp2 = pagemtimes(pagemtimes(FNr,"transpose",Taf_2,"transpose"),FNr*xr); % 子项1.2
        A2_temp2 = pagemtimes(A2_temp2,"none",A2_temp2,"transpose");
        A2_12 = sum(bsxfun(@times, A2_temp2, reshape(omega_alpha_1, 1, 1, [])),3); % 加权求和
        A2_1 = A2_11 + A2_12;
        A2_temp3 = pagemtimes(chi_matrix,"transpose",xr-cr,"none");
        A2_temp3 = pagemtimes(A2_temp3,"none",A2_temp3,"transpose");
        A2_2 = sum(bsxfun(@times, A2_temp3, reshape(rho_1/2, 1, 1, [])),3);
        A2_3 = rho_0/2*eye(2*N,2*N);
        A2 = A2_1 + A2_2 + A2_3;
        BT2_1 = -lambda_0';
        BT2_2 = -rho_0*xr';
        BT2_3_temp1 = pagemtimes(xr-cr,"transpose",chi_matrix,"none");
        BT2_3 = sum(bsxfun(@times,BT2_3_temp1,reshape(lambda_1,1,1,[])),3);
        BT2_4_temp1 = pagemtimes(chi_matrix,"none",xr-cr,"none");
        BT2_4_temp1 = pagemtimes(cr,"transpose",pagemtimes(BT2_4_temp1,"none",BT2_4_temp1,"transpose"),"none");
        BT2_4 = sum(bsxfun(@times, BT2_4_temp1, reshape(-rho_1, 1, 1, [])), 3);
        BT2_5_temp1 = pagemtimes(xr-cr,"transpose",chi_matrix,"none");
        BT2_5 = -(vartheta+h)*sum(bsxfun(@times, BT2_5_temp1, reshape(rho_1, 1, 1, [])), 3);
        % BT2 = BT2_1 + BT2_2 + BT2_3 + BT2_4 + BT2_5;
        BT2 = BT2_2 + BT2_3 + BT2_4 + BT2_5;
    end
    % br = -inv(2*A2)*BT2';

    temp_H = cell(1);temp_k= cell(1);temp_d= cell(1);
    temp_H{1} = eye(2*N)*2;temp_k{1} = zeros(2*N,1);temp_d{1} = -1;
    options = optimoptions(@fmincon,'Algorithm','interior-point',...
        'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
        'HessianFcn',@(x,lambda)quadhess(lambda,A2*2,temp_H),'Display','off');
    [br,~,~,~] = fmincon(@(x) quadobj(x,A2*2,BT2',0), br_in,[],[],[],[],[],[],...
        @(x) quadconstr(x,temp_H,temp_k,temp_d),options);
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