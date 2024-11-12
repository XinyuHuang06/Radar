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
        A1_1 = zeros(2*N, 2*N);
        % 如果GPU可用，切换为GPU计算模式，否则使用parfor并行计算
        flag_GPU = canUseGPU();
        if flag_GPU
            FNr = gpuArray(FNr);
            br  = gpuArray(br);
            Taf_1 = cellfun(@gpuArray, Taf_1, 'UniformOutput', false);
            % Pleasr don't use parfor, it will be more slower
            for k = 1:(N+1)*N/2
                XA1_1_1 = FNr'*(Taf_1{k})'*(FNr*br);
                XA1_1_2 = FNr'*(Taf_2{k})'*(FNr*br);
                A1_1 = A1_1 + omega_alpha_1(k)*(XA1_1_1*XA1_1_1' + XA1_1_2*XA1_1_2');
            end
            A1_1 = gather(A1_1);
            XA1_2 = cellfun( @(chi,rho) rho/2*chi*(br-cr)*(br-cr)'*chi', chi_matrix, num2cell(rho_1), "UniformOutput",false);
            A1_2  = zeros(2*N, 2*N);
            for i = 1:N
                A1_2 = A1_2 + XA1_2{i};
            end
            A1_3 = rho_0/2*eye(2*N,2*N);
            A1 = A1_1 + A1_2 + A1_3;
            A1 = gather(A1);
            BT1_1 = lambda_0'; % Caculate the matrix B_1^T
            BT1_2 = -rho_0*br';
            XBT1_3 = cell(N,1);
            XBT1_4 = cell(N,1);
            for i = 1:N
                XBT1_3{i} = lambda_1(i)*chi_matrix{i}*(br-cr);
                XBT1_4{i} = -rho_1(i)*chi_matrix{i}*(br-cr)*(br-cr)'*(chi_matrix{i})'*cr;
            end
            BT1_3 = zeros(1,2*N);
            BT1_4 = zeros(1,2*N);
            for i = 1:N
                BT1_3 = BT1_3 + (XBT1_3{i})';            
                BT1_4 = BT1_4 + (XBT1_4{i})';
            end
            BT1 = BT1_1 + BT1_2 + BT1_3 + BT1_4;
        else
            for k = 1:(N+1)*N/2
                XA1_1_1 = FNr'*(Taf_1{k})'*(FNr*br);
                XA1_1_2 = FNr'*(Taf_2{k})'*(FNr*br);
                A1_1 = A1_1 + omega_alpha_1(k)*(XA1_1_1*XA1_1_1' + XA1_1_2*XA1_1_2');
            end
            XA1_2 = cellfun( @(chi,rho) rho/2*chi*(br-cr)*(br-cr)'*chi', chi_matrix, num2cell(rho_1), "UniformOutput",false);
            A1_2  = zeros(2*N, 2*N);
            for i = 1:N
                A1_2 = A1_2 + XA1_2{i};
            end
            A1_3 = rho_0/2*eye(2*N,2*N);
            A1 = A1_1 + A1_2 + A1_3;
            % paper definition
            BT1_1 = lambda_0'; % Caculate the matrix B_1^T
            BT1_2 = -rho_0*br';
            XBT1_3 = cellfun( @(rho,chi) rho*chi*(br-cr), num2cell(lambda_1), chi_matrix,"UniformOutput",false);
            XBT1_4 = cellfun( @(rho,chi) -rho/2*chi*(br-cr)*(br-cr)'*chi'*cr, num2cell(rho_1), chi_matrix ,"UniformOutput",false);
            BT1_3 = zeros(1,2*N);
            BT1_4 = zeros(1,2*N);
            for i = 1:N
                BT1_3 = BT1_3 + (XBT1_3{i})';            
                BT1_4 = BT1_4 + (XBT1_4{i})';
            end
            BT1 = BT1_1 + BT1_2 + BT1_3 + BT1_4;
        end
        BT1 = gather(BT1);
        % % % B1 % my definition
        % BT1_1 = lambda_0'; % Caculate the matrix B_1^T
        % BT1_2 = -rho_0*br';
        % XBT1_3 = cellfun(@(lambda,chi) lambda*(br-cr)'*chi', num2cell(lambda_1), chi_matrix, "UniformOutput", false);
        % XBT1_4 = cellfun( @(rho,chi) -rho*cr'*(chi*(br-cr)*(br-cr)'*chi'), num2cell(rho_1), chi_matrix, "UniformOutput", false);
        % XBT1_5 = cellfun(@(rho,chi) -(vartheta+h)*rho*(br-cr)'*chi', num2cell(rho_1), chi_matrix,"UniformOutput",false);
        % BT1_3 = zeros(1,2*N);
        % BT1_4 = zeros(1,2*N);
        % BT1_5 = zeros(1,2*N);
        % for i = 1:N
        %     BT1_3 = BT1_3 + XBT1_3{i};            
        %     BT1_4 = BT1_4 + XBT1_4{i};
        %     BT1_5 = BT1_5 + XBT1_5{i};
        % end
        % BT1 = BT1_1 + BT1_2 + BT1_3 + BT1_4 + BT1_5;

        % % Target Value Test
        % C1
        % C1_1 = -lambda_0'*br;
        % C1_2 = rho_0/2*(br'*br);
        % C1_3 = sum(cellfun(@(lambda,chi) lambda*(-cr'*chi*(br-cr)-vartheta-h),num2cell(lambda_1),chi_matrix));
        % C1_4 = sum(cellfun(@(rho,chi) rho/2*((vartheta+h)^2+cr'*chi*(br-cr)*(br-cr)'*chi'*cr),num2cell(rho_1),chi_matrix));
        % C1 = C1_1 + C1_2 + C1_3 + C1_4;
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