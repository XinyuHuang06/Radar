function [xr] = Update_xr(A1,BT1,xr_in)
    % Constant Modules Constraints 1/2*xT*H*x + k*x + d = 0 => xT*(2*I)*x + 0*x -1 = 0
    N = length(xr_in)/2;
    temp_H = cell(1);temp_k= cell(1);temp_d= cell(1);
    temp_H{1} = eye(2*N)*2;temp_k{1} = zeros(2*N,1);temp_d{1} = -1;
    options = optimoptions(@fmincon,'Algorithm','interior-point',...
        'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
        'HessianFcn',@(x,lambda)quadhess(lambda,A1*2,temp_H),'Display','off');
    [xr,~,~,~] = fmincon(@(x) quadobj(x,A1*2,BT1',0), xr_in,[],[],[],[],[],[],...
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