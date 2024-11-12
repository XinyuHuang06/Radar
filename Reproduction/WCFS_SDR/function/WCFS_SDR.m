function varargout = WCFS_SDR(N,threshold,c,epsilon)
    % % 2. 相似度约束矩阵
    R0 = c*c';
    % % 4. 循环谱约束矩阵
    S = GenerateMatrixS(N, 4, threshold);
    S = (S + S')/2;
    x = exp(1j*2*pi*rand(N,1)); % 随机初始化
    % x = c;
    % X0 = x*x';
    %% SDP问题求解
    delta = (1-epsilon/2)^2;
    delta = N^2*delta;

    cvx_begin sdp quiet
        cvx_solver SDPT3
            variable X(N,N) complex hermitian
            minimize(real(trace(S*X)));
            subject to
                X >= 0;
                diag(X) == 1;
                (trace(R0*X)) >= delta;
    cvx_end
    % % % 向量求解
    [eigenvector,D] = eig(X);
    eigenvalue = diag(D);
    % eigenvalue_threshold = max(eigenvalue)*1e-4;
    % rank_X = sum(eigenvalue > eigenvalue_threshold);
    [value,num] = max(eigenvalue);
    x1 = sqrt(value)*eigenvector(:,num);
    varargout{1} = x1;
end