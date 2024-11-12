X = [1, 0.5 + 0.5i; 0.5 - 0.5i, 1];
NumL = 1000; % 生成Num_L组随机向量
randCN = zeros(NumL, size(X,1));
randCN_out = zeros(NumL, 1);
L = chol(X, 'lower');
for i = 1:NumL
    z = (randn(size(X, 1), 1) + 1i * randn(size(X, 1), 1)) / sqrt(2);
    temp_randCN = L * z;
    randCN(i,:) = temp_randCN;
    randCN_out(i) = temp_randCN*C*temp_randCN';
end
[~, index] = min(randCN_out);
x = randCN(index, :);