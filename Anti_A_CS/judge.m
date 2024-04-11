% 判断两种计算方法是否等价
clc
clear
N = 64;
M = 4;
x = exp(1j*2*pi*rand(N,1));
fs = 1e9;
FN = dftmtx(N);
Fshift = zeros(N);
Fshift(1:N/2,N/2+1:end) = eye(N/2);
Fshift(N/2+1:end,1:N/2) = eye(N/2);
FN = Fshift*FN;
i_temp = 0;
x_dft = FN*x;
% x_dft = fftshift(x_dft);
taf_matrix = cell(N*(N+1)/2,1);
CS_1 = zeros(N,N);
m = -M/2+1:M/2;
df = fs/10;
% 方法一：使用矩阵截取进行计算
for iter_1 = 1:N
    for iter_2 = iter_1:N
        i_temp = i_temp + 1;
        B = max(1-iter_1, 1-iter_2);
        A = min (N-iter_1, N-iter_2);
        n = m((m<=A) & (m>=B));
        if isempty(n)
            temp_taf = zeros(N,N);
        else
            temp_p = iter_1 + n;
            temp_q = iter_2 + n;
            temp_taf = zeros(N,N);
            temp_taf(temp_p(1):temp_p(end),temp_q(1):temp_q(end)) = eye(length(temp_q));
        end
        taf_matrix{i_temp} = temp_taf;
        CS_1(iter_1,iter_2) = transpose(x_dft)*temp_taf*(conj(x_dft));
        if (iter_2 > iter_1) && iter_1 < N
            if CS_1(iter_2, iter_1) == 0
                CS_1(iter_2, iter_1) = transpose(x_dft)*temp_taf*(conj(x_dft));
            else
                fprintf('error');
            end
        end
    end
end
% % 方法2：使用DFSM方法进行计算
CS_1 = real(CS_1./M);
CS_1 = abs(CS_1./max(CS_1(:)));
[CS_2, f, alfa] = Analysis_CS_DFSM(fs, x, df, M, 'bool_draw', 0);
figure
contour(f,alfa,CS_1);
figure
contour(f,alfa,CS_2);