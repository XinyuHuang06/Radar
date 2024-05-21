N = 512;
M = 4;
fs = 4e6;
T = N/fs;
fc = 1e6;
temp_B = 1e6;
s = generator_LFM(fs,fc,temp_B,T);
s = exp(1j*2*pi*randn(N,1));
p = 0:1:(N-M);
q = 0:2:N-M;
s = (s);
x = fftshift(fft(s,N));
xc = conj(x);
W1 = zeros(length(p),length(q));
indexp = zeros(length(p),length(q));
indexq = zeros(length(p),length(q));
for i_p = 1:length(p)
    for i_q = 1:length(q)
        tp = p(i_p);
        tq = q(i_q);
        m = -M/2:1:M/2-1;
        temp_1 = tp + tq/2;
        temp_2 = tp - tq/2;
        index_1 = m + temp_1;
        index_2 = m + temp_2;
        index = index_1>0 & index_2>0 & index_1<=N & index_2<=N;
        a1 = index_1(index);
        a2 = index_2(index);
        if isempty(index)
            W1(i_p,i_q) = -1;
        else
            W1(i_p,i_q) = abs(sum(x(a1).*xc(a2)));
            indexp(i_p,i_q) = i_p;
            indexq(i_p,i_q) = i_q;
        end
    end
end
% [S_xr, ~, ~] = Analysis_CS_DFSM(fs,real(s),fs/256,M,'bool_draw',1);
W2 = fliplr(W1(:,1:end-1));
W = [W2,W1];
plot(W(450,:))
figure 
contour(indexp,indexq,W1)