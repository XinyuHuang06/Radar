function [W1,i_p,i_q]=WCS_test(s,N,M)
p = 0:1:(N);
q = 0:2:N;
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
W1 = W1/max(W1(:));
end