function [CFS, alfa, f] = WCFS(signal, M)
    signal = signal(:);
    N = length(signal);
    alfa = -(N-M)/2:(N-M)/2;
    f = 0:N-M;

    % % FFT 
    x = fftshift(fft(signal,N));
    xc = conj(x);
    % % 
    m = -M/2:1:M/2-1;
    for i_alfa = 1:length(alfa)
        for i_f = 1:length(f)
            t_alfa = alfa(i_alfa);
            t_f = alfa(i_f);
            

        end
    end
end