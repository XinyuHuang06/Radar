fs = 20*1e6;
fc = 1*1e6;
Rb = 1*1e6;
code_seq = round(rand(1,16));
[signal_BPSK,t] = generator_PSK(fs,fc,Rb,code_seq,2);



backer_seq = [0, +1, 0, +1, 0, +1, 0, +1, 0, +1, 0, +1];

figure
subplot(121);
[signal_BPSK,t] = generator_PSK(fs,0,Rb,backer_seq,2);
plot(real(signal_BPSK));
title('基带波形（13位巴克码）');
[signal_BPSK,t] = generator_PSK(fs,fc,Rb,backer_seq,2);
subplot(122);
plot(real(signal_BPSK));
title('调制后波形（13位巴克码）');
exportgraphics(gcf, './output/BPSK_theory.jpg', 'Resolution', 300);
