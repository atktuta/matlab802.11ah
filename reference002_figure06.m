% coba reproduce Ferrand 2013
% figure 6, 
% menggunakan block rayleigh dan block rice, BPSK, N = 312 bits

clear all;
clc;

% N adalah number of symbol per packet.
N = 312/8;  % 312 bits diubah ke bytes

SER_SNR_threshold = 1 - 10^(-0.3/N);

% SER = Q(sqrt(k*SNR)), BPSK k = 1
akar_SNR_threshold = qfuncinv(SER_SNR_threshold);
SNR_threshold = akar_SNR_threshold^2;

SNR = 0:40;

PER_Rayleigh = 1 - exp(-SNR_threshold./(10.^(SNR./10)));
%PER_Rayleigh2 = 1 - exp(-SNR_threshold/(10^(SNR(11)/10)));
K = 4;
PER_Rice = 1 - marcumq(sqrt(2*K),sqrt(2*(K+1)*SNR_threshold./(10.^(SNR./10))));

figure
semilogy(SNR, PER_Rayleigh, 'r');
grid on
axis([0 40 1e-3 1])
hold on
semilogy(SNR, PER_Rice, 'b');

xlabel('Average SNR (dB)')
ylabel('PER')
legend('Rayleigh',...
    'Rice (K=4)'); 

