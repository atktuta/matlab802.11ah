
% Khusus BPSK (berlaku juga di OFDM-BPSK) untuk AWGN dan Rayleigh

clc;
clear all;

Eb_N0_dB = -10:8;

theoryBerAWGN_pake_erfc = 0.5*erfc(sqrt(10.^(Eb_N0_dB/10))); % BPSK
theoryBerAWGN_pake_qfunc = qfunc(sqrt(2*10.^(Eb_N0_dB/10)));

% rumus BER rayleigh BPSK AWGN
EbN0Lin = 10.^(Eb_N0_dB/10);
theoryBerRayleigh = 0.5.*(1-sqrt(EbN0Lin./(EbN0Lin+1)));

% pakai fungsi berfading, uncoded rayleigh
theoryBerRayleigh_pake_berfading = berfading(Eb_N0_dB, 'psk', 2, 1);

% hasilnya sama, jadi terserah mau pakai yang mana.

figure
semilogy(Eb_N0_dB,theoryBerAWGN_pake_erfc,'b-o','LineWidth',1);
hold on
semilogy(Eb_N0_dB,theoryBerAWGN_pake_qfunc,'rp-','LineWidth',1);
semilogy(Eb_N0_dB,theoryBerRayleigh,'g-o','LineWidth',1);
semilogy(Eb_N0_dB,theoryBerRayleigh_pake_berfading,'kp-','LineWidth',1);
grid on
legend("theoryBerAWGN erfc", "theoryBerAWGN qfunc", ...
    "theoryBerRayleigh rumus", "theoryBerRayleigh berfading" );
ylabel("Bit Error Rate");
xlabel("Eb/N0 in dB");
axis([Eb_N0_dB(1) Eb_N0_dB(length(Eb_N0_dB)) 1e-3 1])

%% sedangkan kalau pake SNR
SNRdB=-10:8;
SNRLin=10.^(SNRdB/10);
theory_snr_pake_qfunc = qfunc(sqrt(SNRLin));
theory_snr_pake_erfc = 0.5*erfc(sqrt(SNRLin/2));
figure
semilogy(SNRdB,theory_snr_pake_qfunc,'b-o','LineWidth',1);
hold on
semilogy(SNRdB,theory_snr_pake_erfc,'rp-','LineWidth',1);
grid on
legend("theoryBerAWGN qfunc", "theoryBerAWGN erfc" );
ylabel("Bit Error Rate");
xlabel("SNR in dB");
axis([SNRdB(1) SNRdB(length(SNRdB)) 1e-3 1])
