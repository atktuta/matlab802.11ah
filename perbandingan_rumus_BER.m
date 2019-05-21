
% Khusus BPSK (berlaku juga di OFDM-BPSK untuk rumus) untuk AWGN dan
% Rayleigh.
% kalau simulasi OFDM, Eb_N0 harus diubah ke Es_N0
clc;
clear;

Eb_N0_dB = -10:8;
EbN0Lin = 10.^(Eb_N0_dB/10);

% rumus BER AWGN BPSK
theoryBerAWGN_pake_erfc = 0.5*erfc(sqrt(EbN0Lin)); % BPSK
theoryBerAWGN_pake_qfunc = qfunc(sqrt(2*EbN0Lin));

% rumus BER rayleigh BPSK
theoryBerRayleigh = 0.5.*(1-sqrt(EbN0Lin./(EbN0Lin+1)));
% pakai fungsi berfading, uncoded rayleigh
theoryBerRayleigh_pake_berfading = berfading(Eb_N0_dB, 'psk', 2, 1);

% hasilnya sama, jadi terserah mau pakai yang mana.

% simulation
% ada stream data
random_binary = randi([0, 1], 1, 10000);
% dimodulasi bpsk
hasil_modulasi = qammod(random_binary, 2,'UnitAveragePower',true, ...
    'PlotConstellation',true);

rasio_error_awgn_snr_sendiri = zeros(1,length(Eb_N0_dB));
rasio_error_rayleigh_snr_sendiri = zeros(1,length(Eb_N0_dB));
% lewatkan channel
for i=1:length(Eb_N0_dB)
    noise = 1/sqrt(2)*(randn(1,length(random_binary)) + ...
        1i*randn(1,length(random_binary))); % white gaussian noise, 0dB variance 
    h = 1/sqrt(2)*(randn(1,length(random_binary)) + ...
        1i*randn(1,length(random_binary)));  % Rayleigh channel

    % Channel and noise Noise addition
    y1 = hasil_modulasi + 10^(-Eb_N0_dB(i)/20)*noise; 
    y2 = h.*hasil_modulasi + 10^(-Eb_N0_dB(i)/20)*noise; 

    % equalization
    yHat = y2./h;
    % demodulasi bpsk
    hasil_awgn_sendiri_demod = qamdemod(y1,2);
    hasil_rayleigh_sendiri_demod = qamdemod(yHat,2);
    % lihat BER-nya  
    [num_error2, rasio_error2] = biterr(random_binary,hasil_awgn_sendiri_demod);   
    [num_error3, rasio_error3] = biterr(random_binary,hasil_rayleigh_sendiri_demod);   
    rasio_error_awgn_snr_sendiri(i) = rasio_error2;
    rasio_error_rayleigh_snr_sendiri(i) = rasio_error3;    
end

figure
semilogy(Eb_N0_dB,theoryBerAWGN_pake_erfc,'b-o','LineWidth',1);
hold on
semilogy(Eb_N0_dB,theoryBerAWGN_pake_qfunc,'rp-','LineWidth',1);
semilogy(Eb_N0_dB,theoryBerRayleigh,'g-o','LineWidth',1);
semilogy(Eb_N0_dB,theoryBerRayleigh_pake_berfading,'kp-','LineWidth',1);
semilogy(Eb_N0_dB,rasio_error_awgn_snr_sendiri,'y<-','LineWidth',1);
semilogy(Eb_N0_dB,rasio_error_rayleigh_snr_sendiri,'c>-','LineWidth',1);
grid on
legend("theoryBerAWGN erfc", "theoryBerAWGN qfunc", ...
    "theoryBerRayleigh rumus", "theoryBerRayleigh berfading", ...
    "Sim AWGN", "Sim Rayleigh");
ylabel("Bit Error Rate");
xlabel("Eb/N0 in dB");
axis([Eb_N0_dB(1) Eb_N0_dB(length(Eb_N0_dB)) 1e-3 1])
title("Eb/N0");

%% sedangkan kalau pake SNR
SNRdB=-10:8;
SNRLin=10.^(SNRdB/10);
theory_snr_pake_qfunc = qfunc(sqrt(SNRLin));
theory_snr_pake_erfc = 0.5*erfc(sqrt(SNRLin/2)); % ??? ngasal
%theory_snr_pake_erfc = 0.5*erfc(sqrt(2*SNRLin)/sqrt(2));
theory_snr_rayleigh = 0.5.*(1-sqrt(SNRLin./(SNRLin+1))); % ??? ngasal

% simulasi
% ada stream data
random_binary = randi([0, 1], 1, 10000);
% dimodulasi bpsk
hasil_modulasi = qammod(random_binary, 2,'UnitAveragePower',true, ...
    'PlotConstellation',false);

% preallocation for speed
rasio_error_awgn_snr_fungsi = zeros(1,length(SNRdB));
rasio_error_awgn_snr_sendiri = zeros(1,length(SNRdB));
rasio_error_rayleigh_snr_sendiri = zeros(1,length(SNRdB));
% lewatkan awgn
for i=1:length(SNRdB)
    hasil_awgn_pake_fungsi = awgn(hasil_modulasi,SNRdB(i),'measured');
    noise=1/sqrt(2)*(randn(1,length(random_binary))+1i*randn(1,length(random_binary)));
    hasil_awgn_pake_sendiri=sqrt(SNRLin(i))*hasil_modulasi+noise;
    % hasil_rayleigh = ??;
    taps =1/sqrt(2)*(randn(1,length(random_binary))+1i*randn(1,length(random_binary))); 
    % rumus pertama
    hasil_rayleigh_pake_sendiri=sqrt(SNRLin(i))*abs(taps).*hasil_modulasi+noise;
    % rumus kedua pakai equalization
    tesmp = taps.*hasil_modulasi + 10^(-SNRdB(i)/20)*noise; 
    %hasil_rayleigh_pake_sendiri = tesmp./taps; % equalization
    % hasil rumus pertama dan kedua sama

    % demodulasi bpsk
    hasil_awgn_fungsi_demod = qamdemod(hasil_awgn_pake_fungsi,2);
    hasil_awgn_sendiri_demod = qamdemod(hasil_awgn_pake_fungsi,2);
    hasil_rayleigh_sendiri_demod = qamdemod(hasil_rayleigh_pake_sendiri,2);
    % lihat BER-nya
    [num_error, rasio_error] = biterr(random_binary,hasil_awgn_fungsi_demod);   
    [num_error2, rasio_error2] = biterr(random_binary,hasil_awgn_sendiri_demod);   
    [num_error3, rasio_error3] = biterr(random_binary,hasil_rayleigh_sendiri_demod);   
    rasio_error_awgn_snr_fungsi(i) = rasio_error;
    rasio_error_awgn_snr_sendiri(i) = rasio_error2;
    rasio_error_rayleigh_snr_sendiri(i) = rasio_error3;
end

figure
semilogy(SNRdB,theory_snr_pake_qfunc,'b-o','LineWidth',1);
hold on
semilogy(SNRdB,theory_snr_pake_erfc,'rp-','LineWidth',1);
semilogy(SNRdB,rasio_error_awgn_snr_fungsi,'gp-','LineWidth',1);
semilogy(SNRdB,rasio_error_awgn_snr_sendiri,'k>-','LineWidth',1);
semilogy(SNRdB,theory_snr_rayleigh,'y<-','LineWidth',1);
semilogy(SNRdB,rasio_error_rayleigh_snr_sendiri,'bx-','LineWidth',1);
grid on
legend("theoryBerAWGN qfunc", "theoryBerAWGN erfc", "sim AWGN fungsi", ...
    "sim AWGN sendiri", "theoryRayleigh", "sim sendiri Rayleigh" );
ylabel("Bit Error Rate");
xlabel("SNR in dB");
axis([SNRdB(1) SNRdB(length(SNRdB)) 1e-3 1])
title("SNR");



%% sedangkan kalau pake jarak
% sebenarnya tinggal ganti snr dengan jarak
% berdasarkan relasi SNR = power received - noise
% dimana power received = power transmit - pathloss + gain
% dan noise = noise_temp + noise_figure

jarak = 1:100:801;
power_transmit = -20; % 10 mW = 10 dBm = -20 dBW
path_loss = path_loss_rooftop_macro_deployment(jarak); % dB
gain = 0; % dB
power_received = power_transmit - path_loss + gain;
noise_temp = -145.22; % dB
noise_figure = 3; % dB
SNRdB = power_received - noise_temp - noise_figure;

%SNRdB=-10:8;
SNRLin=10.^(SNRdB/10);
theory_snr_pake_qfunc = qfunc(sqrt(SNRLin));
theory_snr_pake_erfc = 0.5*erfc(sqrt(SNRLin/2)); % ??? ngasal
%theory_snr_pake_erfc = 0.5*erfc(sqrt(2*SNRLin)/sqrt(2));
theory_snr_rayleigh = 0.5.*(1-sqrt(SNRLin./(SNRLin+1))); % ??? ngasal

% simulasi
% ada stream data
random_binary = randi([0, 1], 1, 10000);
% dimodulasi bpsk
hasil_modulasi = qammod(random_binary, 2,'UnitAveragePower',true, ...
    'PlotConstellation',false);

% preallocation for speed
rasio_error_awgn_snr_fungsi = zeros(1,length(SNRdB));
rasio_error_awgn_snr_sendiri = zeros(1,length(SNRdB));
rasio_error_rayleigh_snr_sendiri = zeros(1,length(SNRdB));
% lewatkan awgn
for i=1:length(SNRdB)
    hasil_awgn_pake_fungsi = awgn(hasil_modulasi,SNRdB(i),'measured');
    noise=1/sqrt(2)*(randn(1,length(random_binary))+1i*randn(1,length(random_binary)));
    hasil_awgn_pake_sendiri=sqrt(SNRLin(i))*hasil_modulasi+noise;
    % hasil_rayleigh = ??;
    taps =1/sqrt(2)*(randn(1,length(random_binary))+1i*randn(1,length(random_binary))); 
    % rumus pertama
    hasil_rayleigh_pake_sendiri=sqrt(SNRLin(i))*abs(taps).*hasil_modulasi+noise;
    % rumus kedua pakai equalization
    tesmp = taps.*hasil_modulasi + 10^(-SNRdB(i)/20)*noise; 
    %hasil_rayleigh_pake_sendiri = tesmp./taps; % equalization
    % hasil rumus pertama dan kedua sama

    % demodulasi bpsk
    hasil_awgn_fungsi_demod = qamdemod(hasil_awgn_pake_fungsi,2);
    hasil_awgn_sendiri_demod = qamdemod(hasil_awgn_pake_fungsi,2);
    hasil_rayleigh_sendiri_demod = qamdemod(hasil_rayleigh_pake_sendiri,2);
    % lihat BER-nya
    [num_error, rasio_error] = biterr(random_binary,hasil_awgn_fungsi_demod);   
    [num_error2, rasio_error2] = biterr(random_binary,hasil_awgn_sendiri_demod);   
    [num_error3, rasio_error3] = biterr(random_binary,hasil_rayleigh_sendiri_demod);   
    rasio_error_awgn_snr_fungsi(i) = rasio_error;
    rasio_error_awgn_snr_sendiri(i) = rasio_error2;
    rasio_error_rayleigh_snr_sendiri(i) = rasio_error3;
end

N = 20;
PER = 1 - (1 - theory_snr_rayleigh).^N;
PER_sim = 1 - (1 - rasio_error_rayleigh_snr_sendiri).^N;
figure
semilogy(jarak,theory_snr_pake_qfunc,'b-o','LineWidth',1);
hold on
semilogy(jarak,theory_snr_pake_erfc,'rp-','LineWidth',1);
semilogy(jarak,rasio_error_awgn_snr_fungsi,'gp-','LineWidth',1);
semilogy(jarak,rasio_error_awgn_snr_sendiri,'k>-','LineWidth',1);
semilogy(jarak,theory_snr_rayleigh,'b<-','LineWidth',1);
semilogy(jarak,rasio_error_rayleigh_snr_sendiri,'yx-','LineWidth',1);
semilogy(jarak,PER,'co-','LineWidth',1);
semilogy(jarak,PER_sim,'mx-','LineWidth',1);
grid on
legend("theoryBerAWGN qfunc", "theoryBerAWGN erfc", "sim AWGN fungsi", ...
    "sim AWGN sendiri", "theoryRayleigh", "sim sendiri Rayleigh",...
    "PER theory", "PER sim");
ylabel("Bit and Packet Error Rate");
xlabel("Distance in meter");
axis([jarak(1) jarak(length(jarak)) 1e-3 1])
str = sprintf("Jarak dengan N = %d", N);
title(str);

