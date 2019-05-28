
% Khusus BPSK (berlaku juga di OFDM-BPSK untuk rumus) untuk Rayleigh.
% untuk mencari PER, kita menggunakan simulasi/theory utk memperoleh BER,
% dan memasukkannya ke rumus PER = 1 - (1-BER)^N

clc;
clear;

%% sedangkan kalau pake jarak
% sebenarnya tinggal ganti snr dengan jarak
% berdasarkan relasi SNR = power received - noise
% dimana power received = power transmit - pathloss + gain
% dan noise = noise_temp + noise_figure

jumlah_bit_dikirim = 1000;
jarak = 100:50:801;
power_transmit = -20; % 10 mW = 10 dBm = -20 dBW
gain = 3; % dB
noise_temperature = -145.22; % dB
noise_figure = 5; % dB

%ber_theory_snr_rayleigh = 0.5.*(1-sqrt(SNRLin./(SNRLin+1)));

%% hitung SNR threshold
% SER = Q(sqrt(k*SNR)), BPSK k = 1
NuntukThreshold = 100; % number of symbol per packet
SER_SNR_threshold = 1 - 10^(-0.3/NuntukThreshold); % -0.3 berasal dari log10(1/2)
akar_SNR_threshold = qfuncinv(SER_SNR_threshold);
SNR_threshold_N100 = akar_SNR_threshold^2;

%% hitung PER body shadow 
% rumus PER ini diperoleh dari approximation Ferrand, 2013
% untuk BPSK rayleigh
SNRdB = power_transmit - ...
    path_loss_body_shadow(jarak) - ...
    noise_temperature - noise_figure + gain;
SNRLin_body = 10.^(SNRdB/10);
PER_ferrand_body = 1-exp(-SNR_threshold_N100./SNRLin_body);

%% hitung PER normal = PER cross layer

% rumus PER ini diperoleh dari approximation Ferrand, 2013
% untuk BPSK rayleigh
SNRdB = power_transmit - ...
    path_loss_rooftop_macro_deployment(jarak) - ...
    noise_temperature - noise_figure + gain;
SNRLin_normal = 10.^(SNRdB/10);
PER_ferrand_normal = 1-exp(-SNR_threshold_N100./SNRLin_normal);

% curve perbandingan antara jarak dengan PER
% figure
% semilogy(jarak, PER_ferrand_normal,'b--')
% grid on
% hold on
% semilogy(jarak, PER_ferrand_body,'go-')
% xlabel('Distance (m)')
% ylabel('Average PER') 
% legend("Theory PER normal", "Theory PER bodyloss");
% ylabel("Packet Error Rate");
% xlabel("Distance in meter");
% axis([jarak(1) jarak(length(jarak)) 1e-3 1])
% str = sprintf("Jarak dengan N = %d", NuntukThreshold);
% title(str);


%% hitung throughput normal
Throughput_normal = hitung_throughput_80211ah(PER_ferrand_normal);

%% hitung throughput body pathloss 
Throughput_body = hitung_throughput_80211ah(PER_ferrand_body);

%% curve untuk throughput
% figure
% semilogy(jarak, Throughput_normal,'k--')
% hold on
% semilogy(jarak, Throughput_body,'r-o')
% 
% grid on
% xlabel('Distance AP-ST (m)')
% ylabel('Average Throughput') 
% legend('802.11ah Pathloss','Body Pathloss'); 
% axis([jarak(1) jarak(length(jarak)) 1e1 1e6])
% str = sprintf("Throughput");
% title(str);

%% simulasi
% % ada stream data
% random_binary = randi([0, 1], 1, jumlah_bit_dikirim);
% % dimodulasi bpsk
% hasil_modulasi = qammod(random_binary, 2,'UnitAveragePower',true, ...
%     'PlotConstellation',false);
% 
% % preallocation for speed
% rasio_error_rayleigh_snr_sendiri = zeros(1,length(SNRdB));
% rasio_error_rayleigh_snr_sendiri_body = zeros(1,length(SNRdB));
% % lewatkan channel
% for i=1:length(SNRdB)
%     noise=1/sqrt(2)*(randn(1,length(random_binary))+1i*randn(1,length(random_binary)));
%     taps =1/sqrt(2)*(randn(1,length(random_binary))+1i*randn(1,length(random_binary))); 
% 
%     hasil_rayleigh_pake_sendiri=sqrt(SNRLin_normal(i))*abs(taps).*hasil_modulasi+noise;
%     hasil_rayleigh_pake_sendiri_body=sqrt(SNRLin_body(i))*abs(taps).*hasil_modulasi+noise;
% 
%     % demodulasi bpsk
%     hasil_rayleigh_sendiri_demod = qamdemod(hasil_rayleigh_pake_sendiri,2);
%     hasil_rayleigh_sendiri_demod_body = qamdemod(hasil_rayleigh_pake_sendiri_body,2);
%     % lihat BER-nya 
%     [num_error2, rasio_error2] = biterr(random_binary,hasil_rayleigh_sendiri_demod);   
%     [num_error3, rasio_error3] = biterr(random_binary,hasil_rayleigh_sendiri_demod_body);   
%     rasio_error_rayleigh_snr_sendiri(i) = rasio_error2;
%     rasio_error_rayleigh_snr_sendiri_body(i) = rasio_error3;
% end
% 
% N = 20;
% %PER_theory = 1 - (1 - theory_snr_rayleigh).^N;
% PER_sim = 1 - (1 - rasio_error_rayleigh_snr_sendiri).^N;
% PER_sim_body = 1 - (1 - rasio_error_rayleigh_snr_sendiri_body).^N;
% figure
% semilogy(jarak,PER_ferrand_normal,'b<-','LineWidth',1);
% hold on
% semilogy(jarak,PER_sim,'c>-','LineWidth',1);
% semilogy(jarak,PER_ferrand_body,'go-','LineWidth',1);
% semilogy(jarak,PER_sim_body,'rx-','LineWidth',1);
% grid on
% legend("Ferrand's PER Rayleigh normalloss", "sim PER Rayleigh normalloss", ...
%    "Ferrand's PER Rayleigh bodyloss", "sim PER Rayleigh bodyloss" );
% ylabel("Packet Error Rate");
% xlabel("Distance in meter");
% axis([jarak(1) jarak(length(jarak)) 1e-3 1])
% str = sprintf("PER dengan N = %d", N);
% title(str);

%% untuk cross layer, tambahkan peluang terjadinya body pathloss

% preallocation for speed
rasio_error_rayleigh_snr_sendiri = zeros(1,length(SNRdB));
rasio_error_rayleigh_snr_sendiri_body = zeros(1,length(SNRdB));
num_error_without_XL = zeros(1,length(SNRdB));
num_error_with_XL = zeros(1,length(SNRdB));

simulation_time = 20; % 10 seconds
beacon_interval = 0.01; % 10 ms

% if body shadow. untuk ECG probability-nya 7%
% if body shadow. untuk Glucose Monitor probability-nya 5%
% if body shadow. untuk Blood Pressure probability-nya 10%  
Prob_shadow = 0.07;

% lewatkan channel
for i=1:length(SNRdB)
    pengali_bit_tidak_dikirim = 0;
    for detik=beacon_interval:beacon_interval:simulation_time
        % ada stream data
        random_binary = randi([0, 1], 1, jumlah_bit_dikirim);
        % dimodulasi bpsk
        hasil_modulasi = qammod(random_binary, 2,'UnitAveragePower',true, ...
        'PlotConstellation',false);
        % generate r acak, untuk menghasilkan peluang body pathloss
        r = rand;
        noise=1/sqrt(2)*(randn(1,length(random_binary))+1i*randn(1,length(random_binary)));
        taps =1/sqrt(2)*(randn(1,length(random_binary))+1i*randn(1,length(random_binary))); 
        % if body shadow. untuk ECG probability-nya 7%
        % if body shadow. untuk Glucose Monitor probability-nya 5%
        % if body shadow. untuk Blood Pressure probability-nya 10%        
        if (0<=r) && (r<Prob_shadow)
            % BER without cross layer.
            % karena ini situasi ketika bodyloss, maka SNR yang dipakai adalah
            % SNR bodyloss.

            % lewatkan channel
            hasil_rayleigh_pake_sendiri_body=sqrt(SNRLin_body(i))*abs(taps).*hasil_modulasi+noise;
            % demod
            hasil_rayleigh_sendiri_demod_body = qamdemod(hasil_rayleigh_pake_sendiri_body,2);
            % lihat BER-nya
            [num_error3, rasio_error3] = biterr(random_binary,hasil_rayleigh_sendiri_demod_body);
            num_error_without_XL(i) = num_error_without_XL(i) + num_error3;
            
            % BER with cross layer
            % didn't send, so BER and throughput zero
            % num_error_with_XL(i) = num_error_with_XL(i);
            pengali_bit_tidak_dikirim = pengali_bit_tidak_dikirim + 1;

        else
            % karena ini situasi normal, tanpa body pathloss, maka 
            % SNR yang dipakai adalah SNR normal 
            % lewatkan channel
            hasil_rayleigh_pake_sendiri=sqrt(SNRLin_normal(i))*abs(taps).*hasil_modulasi+noise;
            % demod
            hasil_rayleigh_sendiri_demod = qamdemod(hasil_rayleigh_pake_sendiri,2);
            % lihat BER-nya
            [num_error2, rasio_error2] = biterr(random_binary,hasil_rayleigh_sendiri_demod);
            num_error_without_XL(i) = num_error_without_XL(i) + num_error2;
            
            % BER with cross layer
            % sama dengan BER without cross layer
            num_error_with_XL(i) = num_error_with_XL(i) + num_error2;
        end
            
    end

end

BER_without_XL = num_error_without_XL/(jumlah_bit_dikirim*simulation_time/beacon_interval);
jumlah_bit_tidak_dikirim_with_XL = pengali_bit_tidak_dikirim *jumlah_bit_dikirim;
BER_with_XL = num_error_with_XL/((jumlah_bit_dikirim*simulation_time/beacon_interval)-...
    jumlah_bit_tidak_dikirim_with_XL);

%% coba tampilkan BER
% BER teori with XL
theory_snr_rayleigh = 0.5.*(1-sqrt(SNRLin_normal./(SNRLin_normal+1)));
theory_snr_rayleigh_body = 0.5.*(1-sqrt(SNRLin_body./(SNRLin_body+1)));

theory_tanpa_XL = (Prob_shadow*theory_snr_rayleigh_body) + ((1-Prob_shadow)*theory_snr_rayleigh);

figure
semilogy(jarak,BER_without_XL,'b<-','LineWidth',1);
hold on
semilogy(jarak,theory_tanpa_XL,'ro-','LineWidth',1);
semilogy(jarak,BER_with_XL,'c>-','LineWidth',1);
semilogy(jarak,theory_snr_rayleigh,'m<-','LineWidth',1);
%semilogy(jarak,theory_snr_rayleigh_body,'bo-','LineWidth',1);

grid on
legend("simulation w/o XL", ...
   "theory w/o XL", "simulation with XL",  ...
   "theory with XL");
ylabel("Bit Error Rate");
xlabel("Distance AP-ST (m)");
axis([jarak(1) jarak(length(jarak)) 1e-3 1])
str = sprintf("BER Prob shadow = %f", Prob_shadow);
title(str);

%% kalau PER

N = 23;

PER_sim_without_XL = 1 - (1 - BER_without_XL).^N;
PER_teori_without_XL = 1 - (1 - theory_tanpa_XL).^N;
PER_sim_with_XL = 1 - (1 - BER_with_XL).^N;
PER_teori_with_XL = 1 - (1 - theory_snr_rayleigh).^N;

PER_teori_without_XL = PER_ferrand_body * (Prob_shadow) + ...
         PER_ferrand_normal * (1-Prob_shadow);
PER_teori_with_XL = PER_ferrand_normal;

figure
semilogy(jarak,PER_sim_without_XL,'b<-','LineWidth',1);
hold on
semilogy(jarak,PER_teori_without_XL,'c>-','LineWidth',1);
semilogy(jarak,PER_sim_with_XL,'go-','LineWidth',1);
semilogy(jarak,PER_teori_with_XL,'rx-','LineWidth',1);
grid on
legend("PER simulation w/o XL", "PER analysis w/o XL", ...
   "PER simulation with XL", "PER analysis with XL" );
ylabel("Packet Error Rate");
xlabel("Distance APS-T (m)");
axis([jarak(1) jarak(length(jarak)) 1e-2 1])
%str = sprintf("PER dengan Prob. body pathloss happen = %f", Prob_shadow);
str = sprintf("PER N = %d, Prob shadow = %f", N, Prob_shadow);
%title(str);


%% hitung throughput cross layer
throughput_teori_without_XL = hitung_throughput_80211ah(PER_teori_without_XL);
throughput_sim_without_XL = hitung_throughput_80211ah(PER_sim_without_XL);
throughput_teori_with_XL = hitung_throughput_80211ah(PER_teori_with_XL);
throughput_sim_with_XL = hitung_throughput_80211ah(PER_sim_with_XL);

% throughput_teori_with_XL = hitung_throughput_80211ah(PER_ferrand_normal);
% throughput_teori_without_XL = hitung_throughput_80211ah(PER_ferrand_body) * (Prob_shadow) + ...
%          hitung_throughput_80211ah(PER_ferrand_normal) * (1-Prob_shadow);

figure
semilogy(jarak,throughput_teori_without_XL,'b<-','LineWidth',1);
hold on
semilogy(jarak,throughput_sim_without_XL,'c>-','LineWidth',1);
semilogy(jarak,throughput_teori_with_XL,'go-','LineWidth',1);
semilogy(jarak,throughput_sim_with_XL,'rx-','LineWidth',1);
grid on
legend("throughput analysis w/o XL", "throughput simulation w/o XL", ...
   "throughput analysis with XL", "throughput simulation with XL" );
ylabel("Average Throughput (bps)");
xlabel("Distance AP-ST (m)");
axis([jarak(1) jarak(length(jarak)) 1e4 0.5*1e6])
str = sprintf("Throughput N = %d, Prob shadow = %f", N, Prob_shadow);
%title(str);

%% coba bandingkan throughput 5 %
% throughput normal * (1-5%) * PER_normal untuk XL
% throughput PER 5%

