
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

N = 27;

%% simulasi
for i=1:length(SNRdB)
    PER_total = 0;
    throughput_total = 0;
    k = 1;
    % reset nilai variabel temporary
    PER_tanpa_XL = zeros(1,(simulation_time/beacon_interval));
    throughput_tanpa_XL = zeros(1,(simulation_time/beacon_interval));
    PER_dengan_XL = zeros(1,(simulation_time/beacon_interval));
    throughput_dengan_XL = zeros(1,(simulation_time/beacon_interval));
    
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
            % karena ini situasi ketika bodyloss, maka SNR yang dipakai adalah
            % SNR bodyloss.
            % lewatkan channel
            hasil_rayleigh_pake_sendiri_body=sqrt(SNRLin_body(i))*abs(taps).*hasil_modulasi+noise;
            % demod
            hasil_rayleigh_sendiri_demod_body = qamdemod(hasil_rayleigh_pake_sendiri_body,2);
            % hitung jumlah bit error
            [num_error3, rasio_error3] = biterr(random_binary,hasil_rayleigh_sendiri_demod_body);
            
            % -------------BER without cross layer.            
            BER_tanpa_XL = num_error3/jumlah_bit_dikirim;
            % hitung PER-nya
            PER_tanpa_XL(k) = 1 - (1 - BER_tanpa_XL).^N;
            throughput_tanpa_XL(k) = hitung_throughput_80211ah(PER_tanpa_XL(k));
                        
            % ------------BER with cross layer
            % didn't send, so BER and throughput zero
            PER_dengan_XL(k) = 0;
            % throughput 
            throughput_dengan_XL(k) = 0;
            

        else
            % karena ini situasi normal, tanpa body pathloss, maka 
            % SNR yang dipakai adalah SNR normal 
            % lewatkan channel
            hasil_rayleigh_pake_sendiri=sqrt(SNRLin_normal(i))*abs(taps).*hasil_modulasi+noise;
            % demod
            hasil_rayleigh_sendiri_demod = qamdemod(hasil_rayleigh_pake_sendiri,2);
            % % hitung jumlah bit error
            [num_error2, rasio_error2] = biterr(random_binary,hasil_rayleigh_sendiri_demod);
            
            % -------------BER without cross layer.
            BER_tanpa_XL = num_error2/jumlah_bit_dikirim;
            % hitung PER-nya
            PER_tanpa_XL(k) = 1 - (1 - BER_tanpa_XL).^N;
            throughput_tanpa_XL(k) = hitung_throughput_80211ah(PER_tanpa_XL(k));            
            
            % -------------BER with cross layer
            % BER nya sama dengan yang without cross layer
            % PER nya sama dengan yang without cross layer
            PER_dengan_XL(k) = PER_tanpa_XL(k);
            % throughput 
            throughput_dengan_XL(k) = throughput_tanpa_XL(k);            
        end
        k = k + 1;    
    end
    
    % hitung rata-rata untuk seluruh percobaan
    PER_tanpa_XL_rata_rata(i) = sum(PER_tanpa_XL)/length(PER_tanpa_XL);
    throughput_tanpa_XL_rata_rata(i) = sum(throughput_tanpa_XL)/length(throughput_tanpa_XL);
    PER_dengan_XL_rata_rata(i) = sum(PER_dengan_XL)/length(PER_dengan_XL);
    throughput_dengan_XL_rata_rata(i) = sum(throughput_dengan_XL)/length(throughput_dengan_XL);

end

%% teori PER
PER_teori_tanpa_XL = PER_ferrand_body * (Prob_shadow) + ...
         PER_ferrand_normal * (1-Prob_shadow);
PER_teori_dengan_XL = PER_ferrand_normal;

%% curve PER
figure
semilogy(jarak,PER_tanpa_XL_rata_rata,'b<-','LineWidth',1);
hold on
semilogy(jarak,PER_dengan_XL_rata_rata,'c>-','LineWidth',1);
semilogy(jarak,PER_teori_tanpa_XL,'go-','LineWidth',1);
semilogy(jarak,PER_teori_dengan_XL,'rx-','LineWidth',1);
grid on
legend("PER simulation w/o XL", "PER simulation with XL", ...
    "PER teori w/o XL", "PER teori with XL");
ylabel("Packet Error Rate");
xlabel("Distance APS-T (m)");
axis([jarak(1) jarak(length(jarak)) 1e-2 1])
str = sprintf("PER N = %d, Prob shadow = %f", N, Prob_shadow);
%title(str);

%% teori throughput
throughput_teori_tanpa_XL = hitung_throughput_80211ah(PER_teori_tanpa_XL);
throughput_teori_dengan_XL = hitung_throughput_80211ah(PER_teori_dengan_XL)*(1-Prob_shadow);

%% curve throughput
figure
semilogy(jarak,throughput_tanpa_XL_rata_rata,'b<-','LineWidth',1);
hold on
semilogy(jarak,throughput_dengan_XL_rata_rata,'r>-','LineWidth',1);
semilogy(jarak,throughput_teori_tanpa_XL,'bo--','LineWidth',1);
semilogy(jarak,throughput_teori_dengan_XL,'rx--','LineWidth',1);
grid on
legend("throughput simulation w/o XL", "throughput simulation with XL",...
    "throughput teori w/o XL", "throughput teori with XL");
ylabel("Average Throughput (bps)");
xlabel("Distance AP-ST (m)");
axis([jarak(1) jarak(length(jarak)) 1e4 0.5*1e6])
str = sprintf("Throughput N = %d, Prob shadow = %f", N, Prob_shadow);
%title(str);


