clear
clc

transmission_power_10mW = -20;
distance = 100:800;

% inisiasi untuk perhitungan throughput
m = 6; % maximum number of backoff stages and corresponds to 6 (i.e. CW_max = 2^6 CW_min )
mpdu_size = 475; % 475 Bytes

sigma_11ah = 1561/physconst('LightSpeed'); % jarak maks/kecepatan cahaya

DIFS_11ah = 264e-6; % 264 us untuk 11ah
SIFS_11ah = 160e-6; % 160 us untuk 11ah
CW_min = 15;
CW_max = 1023;
T_SLOT_11ah = 52e-6; % 52 us untuk 11a
T_BACKOFF_11ah = CW_min / 2 * T_SLOT_11ah; % apabila tidak ada collision dll.

T_preamble_header_11ah = 560e-6; % 560 us untuk 11ah
T_Sym_11ah = 40e-6; % 40 us untuk 11ah

L_header_plus_L_data_ack = 14; % 14 bytes untuk ACK
N_Sym_ack_11ah = (14 + L_header_plus_L_data_ack*8)/6;
T_ACK_11ah = T_preamble_header_11ah + (T_Sym_11ah * N_Sym_ack_11ah);

L_header_11ah = 26; % short header, kalau long header = 36
N_Sym_11ah = (14 + (L_header_11ah + mpdu_size) * 8) / 6;
T_DATA_11ah = T_preamble_header_11ah + (T_Sym_11ah * N_Sym_11ah);

T_message_wo_T_BACKOFF = DIFS_11ah + T_DATA_11ah + SIFS_11ah + ...
    T_ACK_11ah + 2 * sigma_11ah;

batas_atas = 60; % sebetulnya batas atas adalah tak hingga
                 % namun di atas 60 sudah convergen utk PER 0.9

%% hitung noise
noise_AWGN = -145.22; % dB
noise_figure = 5; % dB
antenna_gain = 3; % dB
noise_AWGN_figure_temp__fading__gain_dll = noise_AWGN + noise_figure -...
    antenna_gain; % dB

%% hitung SNR threshold
% SER = Q(sqrt(k*SNR)), BPSK k = 1
N = 100; % number of symbol per packet
SER_SNR_threshold = 1 - 10^(-0.3/N); % -0.3 berasal dari log10(1/2)
akar_SNR_threshold = qfuncinv(SER_SNR_threshold);
SNR_threshold_N100 = akar_SNR_threshold^2;

%% hitung PER body shadow 
C_per_N_10mW = transmission_power_10mW - ...
    path_loss_body_shadow(distance) - ...
    noise_AWGN_figure_temp__fading__gain_dll;
C_per_NLin = 10.^(C_per_N_10mW/10);
PER_body_shadow = 1-exp(-SNR_threshold_N100./C_per_NLin);

%% hitung PER normal = PER cross layer
% rumus PER ini diperoleh dari approximation Ferrand, 2013
% untuk BPSK rayleigh
C_per_N_10mW = transmission_power_10mW - ...
    path_loss_normal(distance) - ...
    noise_AWGN_figure_temp__fading__gain_dll;
C_per_NLin = 10.^(C_per_N_10mW/10);
PER_normal = 1-exp(-SNR_threshold_N100./C_per_NLin);

% curve perbandingan antara SNR dengan PER
% figure
% semilogy(C_per_N_10mW, PER_normal,'k--')
% grid on
% xlabel('SNR (dB)')
% ylabel('Average PER') 
% title('SNR vs PER BPSK Rayleigh'); 
% 
% % curve perbandingan antara jarak dengan PER
% figure
% semilogy(distance, PER_normal,'k--')
% grid on
% xlabel('Distance (m)')
% ylabel('Average PER') 
% title('Distance vs PER BPSK Rayleigh'); 

%% hitung throughput normal
T_BACKOFF = 0;
for j=1:batas_atas
    T_BACKOFF = T_BACKOFF + PDR_array(PER_normal,j)*T_backoff(CW_min, ...
        CW_max, T_SLOT_11ah, j, m);
end
T_message = T_message_wo_T_BACKOFF + T_BACKOFF;      
Throughput_normal = (mpdu_size * 8)./T_message .* (1-PER_normal);

%% hitung PER 5%
PER_5_persen_body_shadow = zeros(1,length(distance));
Prob_shadow = 0.05;
for i=1:length(distance)
    PER_5_persen_body_shadow(i) = PER_body_shadow(i) * (Prob_shadow) + ...
        PER_normal(i) * (1-Prob_shadow);
    T_BACKOFF = 0;
    for j=1:batas_atas
        T_BACKOFF = T_BACKOFF + PDR(PER_5_persen_body_shadow(i),j)*T_backoff(CW_min, ...
            CW_max, T_SLOT_11ah, j, m);
    end
    T_message = T_message_wo_T_BACKOFF + T_BACKOFF;      
    S = (mpdu_size * 8)/T_message * (1-PER_5_persen_body_shadow(i));
    Throughput_5_persen_body_shadow(i) = S;      
end

%% hitung PER 7%
PER_7_persen_body_shadow = zeros(1,length(distance));
Prob_shadow = 0.07;
for i=1:length(distance)
    PER_7_persen_body_shadow(i) = PER_body_shadow(i) * (Prob_shadow) + ...
        PER_normal(i) * (1-Prob_shadow);
    T_BACKOFF = 0;
    for j=1:batas_atas
        T_BACKOFF = T_BACKOFF + PDR(PER_7_persen_body_shadow(i),j)*T_backoff(CW_min, ...
            CW_max, T_SLOT_11ah, j, m);
    end
    T_message = T_message_wo_T_BACKOFF + T_BACKOFF;      
    S = (mpdu_size * 8)/T_message * (1-PER_7_persen_body_shadow(i));
    Throughput_7_persen_body_shadow(i) = S;     
end

%% hitung PER 10%
PER_10_persen_body_shadow = zeros(1,length(distance));
Prob_shadow = 0.1;
for i=1:length(distance)
    PER_10_persen_body_shadow(i) = PER_body_shadow(i) * (Prob_shadow) + ...
        PER_normal(i) * (1-Prob_shadow);
    T_BACKOFF = 0;    
    for j=1:batas_atas
        T_BACKOFF = T_BACKOFF + PDR(PER_10_persen_body_shadow(i),j)*T_backoff(CW_min, ...
            CW_max, T_SLOT_11ah, j, m);
    end
    T_message = T_message_wo_T_BACKOFF + T_BACKOFF;      
    S = (mpdu_size * 8)/T_message * (1-PER_10_persen_body_shadow(i));
    Throughput_10_persen_body_shadow(i) = S;      
end

%% curve
figure
semilogy(distance, PER_normal,'k--')
hold on
semilogy(distance, PER_5_persen_body_shadow,'r-o','MarkerIndices',1:40:length(distance))
semilogy(distance, PER_7_persen_body_shadow,'g-*','MarkerIndices',1:50:length(distance))
semilogy(distance, PER_10_persen_body_shadow,'b-+','MarkerIndices',1:60:length(distance))
grid on
xlabel('Distance AP-ST (m)')
ylabel('Average PER') 
legend('Tx=10mW, PER XL','Tx=10mW, PER 5%',...
'Tx=10mW, PER 7%', 'Tx=10mW, PER 10%'); 
axis([distance(1) distance(length(distance)) 1*1e-2 1])

% str = sprintf("SNR dengan N = %d", N);
% title(str);

%% curve throughput
figure
semilogy(distance, Throughput_normal,'k--')
hold on
semilogy(distance, Throughput_5_persen_body_shadow,'r-o','MarkerIndices',1:40:length(distance))
semilogy(distance, Throughput_7_persen_body_shadow,'g-*','MarkerIndices',1:50:length(distance))
semilogy(distance, Throughput_10_persen_body_shadow,'b-+','MarkerIndices',1:60:length(distance))

grid on
xlabel('Distance AP-ST (m)')
ylabel('Average Throughput (bps)') 
legend('Tx=10mW, Throughput XL','Tx=10mW, PER 5%',...
'Tx=10mW, PER 7%', 'Tx=10mW, PER 10%'); 
axis([distance(1) distance(length(distance)) 1*1e4 0.5*1e6])


%% coba bandingkan throughput 5 %
% throughput normal * (1-5%) * PER_normal untuk XL
% throughput PER 5%
throughput_XL_05 = (1-0.05)*Throughput_normal;
throughput_XL_07 = (1-0.07)*Throughput_normal;
throughput_XL_10 = (1-0.10)*Throughput_normal;
figure
semilogy(distance, throughput_XL_05,'r--')
hold on
semilogy(distance, throughput_XL_07,'b--','MarkerIndices',1:50:length(distance))
semilogy(distance, throughput_XL_10,'g--','MarkerIndices',1:60:length(distance))
semilogy(distance, Throughput_5_persen_body_shadow,'r-o','MarkerIndices',1:40:length(distance))
semilogy(distance, Throughput_7_persen_body_shadow,'b-o','MarkerIndices',1:40:length(distance))
semilogy(distance, Throughput_10_persen_body_shadow,'g-o','MarkerIndices',1:40:length(distance))

grid on
xlabel('Distance AP-ST (m)')
ylabel('Average Throughput (bps)') 
legend('Tx=10mW, Throughput XL 5%','XL 7%',...
    'XL 10%','without XL 5%',...
    'without XL 7%','without XL 10%'); 
axis([distance(1) distance(length(distance)) 0.5*1e5 2*1e5])
%axis([distance(1) distance(length(distance)) 1*1e4 0.5*1e6])
%title("ini semuanya mengalami body shadow");

%% fungsi pathloss untuk rooftop
function hasil = path_loss_normal(distance)
    % rumus path loss 802.11ah 15m rooftop
    hasil = 8 + 37.6 * log10(distance);
end

%% fungsi pathloss untuk body shadow
function hasil = path_loss_body_shadow(distance)
    % rumus path loss bedside
    d1 = 0.5; % 50 cm. jarak dari sensor ke body edge
    S = 8.68;
    C = 64.7;
    n = 2; % free space
    hasil = 10*(n-2)*log(d1)+20*log10(distance)+S+C+20*log10(900/2400);
end

%% fungsi PDR
function hasil = PDR(PER,i)
   hasil = (1-PER) * PER^(i-1);
end

%% fungsi PDR
function hasil = PDR_array(PER,i)
   hasil = (1-PER) .* PER.^(i-1);
end

%% fungsi  T_backoff
function hasil = T_backoff(CW_min, CW_max, T_Slot, i, m)
   if i < m
       hasil = (2^(i-1)*(CW_min + 1) - 1) / 2 * T_Slot;
   else
       hasil = CW_max / 2 * T_Slot;
   end
end