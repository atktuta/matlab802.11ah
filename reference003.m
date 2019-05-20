% Berdasarkan script dari dsplog, coba membuat
% analisis dan simulasi BER vsEbN0 untuk OFDM-BPSK
% di channel AWGN dan Rayleigh.
% OFDM mengikuti standar 802.11

clear all
clc

% hitung SNR threshold. Dari Ferrand, 2013
% SER = Q(sqrt(k*SNR)), BPSK k = 1
N = 100; % number of symbol per packet
SER_SNR_threshold = 1 - 10^(-0.3/N); % -0.3 berasal dari log10(1/2)
akar_SNR_threshold = qfuncinv(SER_SNR_threshold);
SNR_threshold = akar_SNR_threshold^2;

% berikut standar 802.11
nFFT = 64; % fft size
nDSC = 52; % number of data subcarriers
nBitPerSym = 52; % number of bits per OFDM symbol (same as the number of subcarriers for BPSK)

% symbol yang dimaksud di sini adalah symbol OFDM
nSym = 10^4; % number of symbols

% hitung SNR dari jarak
jarak = 1:100:801;
transmission_power_10mW = -20; % 10 mW = 10 dBm = -20 dBW

% total noise, dari Domazetovic, 2017
noise_AWGN = -145.22; % dB
noise_figure = 5; % dB
antenna_gain = 3; % dB
noise_AWGN_figure_temp_fading_gain_dll = noise_AWGN + noise_figure -...
    antenna_gain; % dB

SNR_normal = transmission_power_10mW - ...   % dB
    path_loss_rooftop_macro_deployment(jarak) - ...
    noise_AWGN_figure_temp_fading_gain_dll;
SNR_normalLin = 10.^(SNR_normal./10);

bandwidth_MCS0 = 2; % 1 MHz
bit_rate_MCS0 = 0.65; % coded 1/2
bit_rate_MCS0 = bit_rate_MCS0*2; % uncoded
EbN0dB = SNR_normal-10 * log10(bit_rate_MCS0 / bandwidth_MCS0);

EbN0dB = [0:10]; % bit to noise ratio
EbN0dB = [-10:8]; % bit to noise ratio
% untuk simulasi menggunakan EsN0
% karena ada energi tambahan di OFDM yang tidak ada di BPSK biasa.
% yaitu untuk bit SC yang tidak terpakai (nDSC/nFFT), dan cyclic prefix(64/80)
EsN0dB = EbN0dB + 10*log10(nDSC/nFFT) + 10*log10(64/80); % converting to symbol to noise ratio

for ii = 1:length(EbN0dB)

   % Transmitter
   ipBit = rand(1,nBitPerSym*nSym) > 0.5; % random 1's and 0's
   ipMod = 2*ipBit-1; % BPSK modulation 0 --> -1, 1 --> +1
   % serial to parallel
   ipMod = reshape(ipMod,nBitPerSym,nSym).'; % grouping into multiple symbols

   % Assigning modulated symbols to subcarriers from [-26 to -1, +1 to +26]
   % dimasukkan ke bins, ada 64 titik.
   xF = [zeros(nSym,6) ipMod(:,[1:nBitPerSym/2]) zeros(nSym,1) ...
       ipMod(:,[nBitPerSym/2+1:nBitPerSym]) zeros(nSym,5)] ;
   
   % Taking FFT, the term (nFFT/sqrt(nDSC)) is for normalizing the power of transmit symbol to 1 
   xt = (nFFT/sqrt(nDSC))*ifft(fftshift(xF.')).';

   % Appending cylic prefix
   % sehingga ukurannya jadi 64+16 = 80. 
   % CP sebesar 16. seperempat dari 64 (standar 802.11)   
   xt = [xt(:,[49:64]) xt];
   
   % Rayleigh channel
   % kalau 1 tap tinggal dikali nantinya.
   % kalau multitap harus konvolusi
   % multipath channel
   nTap = 10;
   ht = 1/sqrt(2)*1/sqrt(nTap)*(randn(nSym,nTap) + 1i*randn(nSym,nTap));
   % computing and storing the frequency response of the channel, for use at recevier
   hF = fftshift(fft(ht,64,2));
   % convolution of each symbol with the random channel
   parfor jj = 1:nSym
      xht(jj,:) = conv(ht(jj,:),xt(jj,:));
   end
   xt2 = xht;
   
   % Concatenating multiple symbols to form a long vector
   % parallel to serial
   xt = reshape(xt.',1,nSym*80);
   xt2 = reshape(xt2.',1,nSym*(80+nTap-1));
   % Gaussian noise of unit variance, 0 mean
   % sehingga harus dibagi sqrt(2)
   nt = 1/sqrt(2)*(randn(1,nSym*80) + 1i*randn(1,nSym*80));
   nt2 = 1/sqrt(2)*[randn(1,nSym*(80+nTap-1)) + 1i*randn(1,nSym*(80+nTap-1))];
   
   % Adding noise, the term sqrt(80/64) is to account for the wasted energy due to cyclic prefix
   yt = sqrt(80/64)*xt + 10^(-EsN0dB(ii)/20)*nt;
   yt2 = sqrt(80/64)*xt2 + 10^(-EsN0dB(ii)/20)*nt2;
   
   % Receiver
   % serial to parallel
   yt = reshape(yt.',80,nSym).'; % formatting the received vector into symbols
   yt2 = reshape(yt2.',80+nTap-1,nSym).'; % formatting the received vector into symbols
   yt = yt(:,[17:80]); % removing cyclic prefix
   yt2 = yt2(:,[17:80]); % removing cyclic prefix

   % converting to frequency domain
   yF = (sqrt(nDSC)/nFFT)*fftshift(fft(yt.')).';
   yF2 = (sqrt(nDSC)/nFFT)*fftshift(fft(yt2.')).';
   % khusus untuk Rayleigh, di equalization dulu
   % asumsi equalizer nya sudah diketahui
   yF2 = yF2./hF;   
   % buang yang diluar [-26 to -1, +1 to +26] 
   yMod = yF(:,[6+[1:nBitPerSym/2] 7+[nBitPerSym/2+1:nBitPerSym] ]); 
   yMod2 = yF2(:,[6+[1:nBitPerSym/2] 7+[nBitPerSym/2+1:nBitPerSym] ]); 

   % BPSK demodulation
   % +ve value --> 1, -ve value --> -1
   ipModHat = 2*floor(real(yMod/2)) + 1;
   ipModHat2 = 2*floor(real(yMod2/2)) + 1;
   ipModHat(find(ipModHat>1)) = +1;
   ipModHat2(find(ipModHat2>1)) = +1;
   ipModHat(find(ipModHat<-1)) = -1;
   ipModHat2(find(ipModHat2<-1)) = -1;

   % converting modulated values into bits
   ipBitHat = (ipModHat+1)/2;
   ipBitHat2 = (ipModHat2+1)/2;
   % parallel to serial
   ipBitHat = reshape(ipBitHat.',nBitPerSym*nSym,1).';
   ipBitHat2 = reshape(ipBitHat2.',nBitPerSym*nSym,1).';

   % counting the errors
   nErr(ii) = size(find(ipBitHat - ipBit),2);
   nErr2(ii) = size(find(ipBitHat2 - ipBit),2);

end

simBer = nErr/(nSym*nBitPerSym);
simBer2 = nErr2/(nSym*nBitPerSym);
theoryBer = (1/2)*erfc(sqrt(10.^(EbN0dB/10)));
theoryBerRayleigh = 0.5.*(1-sqrt(10.^(EbN0dB/10)./(10.^(EbN0dB/10)+1)));

figure
semilogy(EbN0dB,theoryBer,'bs-','LineWidth',1);
hold on
semilogy(EbN0dB,simBer,'mx-','LineWidth',1);
semilogy(EbN0dB,theoryBerRayleigh,'gs-','LineWidth',1);
semilogy(EbN0dB,simBer2,'rx-','LineWidth',1);
axis([EbN0dB(1) EbN0dB(length(EbN0dB)) 10^-3 1])
grid on
legend('theory AWGN', 'simulation AWGN', 'theory Rayleigh',...
    'simulation Rayleigh');
xlabel('Eb/No, dB')
ylabel('Bit Error Rate')
title('BER for BPSK using OFDM')

% Theory_PER = 1 - (1-simBer2).^4;
% figure
% semilogy(jarak,theoryBer,'bs-','LineWidth',1);
% hold on
% semilogy(jarak,simBer,'mx-','LineWidth',1);
% semilogy(jarak,theoryBerRayleigh,'gs-','LineWidth',1);
% semilogy(jarak,simBer2,'rx-','LineWidth',1);
% semilogy(jarak,Theory_PER,'kx-','LineWidth',1);
% axis([0 10 10^-2 1])
% grid on
% legend('theory AWGN', 'simulation AWGN', 'theory Rayleigh',...
%     'simulation Rayleigh');
% xlabel('Jarak (m)')
% ylabel('Bit Error Rate')
% title('BER for BPSK using OFDM')

