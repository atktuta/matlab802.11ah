% Berdasarkan script dari dsplog, coba membuat
% analisis dan simulasi BER vsEbN0 untuk OFDM-BPSK
% di channel AWGN dan Rayleigh.
% OFDM mengikuti standar 802.11

clear all

% berikut standar 802.11
nFFT = 64; % fft size
nDSC = 52; % number of data subcarriers
nBitPerSym = 52; % number of bits per OFDM symbol (same as the number of subcarriers for BPSK)

% symbol yang dimaksud di sini adalah symbol OFDM
nSym = 10^4; % number of symbols

EbN0dB = [0:10]; % bit to noise ratio
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

   % Concatenating multiple symbols to form a long vector
   % parallel to serial
   xt = reshape(xt.',1,nSym*80);

   % Gaussian noise of unit variance, 0 mean
   nt = 1/sqrt(2)*[randn(1,nSym*80) + 1i*randn(1,nSym*80)];

   % Adding noise, the term sqrt(80/64) is to account for the wasted energy due to cyclic prefix
   yt = sqrt(80/64)*xt + 10^(-EsN0dB(ii)/20)*nt;

   % Receiver
   % serial to parallel
   yt = reshape(yt.',80,nSym).'; % formatting the received vector into symbols
   yt = yt(:,[17:80]); % removing cyclic prefix

   % converting to frequency domain
   yF = (sqrt(nDSC)/nFFT)*fftshift(fft(yt.')).';
   % buang yang diluar [-26 to -1, +1 to +26] 
   yMod = yF(:,[6+[1:nBitPerSym/2] 7+[nBitPerSym/2+1:nBitPerSym] ]); 

   % BPSK demodulation
   % +ve value --> 1, -ve value --> -1
   ipModHat = 2*floor(real(yMod/2)) + 1;
   ipModHat(find(ipModHat>1)) = +1;
   ipModHat(find(ipModHat<-1)) = -1;

   % converting modulated values into bits
   ipBitHat = (ipModHat+1)/2;
   % parallel to serial
   ipBitHat2 = reshape(ipBitHat.',nBitPerSym*nSym,1).';

   % counting the errors
   nErr(ii) = size(find(ipBitHat2 - ipBit),2);

end

simBer = nErr/(nSym*nBitPerSym);
theoryBer = (1/2)*erfc(sqrt(10.^(EbN0dB/10)));

figure
semilogy(EbN0dB,theoryBer,'bs-','LineWidth',2);
hold on
semilogy(EbN0dB,simBer,'mx-','LineWidth',2);
axis([0 10 10^-5 1])
grid on
legend('theory AWGN', 'simulation AWGN');
xlabel('Eb/No, dB')
ylabel('Bit Error Rate')
title('BER for BPSK using OFDM')