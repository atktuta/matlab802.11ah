% coba duplikasi OFDM sederhana, dengan parameter dari 802.11a
% namun tanpa SIGNAL dan tanpa PREAMBLE. hanya OFDM simbol.
% sumber dari Mathuranathan, simulation of digital comm system matlab
% dengan beberapa tambahan modifikasi.
% sinyal OFDM dilewatkan AWGN dan rayleigh.
clear all
clc

% untuk tiap Eb_No dalam dB
eb_no_db = 0:2:10;

% ukuran paket
N = 96; % bytes. usahakan kelipatan 6
N_bits = N * 8; % ukuran paket dalam bit

% untuk BPSK, satu simbol OFDM membawa 48 bit (sama dengan jumlah sc)
% sehingga butuh berapa simbol untuk membawa N_bits?
jumlah_simbol_per_paket = N_bits/48;

% misal untuk simulasi kita ingin xx paket untuk tiap Eb_No
jumlah_paket = 100;

% maka total simbol yang mau disimulasikan untuk tiap Eb_No
% total_symbol = jumlah_paket * jumlah_simbol_per_paket;

% untuk perhitungan noise di simulasi, pakainya Es_No
% 52 adalah sub carrier terpakai, 64 adalah total sub carrier.
% 80 adalah total sub carrier ditambah cyclic prefix
es_no_db = eb_no_db + 10*log10(52/64) + 10*log10(64/80);

% inisialisasi error
total_error_bit_tiap_eb_no = zeros(1,length(eb_no_db));
total_error_packet_tiap_eb_no = zeros(1,length(eb_no_db));
total_error = zeros(length(eb_no_db), jumlah_paket, jumlah_simbol_per_paket);
error_buat_paket = 0;
% untuk tiap Eb_No
for i=1:length(eb_no_db)
    % untuk tiap paket
    for j=1:jumlah_paket
    % untuk tiap symbol OFDM
        for k=1:jumlah_simbol_per_paket
            % ada stream data, misal 48 bit. biar pas untuk BPSK.
            % kalau ngga pas bisa di-padding sih
            randomBinaryNumbers = randi([0, 1], 1, 48); % acak 0 atau 1, matrix 1x48

            % kemudian dipecah menggunakan serial to parallel
            % BPSK berarti dipecah satu-satu
            % sekaligus sisipkan pilot.
            % Assigning subcarriers from 1 to 26 (mapped to 1-26 of IFFT input)
            % and -26 to -1 (mapped to 38 to 63 of IFFT input); Nulls from 27 to 37
            % and at 0 position.
            % pilot ada di posisi -21, -7, 7 dan 21

            serialtoparallel=[zeros(1,1) randomBinaryNumbers(1:6) ones(1,1) ... % tambah pilot
                randomBinaryNumbers(7:19) ones(1,1) ... % tambah pilot
                randomBinaryNumbers(20:24) zeros(1,11) ...
                randomBinaryNumbers(25:29) ones(1,1) ... % tambah pilot
                randomBinaryNumbers(30:42) ones(1,1) ... % tambah pilot
                randomBinaryNumbers(43:48)];
            serialtoparallel = serialtoparallel.'; % transpose tanpa complex conjugate

            % kemudian di mapping
            % kalau BPSK gampang, bit 0 adalah 1 + j0, bit 1 adalah -1 + j0
            % terus tidak perlu dinormalisasi
            M = modem.pskmod(2);
            freq_domain_transmitter = modulate(M,serialtoparallel);

            % ubah ke domain waktu dan 
            % normalisasi (total subcarrier / akar subcarrier  yg digunakan)
            % tujuannya agar power transmit symbol = 1
            time_domain_transmitter = 64/sqrt(52)*ifft(freq_domain_transmitter); 

            % tambahkan cyclic prefix dari 16 titik fft terakhir ke depan
            % jadilah sebuah OFDM symbol
            ofdm_symbol=[time_domain_transmitter(49:64); time_domain_transmitter]; 

            % parallel to serial
            ofdm_symbol = ofdm_symbol.'; % transpose tanpa complex conjugate

            %kemudian digabung dengan simbol berikutnya,misal ada 4 simbol

            % kalau tanpa channel
            % received_ofdm_symbol = ofdm_symbol;

            % lewatkan channel
            noise = 1/sqrt(2)*(randn(1,length(ofdm_symbol))+1i*randn(1,length(ofdm_symbol)));

            %--------------Channel Modeling AWGN ----------------
            Ebs_No = es_no_db(i);
            received_ofdm_symbol = sqrt(80/64)*ofdm_symbol + 10^(-Ebs_No/20)*noise;

            % rayleigh fading
            % y = hx + n
            % h berupa taps
            % first we create taps for each symbol. single tap
            taps = 1/sqrt(2)*(randn(1,length(ofdm_symbol))+1i*randn(1,length(ofdm_symbol)));
            received_rayleigh_ofdm_symbol = abs(taps).*ofdm_symbol + 10^(-Ebs_No/20)*noise;

            % serial to parallel lagi
            received_ofdm_symbol = received_rayleigh_ofdm_symbol.'; 

            % hapus cyclic prefix, yaitu 16 titik fft yang di depan
            time_domain_receiver = received_ofdm_symbol(17:80);

            % kembalikan ke domain frequency
            freq_domain_receiver = sqrt(52)/64*fft(time_domain_receiver);

            % demapping
            Z=modem.pskdemod(2); %demodulation object
            demod_data = demodulate(Z,freq_domain_receiver); 

            % ambil datanya, yang sebanyak 48 bit. Pilot dianggap error-free
            received_binary_number = [demod_data(2:7); demod_data(9:21); demod_data(23:27);...
                demod_data(39:43); demod_data(45:57); demod_data(59:64)];

            received_binary_number = received_binary_number.';

            [num_error, rasio_error] = biterr(randomBinaryNumbers,received_binary_number);
            % ini adalah error tiap simbol, di dalam paket ke j, 
            % untuk ebno ke i
            total_error(i, j, k) = num_error; 
            
            % PER untuk packet size 96 bytes.
            % jadi tiap 16 simbol = 1 packet. cek ada error ngga
            % jadi tiap 800 simbol ada 50 packet.        
        end
        
    end

     bit_error_rate_the(i)=(1/2)*erfc(sqrt(10.^(eb_no_db(i)/10))); %Same as BER for BPSK over AWGN
%     %bit_error_rate_the(i) = qfunc(sqrt(2*10.^(eb_no_db(i)/10))); %Same as BER for BPSK over AWGN
     EbN0Lin = 10.^(eb_no_db(i)/10);
     bit_error_rate_the_rayleigh(i) = 0.5.*(1-sqrt(EbN0Lin./(EbN0Lin+1)));
%     % jadi tiap 800 simbol ada 50 packet.
%     packet_error_rate_sim(i) = total_error_packet_tiap_eb_no(i) / jumlah_paket;
    a(i) = sum(sum(total_error(i,:,:)));
end

bit_error_rate_sim = a/(jumlah_paket*jumlah_simbol_per_paket*48);
packet_error_rate_sim = 1 - (1-bit_error_rate_sim).^N;

% figure
% plot(eb_no_db,log10(bit_error_rate_sim),'r-o');
% grid on;
% hold on;
% plot(eb_no_db,log10(bit_error_rate_the),'k*');

figure
semilogy(eb_no_db, bit_error_rate_sim,'r-o');
grid on;
hold on;
semilogy(eb_no_db, bit_error_rate_the,'g-*');
semilogy(eb_no_db, bit_error_rate_the_rayleigh,'b-x');
title('BER Vs EbNodB for OFDM with BPSK modulation over AWGN-Rayleigh');
xlabel('Eb/N0 (dB)');ylabel('BER');legend('simulated','theoretical AWGN', ...
    'theoretical Rayleigh');

figure
semilogy(eb_no_db, packet_error_rate_sim,'r-o');
grid on;
title('PER Vs EbNodB for OFDM with BPSK modulation over AWGN-Rayleigh');
xlabel('Eb/N0 (dB)');ylabel('PER');legend('simulated');
axis([0 20 1*1e-2 1]);