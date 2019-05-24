%% fungsi hitung throughput dari PER
function hasil = hitung_throughput_80211ah(PER)
    m = 6; % maximum number of backoff stages and corresponds to 6 (i.e. CW_max = 2^6 CW_min )
    mpdu_size = 475; % 475 Bytes

    sigma_11ah = 1561/physconst('LightSpeed'); % jarak maks/kecepatan cahaya

    DIFS_11ah = 264e-6; % 264 us untuk 11ah
    SIFS_11ah = 160e-6; % 160 us untuk 11ah
    CW_min = 15;
    CW_max = 1023;
    T_SLOT_11ah = 52e-6; % 52 us untuk 11a

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
    T_BACKOFF = 0;
    for j=1:batas_atas
        T_BACKOFF = T_BACKOFF + PDR_array(PER,j)*T_backoff(CW_min, ...
            CW_max, T_SLOT_11ah, j, m);
    end
    T_message = T_message_wo_T_BACKOFF + T_BACKOFF;      
    hasil = (mpdu_size * 8)./T_message .* (1-PER);   
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