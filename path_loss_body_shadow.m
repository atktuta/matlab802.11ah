%% fungsi pathloss untuk body shadow
function hasil = path_loss_body_shadow(distance)
    % rumus path loss bedside
    d1 = 0.5; % 50 cm. jarak dari sensor ke body edge
    S = 8.68;
    C = 64.7;
    n = 2; % free space
    hasil = 10*(n-2)*log(d1)+20*log10(distance)+S+C+20*log10(900/2400);
end