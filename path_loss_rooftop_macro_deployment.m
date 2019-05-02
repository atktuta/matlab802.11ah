%% fungsi pathloss untuk rooftop
function hasil = path_loss_rooftop_macro_deployment(distance)
    % rumus path loss 802.11ah 15m rooftop
    hasil = 8 + 37.6 * log10(distance);
end