clear;
clc;

transmission_power_10mW =  10; % 10 mW = 10 dBm = -20 dBW
 
distance = 1:10:600; % meter

minimum_sensitivity_BPSK_MCS10 =  ones(1,length(distance))*-98; % standard page 302

path_loss = path_loss_rooftop_macro_deployment(distance);

%signal power level yang diterima (Power Rx)
signal_power_rx = transmission_power_10mW - path_loss;

figure
plot(distance,signal_power_rx)
hold on
plot(distance,minimum_sensitivity_BPSK_MCS10)
grid on
axis([0 600 -120 -20])
xlabel('Distance between Tx and Rx (m)')
ylabel('Signal Power (dBm)')
legend('Power vs distance for transmission power=10mW path loss',...
    'Receiver min input sensitivity for MCS10');
