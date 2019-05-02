clear;
clc;

% transmission_power_10mW =  10; % 10 mW = 10 dBm = -20 dBW
% transmission_power_200mW =  23.010299957; % 200 mW = 23 dBm
% transmission_power_1000mW =  30; % 1000 mW = 30 dBm

transmission_power = [10 23.013 30];

distance = 1:10:1600; % meter

minimum_sensitivity_BPSK_MCS10 =  ones(1,length(distance))*-98; % standard page 302
minimum_sensitivity_BPSK_MCS0_BW1 =  ones(1,length(distance))*-95; % standard page 302
minimum_sensitivity_BPSK_MCS0_BW2 =  ones(1,length(distance))*-92; % standard page 302

path_loss = path_loss_rooftop_macro_deployment(distance);

%signal power level yang diterima (Power Rx)
signal_power_rx = zeros(length(transmission_power),length(path_loss));
for i=1:length(transmission_power)
    signal_power_rx(i,:) = transmission_power(i) - path_loss;
end

figure
plot(distance,signal_power_rx(1,:))
hold on
plot(distance,signal_power_rx(2,:))
plot(distance,signal_power_rx(3,:))
plot(distance,minimum_sensitivity_BPSK_MCS10)
plot(distance,minimum_sensitivity_BPSK_MCS0_BW1)
plot(distance,minimum_sensitivity_BPSK_MCS0_BW2)
grid on
axis([0 1600 -120 0])
xlabel('Distance between Tx and Rx (m)')
ylabel('Signal Power (dBm)')
legend('Power vs distance for transmission power=10mW','Power vs distance for transmission power=200mW',...
    'Power vs distance for transmission power=1000mW','Receiver min input sensitivity for MCS10',...
    'Receiver min input sensitivity for MCS0(1 MHz BW)','Receiver min input sensitivity for MCS0(2 MHz BW)');
