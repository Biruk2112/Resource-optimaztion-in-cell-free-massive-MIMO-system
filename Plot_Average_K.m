%% Main Channels generation
clear all
close all
clc

%% Downlink

figure
load('Pilot_Assignment_Orthogonal_64');
plot(Pt_dB_DL,R_GUEs_DL_MMSE_CF_Perfect_CSI,'r-','LineWidth',2);
hold on
load('Pilot_Assignment_Orthogonal_32');
plot(Pt_dB_DL,R_GUEs_DL_MMSE_CF_WFPC,'r-','LineWidth',2);
grid
legend('Our Rate, K=20',...
    'Ngo Rate, K=20',...
    'Our Rate, K=40',...
    'Ngo Rate, K=40',...
    'Our Rate, K=60',...
    'Ngo Rate, K=60',...
    'Location', 'NorthWest');
xlabel('P_T [dBW]');
ylabel('Average Rate per user [Mbit/s]');
title('Downlink')

