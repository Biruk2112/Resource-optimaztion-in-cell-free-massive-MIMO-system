clear all
close all
clc
load('Pilot_Assignment_32')
figure(1)
subplot(2,1,1), plot(x_GUEs_DL_MMSE_CF,y_GUEs_DL_MMSE_CF,'b-','LineWidth',1.2);
hold on
subplot(2,1,1), plot(x_GUEs_DL_MMSE_CF_Perfect_CSI,y_GUEs_DL_MMSE_CF_Perfect_CSI,'b--','LineWidth',1.2);
hold on
subplot(2,1,1), plot(x_GUEs_DL_MMSE_CF_UB,y_GUEs_DL_MMSE_CF_UB,'b-.','LineWidth',1.2);

load('Pilot_Assignment_32_random')
hold on
subplot(2,1,1), plot(x_GUEs_DL_MMSE_CF,y_GUEs_DL_MMSE_CF,'r-','LineWidth',1.2);
hold on
subplot(2,1,1), plot(x_GUEs_DL_MMSE_CF_Perfect_CSI,y_GUEs_DL_MMSE_CF_Perfect_CSI,'r--','LineWidth',1.2);
hold on
subplot(2,1,1), plot(x_GUEs_DL_MMSE_CF_UB,y_GUEs_DL_MMSE_CF_UB,'r-.','LineWidth',1.2);

grid
legend('UPA, ICSI',...
    'UPA, PCSI',...
     'UPA, UB',...
    'Location', 'SouthEastOutside');
xlabel('Rate per user [Mbit/s]');
ylabel('CDF');
title('Non-Orthogonal VS Random Pilot Assignment ')
xlim([1 60]);

load('Pilot_Assignment_Orthogonal_32')
subplot(2,1,2), plot(x_GUEs_DL_MMSE_CF,y_GUEs_DL_MMSE_CF,'b-','LineWidth',1.2);
hold on
subplot(2,1,2), plot(x_GUEs_DL_MMSE_CF_Perfect_CSI,y_GUEs_DL_MMSE_CF_Perfect_CSI,'b--','LineWidth',1.2);
hold on
subplot(2,1,2), plot(x_GUEs_DL_MMSE_CF_UB,y_GUEs_DL_MMSE_CF_UB,'b-.','LineWidth',1.2);
hold on

load('Pilot_Assignment_Orthogonal_32_random')
subplot(2,1,2), plot(x_GUEs_DL_MMSE_CF,y_GUEs_DL_MMSE_CF,'r-','LineWidth',1.2);
hold on
subplot(2,1,2), plot(x_GUEs_DL_MMSE_CF_Perfect_CSI,y_GUEs_DL_MMSE_CF_Perfect_CSI,'r--','LineWidth',1.2);
hold on
subplot(2,1,2), plot(x_GUEs_DL_MMSE_CF_UB,y_GUEs_DL_MMSE_CF_UB,'r-.','LineWidth',1.2);
grid on
legend('UPA, ICSI',...
    'UPA, PCSI',...
     'UPA, UB',...
    'Location', 'SouthEastOutside');
xlabel('Rate per user [Mbit/s]');
ylabel('CDF');
title('Random Orthogonal Pilot Assignment VS Orthogonal Pilot Assignment ')
xlim([0 60]);

% hold on
% %%%WFPC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  load('Results_Rician_DL_WFPC_UC_CF_mM_noUAVs_partial_16Pilot')
% figure(2)
% subplot(2,1,2),plot(x_GUEs_DL_MMSE_CF_WFPC,y_GUEs_DL_MMSE_CF_WFPC,'b-','LineWidth',1.2);
% hold on
% subplot(2,1,2),plot(x_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC,y_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC,'b--','LineWidth',1.2);
% hold on
% subplot(2,1,2),plot(x_GUEs_DL_MMSE_CF_UB_WFPC,y_GUEs_DL_MMSE_CF_UB_WFPC,'b-.','LineWidth',1.2);
% grid on
% legend('WFPC, ICSI',...
%     'WFPC, PCSI',...
%      'WFPC, UB',...
%     'Location', 'SouthEastOutside');
% xlabel('Rate per user [Mbit/s]');
% ylabel('CDF');
% title('32 Orthogonal Pilots ')
% xlim([0 60]);
% 
% 
% 
%  load('Results_Rician_DL_WFPC_UC_CF_mM_noUAVs_partial_64Pilot')
% subplot(2,1,2),plot(x_GUEs_DL_MMSE_CF_WFPC,y_GUEs_DL_MMSE_CF_WFPC,'b-','LineWidth',1.2);
% hold on
% subplot(2,1,2),plot(x_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC,y_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC,'b--','LineWidth',1.2);
% hold on
% subplot(2,1,2),plot(x_GUEs_DL_MMSE_CF_UB_WFPC,y_GUEs_DL_MMSE_CF_UB_WFPC,'b-.','LineWidth',1.2);
% grid on
% legend('WFPC, ICSI',...
%     'WFPC, PCSI',...
%      'WFPC, UB',...
%     'Location', 'SouthEastOutside');
% xlabel('Rate per user [Mbit/s]');
% ylabel('CDF');
% title('64 Orthogonal Pilots ')
% xlim([0 60]);
% 
