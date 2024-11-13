clear all
close all
clc
% load('Pilot_Assignment_Orthogonal_32')
% figure(1)
% subplot(2,1,1), plot(x_GUEs_DL_MMSE_CF,y_GUEs_DL_MMSE_CF,'r-','LineWidth',1.2);
% hold on
% subplot(2,1,1), plot(x_GUEs_DL_MMSE_CF_Perfect_CSI,y_GUEs_DL_MMSE_CF_Perfect_CSI,'r--','LineWidth',1.2);
% hold on
% subplot(2,1,1), plot(x_GUEs_DL_MMSE_CF_UB,y_GUEs_DL_MMSE_CF_UB,'r-.','LineWidth',1.2);
% grid
% legend('UPA, ICSI',...
%     'UPA, PCSI',...
%      'UPA, UB',...
%     'Location', 'SouthEastOutside');
% xlabel('Rate per user [Mbit/s]');
% ylabel('CDF');
% title('32 Orthogonal Pilots')
% xlim([1 60]);
% 
% load('Pilot_Assignment_Orthogonal_64')
% subplot(2,1,2), plot(x_GUEs_DL_MMSE_CF,y_GUEs_DL_MMSE_CF,'b-','LineWidth',1.2);
% hold on
% subplot(2,1,2), plot(x_GUEs_DL_MMSE_CF_Perfect_CSI,y_GUEs_DL_MMSE_CF_Perfect_CSI,'b--','LineWidth',1.2);
% hold on
% subplot(2,1,2), plot(x_GUEs_DL_MMSE_CF_UB,y_GUEs_DL_MMSE_CF_UB,'b-.','LineWidth',1.2);
% grid on
% legend('UPA, ICSI',...
%     'UPA, PCSI',...
%      'UPA, UB',...
%     'Location', 'SouthEastOutside');
% xlabel('Rate per user [Mbit/s]');
% ylabel('CDF');
% title('64 Orthogonal Pilots ')
% xlim([0 60]);



%%%WFPC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 load('Pilot_Assignment_random_Conta_8')
figure
hold on; box on; 
set(gca,'fontsize',11);
p=plot(x_GUEs_DL_MMSE_CF_WFPC,y_GUEs_DL_MMSE_CF_WFPC,'LineWidth',2);
p.Color='#7E2F8E';
p.LineStyle='--';
hold on

p=plot(x_GUEs_DL_MMSE_CF_UB_WFPC,y_GUEs_DL_MMSE_CF_UB_WFPC,'LineWidth',2);
p.Color='#7E2F8E';
p.LineStyle='-';

hold on
p=plot(x_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC,y_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC,'LineWidth',2);
p.Color='#7E2F8E';
p.LineStyle='-.';
hold on


 load('Pilot_Assignment_random_ContaEf_20')
p=plot(x_GUEs_DL_MMSE_CF_WFPC,y_GUEs_DL_MMSE_CF_WFPC,'LineWidth',2);
p.Color='#4DBEEE';
p.LineStyle='--';
hold on

p=plot(x_GUEs_DL_MMSE_CF_UB_WFPC,y_GUEs_DL_MMSE_CF_UB_WFPC,'LineWidth',2);
p.Color='#4DBEEE';
p.LineStyle='-';

p=plot(x_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC,y_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC,'LineWidth',2);
p.Color='#4DBEEE';
p.LineStyle='-.';
grid on
legend('WFPC, ICSI, $\tau =8$',...
      'WFPC, UB , $\tau =8 $',...
    'WFPC, PCSI, $\tau =8 $',...
     'WFPC, ICSI, $\tau =20 $',...
      'WFPC, UB , $\tau =20 $',...
    'WFPC, PCSI, $\tau =20 $',...
    'Interpreter','latex', 'Location','SouthEast');
xlabel('Rate per user [Mbit/s]');
ylabel('CDF');
 title('Pilot contamination effect ')
 xlim([0 100]);


% m = [ x_GUEs_DL_MMSE_CF_WFPC x_GUEs_DL_MMSE_CF_UB_WFPC x_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC ;y_GUEs_DL_MMSE_CF_WFPC y_GUEs_DL_MMSE_CF_UB_WFPC x_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC;735 1 5;264 2 7;346 9 7] 
%  writematrix(y_GUEs_DL_MMSE_CF_WFPC,'Export_Casey.csv')

% writematrix(x_GUEs_DL_MMSE_CF_WFPC,'Export_Casey.csv','WriteMode','append')

% X=get(x_GUEs_DL_MMSE_CF_WFPC, 'XData');
% Y=get(y_GUEs_DL_MMSE_CF_WFPC, 'YData');
% 
% m= cell2mat(Y');
% maxNumCol = max
% writematrix(m,'Export_Casey00.csv')

% 
% grid on
% legend('WFPC-, ICSI',...
%     'WFPC-, PCSI',...
%      'WFPC-, UB',...
%     'Location', 'SouthEastOutside');
% xlabel('Rate per user [Mbit/s]');
% ylabel('CDF');
% title('32 Orthogonal Pilots ')
% xlim([0 60]);