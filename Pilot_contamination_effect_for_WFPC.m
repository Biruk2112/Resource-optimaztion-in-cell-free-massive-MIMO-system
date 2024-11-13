clear all
close all
clc
%%%WFPC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
load('Pilot_Assignment_Orthogonal_32')
hold on; box on; 
set(gca,'fontsize',11);
p=plot(x_GUEs_DL_MMSE_CF_WFPC,y_GUEs_DL_MMSE_CF_WFPC,'LineWidth',2);
p.Color='#4DBEEE';
p.LineStyle='-';

hold on
load('Pilot_Assignment_Orthogonal_32_random')
p=plot(x_GUEs_DL_MMSE_CF_WFPC,y_GUEs_DL_MMSE_CF_WFPC,'LineWidth',2);
p.Color='#4DBEEE';
p.LineStyle='-.';
hold on

 load('Pilot_Assignment_random')
p=plot(x_GUEs_DL_MMSE_CF_WFPC,y_GUEs_DL_MMSE_CF_WFPC,'LineWidth',2);
p.Color='#7E2F8E';
p.LineStyle='-.';

hold on
load('Pilot_Assignment_32')
p=plot(x_GUEs_DL_MMSE_CF_WFPC,y_GUEs_DL_MMSE_CF_WFPC,'LineWidth',2);
p.Color='#7E2F8E';
p.LineStyle='-';
hold on
hold on
% grid on
% %%%%%%%%%
% figure(2) 
load('Pilot_Assignment_Orthogonal_32')
% hold on; box on; 
% set(gca,'fontsize',11);
p=plot(x_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC,y_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC,'LineWidth',2);
p.Color='#4DBEEE';
p.LineStyle='-';

hold on
load('Pilot_Assignment_Orthogonal_32_random')
p=plot(x_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC,y_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC,'LineWidth',2);
p.Color='#4DBEEE';
p.LineStyle='-.';
hold on

 load('Pilot_Assignment_random')
p=plot(x_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC,y_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC,'LineWidth',2);
p.Color='#7E2F8E';
p.LineStyle='-.';

hold on
load('Pilot_Assignment_32')
p=plot(x_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC,y_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC,'LineWidth',2);
p.Color='#7E2F8E';
p.LineStyle='-';
hold on

grid on
legend('Orthogonal pilot assignment ',...
      'Greedy random pilot assignment',...
      'Random pilot assignment',...
    'Non-orthogonal pilot assignment',...
    'Interpreter','latex', 'Location','SouthEast');
xlabel('Rate per user [Mbit/s]');
ylabel('CDF of WFPC');
title('Pilot Assignment ')
 xlim([0 60]);

