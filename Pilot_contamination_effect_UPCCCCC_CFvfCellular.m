clear all
close all
clc

%%%Imperfect%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 load('Pilot_Assignment_random_Conta_8')
figure
hold on; box on; 
set(gca,'fontsize',11);

p=plot(x_GUEs_DL_MMSE_CF_Perfect_CSI,y_GUEs_DL_MMSE_CF_Perfect_CSI,'Color','#4DBEEE','LineWidth',2);
p.LineStyle='-.';

 hold on
p=plot(x_GUEs_DL_MMSE_mM,y_GUEs_DL_MMSE_mM,'LineWidth',2);
p.Color='#7E2F8E';
p.LineStyle='-.';
hold on


load('Pilot_Assignment_random_ContaEf_20')

p=plot(x_GUEs_DL_MMSE_CF_Perfect_CSI,y_GUEs_DL_MMSE_CF_Perfect_CSI,'Color','#4DBEEE','LineWidth',2);
p.LineStyle='-';
hold on
p=plot(x_GUEs_DL_MMSE_mM_Perfect_CSI,y_GUEs_DL_MMSE_mM_Perfect_CSI,'LineWidth',2);
p.Color='#7E2F8E';
p.LineStyle='-';

grid on
legend('Cell-free mMIMO, PCSI, $\tau =8$',...
     'massive MIMO, PCSI, $\tau =8 $',...
     'Cell-free mMIMO, PCSI, $\tau =20$',...
     'massive MIMO, PCSI, $\tau = 20$',...
    'Interpreter','latex', 'Location','SouthEast');
xlabel('Rate per user [Mbit/s]');
ylabel('CDF');
% title('Pilot contamination effect ')
 xlim([0 60]);



% 
% p=plot(x_GUEs_DL_MMSE_CF_Perfect_CSI,y_GUEs_DL_MMSE_CF_Perfect_CSI,'Color','#4DBEEE','LineWidth',2);
% p.LineStyle='-.';
% % p.Marker='d';
% 
% p=plot(x_GUEs_DL_MMSE_mM_Perfect_CSI,y_GUEs_DL_MMSE_mM_Perfect_CSI,'LineWidth',2);
% p.Color='#7E2F8E';
% p.LineStyle='-.';
% hold on
% 
%  
% % p.Marker='o';
% % p.MarkerSize = 5;
% % p.MarkerIndices = 1:5:length(y_GUEs_DL_MMSE_CF);
%  hold on
% p=plot(x_GUEs_DL_MMSE_CF_Perfect_CSI,y_GUEs_DL_MMSE_CF_Perfect_CSI,'Color','#4DBEEE','LineWidth',2);
% p.LineStyle='-.';
% % p.Marker='d';
%  
% % p.Marker='o';
% % p.MarkerSize = 5;
% % p.MarkerIndices = 1:5:length(5);
% % hold on
% p=plot(x_GUEs_DL_MMSE_mM_Perfect_CSI,y_GUEs_DL_MMSE_mM_Perfect_CSI,'LineWidth',2);
% p.Color='#7E2F8E';
% p.LineStyle='-.';
% 
% 
% 
% %  writematrix(y_GUEs_DL_MMSE_CF_WFPC,'Export_Casey.csv')
% %  writematrix(y_GUEs_DL_MMSE_CF_UB_WFPC,'Export_CaseUB.csv')
% %  writematrix(y_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC,'Export_CasePerf.csv')
% 
% %   writematrix(x_GUEs_DL_MMSE_CF_WFPC,'EExport_Casey.csv','WriteMode','append')
% %   writematrix(x_GUEs_DL_MMSE_CF_UB_WFPC,'EExport_CaseUB.csv','WriteMode','append')
% %   writematrix(x_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC,'EExport_CasePerf.csv','WriteMode','append')
% 
% 
% % 
% % grid on
% % legend('WFPC-, ICSI',...
% %     'WFPC-, PCSI',...
% %      'WFPC-, UB',...
% %     'Location', 'SouthEastOutside');
% % xlabel('Rate per user [Mbit/s]');
% % ylabel('CDF');
% % title('32 Orthogonal Pilots ')
% % xlim([0 60]);

%%%%________________________________________________________________________________
%%%%%%%%%%%%%%%%_____________________________________________________________________________
%%%%%%%%%%%%%%%%%%%%%__________________________________________________________________________________
clear all
close all
clc

%%%Imperfect%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  load('Pilot_Assignment_random_Conta_8')
% figure
% hold on; box on; 
% set(gca,'fontsize',11);
% 
% p=plot(x_GUEs_DL_MMSE_CF_Perfect_CSI,y_GUEs_DL_MMSE_CF_Perfect_CSI,'Color','#4DBEEE','LineWidth',2);
% p.LineStyle='-.';
% 
%  hold on
% p=plot(x_GUEs_DL_MMSE_mM,y_GUEs_DL_MMSE_mM,'LineWidth',2);
% p.Color='#7E2F8E';
% p.LineStyle='-.';
% hold on
% 
% 
% load('Pilot_Assignment_random_ContaEf_20')
% 
% p=plot(x_GUEs_DL_MMSE_CF_Perfect_CSI,y_GUEs_DL_MMSE_CF_Perfect_CSI,'Color','#4DBEEE','LineWidth',2);
% p.LineStyle='-';
% hold on
% p=plot(x_GUEs_DL_MMSE_mM_Perfect_CSI,y_GUEs_DL_MMSE_mM_Perfect_CSI,'LineWidth',2);
% p.Color='#7E2F8E';
% p.LineStyle='-';
% 
% grid on
% legend('Cell-free mMIMO, PCSI, $\tau =8$',...
%      'massive MIMO, PCSI, $\tau =8 $',...
%      'Cell-free mMIMO, PCSI, $\tau =20$',...
%      'massive MIMO, PCSI, $\tau = 20$',...
%     'Interpreter','latex', 'Location','SouthEast');
% xlabel('Rate per user [Mbit/s]');
% ylabel('CDF');
% % title('Pilot contamination effect ')
%  xlim([0 60]);





