clear all
close all
clc

load('Pilot_Assignment_Orthogonal_32')
figure(1)
hold on; box on;
set(gca,'fontsize',11);
% Uniform Power allocatoin in UB
p=plot(x_GUEs_DL_MMSE_CF_UB,y_GUEs_DL_MMSE_CF_UB,'Color','#4DBEEE','LineWidth',2);
p.LineStyle='-';
hold on
% Uniform Power allocatoin in Perfect CSI
p=plot(x_GUEs_DL_MMSE_CF_Perfect_CSI,y_GUEs_DL_MMSE_CF_Perfect_CSI,'Color','#4DBEEE','LineWidth',2);
p.LineStyle='-.';
hold on
% Uniform Power allocatoin in Imperfect CSI
p=plot(x_GUEs_DL_MMSE_CF,y_GUEs_DL_MMSE_CF,'Color','#4DBEEE','LineWidth',2);
p.LineStyle='--';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% WFPA Power allocatoin in UB
p=plot(x_GUEs_DL_MMSE_CF_UB_WFPC,y_GUEs_DL_MMSE_CF_UB_WFPC,'LineWidth',2);
p.Color='#7E2F8E';
p.LineStyle='-';

% WFPA Power allocatoin in Perfect CSI
p=plot(x_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC,y_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC,'LineWidth',2);
p.Color='#7E2F8E';
p.LineStyle='-.';
hold on
% WFPA Power allocatoin in Imperfect CSI
p=plot(x_GUEs_DL_MMSE_CF_WFPC,y_GUEs_DL_MMSE_CF_WFPC,'LineWidth',2);
p.Color='#7E2F8E';
p.LineStyle='--';
grid on
legend('Uniform Power control, UB', 'Uniform Power control, PCSI', 'Uniform Power control, ICSI',...
     'WF Power control, UB', 'WF Power control, PCSI [32]','WF Power control, ICSI',...
    'Location', 'SouthEast');
xlabel('Rate Per-UE [Mbit/s]');
ylabel('CDF');
xlim([0 100]);