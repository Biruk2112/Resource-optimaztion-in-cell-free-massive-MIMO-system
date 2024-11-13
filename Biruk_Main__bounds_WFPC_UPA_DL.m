%% Main Channels generation
tic
clear all
close all
clc

% M=100;
% N_AP=4;
% N_AP_mM=100;
% M_mM=4;
% K_GUE=60; % total number of users

M=40; %Number of APs
N_AP=2; %Number of antennas per AP
N_AP_mM=40;
M_mM=2;
K_GUE=64; % total number of users
% K_GUE NN=[0:10:K_GUE, 200];
 NN=0:K_GUE;
dim_N=length(NN);
K=K_GUE;
f=1.9e9;
lambda=3e8/f;
antenna_spacing=lambda/2;
W=20e6;
h_AP=15;
h_MS=1.65;
D=1000; %

% Noise variance
noise_figure=9; % noise figure in dB
N0=-174;% PSD noise in dBm/Hz
 noise_variance=W*10^(0.1*noise_figure)*10^(0.1*N0)*10^-3; % F*N0*B
%  noise_variance=1e-18;

P_DL=200e-3;
P_UL=100e-3;

Pt_dB_DL=10*log10(P_DL);
Pt_dB_DL_mM=10*log10(P_DL*M/M_mM);
Pt_dB_UL=10*log10(P_UL);

% Pt_dB_DL=10*log10(P_DL/N_AP);
% Pt_dB_DL_mM=10*log10(P_DL/N_AP_mM*M/M_mM);
% Pt_dB_UL=10*log10(P_UL);
% save('Results_Multiple_Antenna_CDF_UL_UC_SIC_N2');
toc
%% Pilot_Assignment
tau_p=8;
%    Phi=Pilot_Assignment_Orthogonal_random(K,tau_p);%Orthogonal Pilot Assignment
%   Phi=Pilot_Assignment_Orthogonal(K,tau_p); %random Orthogonal Pilot
%  Assignment [I've finishied this]
%   Phi=Pilot_Assignment_Random(K,tau_p); %Random Pilot Assignment [I've finishied this]
  Phi=Pilot_Assignment(K,tau_p); % Non-Orth Pilot Assignment


P_training=100e-3;
P_max_dBm_training=10*log10(P_training/1e-3);
P_max_dBm=10*log10(P_UL/1e-3);
P0_dBm_GUE=-45;
P0_dBm_UAV=-45;
alpha=0.5;
%FPC
% P0_dBm_GUE=-10;
 
%% Rate 
tau_c=200; % length of coherence time in samples
N_Channels=80;

alpha_power1=0.1;
alpha_power2=0.2;


%%% Uniform power allocation in Cell free imperfect CSI
R_GUEs_DL_MMSE_CF=[];

%%% Uniform power allocation in Cell free Perfect CSI
R_GUEs_DL_MMSE_CF_Perfect_CSI=[];

R_GUEs_DL_MMSE_CF_UB=[];

% WFPC Power allocatoin in Imperfect CSI
R_GUEs_DL_MMSE_CF_WFPC=[];

% WFPC Power allocatoin in Perfect CSI
R_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC=[];

% WFPC Power allocatoin in UB
R_GUEs_DL_MMSE_CF_UB_WFPC=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WFPC
% R_GUEs_DL_MMSE_CF_WFPC_Alpha1=zeros(K,dim_N);
R_GUEs_DL_MMSE_CF_Perfect_CSI_WFPCC= zeros(K_GUE,dim_N);
R_GUEs_DL_MMSE_CF_WFPCC=zeros(K_GUE,dim_N);

% media_R_k_zf_upa_ns = zeros(1,length(snr_db));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Imperfect massive MIMO
R_GUEs_DL_MMSE_mM=[];

%Perfect massive MIMO
R_GUEs_DL_MMSE_mM_Perfect_CSI=[];

%UB massive MIMO
R_GUEs_DL_MMSE_mM_UB=[];

N_slow_fading=1;

for nn=1:N_Channels
    
    %% Scenario cell FREE
     AP_pos=[D*rand(M,2) h_AP*ones(M,1)];
    AP_orientations=pi/2*ones(M,1);
    
    MS_pos=[D*rand(K_GUE,2) h_MS*ones(K_GUE,1)];
    
    
    K_factors_GUE=zeros(M,K_GUE); % the channels of the GUEs are purely Rayleigh
    
    K_factors=K_factors_GUE;

    UAV_GUE_pos=MS_pos;

    %% Calculate the Betas
    %% %% Calculate the Betas
    [~,Beta_GUE]=microwaveChannel_GUE_wrapped_area(f,AP_pos,MS_pos,D,N_AP); % I only calculate the Beta of the GUEs
    Beta=Beta_GUE;
    %% Channels UAVs GUEs
    [G,A_vectors,G_matrices]=RicianChannel_wrapped_area_GUE_UAV(Beta,AP_pos,UAV_GUE_pos,D,N_AP,AP_orientations,antenna_spacing,K_factors,lambda);
   
    %% Uniform power allocation
    
    % Channel estimation
   Eta_pilot=tau_p*P_training*ones(K,1);
    [Gamma,D_matrices,G_estimate]=Channel_Estimation_MMSE_Multiantenna_Rice(G,Phi,Eta_pilot,noise_variance,G_matrices);
    
    % Cell_Free downlink lower bound
    % Imperfect CSI uniform power allocation
    K_ASS_CF=cell(M,1);
    
    for m=1:M
        K_ASS_CF{m,1}=1:K;
    end
    
    M_ASS_CF=cell(K,1);
    for k=1:K
        M_ASS_CF{k,1}=1:M;
    end
    
     Eta_DL_CF=Eta_DL_Multiantenna(Gamma,Pt_dB_DL,K_ASS_CF);
    [SE_DL_CF]=Rate_Downlink_Multiantenna_Rice(Eta_DL_CF,Eta_pilot,Beta,Gamma,G_matrices,noise_variance,Phi,M_ASS_CF,A_vectors,D_matrices,K_factors);
    
    R_GUEs_DL_MMSE_CF=[R_GUEs_DL_MMSE_CF; SE_DL_CF(1:K_GUE)];
    % Cell_Free downlink upper bound    
    Q_DL_norm=Precoding_Downlink(G_estimate);
    Eta_DL_norm=Eta_Q_DL_norm(K_ASS_CF,Q_DL_norm,Pt_dB_DL);
    SE_DL_CF=Rate_Downlink_Multiantenna_Upper_Bound2(G,Eta_DL_norm,Q_DL_norm,noise_variance);

    R_GUEs_DL_MMSE_CF_UB=[R_GUEs_DL_MMSE_CF_UB; SE_DL_CF(1:K_GUE)];
    
    % Cell Free Perfect CSI
    % Perfect CSI uniform power allocation
    Q_DL_norm=Precoding_Downlink(G);
    Eta_DL_norm=Eta_Q_DL_norm(K_ASS_CF,Q_DL_norm,Pt_dB_DL);
    [SE_DL_CF]=Rate_Downlink_Multiantenna_Upper_Bound2(G,Eta_DL_norm,Q_DL_norm,noise_variance);

    R_GUEs_DL_MMSE_CF_Perfect_CSI=[R_GUEs_DL_MMSE_CF_Perfect_CSI; SE_DL_CF(1:K_GUE)];
    
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Waterfilling power allocation
    %Imperfect WFPC
    Eta_DL_CF_WFPC=Waterfilling_power_control_generic_Gamma(Gamma,Pt_dB_DL,K_ASS_CF,noise_variance);
    [SE_DL_CF]=Rate_Downlink_Multiantenna_Rice(Eta_DL_CF_WFPC,Eta_pilot,Beta,Gamma,G_matrices,noise_variance,Phi,M_ASS_CF,A_vectors,D_matrices,K_factors);
    [R_CF_DL_WFPC_temp]=Rate_Downlink_Multiantenna_Rice(Eta_DL_CF_WFPC,Eta_pilot,Beta,Gamma,G_matrices,noise_variance,Phi,M_ASS_CF,A_vectors,D_matrices,K_factors);
    
    R_GUEs_DL_MMSE_CF_WFPC=[R_GUEs_DL_MMSE_CF_WFPC; SE_DL_CF(1:K_GUE)];
%     R_GUEs_DL_MMSE_CF_WFPCC=[R_GUEs_DL_MMSE_CF_WFPCC; R_CF_DL_WFPC_temp(1:NN)];
    
    App3=zeros(K_GUE,dim_N);
    for nn=1:dim_N
%         App1(:,nn)=R_CF_DL_temp;
%         App2(:,nn)=R_CF_DL_SR_temp;
        App3(:,nn)=R_CF_DL_WFPC_temp;
    end
%      R_CF_DL=R_CF_DL+App1;
%     R_CF_DL_SR=R_CF_DL_SR+App2;
    R_GUEs_DL_MMSE_CF_WFPCC=R_GUEs_DL_MMSE_CF_WFPCC+App3;
    
    % Cell_Free downlink upper bound    
    Q_DL_norm=Precoding_Downlink(G_estimate);
    SE_DL_CF=Rate_Downlink_Multiantenna_Upper_Bound2(G,Eta_DL_CF_WFPC,Q_DL_norm,noise_variance);

    R_GUEs_DL_MMSE_CF_UB_WFPC=[R_GUEs_DL_MMSE_CF_UB_WFPC; SE_DL_CF(1:K_GUE)];
    
    % Cell Free Perfect CSI
    % WFPA Power allocatoin in Perfect CSI
    Q_DL_norm=Precoding_Downlink(G);
    trace_G_matrices=zeros(M,K); % Gamma in the case of Perfect CSI
    for mm=1:M
        for kk=1:K
            trace_G_matrices(mm,kk)=trace(G_matrices(:,:,mm,kk));
        end
    end
    
    Eta_DL_CF_WFPC_PCSI=Waterfilling_power_control_generic_Gamma(trace_G_matrices,Pt_dB_DL,K_ASS_CF,noise_variance);
    [SE_DL_CF]=Rate_Downlink_Multiantenna_Upper_Bound2(G,Eta_DL_CF_WFPC_PCSI,Q_DL_norm,noise_variance);
    [R_CF_Perfect_DL_WFPA_temp]=Rate_Downlink_Multiantenna_Upper_Bound2(G,Eta_DL_CF_WFPC_PCSI,Q_DL_norm,noise_variance);

    R_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC=[R_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC; SE_DL_CF(1:K_GUE)];
%     R_GUEs_DL_MMSE_CF_Perfect_CSI_WFPCC=[R_GUEs_DL_MMSE_CF_Perfect_CSI_WFPCC; R_CF_Perfect_DL_WFPA_temp(1:NN)];
    
%     App1=zeros(K_GUE,dim_N);
%     App2=zeros(K_GUE,dim_N);
    App3=zeros(K_GUE,dim_N);
    for nn=1:dim_N
%        App1(:,nn)=R_CF_Perfect_DL_temp;
%        App2(:,nn)=R_CF_Perfect_DL_SR_temp;
%         App2(:,nn)=R_CF_Perfect_DL_SR_temp;
       App3(:,nn)= R_CF_Perfect_DL_WFPA_temp;
       
    end

% R_CF_Perfect_DL=R_CF_Perfect_DL+App1;
%     R_CF_Perfect_DL_SR=R_CF_Perfect_DL_SR+App2;
    R_GUEs_DL_MMSE_CF_Perfect_CSI_WFPCC=R_GUEs_DL_MMSE_CF_Perfect_CSI_WFPCC+App3;

   
    
    
    %% Scenario Massive MIMO
    %I put 4 fixed BS with N_AP antenna each (very large N_AP)
    AP_pos_mM=zeros(M_mM,3);
    AP_pos_mM(1,:)=[D/4 D/4 h_AP];
    AP_pos_mM(2,:)=[3*D/4 D/4 h_AP];
    AP_pos_mM(3,:)=[3*D/4 3*D/4 h_AP];
    AP_pos_mM(4,:)=[D/4 3*D/4 h_AP];
    %% Calculate the Betas
    [~,Beta_GUE_mM]=microwaveChannel_GUE_wrapped_area(f,AP_pos_mM,MS_pos,D,N_AP_mM); % I only calculate the Betas of the GUEs
    Beta_mM=Beta_GUE_mM;
    %% Channels UAVs GUEs
    [G_mM,A_vectors_mM,G_matrices_mM]=RicianChannel_wrapped_area_GUE_UAV(Beta_mM,AP_pos_mM,UAV_GUE_pos,D,N_AP_mM,AP_orientations,antenna_spacing,K_factors,lambda);
    
    
    %% Uniform power allocation
    
    % Channel estimation
    Eta_pilot=P_training*ones(K,1);
    [Gamma_mM,D_matrices_mM,G_estimate_mM]=Channel_Estimation_MMSE_Multiantenna_Rice(G_mM,Phi,Eta_pilot,noise_variance,G_matrices_mM);
    % Massive MIMO downlink
    [K_ASS_mM,M_ASS_mM]=User_Association_N_User_centric(Beta_mM,1);
 
    Eta_DL_mM=Eta_DL_Multiantenna_Uniformly(Gamma_mM,Pt_dB_DL_mM,K_ASS_mM);
    [SE_DL_mM]=Rate_Downlink_Multiantenna_Rice(Eta_DL_mM,Eta_pilot,Beta_mM,Gamma_mM,G_matrices_mM,noise_variance,Phi,M_ASS_mM,A_vectors_mM,D_matrices_mM,K_factors);
    
    R_GUEs_DL_MMSE_mM=[R_GUEs_DL_MMSE_mM; SE_DL_mM(1:K_GUE)];
    %upper bound    
    Q_DL_norm=Precoding_Downlink(G_estimate_mM);
    Eta_DL_norm=Eta_Q_DL_norm_Uniformly(K_ASS_mM,Q_DL_norm,Pt_dB_DL_mM);
    SE_DL_mM=Rate_Downlink_Multiantenna_Upper_Bound2(G_mM,Eta_DL_norm,Q_DL_norm,noise_variance);

    R_GUEs_DL_MMSE_mM_UB=[R_GUEs_DL_MMSE_mM_UB; SE_DL_mM(1:K_GUE)];
    
    % Perfect CSI
    
    Q_DL_norm=Precoding_Downlink(G_mM);
    Eta_DL_norm=Eta_Q_DL_norm_Uniformly(K_ASS_mM,Q_DL_norm,Pt_dB_DL_mM);
    [SE_DL_mM]=Rate_Downlink_Multiantenna_Upper_Bound2(G_mM,Eta_DL_norm,Q_DL_norm,noise_variance);

    R_GUEs_DL_MMSE_mM_Perfect_CSI=[R_GUEs_DL_MMSE_mM_Perfect_CSI; SE_DL_mM(1:K_GUE)];
    nn,
    
%   
end

% Imperfect CSI uniform power allocation
R_GUEs_DL_MMSE_CF=W*(1-tau_p/tau_c)/2*R_GUEs_DL_MMSE_CF;

% sum_R_CF_Perfect_DL=W*sum(R_CF_Perfect_DL,1)/N_Channels/(1e6);

% Perfect CSI uniform power allocation
R_GUEs_DL_MMSE_CF_Perfect_CSI=W*1/2*R_GUEs_DL_MMSE_CF_Perfect_CSI;

R_GUEs_DL_MMSE_CF_UB=W*(1-tau_p/tau_c)/2*R_GUEs_DL_MMSE_CF_UB;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% WFPA Power allocatoin in Imperfect CSI
R_GUEs_DL_MMSE_CF_WFPC=W*(1-tau_p/tau_c)/2*R_GUEs_DL_MMSE_CF_WFPC;
R_GUEs_DL_MMSE_CF_WFPCC=W*(1-tau_p/tau_c)/2*R_GUEs_DL_MMSE_CF_WFPCC;

% WFPA Power allocatoin in Perfect CSI
R_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC=W*1/2*R_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC;
R_GUEs_DL_MMSE_CF_Perfect_CSI_WFPCC=W*1/2*R_GUEs_DL_MMSE_CF_Perfect_CSI_WFPCC;

R_GUEs_DL_MMSE_CF_UB_WFPC=W*(1-tau_p/tau_c)/2*R_GUEs_DL_MMSE_CF_UB_WFPC;

R_GUEs_DL_MMSE_mM=W*(1-tau_p/tau_c)/2*R_GUEs_DL_MMSE_mM;

% Perfect CSI uniform power allocation for massive MIMO
R_GUEs_DL_MMSE_mM_Perfect_CSI=W*1/2*R_GUEs_DL_MMSE_mM_Perfect_CSI;

% Imperfect CSI uniform power allocation for massive MIMO
R_GUEs_DL_MMSE_mM_UB=W*(1-tau_p/tau_c)/2*R_GUEs_DL_MMSE_mM_UB;

% imperfect CSI uniform power allocation
[ x_GUEs_DL_MMSE_CF,y_GUEs_DL_MMSE_CF] = Empirical_CDF(R_GUEs_DL_MMSE_CF./(1e6));

% Perfect CSI uniform power allocation
[ x_GUEs_DL_MMSE_CF_Perfect_CSI,y_GUEs_DL_MMSE_CF_Perfect_CSI] = Empirical_CDF(R_GUEs_DL_MMSE_CF_Perfect_CSI./(1e6));

% UB uniform power allocation
[ x_GUEs_DL_MMSE_CF_UB,y_GUEs_DL_MMSE_CF_UB] = Empirical_CDF(R_GUEs_DL_MMSE_CF_UB./(1e6));

%
% % WFPA Power allocatoin in Imperfect CSI
[ x_GUEs_DL_MMSE_CF_WFPC,y_GUEs_DL_MMSE_CF_WFPC] = Empirical_CDF(R_GUEs_DL_MMSE_CF_WFPC./(1e6));
[ x_GUEs_DL_MMSE_CF_WFPCC,y_GUEs_DL_MMSE_CF_WFPCC] = Empirical_CDF(R_GUEs_DL_MMSE_CF_WFPCC./(1e6));

% WFPA Power allocatoin in Perfect CSI
[ x_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC,y_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC] = Empirical_CDF(R_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC./(1e6));
[ x_GUEs_DL_MMSE_CF_Perfect_CSI_WFPCC,y_GUEs_DL_MMSE_CF_Perfect_CSI_WFPCC] = Empirical_CDF(R_GUEs_DL_MMSE_CF_Perfect_CSI_WFPCC./(1e6));

[ x_GUEs_DL_MMSE_CF_UB_WFPC,y_GUEs_DL_MMSE_CF_UB_WFPC] = Empirical_CDF(R_GUEs_DL_MMSE_CF_UB_WFPC./(1e6));

[ x_GUEs_DL_MMSE_mM,y_GUEs_DL_MMSE_mM] = Empirical_CDF(R_GUEs_DL_MMSE_mM./(1e6));

%Perfect CSI for massive MIMO
[ x_GUEs_DL_MMSE_mM_Perfect_CSI,y_GUEs_DL_MMSE_mM_Perfect_CSI] = Empirical_CDF(R_GUEs_DL_MMSE_mM_Perfect_CSI./(1e6));

%Imperfect CSI for massive MIMO
[ x_GUEs_DL_MMSE_mM_UB,y_GUEs_DL_MMSE_mM_UB] = Empirical_CDF(R_GUEs_DL_MMSE_mM_UB./(1e6));

% R_GUEs_DL_MMSE_CF=R_GUEs_DL_MMSE_CF./(1e6);
% R_GUEs_DL_MMSE_CF_Perfect_CS=R_GUEs_DL_MMSE_CF_Perfect_CSI./(1e6);
% R_GUEs_DL_MMSE_CF_WFPC=R_GUEs_DL_MMSE_CF_WFPC./(1e6);
% R_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC=R_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC./(1e6);
% WFPA maximization
% sum_R_CF_Perfect_DL_WFPC=R_GUEs_DL_MMSE_CF_Perfect_CSI_WFPCC./(1e6);
% sum_R_CF_DL_WFPC=R_GUEs_DL_MMSE_CF_WFPCC./(1e6);

sum_R_CF_Perfect_DL_WFPC=sum(R_GUEs_DL_MMSE_CF_Perfect_CSI_WFPCC,1)/N_Channels/(1e6);
sum_R_CF_DL_WFPC=sum(R_GUEs_DL_MMSE_CF_WFPCC,1)/N_Channels/(1e6);

%%%%% Pilot Assignment %%%%%%%%%%
 %%%%% Pilot Assignment %%%%%%%%%%
% %      save('Pilot_Assignment_Orthogonal_random'); %I've finishied this
     save('Pilot_Assignment_Orthogonal_60AP');  %I've finishied this
% %       save('Pilot_Assignment_Orthogonal_32_random'); %I've finishied this
% %     save('Pilot_Assignment_Random');
  save('Pilot_Assignment_random_Conta_8'); %
% %       save('Pilot_Assignment_Orthogonal_644'); % this a 64 PA
%     save('Pilot_Assignment_322X_random');
 

toc

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
     'WF Power control, UB', 'WF Power control, PCSI [49]','WF Power control, ICSI',...
    'Location', 'SouthEast');
xlabel('Rate Per-UE [Mbit/s]');
ylabel('CDF');


figure(2)
hold on; box on;
set(gca,'fontsize',10);
grid on
p=plot(x_GUEs_DL_MMSE_CF,y_GUEs_DL_MMSE_CF,'Color','#4DBEEE','LineWidth',2);
p.LineStyle='-';
% p.Marker='o';
% p.MarkerSize = 5;
% p.MarkerIndices = 1:5:length(y_GUEs_DL_MMSE_CF);
 hold on
p=plot(x_GUEs_DL_MMSE_CF_Perfect_CSI,y_GUEs_DL_MMSE_CF_Perfect_CSI,'Color','#4DBEEE','LineWidth',2);
p.LineStyle='-.';
% p.Marker='d';
 hold on
p=plot(x_GUEs_DL_MMSE_mM,y_GUEs_DL_MMSE_mM,'LineWidth',2);
p.Color='#7E2F8E';
p.LineStyle='-';
% p.Marker='o';
% p.MarkerSize = 5;
% p.MarkerIndices = 1:5:length(5);
% hold on
p=plot(x_GUEs_DL_MMSE_mM_Perfect_CSI,y_GUEs_DL_MMSE_mM_Perfect_CSI,'LineWidth',2);
p.Color='#7E2F8E';
p.LineStyle='-.';
grid on
legend('UPA Cell-free mMIMO, ICSI', 'UPA Cell-free mMIMO, PCSI', 'UPA massive MIMO, ICSI','UPA massive MIMO, PCSI',...
    'Location', 'SouthEast');
xlabel('Rate Per-UE (Mbit/s)');
ylabel('CDF');
title('CF mMIMO VS Massive MIMO in UPA')

% xlim([0 60]);


 %%%%%       $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 figure(3)
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
grid on
legend('Uniform Power control, UB, $ tau = 8$', 'Uniform Power control, PCSI, $ tau = 8$', 'Uniform Power control, ICSI, $ tau = 8$',...
    'Location','SouthEast','TextColor','blue');
xlabel('Rate Per-UE [Mbit/s]');
ylabel('CDF');

figure(4)
hold on; box on;
set(gca,'fontsize',11);
% WFPA Power allocatoin in UB
p=plot(x_GUEs_DL_MMSE_CF_UB_WFPC,y_GUEs_DL_MMSE_CF_UB_WFPC,'LineWidth',2);
p.Color='#4DBEEE';
p.LineStyle='-';

% WFPA Power allocatoin in Perfect CSI
p=plot(x_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC,y_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC,'LineWidth',2);
p.Color='#4DBEEE';
p.LineStyle='-.';
hold on
% WFPA Power allocatoin in Imperfect CSI
p=plot(x_GUEs_DL_MMSE_CF_WFPC,y_GUEs_DL_MMSE_CF_WFPC,'LineWidth',2);
p.Color='#4DBEEE';
p.LineStyle='--';
grid on
legend('WF Power control, UB, $ tau = 8$', 'WF Power control, PCSI, $ tau = 8$','WF Power control, ICSI, $ tau = 8$',...
    'Location','SouthEast','TextColor','blue');
xlabel('Rate Per-UE [Mbit/s]');
ylabel('CDF');
 
figure(5)
plot(Pt_dB_DL,R_GUEs_DL_MMSE_CF,'b','LineWidth',2);
hold on
plot(Pt_dB_DL,R_GUEs_DL_MMSE_CF_Perfect_CSI,'b--','LineWidth',2);
hold on
plot(Pt_dB_DL,R_GUEs_DL_MMSE_CF_WFPC,'r','LineWidth',2);
hold on
plot(Pt_dB_DL,R_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC,'r--','LineWidth',2);
grid on
AX=legend('UPA, ICSI','UPA, PCSI','WFPC, ICSI','WFPC, PCSI','Location', 'SouthEast');
LEG = findobj(AX,'type','text');
set(LEG,'FontSize',12)
xlabel('P_T [dBW]');
ylabel('Sum Rate [Mbit/s]');

 figure(6) 
 hold on; box on;
set(gca,'fontsize',11);
plot(NN,sum_R_CF_Perfect_DL_WFPC,'b','LineWidth',2);
% hold on
plot(NN,sum_R_CF_DL_WFPC,'r','LineWidth',2);
% hold on
grid on

AX=legend('WFPC, PCSI','WFPC, ICSI','Location', 'SouthEast');
LEG = findobj(AX,'type','text');
set(LEG,'FontSize',12)

xlabel('K');
ylabel('Sum Rate [Mbit/s]');
title('SE vs Number of users ($K$)','Interpreter','latex');
xlabel('Number of User per APs, ($K$)','Interpreter','latex');

           %%%%%%%%%%%%%%%%%%%%%_____________________________________%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 

hold on


M=40; %Number of APs
N_AP=2; %Number of antennas per AP
N_AP_mM=40;
M_mM=2;
K_GUE=64; % total number of users
% K_GUE NN=[0:10:K_GUE, 200];
 NN=0:K_GUE;
dim_N=length(NN);
K=K_GUE;
f=1.9e9;
lambda=3e8/f;
antenna_spacing=lambda/2;
W=20e6;
h_AP=15;
h_MS=1.65;
D=1000; %

% Noise variance
noise_figure=9; % noise figure in dB
N0=-174;% PSD noise in dBm/Hz
 noise_variance=W*10^(0.1*noise_figure)*10^(0.1*N0)*10^-3; % F*N0*B
%  noise_variance=1e-18;

P_DL=200e-3;
P_UL=100e-3;

Pt_dB_DL=10*log10(P_DL);
Pt_dB_DL_mM=10*log10(P_DL*M/M_mM);
Pt_dB_UL=10*log10(P_UL);

% Pt_dB_DL=10*log10(P_DL/N_AP);
% Pt_dB_DL_mM=10*log10(P_DL/N_AP_mM*M/M_mM);
% Pt_dB_UL=10*log10(P_UL);
% save('Results_Multiple_Antenna_CDF_UL_UC_SIC_N2');
toc
%% Pilot_Assignment
tau_p=20;
%    Phi=Pilot_Assignment_Orthogonal_random(K,tau_p);%Orthogonal Pilot Assignment
%   Phi=Pilot_Assignment_Orthogonal(K,tau_p); %random Orthogonal Pilot
%  Assignment [I've finishied this]
%    Phi=Pilot_Assignment_Random(K,tau_p); %Random Pilot Assignment [I've finishied this]
   Phi=Pilot_Assignment(K,tau_p); % Non-Orth Pilot Assignment


P_training=100e-3;
P_max_dBm_training=10*log10(P_training/1e-3);
P_max_dBm=10*log10(P_UL/1e-3);
P0_dBm_GUE=-45;
P0_dBm_UAV=-45;
alpha=0.5;
%FPC
% P0_dBm_GUE=-10;
 
%% Rate 
tau_c=200; % length of coherence time in samples
N_Channels=80;

alpha_power1=0.1;
alpha_power2=0.2;


%%% Uniform power allocation in Cell free imperfect CSI
R_GUEs_DL_MMSE_CF=[];

%%% Uniform power allocation in Cell free Perfect CSI
R_GUEs_DL_MMSE_CF_Perfect_CSI=[];

R_GUEs_DL_MMSE_CF_UB=[];

% WFPC Power allocatoin in Imperfect CSI
R_GUEs_DL_MMSE_CF_WFPC=[];

% WFPC Power allocatoin in Perfect CSI
R_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC=[];

% WFPC Power allocatoin in UB
R_GUEs_DL_MMSE_CF_UB_WFPC=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WFPC
% R_GUEs_DL_MMSE_CF_WFPC_Alpha1=zeros(K,dim_N);
R_GUEs_DL_MMSE_CF_Perfect_CSI_WFPCC= zeros(K_GUE,dim_N);
R_GUEs_DL_MMSE_CF_WFPCC=zeros(K_GUE,dim_N);

% media_R_k_zf_upa_ns = zeros(1,length(snr_db));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Imperfect massive MIMO
R_GUEs_DL_MMSE_mM=[];

%Perfect massive MIMO
R_GUEs_DL_MMSE_mM_Perfect_CSI=[];

%UB massive MIMO
R_GUEs_DL_MMSE_mM_UB=[];

N_slow_fading=1;

for nn=1:N_Channels
    
    %% Scenario cell FREE
     AP_pos=[D*rand(M,2) h_AP*ones(M,1)];
    AP_orientations=pi/2*ones(M,1);
    
    MS_pos=[D*rand(K_GUE,2) h_MS*ones(K_GUE,1)];
    
    
    K_factors_GUE=zeros(M,K_GUE); % the channels of the GUEs are purely Rayleigh
    
    K_factors=K_factors_GUE;

    UAV_GUE_pos=MS_pos;

    %% Calculate the Betas
    %% %% Calculate the Betas
    [~,Beta_GUE]=microwaveChannel_GUE_wrapped_area(f,AP_pos,MS_pos,D,N_AP); % I only calculate the Beta of the GUEs
    Beta=Beta_GUE;
    %% Channels UAVs GUEs
    [G,A_vectors,G_matrices]=RicianChannel_wrapped_area_GUE_UAV(Beta,AP_pos,UAV_GUE_pos,D,N_AP,AP_orientations,antenna_spacing,K_factors,lambda);
   
    %% Uniform power allocation
    
    % Channel estimation
   Eta_pilot=tau_p*P_training*ones(K,1);
    [Gamma,D_matrices,G_estimate]=Channel_Estimation_MMSE_Multiantenna_Rice(G,Phi,Eta_pilot,noise_variance,G_matrices);
    
    % Cell_Free downlink lower bound
    % Imperfect CSI uniform power allocation
    K_ASS_CF=cell(M,1);
    
    for m=1:M
        K_ASS_CF{m,1}=1:K;
    end
    
    M_ASS_CF=cell(K,1);
    for k=1:K
        M_ASS_CF{k,1}=1:M;
    end
    
     Eta_DL_CF=Eta_DL_Multiantenna(Gamma,Pt_dB_DL,K_ASS_CF);
    [SE_DL_CF]=Rate_Downlink_Multiantenna_Rice(Eta_DL_CF,Eta_pilot,Beta,Gamma,G_matrices,noise_variance,Phi,M_ASS_CF,A_vectors,D_matrices,K_factors);
    
    R_GUEs_DL_MMSE_CF=[R_GUEs_DL_MMSE_CF; SE_DL_CF(1:K_GUE)];
    % Cell_Free downlink upper bound    
    Q_DL_norm=Precoding_Downlink(G_estimate);
    Eta_DL_norm=Eta_Q_DL_norm(K_ASS_CF,Q_DL_norm,Pt_dB_DL);
    SE_DL_CF=Rate_Downlink_Multiantenna_Upper_Bound2(G,Eta_DL_norm,Q_DL_norm,noise_variance);

    R_GUEs_DL_MMSE_CF_UB=[R_GUEs_DL_MMSE_CF_UB; SE_DL_CF(1:K_GUE)];
    
    % Cell Free Perfect CSI
    % Perfect CSI uniform power allocation
    Q_DL_norm=Precoding_Downlink(G);
    Eta_DL_norm=Eta_Q_DL_norm(K_ASS_CF,Q_DL_norm,Pt_dB_DL);
    [SE_DL_CF]=Rate_Downlink_Multiantenna_Upper_Bound2(G,Eta_DL_norm,Q_DL_norm,noise_variance);

    R_GUEs_DL_MMSE_CF_Perfect_CSI=[R_GUEs_DL_MMSE_CF_Perfect_CSI; SE_DL_CF(1:K_GUE)];
    
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Waterfilling power allocation
    %Imperfect WFPC
    Eta_DL_CF_WFPC=Waterfilling_power_control_generic_Gamma(Gamma,Pt_dB_DL,K_ASS_CF,noise_variance);
    [SE_DL_CF]=Rate_Downlink_Multiantenna_Rice(Eta_DL_CF_WFPC,Eta_pilot,Beta,Gamma,G_matrices,noise_variance,Phi,M_ASS_CF,A_vectors,D_matrices,K_factors);
    [R_CF_DL_WFPC_temp]=Rate_Downlink_Multiantenna_Rice(Eta_DL_CF_WFPC,Eta_pilot,Beta,Gamma,G_matrices,noise_variance,Phi,M_ASS_CF,A_vectors,D_matrices,K_factors);
    
    R_GUEs_DL_MMSE_CF_WFPC=[R_GUEs_DL_MMSE_CF_WFPC; SE_DL_CF(1:K_GUE)];
%     R_GUEs_DL_MMSE_CF_WFPCC=[R_GUEs_DL_MMSE_CF_WFPCC; R_CF_DL_WFPC_temp(1:NN)];
    
    App3=zeros(K_GUE,dim_N);
    for nn=1:dim_N
%         App1(:,nn)=R_CF_DL_temp;
%         App2(:,nn)=R_CF_DL_SR_temp;
        App3(:,nn)=R_CF_DL_WFPC_temp;
    end
%      R_CF_DL=R_CF_DL+App1;
%     R_CF_DL_SR=R_CF_DL_SR+App2;
    R_GUEs_DL_MMSE_CF_WFPCC=R_GUEs_DL_MMSE_CF_WFPCC+App3;
    
    % Cell_Free downlink upper bound    
    Q_DL_norm=Precoding_Downlink(G_estimate);
    SE_DL_CF=Rate_Downlink_Multiantenna_Upper_Bound2(G,Eta_DL_CF_WFPC,Q_DL_norm,noise_variance);

    R_GUEs_DL_MMSE_CF_UB_WFPC=[R_GUEs_DL_MMSE_CF_UB_WFPC; SE_DL_CF(1:K_GUE)];
    
    % Cell Free Perfect CSI
    % WFPA Power allocatoin in Perfect CSI
    Q_DL_norm=Precoding_Downlink(G);
    trace_G_matrices=zeros(M,K); % Gamma in the case of Perfect CSI
    for mm=1:M
        for kk=1:K
            trace_G_matrices(mm,kk)=trace(G_matrices(:,:,mm,kk));
        end
    end
    
    Eta_DL_CF_WFPC_PCSI=Waterfilling_power_control_generic_Gamma(trace_G_matrices,Pt_dB_DL,K_ASS_CF,noise_variance);
    [SE_DL_CF]=Rate_Downlink_Multiantenna_Upper_Bound2(G,Eta_DL_CF_WFPC_PCSI,Q_DL_norm,noise_variance);
    [R_CF_Perfect_DL_WFPA_temp]=Rate_Downlink_Multiantenna_Upper_Bound2(G,Eta_DL_CF_WFPC_PCSI,Q_DL_norm,noise_variance);

    R_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC=[R_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC; SE_DL_CF(1:K_GUE)];
%     R_GUEs_DL_MMSE_CF_Perfect_CSI_WFPCC=[R_GUEs_DL_MMSE_CF_Perfect_CSI_WFPCC; R_CF_Perfect_DL_WFPA_temp(1:NN)];
    
%     App1=zeros(K_GUE,dim_N);
%     App2=zeros(K_GUE,dim_N);
    App3=zeros(K_GUE,dim_N);
    for nn=1:dim_N
%        App1(:,nn)=R_CF_Perfect_DL_temp;
%        App2(:,nn)=R_CF_Perfect_DL_SR_temp;
%         App2(:,nn)=R_CF_Perfect_DL_SR_temp;
       App3(:,nn)= R_CF_Perfect_DL_WFPA_temp;
       
    end

% R_CF_Perfect_DL=R_CF_Perfect_DL+App1;
%     R_CF_Perfect_DL_SR=R_CF_Perfect_DL_SR+App2;
    R_GUEs_DL_MMSE_CF_Perfect_CSI_WFPCC=R_GUEs_DL_MMSE_CF_Perfect_CSI_WFPCC+App3;

   
    
    
    %% Scenario Massive MIMO
    %I put 4 fixed BS with N_AP antenna each (very large N_AP)
    AP_pos_mM=zeros(M_mM,3);
    AP_pos_mM(1,:)=[D/4 D/4 h_AP];
    AP_pos_mM(2,:)=[3*D/4 D/4 h_AP];
    AP_pos_mM(3,:)=[3*D/4 3*D/4 h_AP];
    AP_pos_mM(4,:)=[D/4 3*D/4 h_AP];
    %% Calculate the Betas
    [~,Beta_GUE_mM]=microwaveChannel_GUE_wrapped_area(f,AP_pos_mM,MS_pos,D,N_AP_mM); % I only calculate the Betas of the GUEs
    Beta_mM=Beta_GUE_mM;
    %% Channels UAVs GUEs
    [G_mM,A_vectors_mM,G_matrices_mM]=RicianChannel_wrapped_area_GUE_UAV(Beta_mM,AP_pos_mM,UAV_GUE_pos,D,N_AP_mM,AP_orientations,antenna_spacing,K_factors,lambda);
    
    
    %% Uniform power allocation
    
    % Channel estimation
    Eta_pilot=P_training*ones(K,1);
    [Gamma_mM,D_matrices_mM,G_estimate_mM]=Channel_Estimation_MMSE_Multiantenna_Rice(G_mM,Phi,Eta_pilot,noise_variance,G_matrices_mM);
    % Massive MIMO downlink
    [K_ASS_mM,M_ASS_mM]=User_Association_N_User_centric(Beta_mM,1);
 
    Eta_DL_mM=Eta_DL_Multiantenna_Uniformly(Gamma_mM,Pt_dB_DL_mM,K_ASS_mM);
    [SE_DL_mM]=Rate_Downlink_Multiantenna_Rice(Eta_DL_mM,Eta_pilot,Beta_mM,Gamma_mM,G_matrices_mM,noise_variance,Phi,M_ASS_mM,A_vectors_mM,D_matrices_mM,K_factors);
    
    R_GUEs_DL_MMSE_mM=[R_GUEs_DL_MMSE_mM; SE_DL_mM(1:K_GUE)];
    %upper bound    
    Q_DL_norm=Precoding_Downlink(G_estimate_mM);
    Eta_DL_norm=Eta_Q_DL_norm_Uniformly(K_ASS_mM,Q_DL_norm,Pt_dB_DL_mM);
    SE_DL_mM=Rate_Downlink_Multiantenna_Upper_Bound2(G_mM,Eta_DL_norm,Q_DL_norm,noise_variance);

    R_GUEs_DL_MMSE_mM_UB=[R_GUEs_DL_MMSE_mM_UB; SE_DL_mM(1:K_GUE)];
    
    % Perfect CSI
    
    Q_DL_norm=Precoding_Downlink(G_mM);
    Eta_DL_norm=Eta_Q_DL_norm_Uniformly(K_ASS_mM,Q_DL_norm,Pt_dB_DL_mM);
    [SE_DL_mM]=Rate_Downlink_Multiantenna_Upper_Bound2(G_mM,Eta_DL_norm,Q_DL_norm,noise_variance);

    R_GUEs_DL_MMSE_mM_Perfect_CSI=[R_GUEs_DL_MMSE_mM_Perfect_CSI; SE_DL_mM(1:K_GUE)];
    nn,
    
%   
end

% Imperfect CSI uniform power allocation
R_GUEs_DL_MMSE_CF=W*(1-tau_p/tau_c)/2*R_GUEs_DL_MMSE_CF;

% sum_R_CF_Perfect_DL=W*sum(R_CF_Perfect_DL,1)/N_Channels/(1e6);

% Perfect CSI uniform power allocation
R_GUEs_DL_MMSE_CF_Perfect_CSI=W*1/2*R_GUEs_DL_MMSE_CF_Perfect_CSI;

R_GUEs_DL_MMSE_CF_UB=W*(1-tau_p/tau_c)/2*R_GUEs_DL_MMSE_CF_UB;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% WFPA Power allocatoin in Imperfect CSI
R_GUEs_DL_MMSE_CF_WFPC=W*(1-tau_p/tau_c)/2*R_GUEs_DL_MMSE_CF_WFPC;
R_GUEs_DL_MMSE_CF_WFPCC=W*(1-tau_p/tau_c)/2*R_GUEs_DL_MMSE_CF_WFPCC;

% WFPA Power allocatoin in Perfect CSI
R_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC=W*1/2*R_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC;
R_GUEs_DL_MMSE_CF_Perfect_CSI_WFPCC=W*1/2*R_GUEs_DL_MMSE_CF_Perfect_CSI_WFPCC;

R_GUEs_DL_MMSE_CF_UB_WFPC=W*(1-tau_p/tau_c)/2*R_GUEs_DL_MMSE_CF_UB_WFPC;

R_GUEs_DL_MMSE_mM=W*(1-tau_p/tau_c)/2*R_GUEs_DL_MMSE_mM;

% Perfect CSI uniform power allocation for massive MIMO
R_GUEs_DL_MMSE_mM_Perfect_CSI=W*1/2*R_GUEs_DL_MMSE_mM_Perfect_CSI;

% Imperfect CSI uniform power allocation for massive MIMO
R_GUEs_DL_MMSE_mM_UB=W*(1-tau_p/tau_c)/2*R_GUEs_DL_MMSE_mM_UB;

% imperfect CSI uniform power allocation
[ x_GUEs_DL_MMSE_CF,y_GUEs_DL_MMSE_CF] = Empirical_CDF(R_GUEs_DL_MMSE_CF./(1e6));

% Perfect CSI uniform power allocation
[ x_GUEs_DL_MMSE_CF_Perfect_CSI,y_GUEs_DL_MMSE_CF_Perfect_CSI] = Empirical_CDF(R_GUEs_DL_MMSE_CF_Perfect_CSI./(1e6));

% UB uniform power allocation
[ x_GUEs_DL_MMSE_CF_UB,y_GUEs_DL_MMSE_CF_UB] = Empirical_CDF(R_GUEs_DL_MMSE_CF_UB./(1e6));

%
% % WFPA Power allocatoin in Imperfect CSI
[ x_GUEs_DL_MMSE_CF_WFPC,y_GUEs_DL_MMSE_CF_WFPC] = Empirical_CDF(R_GUEs_DL_MMSE_CF_WFPC./(1e6));
[ x_GUEs_DL_MMSE_CF_WFPCC,y_GUEs_DL_MMSE_CF_WFPCC] = Empirical_CDF(R_GUEs_DL_MMSE_CF_WFPCC./(1e6));

% WFPA Power allocatoin in Perfect CSI
[ x_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC,y_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC] = Empirical_CDF(R_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC./(1e6));
[ x_GUEs_DL_MMSE_CF_Perfect_CSI_WFPCC,y_GUEs_DL_MMSE_CF_Perfect_CSI_WFPCC] = Empirical_CDF(R_GUEs_DL_MMSE_CF_Perfect_CSI_WFPCC./(1e6));

[ x_GUEs_DL_MMSE_CF_UB_WFPC,y_GUEs_DL_MMSE_CF_UB_WFPC] = Empirical_CDF(R_GUEs_DL_MMSE_CF_UB_WFPC./(1e6));

[ x_GUEs_DL_MMSE_mM,y_GUEs_DL_MMSE_mM] = Empirical_CDF(R_GUEs_DL_MMSE_mM./(1e6));

%Perfect CSI for massive MIMO
[ x_GUEs_DL_MMSE_mM_Perfect_CSI,y_GUEs_DL_MMSE_mM_Perfect_CSI] = Empirical_CDF(R_GUEs_DL_MMSE_mM_Perfect_CSI./(1e6));

%Imperfect CSI for massive MIMO
[ x_GUEs_DL_MMSE_mM_UB,y_GUEs_DL_MMSE_mM_UB] = Empirical_CDF(R_GUEs_DL_MMSE_mM_UB./(1e6));

% R_GUEs_DL_MMSE_CF=R_GUEs_DL_MMSE_CF./(1e6);
% R_GUEs_DL_MMSE_CF_Perfect_CS=R_GUEs_DL_MMSE_CF_Perfect_CSI./(1e6);
% R_GUEs_DL_MMSE_CF_WFPC=R_GUEs_DL_MMSE_CF_WFPC./(1e6);
% R_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC=R_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC./(1e6);
% WFPA maximization
% sum_R_CF_Perfect_DL_WFPC=R_GUEs_DL_MMSE_CF_Perfect_CSI_WFPCC./(1e6);
% sum_R_CF_DL_WFPC=R_GUEs_DL_MMSE_CF_WFPCC./(1e6);

sum_R_CF_Perfect_DL_WFPC=sum(R_GUEs_DL_MMSE_CF_Perfect_CSI_WFPCC,1)/N_Channels/(1e6);
sum_R_CF_DL_WFPC=sum(R_GUEs_DL_MMSE_CF_WFPCC,1)/N_Channels/(1e6);

%%%%% Pilot Assignment %%%%%%%%%%
 %%%%% Pilot Assignment %%%%%%%%%%
% %      save('Pilot_Assignment_Orthogonal_random'); %I've finishied this
  %% save('Pilot_Assignment_Orthogonal_80_AP');  %I've finishied this
% %       save('Pilot_Assignment_Orthogonal_32_random'); %I've finishied this
% %     save('Pilot_Assignment_Random');
   save('Pilot_Assignment_random_ContaEf_20'); %
% %       save('Pilot_Assignment_Orthogonal_644'); % this a 64 PA
%     save('Pilot_Assignment_322X_random');
 

toc

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
     'WF Power control, UB', 'WF Power control, PCSI [49]','WF Power control, ICSI',...
    'Location', 'SouthEast');
xlabel('Rate Per-UE [Mbit/s]');
ylabel('CDF');


figure(2)
hold on; box on;
set(gca,'fontsize',10);
grid on
p=plot(x_GUEs_DL_MMSE_CF,y_GUEs_DL_MMSE_CF,'Color','#4DBEEE','LineWidth',2);
p.LineStyle='-';
% p.Marker='o';
% p.MarkerSize = 5;
% p.MarkerIndices = 1:5:length(y_GUEs_DL_MMSE_CF);
 hold on
p=plot(x_GUEs_DL_MMSE_CF_Perfect_CSI,y_GUEs_DL_MMSE_CF_Perfect_CSI,'Color','#4DBEEE','LineWidth',2);
p.LineStyle='-.';
% p.Marker='d';
 hold on
p=plot(x_GUEs_DL_MMSE_mM,y_GUEs_DL_MMSE_mM,'LineWidth',2);
p.Color='#7E2F8E';
p.LineStyle='-';
% p.Marker='o';
% p.MarkerSize = 5;
% p.MarkerIndices = 1:5:length(5);
% hold on
p=plot(x_GUEs_DL_MMSE_mM_Perfect_CSI,y_GUEs_DL_MMSE_mM_Perfect_CSI,'LineWidth',2);
p.Color='#7E2F8E';
p.LineStyle='-.';
grid on
legend('UPA Cell-free mMIMO, ICSI', 'UPA Cell-free mMIMO, PCSI', 'UPA massive MIMO, ICSI','UPA massive MIMO, PCSI',...
    'Location', 'SouthEast');
xlabel('Rate Per-UE (Mbit/s)');
ylabel('CDF');
title('CF mMIMO VS Massive MIMO in UPA')

% xlim([0 60]);


 %%%%%       $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 figure(3)
hold on; box on;
set(gca,'fontsize',11);
% Uniform Power allocatoin in UB
p=plot(x_GUEs_DL_MMSE_CF_UB,y_GUEs_DL_MMSE_CF_UB,'Color','#7E2F8E','LineWidth',2);
p.LineStyle='-';
hold on
% Uniform Power allocatoin in Perfect CSI
p=plot(x_GUEs_DL_MMSE_CF_Perfect_CSI,y_GUEs_DL_MMSE_CF_Perfect_CSI,'Color','#7E2F8E','LineWidth',2);
p.LineStyle='-.';
hold on
% Uniform Power allocatoin in Imperfect CSI
p=plot(x_GUEs_DL_MMSE_CF,y_GUEs_DL_MMSE_CF,'Color','#7E2F8E','LineWidth',2);
p.LineStyle='--';
grid on
legend('Uniform Power control, UB', 'Uniform Power control, PCSI', 'Uniform Power control, ICSI',...
    'Location','SouthEast','TextColor','blue');
xlabel('Rate Per-UE [Mbit/s]');
ylabel('CDF');

figure(4)
hold on; box on;
set(gca,'fontsize',11);
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

legend('WF Power control, UB, $ tau = 8$', 'WF Power control, PCSI, $ tau = 8$','WF Power control, ICSI, $ tau = 8$',...
    'Location','SouthEast','TextColor','blue');

xlabel('Rate Per-UE [Mbit/s]');
ylabel('CDF');
 
figure(5)
plot(Pt_dB_DL,R_GUEs_DL_MMSE_CF,'b','LineWidth',2);
hold on
plot(Pt_dB_DL,R_GUEs_DL_MMSE_CF_Perfect_CSI,'b--','LineWidth',2);
hold on
plot(Pt_dB_DL,R_GUEs_DL_MMSE_CF_WFPC,'r','LineWidth',2);
hold on
plot(Pt_dB_DL,R_GUEs_DL_MMSE_CF_Perfect_CSI_WFPC,'r--','LineWidth',2);
grid on
AX=legend('UPA, ICSI','UPA, PCSI','WFPC, ICSI','WFPC, PCSI','Location', 'SouthEast');
LEG = findobj(AX,'type','text');
set(LEG,'FontSize',12)
xlabel('P_T [dBW]');
ylabel('Sum Rate [Mbit/s]');

 figure(6) 
 hold on; box on;
set(gca,'fontsize',11);
plot(NN,sum_R_CF_Perfect_DL_WFPC,'b','LineWidth',2);
% hold on
plot(NN,sum_R_CF_DL_WFPC,'r','LineWidth',2);
% hold on
grid on

AX=legend('WFPC, PCSI','WFPC, ICSI','Location', 'SouthEast');
LEG = findobj(AX,'type','text');
set(LEG,'FontSize',12)

xlabel('K');
ylabel('Sum Rate [Mbit/s]');
title('SE vs Number of users ($K$)','Interpreter','latex');
xlabel('Number of User per APs, ($K$)','Interpreter','latex');
