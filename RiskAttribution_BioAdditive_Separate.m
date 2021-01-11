%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% STAKEHOLDER'S RISK ATTRIBUTION & UNCERTAINTY ANALYSIS %%%%%%%%%%
%%%%%%%%%%         ALGAE PLASTIC BLEND ADDITIVE PROCESS          %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Jeehwan Lee, KAIST, stevelee@kaist.ac.kr %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clc;
tic;
format longG;
progress = waitbar(0, 'Running...', 'Name', 'Running Simulations...');
total_steps = 1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% [1] Initialize Parameters and Spaces
% (a) Define # of Model Evaluation (# of MC sampling of Parameter Sets)
n_sim = 200000;
% (b) Automatically generate Bayesian target values, or manually enter them
bayes_auto = 1;
bayes_target = questdlg('Automatically generate Target Output Criteria for Naive Bayesian Classification?',...
        'USER INPUT',...
        'No','Yes','Yes');
switch bayes_target
    case 'Yes'
        bayes_auto = 1;
    case 'No'
        bayes_auto = 0;
end
% (c) Initialize Parameter Values for Key Parameters with kp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kp = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];     
alpha = 0.05;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
% (d) Scalar # of Key Parameters. Count columns of kp array.
n_kp = size(kp, 2);
% (d) Define two identical parameter spaces for Sobol Analysis
%     If n_kp params are simulated n_sim times, space is n_sim x n_kp 
ParSpace = [];                 % Space for the Key Parameter of Interest
c_ParSpace = [];               % Space for Complementary Parameters to kp
% (e) Define Resolution of Kernel Distributions
np_1 = 10000;               % Resolution for Parameter KDs for MC sampling
np_2 = 10000;               % Resolution for Bayesian KDs for computing area difference
%%Progress Bar%%
waitbar(20/total_steps, progress, 'Generating Parameter Spaces...');




%% [2] Populate Parameter Space via MC Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Define Parametric Distributions for Key Params]
ParSpace(:,1) = 0.0016.*randn(n_sim, 1) + 0.0224;            % avg_GPC
ParSpace(:,9) = 40076.5.*randn(n_sim, 1) + 222647;           % C_PBR
ParSpace(:,16) = 272.5.*randn(n_sim, 1) + 1090;              % OP_PE_MAINT
ParSpace(:,17) = 116.3.*randn(n_sim, 1) + 775.4;             % OP_PE_PBR
c_ParSpace(:,1) = 0.0016.*randn(n_sim, 1) + 0.0224;          % avg_GPC
c_ParSpace(:,9) = 40076.5.*randn(n_sim, 1) + 222647;         % C_PBR
c_ParSpace(:,16) = 272.5.*randn(n_sim, 1) + 1090;            % OP_PE_MAINT
c_ParSpace(:,17) = 116.3.*randn(n_sim, 1) + 775.4;           % OP_PE_PBR

% [List Kernel Data Points]
data_x2 = [0.89, 1.07, 1.24, 1.00, 0.4, 0.8, 0.84, 0.81, 1.11, 0.82, 1];                           % KD for Mu_max
data_x3 = [3.6, 2.042, 1.4];               % KD for kLa
data_x4 = [0.1358, 0.1184, 0.1283, 0.116]; % KD for P_CO2
data_x5 = [0.402, 0.424, 0.465, 0.515, 0.562, 0.592, 0.599, 0.581, 0.542, 0.493, 0.445, 0.411];    % KD for Daylight %
data_x6 = [442.7, 598.0, 759.0, 935.3, 983.2, 962.1, 782.0, 828.0, 772.4, 663.1, 448.5, 381.4];    % KD for Solar Irradiance, Korea
data_x7 = [0.95, 0.93, 0.99, 0.97];        % KD for Recovery
data_x8 = [0.1, 0.2, 0.2, 0.2, 0.22];      % KD for slurry solids
data_x10 = [4933, 5513, 5803, 6093, 6964]; % KD for C_FG_Blower
data_x11 = [1552, 1735, 1826, 1918, 2192]; % KD for C_Stock_Pump
data_x12 = [17199, 19223, 20235, 21246, 24282];  % KD for C_GM_Pump
data_x13 = [53203, 59462, 62591, 65721, 75109];  % KD for C_SEP
data_x14 = [37911, 46035, 54159, 64991, 81239];        % KD for C_SS_T
data_x15 = [334679, 406396, 478113, 573736, 717170];   % KD for C_GM_T
data_x18 = [110, 140, 105, 180];                 % KD for OP_RM_NSOURCE, USDA NASS Quick Stats
data_x19 = [800, 1100, 1400, 1600, 1750, 1800];  % KD for OP_RM_PSOURCE, USDA NASS Quick Stats
data_x20 = [80, 115, 150];                   % KD for OP_RM_SSOURCE
data_x21 = [160, 180, 200, 250];             % KD for OP_RM_CaSOURCE
data_x22 = [1000, 1500, 2200, 3000];         % KD for OP_RM_NaSOURCE
data_x23 = [0.068];                          % KD for OP_RM_Water. 0 from NREL, 0.068 from IPON
data_x24 = [13.19, 14.5, 24.10, 15.16, 18.94];         % KD for UT_ELEC
data_x25 = [1.18, 1.86328, 2.29408, 1.43248];                      % KD for GWI_NSOURCE, Winnipeg Data
data_x26 = [5.60358, 6.19445, 6.78533, 6.76642, 7.52259, 6.01024]; % KD for GWI_PSOURCE KH2PO4 and K2HPO4 which is produced from KCl and H3PO4 which are mined
data_x27 = [0.3];                                                  % KD for GWI_SSOURCE, MgSO4 Winnipeg Data
data_x28 = [0.89];                                                 % KD for CaCl2-2H2O Winnipeg Data
data_x29 = [0.053];                                                % KD for EDTA-Na2
data_x30 = [0.00003, 0.00001];                                     % KD for Water, from Winnipeg
data_x31 = [108.800, 109.7, 127.000, 82.0, 125.0, 166];            % KD for GWI_ELEC

% [Define Kernel Distributions]
f_x2 = fitdist(data_x2', 'Kernel', 'Support', [0.01, 10]);      % Mu_max
f_x3 = fitdist(data_x3', 'Kernel', 'Support', [0.01, 100]);     % KLa
f_x4 = fitdist(data_x4', 'Kernel', 'Support', [0.01, 1]);       % P_CO2
f_x5 = fitdist(data_x5', 'Kernel', 'Support', [0.25, 1]);       % Daylight%
f_x6 = fitdist(data_x6', 'Kernel', 'Support', 'Positive');      % Solar Irradiance
f_x7 = fitdist(data_x7', 'Kernel', 'Support', [0.5, 1]);        % Recovery%
f_x8 = fitdist(data_x8', 'Kernel', 'Support', [0.05, 1]);       % SlurrySolids%
f_x10 = fitdist(data_x10', 'Kernel', 'Support', 'Positive');    % C_FG_Blower
f_x11 = fitdist(data_x11', 'Kernel', 'Support', 'Positive');    % C_Stock_Pump
f_x12 = fitdist(data_x12', 'Kernel', 'Support', 'Positive');    % C_GM_Pump
f_x13 = fitdist(data_x13', 'Kernel', 'Support', 'Positive');    % C_SEP
f_x14 = fitdist(data_x14', 'Kernel', 'Support', 'Positive');    % C_SS_T
f_x15 = fitdist(data_x15', 'Kernel', 'Support', 'Positive');    % C_GM_T
f_x18 = fitdist(data_x18', 'Kernel', 'Support', 'Positive');    % OP_RM_NSOURCE
f_x19 = fitdist(data_x19', 'Kernel', 'Support', 'Positive');    % OP_RM_PSOURCE
f_x20 = fitdist(data_x20', 'Kernel', 'Support', 'Positive');    % OP_RM_SSOURCE
f_x21 = fitdist(data_x21', 'Kernel', 'Support', 'Positive');    % OP_RM_CaSOURCE
f_x22 = fitdist(data_x22', 'Kernel', 'Support', 'Positive');    % OP_RM_NaSOURCE
f_x23 = fitdist(data_x23', 'Kernel', 'Support', 'Positive');    % OP_RM_WATER
f_x24 = fitdist(data_x24', 'Kernel', 'Support', 'Positive');    % UT_ELEC
f_x25 = fitdist(data_x25', 'Kernel', 'Support', 'Positive');    % GWI_NSOURCE
f_x26 = fitdist(data_x26', 'Kernel', 'Support', 'Positive');    % GWI_PSOURCE
f_x27 = fitdist(data_x27', 'Kernel', 'Support', 'Positive');    % GWI_SSOURCE
f_x28 = fitdist(data_x28', 'Kernel', 'Support', 'Positive');    % GWI_CaSOURCE
f_x29 = fitdist(data_x29', 'Kernel', 'Support', 'Positive');    % GWI_NaSOURCE
f_x30 = fitdist(data_x30', 'Kernel', 'Support', 'Positive');    % GWI_WATER
f_x31 = fitdist(data_x31', 'Kernel', 'Support', 'Positive');    % GWI_ELEC

% [Sample from Kernel Distributions]
% Random Sample to Populate Parameter Space
ParSpace(:,2) = random(f_x2, n_sim, 1);
ParSpace(:,3) = random(f_x3, n_sim, 1);
ParSpace(:,4) = random(f_x4, n_sim, 1);
ParSpace(:,5) = random(f_x5, n_sim, 1);
ParSpace(:,6) = random(f_x6, n_sim, 1);
ParSpace(:,7) = random(f_x7, n_sim, 1);
ParSpace(:,8) = random(f_x8, n_sim, 1);
ParSpace(:,10) = random(f_x10, n_sim, 1);
ParSpace(:,11) = random(f_x11, n_sim, 1);
ParSpace(:,12) = random(f_x12, n_sim, 1);
ParSpace(:,13) = random(f_x13, n_sim, 1);
ParSpace(:,14) = random(f_x14, n_sim, 1);
ParSpace(:,15) = random(f_x15, n_sim, 1);
ParSpace(:,18) = random(f_x18, n_sim, 1);
ParSpace(:,19) = random(f_x19, n_sim, 1);
ParSpace(:,20) = random(f_x20, n_sim, 1);
ParSpace(:,21) = random(f_x21, n_sim, 1);
ParSpace(:,22) = random(f_x22, n_sim, 1);
ParSpace(:,23) = random(f_x23, n_sim, 1);
ParSpace(:,24) = random(f_x24, n_sim, 1);
ParSpace(:,25) = random(f_x25, n_sim, 1);
ParSpace(:,26) = random(f_x26, n_sim, 1);
ParSpace(:,27) = random(f_x27, n_sim, 1);
ParSpace(:,28) = random(f_x28, n_sim, 1);
ParSpace(:,29) = random(f_x29, n_sim, 1);
ParSpace(:,30) = random(f_x30, n_sim, 1);
ParSpace(:,31) = random(f_x31, n_sim, 1);
% Random Sample to Populate Complementary Parameter Space
c_ParSpace(:,2) = random(f_x2, n_sim, 1);
c_ParSpace(:,3) = random(f_x3, n_sim, 1);
c_ParSpace(:,4) = random(f_x4, n_sim, 1);
c_ParSpace(:,5) = random(f_x5, n_sim, 1);
c_ParSpace(:,6) = random(f_x6, n_sim, 1);
c_ParSpace(:,7) = random(f_x7, n_sim, 1);
c_ParSpace(:,8) = random(f_x8, n_sim, 1);
c_ParSpace(:,10) = random(f_x10, n_sim, 1);
c_ParSpace(:,11) = random(f_x11, n_sim, 1);
c_ParSpace(:,12) = random(f_x12, n_sim, 1);
c_ParSpace(:,13) = random(f_x13, n_sim, 1);
c_ParSpace(:,14) = random(f_x14, n_sim, 1);
c_ParSpace(:,15) = random(f_x15, n_sim, 1);
c_ParSpace(:,18) = random(f_x18, n_sim, 1);
c_ParSpace(:,19) = random(f_x19, n_sim, 1);
c_ParSpace(:,20) = random(f_x20, n_sim, 1);
c_ParSpace(:,21) = random(f_x21, n_sim, 1);
c_ParSpace(:,22) = random(f_x22, n_sim, 1);
c_ParSpace(:,23) = random(f_x23, n_sim, 1);
c_ParSpace(:,24) = random(f_x24, n_sim, 1);
c_ParSpace(:,25) = random(f_x25, n_sim, 1);
c_ParSpace(:,26) = random(f_x26, n_sim, 1);
c_ParSpace(:,27) = random(f_x27, n_sim, 1);
c_ParSpace(:,28) = random(f_x28, n_sim, 1);
c_ParSpace(:,29) = random(f_x29, n_sim, 1);
c_ParSpace(:,30) = random(f_x30, n_sim, 1);
c_ParSpace(:,31) = random(f_x31, n_sim, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Progress Bar%%
waitbar(100/total_steps, progress, 'Initializing Monte Carlo Simulations');
%% Prior Distribution Plots
% Average Gram per Cell
figure(1)
X_1 = linspace(min(ParSpace(:,1)),max(ParSpace(:,1)),n_sim);
Y_1 = normpdf(X_1, 0.0224, 0.0016);
plot(X_1,Y_1,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P1 (Avg. Weight per 10? Cells, g/10?cells)')
saveas(gcf, 'PRIOR_BioAdditive_param1.png')
% Maximum Specific Growth Rate
figure(2)
[Y_2, X_2] = ksdensity(data_x2, 'npoints', np_2);
plot(X_2,Y_2,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P2 (Maximum specific growth rate, 1/hr)')
saveas(gcf, 'PRIOR_BioAdditive_param2.png')
% Volumetric Mass Transfer Coefficient
figure(3)
[Y_3, X_3] = ksdensity(data_x3, 'npoints', np_2);
plot(X_3,Y_3,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P3 (Volumetric Mass Transfer Coefficient, 1/hr)')
saveas(gcf, 'PRIOR_BioAdditive_param3.png')
% CO2 Partial Pressure
figure(4)
[Y_4, X_4] = ksdensity(data_x4, 'npoints', np_2);
plot(X_4,Y_4,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P4 (Partial Pressure of CO2, atm)')
saveas(gcf, 'PRIOR_BioAdditive_param4.png')
% Photosynthesis Daylight Fraction
figure(5)
[Y_5, X_5] = ksdensity(data_x5, 'npoints', np_2);
plot(X_5,Y_5,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P5 (Daily Photosynthesis Fraction, %)')
saveas(gcf, 'PRIOR_BioAdditive_param5.png')
% Photosynthetic Solar Irradiance
figure(6)
[Y_6, X_6] = ksdensity(data_x6, 'npoints', np_2);
plot(X_6,Y_6,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P6 (PAR, ?E/m2s)')
saveas(gcf, 'PRIOR_BioAdditive_param6.png')
% Biomass Recovery (Harvest+Dewater)
figure(7)
[Y_7, X_7] = ksdensity(data_x7, 'npoints', np_2);
plot(X_7,Y_7,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P7 (Harvest/Dewater Recovery, %)')
saveas(gcf, 'PRIOR_BioAdditive_param7.png')
% Biomass Concentration in Slurry
figure(8)
[Y_8, X_8] = ksdensity(data_x8, 'npoints', np_2);
plot(X_8,Y_8,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P8 (Biomass Concentration in Slurry, %)')
saveas(gcf, 'PRIOR_BioAdditive_param8.png')
% Photobioreactor Cost
figure(9)
X_9 = linspace(min(ParSpace(:,9)),max(ParSpace(:,9)),n_sim);
Y_9 = normpdf(X_9, 222647, 40076.5);
plot(X_9,Y_9,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P9 (Cultivation Cost, USD/Acre)')
saveas(gcf, 'PRIOR_BioAdditive_param9.png')
% FG Blower Cost
figure(10)
[Y_10, X_10] = ksdensity(data_x10, 'npoints', np_2);
plot(X_10,Y_10,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P10 (FG Blower Cost, USD/unit)')
saveas(gcf, 'PRIOR_BioAdditive_param10.png')
% Stock Pump Cost
figure(11)
[Y_11, X_11] = ksdensity(data_x11, 'npoints', np_2);
plot(X_11,Y_11,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P11 (Stock Pump Cost, USD/unit)')
saveas(gcf, 'PRIOR_BioAdditive_param11.png')
% GM Pump Cost
figure(12)
[Y_12, X_12] = ksdensity(data_x12, 'npoints', np_2);
plot(X_12,Y_12,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P12 (GM Pump Cost, USD/unit)')
saveas(gcf, 'PRIOR_BioAdditive_param12.png')
% Centrifugal Separator Cost
figure(13)
[Y_13, X_13] = ksdensity(data_x13, 'npoints', np_2);
plot(X_13,Y_13,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P13 (Centrifugal Separator Cost, USD/unit)')
saveas(gcf, 'PRIOR_BioAdditive_param13.png')
% Stock Solution Tank (w/Agitator) Cost
figure(14)
[Y_14, X_14] = ksdensity(data_x15, 'npoints', np_2);
plot(X_14,Y_14,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P14 (Stock Solution Tank Cost, USD/unit)')
saveas(gcf, 'PRIOR_BioAdditive_param14.png')
% Growth Media Tank (Storage) Cost
figure(15)
[Y_15, X_15] = ksdensity(data_x15, 'npoints', np_2);
plot(X_15,Y_15,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P15 (Growth Media Tank Cost, USD/unit)')
saveas(gcf, 'PRIOR_BioAdditive_param15.png')
% Plant: Maintenance and CIP Cost
figure(16)
X_16 = linspace(min(ParSpace(:,16)),max(ParSpace(:,16)),n_sim);
Y_16 = normpdf(X_16, 1090, 272.5);
plot(X_16,Y_16,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P16 (PBR Maintenance Cost, USD/acre)')
saveas(gcf, 'PRIOR_BioAdditive_param16.png')
% PBR LDPE Cost
figure(17)
X_17 = linspace(min(ParSpace(:,17)),max(ParSpace(:,17)),n_sim);
Y_17 = normpdf(X_17, 775.4, 116.3);
plot(X_17,Y_17,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P17 (PBR LDPE Cost, USD/acre)')
saveas(gcf, 'PRIOR_BioAdditive_param17.png')
% NUTRIENT: N-Nutrient Cost
figure(18)
[Y_18, X_18] = ksdensity(data_x18, 'npoints', np_2);
plot(X_18,Y_18,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P18 (N Nutrient Cost, USD/ton)')
saveas(gcf, 'PRIOR_BioAdditive_param18.png')
% NUTRIENT: P-Nutrient Cost
figure(19)
[Y_19, X_19] = ksdensity(data_x19, 'npoints', np_2);
plot(X_19,Y_19,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P19 (P Nutrient Cost, USD/ton)')
saveas(gcf, 'PRIOR_BioAdditive_param19.png')
% NUTRIENT: S-Nutrient Cost
figure(20)
[Y_20, X_20] = ksdensity(data_x20, 'npoints', np_2);
plot(X_20,Y_20,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P20 (S Nutrient Cost, USD/ton)')
saveas(gcf, 'PRIOR_BioAdditive_param20.png')
% NUTRIENT: Ca-Nutrient Cost
figure(21)
[Y_21, X_21] = ksdensity(data_x21, 'npoints', np_2);
plot(X_21,Y_21,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P21 (Ca Nutrient Cost, USD/ton)')
saveas(gcf, 'PRIOR_BioAdditive_param21.png')
% NUTRIENT: Na-Nutrient Cost
figure(22)
[Y_22, X_22] = ksdensity(data_x22, 'npoints', np_2);
plot(X_22,Y_22,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P22 (P Nutrient Cost, USD/ton)')
saveas(gcf, 'PRIOR_BioAdditive_param22.png')
% NUTRIENT: Water Cost
figure(23)
[Y_23, X_23] = ksdensity(data_x23, 'npoints', np_2);
plot(X_23,Y_23,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P23 (Water Cost, USD/ton)')
saveas(gcf, 'PRIOR_BioAdditive_param23.png')
% Korean Grid Mix Electricity Cost
figure(24)
[Y_24, X_24] = ksdensity(data_x24, 'npoints', np_2);
plot(X_24,Y_24,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P24 (Korean Grid Mix Electricity, USD/GJ)')
saveas(gcf, 'PRIOR_BioAdditive_param24.png')
% GWI: N-Nutrient
figure(25)
[Y_25, X_25] = ksdensity(data_x25, 'npoints', np_2);
plot(X_25,Y_25,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P25 (GWI of N Nutrients, kgCO2eq/ton)')
saveas(gcf, 'PRIOR_BioAdditive_param25.png')
% GWI: P-Nutrient
figure(26)
[Y_26, X_26] = ksdensity(data_x26, 'npoints', np_2);
plot(X_26,Y_26,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P26 (GWI of P Nutrients, kgCO2eq/ton)')
saveas(gcf, 'PRIOR_BioAdditive_param26.png')
% GWI: S-Nutrient
figure(27)
[Y_27, X_27] = ksdensity(data_x27, 'npoints', np_2);
plot(X_27,Y_27,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P27 (GWI of S Nutrients, kgCO2eq/ton)')
saveas(gcf, 'PRIOR_BioAdditive_param27.png')
% GWI: Ca-Nutrient
figure(28)
[Y_28, X_28] = ksdensity(data_x28, 'npoints', np_2);
plot(X_28,Y_28,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P28 (GWI of Ca Nutrients, kgCO2eq/ton)')
saveas(gcf, 'PRIOR_BioAdditive_param28.png')
% GWI: Na-Nutrient
figure(29)
[Y_29, X_29] = ksdensity(data_x29, 'npoints', np_2);
plot(X_29,Y_29,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P29 (GWI of Na Nutrients, kgCO2eq/ton)')
saveas(gcf, 'PRIOR_BioAdditive_param29.png')
% GWI: Local Water
figure(30)
[Y_30, X_30] = ksdensity(data_x30, 'npoints', np_2);
plot(X_30,Y_30,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P30 (GWI of Local Water Aquifer, kgCO2eq/ton)')
saveas(gcf, 'PRIOR_BioAdditive_param30.png')
% GWI of Grid Mix Electricity
figure(31)
[Y_31, X_31] = ksdensity(data_x31, 'npoints', np_2);
plot(X_31,Y_31,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P31 (GWI of Korean Grid Mix, kgCO2eq/GJ)')
saveas(gcf, 'PRIOR_BioAdditive_param31.png')




%% [3] Process Parameter Space for System Evaluation
% (a) Determine the # of System Model Outputs (# of Evaluation Metrics)
n_out = size(MODEL_BioAdditive(kp),2);
% (b) Create containers for Evaluated Model Outputs
% NOTE: eval_P and eval_C are n_sim x n_out x n_kp matrix. 
fx = zeros(n_sim, n_out);  % Evaluated outputs with all Params from ParSpace
fx_P = zeros(n_sim, n_out);% Evaluated outputs with i from ParSpace, ~i from c_ParSpace
fx_C = zeros(n_sim, n_out);% Evaluated outputs with i from c_ParSpacek, ~i from ParSpace
% (c) Evaluate Model from Monte Carlo Sampled Inputs (ParSpace)
parfor i = 1:n_sim
    % Each Parameter Set is a Row in the ParSpace matrix
    Parameter_Set = ParSpace(i,:); 
    fx(i,:) = MODEL_BioAdditive(Parameter_Set);
end
% (d) Generate Function Output Space based on i and ~i
for i = 1:n_sim
    for j = 1:n_kp
        %%%% fx_P = f(x_ik, x'_~ik)
        kp = [c_ParSpace(i,1:j-1), ParSpace(i,j), c_ParSpace(i,j+1:n_kp)];
        fx_P(i,:,j) = MODEL_BioAdditive(kp);
        %%%% fx_C = f(x'_ik, x_~ik)
        kp = [ParSpace(i,1:j-1), c_ParSpace(i,j), ParSpace(i,j+1:n_kp)];
        fx_C(i,:,j) = MODEL_BioAdditive(kp);
    end
    %%Progress Bar%%
    waitbar((100 + ((i/n_sim)*520))/total_steps, progress,...
    sprintf('Processing %d of %d Simulations', i, n_sim));
end
%%Progress Bar%%
waitbar(580/total_steps, progress, 'Computing Sobol Variances');




%% [4-SOBOL] Computing Variances
% (a) Initialize Integrals and Total Variances
% Initialize the Integral of Model Outputs (f0^2) 
f0 = zeros(1, n_out);            
% Initialize Total Variance of Model Outputs 
D = zeros(1, n_out);            
% (b) Compute the Average of Model Outputs, f0
parfor i = 1:n_out 
    for j = 1:n_sim
        f0(i) = f0(i) + fx(j,i);
    end
    f0(i) = f0(i)/n_sim;
end
% (c) Estimating Total Output Variance, D, using MC Model Outputs
parfor i = 1:n_out
    for j = 1:n_sim
        D(i) = D(i) + ((fx(j,i)^2)-(f0(i)^2));
    end
    D(i) = D(i)/n_sim;
end
% (d) Compute Partial Variances (1st Order Sobol)
%D_1st = zeros(n_kp, n_out);
%F_Par = zeros(n_kp, n_out);           % Initialize Partial Variance Factors
%parfor i = 1:n_kp
%    for j = 1:n_out
%        for k = 1:n_sim
%            F_Par(i,j) = F_Par(i,j) + (fx(k,j) - fx_P(k,j,i))^2;
%        end
%        D_1st(i,j) = D(j) - (F_Par(i,j)/(n_sim*2));
%    end
%end
%(e) Compute Total Variances (Total Sobol)
D_Tot = zeros(n_kp, n_out);
F_cPar = zeros(n_kp, n_out);          % Initialize Total Sobol Factors
parfor i = 1:n_kp
    for j = 1:n_out
        for k = 1:n_sim
            F_cPar(i,j) = F_cPar(i,j) + (fx(k,j) - fx_C(k,j,i))^2;
        end
        D_Tot(i,j) = F_cPar(i,j)/(n_sim*2);
    end
end
%%Progress Bar%%
waitbar(630/total_steps, progress, 'Calculating Sobol Indices for each System Output');




%% [4-SOBOL] Determine Sobol Indices
% (a) Compute Sum of Partial and Total Variances
%Sum_Partial = zeros(n_out);
Sum_Total = zeros(n_out);
for i = 1:n_out
    %Sum_Partial(i) = sum(D_1st(:,i));
    Sum_Total(i) = sum(D_Tot(:,i));
end
% (b) Compute Rank Scores based on Partial/Total Variances
%Sobol_1st = zeros(n_kp, n_out);
Sobol_Total = zeros(n_kp, n_out);
for i = 1:n_out
    for j = 1:n_kp
        %Sobol_1st(j,i) = D_1st(j,i)/D(i);
        Sobol_Total(j,i) = D_Tot(j,i)/D(i);
    end
end
% (c) Rank Input Parameters via Sobol Indices (Total)
for i = 1:n_out
    [score2, rank2] = sort(Sobol_Total(:,i), 'descend');
    fprintf('Total Sobol Ranks for Model Evaluation Output %.0f \n', i)
    sprintf('The Total Sobol Indices for the Above Order are: ')
    Sobol_Total(:,i)
    fprintf('Sum of Total Sobol Indices for Model Evaluation Output %.0f \n', i)
    sum(Sobol_Total(:,i))
end
waitbar(640/total_steps, progress, 'Calculating Normalized Total Sobol Indices');




%% [4-SOBOL] Normalized Total Sobol Indices
% (a) Define a vector of weights for each model output
%for i = 1:n_out
%    weight_prompt = sprintf('Enter Stakeholder Weight for Sustainability Criteria %.0f \n', i);
%    weight_prompt_title = 'Enter Stakeholder Weight';
%    dims1 = [1 35];
%    weight_answer = inputdlg(weight_prompt, weight_prompt_title, dims1);
%    weight(i) = str2double(weight_answer{1});
%end
% (b) Check if the sum of the weight vector = 1
%if sum(weight)~= 1
%    error('The sum of the metric weights must equal to 1')
%end
%if size(weight,2)~= n_out
%    error('A numerical weight must be assigned to each system evaluation output. For no weights, enter 0')
%end
% (c) Compute Normalized Total Sobol Indices
%NTS = zeros(1, n_kp);
%for i = 1:n_kp
%    NTS(i) = sum(weight.*Sobol_Total(i,:));
%end
%NTS_Sum = sum(NTS);
%waitbar(650/total_steps, progress, 'Classifying System Model Outputs');




%% [5-GNB Class] Gaussian Naive Bayes Classification of System Model Outputs
% (a) Define Target Values for each System Model Outputs
targets = zeros(1, n_out);  % Array of Target Bayesian Output Metric Values
if bayes_auto == 1
    for i = 1:n_out
        targets(i) = mean(fx(:,i));
    end
else
    for i = 1:n_out
        prompt = sprintf('Enter Stakeholder Decision Criteria for Output Metric %.0f \n', i);
        prompt_title = 'Enter Stakeholder Decision Criteria Value';
        dims = [1 35];
        answer = inputdlg(prompt, prompt_title, dims);
        targets(i) = str2double(answer{1});
    end
end              
% (b) Generate Bayesian Classification Matrix
classmat = zeros(n_sim, n_out);
parfor i = 1:n_sim
    for j = 1:n_out
        if fx(i,j) <= targets(j)
            classmat(i,j) = 1;
        else
            classmat(i,j) = 0;
        end
    end
end
%%Progress Bar%%
waitbar(700/total_steps, progress, 'Classifying Parameter Inputs based on Criteria');




%% [5-GNB Class] Factorize Classification Matrix and Append Values
% (a) Generate Empty Factorized Containers. The maximum dimensions of each
%     Label matrix is n_sim x n_kp x n_out
Success_Mat = zeros(n_sim, n_kp, n_out);
Failure_Mat = zeros(n_sim, n_kp, n_out);
% (b) Loop through "classmat" then append values from ParamSpace to 
%     Success_Mat or Failure_Mat accordingly
parfor i = 1:n_sim
    for j = 1:n_out
        if classmat(i,j) == 1
            Success_Mat(i,:,j) = ParSpace(i,:);
            Failure_Mat(i,:,j) = zeros(1,n_kp,'uint32');
        else
            Success_Mat(i,:,j) = zeros(1,n_kp,'uint32');
            Failure_Mat(i,:,j) = ParSpace(i,:);
        end
    end
end
%%Progress Bar%%
waitbar(750/total_steps, progress, 'Generating Bayesian Kernel Distributions');




%% [6-GNB Risk] Generate Kernel Distributions for GNB Sort Matrix
% (a) Initialize Containers
S_KernelData = zeros(n_sim, n_kp, n_out); % Sampled matrix for "success" 
F_KernelData = zeros(n_sim, n_kp, n_out); % Sampled matrix for "failure" 
Stake_Risk_Score = zeros(n_kp, n_out);    % A n_kp x n_out matrix of ranks
Posterior_Diff = zeros(1,np_2);           
X_LIST = zeros(1,np_2);                   % Normalized X-Axis array
% (b) Populate Stakeholders Risk Score by evaluating the differences
for i = 1:n_out
    for j = 1:n_kp
        % Generate array of sample data for Bayesian Kernel Distributions
        S_KernelData = nonzeros(Success_Mat(:,j,i))' ;   
        F_KernelData = nonzeros(Failure_Mat(:,j,i))';
        % Generate the Kernel Distribution for each Bayesian Outcome
        [F_s, x_s] = ksdensity(S_KernelData, 'npoints', np_2);
        [F_f, x_f] = ksdensity(F_KernelData, 'npoints', np_2);
        % Because Success/Failure are 2 data sets, need to interpolate on a
        % consistent X-Axis. Generate the consistent X-Axis array
        X_APPEND = [x_s, x_f];
        X_MIN = min(X_APPEND);
        X_MAX = max(X_APPEND);
        X_LIST = linspace(X_MIN,X_MAX,np_2);
        % Interpolate the Success/Failure probabilities with above array
        Y_s = interp1(x_s, F_s, X_LIST);
        Y_f = interp1(x_f, F_f, X_LIST);
        % Calculate the Posterior Probability differences as Risk Scores
        Posterior_Diff = abs(Y_s-Y_f);
        Posterior_Diff(isnan(Posterior_Diff))=0;
        Stake_Risk_Score(j,i)  = trapz(X_LIST, Posterior_Diff);
        % Plot the Posterior Probability Differences
        plot(X_LIST,Y_s,'b','LineWidth',2);
        hold on
        plot(X_LIST,Y_f,'r','LineWidth',2);
        plot(X_LIST, abs(Y_s-Y_f), '--k', 'LineWidth', 1);
        if i == 1
            title(['For Parameter ', num2str(j), ' for criteria COGM'])
        elseif i == 2
            title(['For Parameter ', num2str(j), ' for criteria CO2 GWI'])
        end
        legend('Success Posterior', 'Failure Posterior', 'Posterior Prob. Diff.')
        hold off
        saveas(gcf, ['POSTERIOR_BioAdditive_output_',num2str(i),'_parameter_',num2str(j),'.png'])
    end
end
% (c) Generate a Normalized Rank Score for Posterior Probability Differences
%Norm_PosProb_Score = zeros(1, n_kp);
%for i = 1:n_kp
%    Norm_PosProb_Score(i) = sum(weight.*Stake_Risk_Score(i,:));
%end
%Sum_NPPS = sum(Norm_PosProb_Score);
%%Progress Bar%%
waitbar(970/total_steps, progress, 'Sorting Parameters by Stakeholders Risk Score');




%% [7] Rank Parameters 
% (a) Normalized Total Sobol Indices for Current Model Input Parameters
%sprintf('Normalized Total Sobol Indices for the Current Parameter Set')
%NTS'
%fprintf('Sum of the Normalized Total Sobol Indices are:')
%NTS_Sum
% (b) Rank Parameters by Normalized Posterior Probability Difference Score
%sprintf('Normalized Posterior Probability Difference Scores for the Current Parameter Set')
%Norm_PosProb_Score'
%sprintf('Sum of the Normalized Posterior Probability Difference Scores are:')
%Sum_NPPS
% (c) Display Posterior Probability Difference Scores for Each System Output
parfor i = 1:n_out
    fprintf('Stakeholders Risk Score for Sustainability Criteria %.0f \n',i)
    Stake_Risk_Score(:,i)
end
waitbar(990/total_steps, progress, 'Plotting Distributions of System Outputs');




%% [8] Plot and Export Results
% (a) Plot Distributions of Model Outputs
%nbins = 50;             % Define resolution of output metric's historgram
% (b) Generate plots (Example for via For loop below)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %figure('Name', 'Variance of Cost of Products Manufactured')
    %histogram(fx(:,1), nbins, 'facecolor', [0, 0, 0])
    %xlabel('Net Unit Production Cost, $/kg')
    %ylabel('Frequency')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Save Results
% Define Vars
Sort_COGM = sort(fx(:,1));
Sort_GWI = sort(fx(:,2));
INTQUART = zeros(2,n_out);
TAILS = zeros(2,n_out);
VAR = zeros(1,n_out);
CVAR = zeros(1,n_out);
% Inter-quartile Range
INTQUART(1,1) = Sort_COGM(n_sim*0.25);
INTQUART(2,1) = Sort_GWI(n_sim*0.25);
INTQUART(1,2) = Sort_COGM(n_sim*0.75);
INTQUART(2,2) = Sort_GWI(n_sim*0.75);
% 5% - 95% range
TAILS(1,1) = Sort_COGM(n_sim*0.05);
TAILS(2,1) = Sort_GWI(n_sim*0.05);
TAILS(1,2) = Sort_COGM(n_sim*0.95);
TAILS(2,2) = Sort_GWI(n_sim*0.95);
% Value-at-Risk for alpha=0.05
VAR(1) = Sort_COGM((1-alpha)*n_sim);
VAR(2) = Sort_GWI((1-alpha)*n_sim);
% Conditional Value-at-Risk for alpha = 0.05
CVAR(1) = mean(Sort_COGM((1-alpha)*n_sim:n_sim));
CVAR(2) = mean(Sort_GWI((1-alpha)*n_sim:n_sim));
%%Progress Bar%%
waitbar(1000/total_steps, progress, 'Complete');
delete(progress)
toc




%% EXPORT RESULTS
RESULTS = struct();
RESULTS.SobolTot = Sobol_Total;
RESULTS.RiskScore = Stake_Risk_Score;
RESULTS.Output = fx;
RESULTS.IQ = INTQUART;
RESULTS.TAILS = TAILS;
RESULTS.VAR = VAR;
RESULTS.CVAR = CVAR;
% Export to CSV
BioAdditive_Res = struct2cell(RESULTS);
writecell( BioAdditive_Res, 'BioAdditive_Results_Separate.csv');
