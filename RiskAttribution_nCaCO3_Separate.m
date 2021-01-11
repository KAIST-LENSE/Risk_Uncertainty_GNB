%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% STAKEHOLDER'S RISK ATTRIBUTION & UNCERTAINTY ANALYSIS %%%%%%%%%%
%%%%%%%%%%%%%         NANO-CALCIUM CARBONATE PROCESS          %%%%%%%%%%%%%
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
kp = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];    
alpha = 0.05;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
% (d) Scalar # of Key Parameters. Count columns of kp array.
n_kp = size(kp, 2);
% (e) Define two identical parameter spaces for Sobol Analysis
%     If n_kp params are simulated n_sim times, space is n_sim x n_kp 
ParSpace = [];                 % Space for the Key Parameter of Interest
c_ParSpace = [];               % Space for Complementary Parameters to kp
% (f) Define Resolution of Kernel Distributions
np_1 = 10000;               % Resolution for Parameter KDs for MC sampling
np_2 = 10000;               % Resolution for Bayesian KDs for computing area difference
%%Progress Bar%%
waitbar(20/total_steps, progress, 'Generating Parameter Spaces...');




%% [2] Populate Parameter Space via MC Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Define Parametric Distributions for Key Params]
% [List Kernel Data Points]
data_x1 = [0.90, 0.92, 0.88, 0.87, 0.85, 0.90, 0.91, 0.90]; % Faradaic Efficiency
data_x2 = [1224, 183.6];                                    % Electrolysis Energy Requirement kWh/ton nCaCO3
data_x3 = [12700, 10936, 11025, 10780];                     % Equipment: TP
data_x4 = [1000, 980, 1500, 600];                           % Equipment: Electrolysis Cell
data_x5 = [0.25, 0.18, 0.2, 0.27, 0.28];                    % Equipment: Fraction of Electrolyzer Cost
data_x6 = [163200, 140108, 117133, 111831, 215119];         % Equipment: Heat Exchanger, CS&SUS304
data_x7 = [4800, 7863, 8300, 19524];                        % Equipment: Heater, CS
data_x8 = [75000, 139767, 65629];                           % Equipment: Reactor
data_x9 = [29800, 29000, 23221];                            % Equipment: Filter
data_x10 = [5400, 5944, 7395];                              % Equipment: Blower
data_x11 = [100000, 154340, 156570];                        % Equipment: Drum
data_x12 = [118400, 78095, 100601];                         % Equipment: Centrifuge
data_x13 = [74800, 71268, 96012, 62713];                    % Equipment: Dryer
data_x14 = [40000, 50942, 35687];                           % Equipment: Crusher
data_x15 = [88300, 81905, 78836, 80905];                    % Equipment: Heat Exchanger, CS&CS
data_x16 = [11200, 10675, 8500];                            % Equipment: Condenser
data_x17 = [13.19, 14.5, 24.10, 15.16, 18.94, 19.65];       % Electricity Cost, $/GJ
data_x18 = [6.08, 6.35, 5.194, 5.303, 5.074, 3.39];         % Steam Cost, $/GJ
data_x19 = [108.800, 109.7, 127.000, 82.0, 125.0, 166];     % Electricity GWI, kgCO2eq/GJ
data_x20 = [203.61, 270.27, 384.9, 436.09];                 % Steam GWI, kgCO2eq/GJ

% [Define Kernel Distributions]
f_x1 = fitdist(data_x1', 'Kernel', 'Support', [0.75, 0.95]);
f_x2 = fitdist(data_x2', 'Kernel', 'Support', [0.01, 2000]);
f_x3 = fitdist(data_x3', 'Kernel', 'Support', 'Positive');
f_x4 = fitdist(data_x4', 'Kernel', 'Support', [500, 3000]);
f_x5 = fitdist(data_x5', 'Kernel', 'Support', [0.01, 1]);
f_x6 = fitdist(data_x6', 'Kernel', 'Support', 'Positive');
f_x7 = fitdist(data_x7', 'Kernel', 'Support', 'Positive');
f_x8 = fitdist(data_x8', 'Kernel', 'Support', 'Positive');
f_x9 = fitdist(data_x9', 'Kernel', 'Support', 'Positive');
f_x10 = fitdist(data_x10', 'Kernel', 'Support', 'Positive');
f_x11 = fitdist(data_x11', 'Kernel', 'Support', 'Positive');
f_x12 = fitdist(data_x12', 'Kernel', 'Support', 'Positive');
f_x13 = fitdist(data_x13', 'Kernel', 'Support', 'Positive');
f_x14 = fitdist(data_x14', 'Kernel', 'Support', 'Positive');
f_x15 = fitdist(data_x15', 'Kernel', 'Support', 'Positive');
f_x16 = fitdist(data_x16', 'Kernel', 'Support', 'Positive');
f_x17 = fitdist(data_x17', 'Kernel', 'Support', [0, 100]);
f_x18 = fitdist(data_x18', 'Kernel', 'Support', [0, 100]);
f_x19 = fitdist(data_x19', 'Kernel', 'Support', [0, 1000]);
f_x20 = fitdist(data_x20', 'Kernel', 'Support', [0, 1000]);

% [Sample from Kernel Distributions]
% Random Sample to Populate Parameter Space
ParSpace(:,1) = random(f_x1, n_sim, 1);
ParSpace(:,2) = random(f_x2, n_sim, 1);
ParSpace(:,3) = random(f_x3, n_sim, 1);
ParSpace(:,4) = random(f_x4, n_sim, 1);
ParSpace(:,5) = random(f_x5, n_sim, 1);
ParSpace(:,6) = random(f_x6, n_sim, 1);
ParSpace(:,7) = random(f_x7, n_sim, 1);
ParSpace(:,8) = random(f_x8, n_sim, 1);
ParSpace(:,9) = random(f_x9, n_sim, 1);
ParSpace(:,10) = random(f_x10, n_sim, 1);
ParSpace(:,11) = random(f_x11, n_sim, 1);
ParSpace(:,12) = random(f_x12, n_sim, 1);
ParSpace(:,13) = random(f_x13, n_sim, 1);
ParSpace(:,14) = random(f_x14, n_sim, 1);
ParSpace(:,15) = random(f_x15, n_sim, 1);
ParSpace(:,16) = random(f_x16, n_sim, 1);
ParSpace(:,17) = random(f_x17, n_sim, 1);
ParSpace(:,18) = random(f_x18, n_sim, 1);
ParSpace(:,19) = random(f_x19, n_sim, 1);
ParSpace(:,20) = random(f_x20, n_sim, 1);
% Random Sample to Populate Complementary Parameter Space
c_ParSpace(:,1) = random(f_x1, n_sim, 1);
c_ParSpace(:,2) = random(f_x2, n_sim, 1);
c_ParSpace(:,3) = random(f_x3, n_sim, 1);
c_ParSpace(:,4) = random(f_x4, n_sim, 1);
c_ParSpace(:,5) = random(f_x5, n_sim, 1);
c_ParSpace(:,6) = random(f_x6, n_sim, 1);
c_ParSpace(:,7) = random(f_x7, n_sim, 1);
c_ParSpace(:,8) = random(f_x8, n_sim, 1);
c_ParSpace(:,9) = random(f_x9, n_sim, 1);
c_ParSpace(:,10) = random(f_x10, n_sim, 1);
c_ParSpace(:,11) = random(f_x11, n_sim, 1);
c_ParSpace(:,12) = random(f_x12, n_sim, 1);
c_ParSpace(:,13) = random(f_x13, n_sim, 1);
c_ParSpace(:,14) = random(f_x14, n_sim, 1);
c_ParSpace(:,15) = random(f_x15, n_sim, 1);
c_ParSpace(:,16) = random(f_x16, n_sim, 1);
c_ParSpace(:,17) = random(f_x17, n_sim, 1);
c_ParSpace(:,18) = random(f_x18, n_sim, 1);
c_ParSpace(:,19) = random(f_x19, n_sim, 1);
c_ParSpace(:,20) = random(f_x20, n_sim, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prior Distribution Plots
% Faradaic Efficiency
figure(1)
[Y_1, X_1] = ksdensity(data_x1, 'npoints', np_2);
plot(X_1,Y_1,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P1 (Faradaic Efficiency, %)')
saveas(gcf, 'PRIOR_nCaCO3_param1.png')
% Electrolyzer Energy Requirement
figure(2)
[Y_2, X_2] = ksdensity(data_x2, 'npoints', np_2);
plot(X_2,Y_2,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P2(Electrolyzer Energy Requirement, kWh/ton)')
saveas(gcf, 'PRIOR_nCaCO3_param2.png')
% EQ: TP
figure(3)
[Y_3, X_3] = ksdensity(data_x3, 'npoints', np_2);
plot(X_3,Y_3,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P3(Equipment- TP, USD)')
saveas(gcf, 'PRIOR_nCaCO3_param3.png')
% EQ: Electrolyzer Cell
figure(4)
[Y_4, X_4] = ksdensity(data_x4, 'npoints', np_2);
plot(X_4,Y_4,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P4(Equipment- Electrolyzer Cell, USD/kW)')
saveas(gcf, 'PRIOR_nCaCO3_param4.png')
% Stack Cost Fraction
figure(5)
[Y_5, X_5] = ksdensity(data_x5, 'npoints', np_2);
plot(X_5,Y_5,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P5(Fraction of Stack in Electrolyzer Cost)')
saveas(gcf, 'PRIOR_nCaCO3_param5.png')
% EQ: HEATEX1
figure(6)
[Y_6, X_6] = ksdensity(data_x6, 'npoints', np_2);
plot(X_6,Y_6,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P6(Equipment- Heat Exchanger (CS%SUS304), USD)')
saveas(gcf, 'PRIOR_nCaCO3_param6.png')
% EQ: Heater
figure(7)
[Y_7, X_7] = ksdensity(data_x7, 'npoints', np_2);
plot(X_7,Y_7,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P7(Equipment- Heater, USD)')
saveas(gcf, 'PRIOR_nCaCO3_param7.png')
% EQ: Reactor
figure(8)
[Y_8, X_8] = ksdensity(data_x8, 'npoints', np_2);
plot(X_8,Y_8,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P8(Equipment- Reactor, USD)')
saveas(gcf, 'PRIOR_nCaCO3_param8.png')
% EQ: Filter
figure(9)
[Y_9, X_9] = ksdensity(data_x9, 'npoints', np_2);
plot(X_9,Y_9,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P9(Equipment- Filter, USD)')
saveas(gcf, 'PRIOR_nCaCO3_param9.png')
% EQ: Blower
figure(10)
[Y_10, X_10] = ksdensity(data_x10, 'npoints', np_2);
plot(X_10,Y_10,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P10(Equipment- Blower, USD)')
saveas(gcf, 'PRIOR_nCaCO3_param10.png')
% EQ: Drum
figure(11)
[Y_11, X_11] = ksdensity(data_x11, 'npoints', np_2);
plot(X_11,Y_11,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P11(Equipment- Drum, USD)')
saveas(gcf, 'PRIOR_nCaCO3_param11.png')
% EQ: Centrifuge
figure(12)
[Y_12, X_12] = ksdensity(data_x12, 'npoints', np_2);
plot(X_12,Y_12,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P12(Equipment- Centrifuge, USD)')
saveas(gcf, 'PRIOR_nCaCO3_param12.png')
% EQ: Dryer
figure(13)
[Y_13, X_13] = ksdensity(data_x13, 'npoints', np_2);
plot(X_13,Y_13,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P13(Equipment- Dryer, USD)')
saveas(gcf, 'PRIOR_nCaCO3_param13.png')
% EQ: Crusher
figure(14)
[Y_14, X_14] = ksdensity(data_x14, 'npoints', np_2);
plot(X_14,Y_14,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P14(Equipment- Crusher, USD)')
saveas(gcf, 'PRIOR_nCaCO3_param14.png')
% EQ: HEATEX2
figure(15)
[Y_15, X_15] = ksdensity(data_x15, 'npoints', np_2);
plot(X_15,Y_15,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P15(Equipment- Heat Exchanger (CS&CS), USD)')
saveas(gcf, 'PRIOR_nCaCO3_param15.png')
% EQ: Condenser
figure(16)
[Y_16, X_16] = ksdensity(data_x16, 'npoints', np_2);
plot(X_16,Y_16,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P16(Equipment- Condenser, USD)')
saveas(gcf, 'PRIOR_nCaCO3_param16.png')
% Electricity Grid Mix Cost, USD/GJ
figure(17)
[Y_17, X_17] = ksdensity(data_x17, 'npoints', np_2);
plot(X_17,Y_17,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P17(Korean Grid Mix Electricity, USD/GJ)')
saveas(gcf, 'PRIOR_nCaCO3_param17.png')
% Steam Cost, USD/GJ
figure(18)
[Y_18, X_18] = ksdensity(data_x18, 'npoints', np_2);
plot(X_18,Y_18,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P18(Steam Costs, USD/GJ)')
saveas(gcf, 'PRIOR_nCaCO3_param18.png')
% GWI of Grid Mix Electricity
figure(19)
[Y_19, X_19] = ksdensity(data_x19, 'npoints', np_2);
plot(X_19,Y_19,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P19(GWI of Korean Grid Mix, kgCO2eq/GJ)')
saveas(gcf, 'PRIOR_nCaCO3_param19.png')
% GWI of Steam
figure(20)
[Y_20, X_20] = ksdensity(data_x20, 'npoints', np_2);
plot(X_20,Y_20,'k-','LineWidth',2);
title('Prior Uncertainty Distribution, P20(GWI of Steam, kgCO2eq/GJ)')
saveas(gcf, 'PRIOR_nCaCO3_param20.png')
%%Progress Bar%%
waitbar(100/total_steps, progress, 'Initializing Monte Carlo Simulations');




%% [3] Process Parameter Space for System Evaluation
% (a) Determine the # of System Model Outputs (# of Evaluation Metrics)
n_out = size(MODEL_nCaCO3(kp),2);
% (b) Create containers for Evaluated Model Outputs
% NOTE: eval_P and eval_C are n_sim x n_out x n_kp matrix. 
fx = zeros(n_sim, n_out);  % Evaluated outputs with all Params from ParSpace
fx_P = zeros(n_sim, n_out);% Evaluated outputs with i from ParSpace, ~i from c_ParSpace
fx_C = zeros(n_sim, n_out);% Evaluated outputs with i from c_ParSpacek, ~i from ParSpace

% (c) Evaluate Model from Monte Carlo Sampled Inputs (ParSpace)
parfor i = 1:n_sim
    % Each Parameter Set is a Row in the ParSpace matrix
    Parameter_Set = ParSpace(i,:); 
    fx(i,:) = MODEL_nCaCO3(Parameter_Set);
end
% (d) Generate Function Output Space based on i and ~i
for i = 1:n_sim
    for j = 1:n_kp
        %%%% fx_P = f(x_ik, x'_~ik)
        kp = [c_ParSpace(i,1:j-1), ParSpace(i,j), c_ParSpace(i,j+1:n_kp)];
        fx_P(i,:,j) = MODEL_nCaCO3(kp);
        %%%% fx_C = f(x'_ik, x_~ik)
        kp = [ParSpace(i,1:j-1), c_ParSpace(i,j), ParSpace(i,j+1:n_kp)];
        fx_C(i,:,j) = MODEL_nCaCO3(kp);
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
waitbar(750/total_steps, progress, 'Generating Bayesian Parametric Posteriors');




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
        saveas(gcf, ['POSTERIOR_nCaCO3_output_',num2str(i),'_parameter_',num2str(j),'.png'])
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
RESULTS.Input = ParSpace;
RESULTS.SobolTot = Sobol_Total;
RESULTS.RiskScore = Stake_Risk_Score;
RESULTS.Output = fx;
RESULTS.IQ = INTQUART;
RESULTS.TAILS = TAILS;
RESULTS.VAR = VAR;
RESULTS.CVAR = CVAR;
% Export to CSV
nCaCO3_Res = struct2cell(RESULTS);
writecell( nCaCO3_Res, 'nCaCO3_Results_Separate.csv');