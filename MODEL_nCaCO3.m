%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Process+Evaluation Model for Waste Steel Slag nCaCO3 Production %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function EvalMetrics = MODEL_nCaCO3(kp)   
%% Declare Global Variables
format long g
%% KEY PARAMETERS WITH UNCERTAINTY
% [1] NaCl Electrolysis Process
FE = kp(1); % faraday efficiency, 0.9
Spec_electrolysis = kp(2); %kWh/tonne nCaCO3, 1224

% [2] Extraction process
% [3] Purification process
% [4] Carbonation process
% [5] Recycling process
% [6] Crushing process

% [TEA] Techno-Economic Analysis Parameters
%%% [EQUIPMENT] Equipment Sizing Costs
REFCOST_TP = kp(3); % USD, 16.7 L/s based, 12700
REFCOST_cell = kp(4); %USD/kW, 1000
cost_fraction = kp(5); % portion of stack among overall electrolyzer cost, 0.25
REFCOST_HEATEX1 = kp(6); % USD, CS tube, sus304 shell, 163200
REFCOST_HEAT = kp(7); % USD, CS, 4800
REFCOST_reactor = kp(8); % USD, 75000
REFCOST_filter = kp(9); % USD, sus304, 29800
REFCOST_Blower = kp(10); % USD, CS, 1 tonne/hr nCaCO3 production based, 5400
REFCOST_Drum = kp(11); % USD, 100000
REFCOST_Centrifuge = kp(12); % USD per one, 118400
REFCOST_Dryer = kp(13); % USD, 74800
REFCOST_Crusher = kp(14); % USD, 40000
REFCOST_HEATEX2 = kp(15); % USD, 88300, CS & CS
REFCOST_COND = kp(16); % USD, 11200
%%% [OPERATING] Utility Costs
UT_ELEC = kp(17); % USD/GJ, 19.650
UT_LP_STEAM = kp(18); % USD/GJ, 6.08

% [LCA] CO2 Life Cycle Assessment Parameters
%%% [UTILITY] Global Warming Potentials for Utility Consumption
GWI_ELEC_GRID = kp(19); % kgCO2-eq/GJ, 108.800
GWI_STEAM_LP = kp(20); % kgCO2-eq/GJ,  82.894




%% NON VARYING FIXED-PARAMETERS & DESIGN SPECS
% [0] General
Convert_frate = 0.277778;             % cum/hr to L/s
Convert_energy = 0.0036;              % GJ/hr to kW

% [1] NaCl Electrolysis Process
Molar_weight_water = 0.018;           % tonne/kmol
Molar_weight_nCaCO3 = 0.10012;        % tonne/kmol 
Molar_weight_NaCl = 0.047648344;      % tonne/kmol
Viscosity = 0.002;                    % viscosity, Pa*s
Conc_HCl = 0.5;                       % kmol/cum, concentration of HCl product
Conc_NaCl = 1.0;                      % kmol/cum, concentration of NaCl
Conc_NaOH = 1.0;                      % kmol/cum, concentration of NaOH product
V_NaCl_initial_ref = 30.2;            % cum/hr, volumetric flow of feed

% [2] Extraction Process
Slag_mass_inlet_ref = 2.72;           % tonne/hr
Wf_Slag_CaO = 0.471;                  % weight fraction of CaO in steel slag
Wf_Slag_SiO2 = 0.315;                 % weight fraction of SiO2 in steel slag
Wf_Slag_Al2O3 = 0.152;                % weight fraction of Al2O3 in steel slag
Wf_Slag_K2O = 0.005;                  % weight fraction of K2O in steel slag
Wf_Slag_MgO = 0.0039;                 % weight fraction of MgO in steel slag
Wf_Slag_SO3 = 0.0064;                 % weight fraction of SO3 in steel slag
Wf_Slag_TiO2 = 0.0018;                % weight fraction of TiO2 in steel slag
Wf_Slag_MnO = 0.0038;                 % weight fraction of MnO in steel slag
Wf_Slag_FeO = 0.006;                  % weight fraction of FeO in steel slag
Ex_FC_Ca = 0.493865;                  % extraction efficieny of Ca2+ ion
Ex_FC_Mg = 0.4973;                    % extraction efficieny of Mg2+ ion
Ex_FC_Fe = 0.25;                      % extraction efficieny of Fe2+ ion
Ex_FC_Al = 0.077;                     % extraction efficieny of Al3+ ion
L_D = 2.5;                            % ratio of length and diameter
d_D = 0.5;                            % ratio of paddle diameter and diameter
b_D = 0.2 * d_D;                      % ratio of paddle width and diameter
Filter_FR_liq = 0.99999;              % liquid to liquid
Filter_FR_sol = 0.99999;              % Solid to solid

% [3] Impurity Removal and Separation
IMRE_FC_Ca = 0.0304;                  % conversion of Ca2+ ion
Filter_IMRE_liq = 0.99999;            % liquid to liquid
Filter_IMRE_sol = 0.99999;            % Solid to solid

% [4] Carbonation Process
V_NaOH_nCaCO3_ref = 20.0;             % cum/hr, volumetric flow of NaOH which is fed to carbonation process
Wf_flue_H2O = 0.174;                  % weight fraction of H2O in flue gas
Wf_flue_CO2 = 0.186;                  % weight fraction of CO2 in flue gas
Wf_flue_N2 = 0.611;                   % weight fraction of N2 in flue gas
Wf_flue_O2 = 0.029;                   % weight fraction of O2 in flue gas
V_flue_ref = 3043.32;                 % cum/hr, volumetric flow of flue gas

% [5] Recycling Process
% [6] Crushing Processs
Work_index = 30.4;                    % kWh/tonneSlag
W_F80 = 20000;                        % 1E-6 m
W_P80 = 30;                           % 1E-6 m

% [TEA] Techno-Economic Analysis Parameters
%%% Chemical Engineering Plant Cost Indices (CEPCI)
CEPCI_ref = 1000;
CEPCI = 619.2;                 % 2019 based
CEPCI_ref2 = 576.1;            % 2014 based
%%% Nth-Plant Assumption Exponential Constants
NTH_HTX = 0.71;
NTH_reactor = 0.53;
NTH_Drum = 0.81;
NTH_Dryer = 0.4;
NTH_Crusher = 0.33;
%%% Capital Cost Lang Factors
REFfrate_TP = 16.7;            % L/s
Lifetime_stack = 6;
LIFET = 30;
DISCO = 0.07; % annual discount rate
%%% Capital Cost Lang Factors
F_m_sus304_TP = 1.9;
F_m_sus304_reactor = 1.8;
OSBL_OS = 0.4;                        % Offsite Operations, from [23]
OSBL_DE = 0.25;                       % Design and Engineering, from [23]
OSBL_CN = 0.1;                        % Contingency, from [23]
f_IEC = 1.1; % Installation factor
%%% Equipment Costs
REFarea_HEATEX1 = 169.14; % m2
REFarea_HEAT = 10.5; % m2
REFvolume_reactor = 3; % cum;
REFarea_filter = 11.32; % m2
REFfrate_Blower = 845.3673; % L/s
REFfactor_Drum = 20; % height * diameter^1.5
REFNum_Centrifuge = 12; % number of centrifuge
REFdiameter_Centrifuge = 1.28; % m
REFArea_Dryer = 5.78; % m2
REFpower_Crusher = 7.5; % kW
REFarea_HEATEX2 = 271.27; % m2
REFarea_COND = 2.63; % m2
%%% [OPERATING COST] Raw Material Costs
Operating_hour = 8000; % hr/year
UT_slag_treatment = 0.05; % USD/ton
OP_RM_NaCl = 60; % USD/ton
OP_RM_H2O = 0.15; % USD/m3
delta_T_min = 20; % K
heat_capacity = 0.00418; % GJ/ton K
REFCOST_Labor = 970797; % USD/year, 200,000 tonne/yr production based
f_Maintenance = 0.03;
f_admin = 0.2;
f_overhead = 0.6;
f_Laboratory = 0.01;

% [LCA] CO2 Life Cycle Assessment Parameters
%%% [MATERIAL] Global Warming Potentials for Material Consumption
GWI_WATER = 0;
GWI_steel_slag = 0.0589; %kgCO2-eq/kg
GWI_disposal = 6.84E-3; %kgCO2-eq/kg
GWI_CO2_FG = -1; %kgCO2-eq/kg




%% DERIVED PARAMETERS AND DESIGN SPECS
%% PROCESS MASS BALANCE MODEL
% [1] NaCl Electrolysis Process
% (a) Setup NaCl Reaction
target_product_amount = 1;            % ton/hr, hourly nCaCO3 production amount
component1 = 6;                       % molecule component considering in NaCl electrolysis process (water, NaCl, H+, Cl-, Na+, OH-)
Conc_NaCl_ref = 1;                    % kmol/cum, based concentration
Rho_NaCl = 1;                         % ton/cum, concentration of NaCl solution
V_NaCl_initial = V_NaCl_initial_ref * target_product_amount/Conc_NaCl;
C_NaCl_initial = [V_NaCl_initial * Rho_NaCl/Molar_weight_water, V_NaCl_initial * Conc_NaCl]; % kmol/hr
C_NaCl_out = size(component1,1); % kmol/hr
% (b) Electrolysis Reaction
for i = 1:component1
    if i == 1
        C_NaCl_out(i) = C_NaCl_initial(1) - C_NaCl_initial(2) * FE;    
    elseif i == 2
           C_NaCl_out(i) = C_NaCl_initial(2) * (1-FE);
    else
         C_NaCl_out(i) = C_NaCl_initial(2) * FE;
    end
end
% (c) Outlet Stream Flowrate
exit_HCl_1 = size(component1,1);                     % kmol/hr
exit_NaOH_1 = size(component1,1);                    % kmol/hr
exit_NaCl = size(component1,1);                      % kmol/hr
Sep_fraction = [0.9429, 0, 1, 1, 0, 0;
                0.0347, 1, 0, 0, 0.064, 0.064;];     % Separation fraction to HCl and NaCl,respectively
for i = 1:component1
        exit_HCl_1(i) = C_NaCl_out(i) * Sep_fraction(1,i);
        exit_NaCl(i) = C_NaCl_out(i) * Sep_fraction(2,i);
        exit_NaOH_1(i) = C_NaCl_out(i) - exit_HCl_1(i) - exit_NaCl(i);
end
V_sWater_HCl = exit_HCl_1(3) * (1/Conc_HCl) - exit_HCl_1(1) * Molar_weight_water; % supplementary water for meeting the concentration, cum/hr
M_sWater_HCl = V_sWater_HCl * (1/Molar_weight_water);
V_sWater_NaOH = exit_NaOH_1(5) * (1/Conc_NaOH) - exit_NaOH_1(1) * Molar_weight_water; % supplementary water for meeting the concentration, cum/hr
M_sWater_NaOH = V_sWater_NaOH * (1/Molar_weight_water);
for i = 1:component1
    if i == 1
        exit_HCl(i) = exit_HCl_1(i) + M_sWater_HCl;
        exit_NaOH(i) = exit_NaOH_1(i) + M_sWater_NaOH;
    else
        exit_HCl(i) = exit_HCl_1(i);
        exit_NaOH(i) = exit_NaOH_1(i);
    end
end
V_exit_HCl = V_NaCl_initial * Sep_fraction(1) + V_sWater_HCl;
V_exit_NaOH = exit_NaOH(1) * Molar_weight_water;
V_exit_NaCl = exit_NaCl(1) * Molar_weight_water;

% [2] Metal Ion Extraction Process
% (a) Setup MIE Reaction
Slag_mass_inlet = Slag_mass_inlet_ref * target_product_amount; %tonne/hr
component2 = 16;
W_Slag_inlet = size(component2); %tonne/hr
Wf_Slag = [Wf_Slag_CaO, Wf_Slag_SiO2, Wf_Slag_Al2O3, Wf_Slag_K2O, Wf_Slag_MgO, Wf_Slag_SO3, Wf_Slag_TiO2, Wf_Slag_MnO, Wf_Slag_FeO, 0, 0, 0, 0, 0, 0, 0];
Molar_weight_slag = [0.056077, 0.060084, 0.101961, 0.094196, 0.040304, 0.080064, 0.079879, 0.070937, 0.071846, 1, 1, 1, 1, 1, 1, 1]; % tonne/kmol
MW_Ex = [0.056077, 0.060084, 0.101961, 0.094196, 0.040304, 0.080064, 0.079879, 0.070937, 0.071846, 0, 0, 0, 0, 0, 0, Molar_weight_water]; %tonne/kmol
% (b) Metal Ion Extractor Mass Balance
Ex_FC = [Ex_FC_Ca, Ex_FC_Mg, Ex_FC_Fe, Ex_FC_Al];
Ex_Sto = [-1, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 1, 0, 0, 0, 1;
          0, 0, 0, 0, -1, 0, 0, 0, 0, -2, 0, 0, 1, 0, 0, 1;
          0, 0, 0, 0, 0, 0, 0, 0, -1, -2, 0, 0, 0, 0, 1, 1;
          0, 0, -1, 0, 0, 0, 0, 0, 0, -6, 0, 0, 0, 2, 0, 3];
Num_Ex_Rxn = size(Ex_Sto,1);
for i = 1:component2
    if i == 10
        C_Slag_inlet(i) = exit_HCl(3);
        W_Slag_inlet(i) = C_Slag_inlet(i);
    elseif i == 11
        C_Slag_inlet(i) = exit_HCl(4);
        W_Slag_inlet(i) = C_Slag_inlet(i);
    elseif i == 16
        C_Slag_inlet(i) = exit_HCl(1);
        W_Slag_inlet(i) = C_Slag_inlet(i) * Molar_weight_water;
    else
        W_Slag_inlet(i) = Slag_mass_inlet * Wf_Slag(i);
        C_Slag_inlet(i) = W_Slag_inlet(i) / Molar_weight_slag(i); %kmol/hr
    end
end
Ex_key = [C_Slag_inlet(1), C_Slag_inlet(5), C_Slag_inlet(9), C_Slag_inlet(3)];
for i = 1:Num_Ex_Rxn
    if i == 1
        exit_Ex = C_Slag_inlet + Ex_key(i) * Ex_Sto(i,:) * Ex_FC(i);
    else
        exit_Ex = exit_Ex + Ex_key(i) * Ex_Sto(i,:) * Ex_FC(i);
    end
end
% (c) Filter Mass Balance
for i = 1:component2
    if i <= 9
        exit_FT_solid(i) = exit_Ex(i);
    else
        exit_FT_solid(i) = (1-Filter_FR_liq)*exit_Ex(i);
    end
end
for i = 1:component2
    if i == 1
        W_waste = MW_Ex(i)*exit_FT_solid(i);
    else
        W_waste = W_waste + MW_Ex(i)*exit_FT_solid(i);
    end
end
for i = 1:7
    exit_FT_liquid(i) = Filter_FR_liq * exit_Ex(i+9); %electrolyte part (H+,Cl-, Ca2+, Mg2+, Al3+, Fe2+, water)
end
% (d) NaOH Split Mass Balance
Cf_NaOH = V_NaOH_nCaCO3_ref/exit_NaOH(5);
exit_NaOH_impu = (1-Cf_NaOH) * exit_NaOH;
exit_NaOH_Car = Cf_NaOH * exit_NaOH;
V_NaOH_IMRE = V_exit_NaOH * (1 - Cf_NaOH);
V_NaOH_CR = V_exit_NaOH * Cf_NaOH;

% [3] Impurity Removal and Separation
% (a) Impurity Removal
IMRE_FC = [IMRE_FC_Ca, 1, 1, 1, 1];
component3 = 13;
MW_IMRE = [0.0, 0.035453, 0.040077, 0.0, 0.0, 0.0, Molar_weight_water, 0.022989, 0.017008, 0.074092, 0.05832, 0.078003, 0.089862]; %tonne/kmol
for i = 1:component3
    if i < 7
        inlet_IMRE(i) = exit_FT_liquid(i); %H+,Cl-, Ca2+, Mg2+, Al3+, Fe2+, water, Na+, OH-, Ca(OH)2, MG(OH)2, AL(OH)3, FE(OH)2 kmol/hr
    elseif i == 7
        inlet_IMRE(i) = exit_FT_liquid(i) + exit_NaOH_impu(1);
    elseif i == 8
        inlet_IMRE(i) = exit_NaOH_impu(5);
    elseif i == 9
        inlet_IMRE(i) = exit_NaOH_impu(6);
    else
        inlet_IMRE(i) = 0;
    end 
end
IMRE_key = [exit_FT_liquid(3), exit_FT_liquid(4), exit_FT_liquid(5), exit_FT_liquid(6), exit_FT_liquid(1)];
IMRE_Sto = [0, 0, -1, 0, 0, 0, 0, 0, -2, 1, 0, 0, 0;
            0, 0, 0, -1, 0, 0, 0, 0, -2, 0, 1, 0, 0;
            0, 0, 0, 0, -1, 0, 0, 0, -3, 0, 0, 1, 0;
            0, 0, 0, 0, 0, -1, 0, 0, -2, 0, 0, 0, 1;
            -1, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0];
Num_IMRE_Rxn = size(IMRE_Sto,1);
for i = 1:Num_IMRE_Rxn
    if i == 1
        exit_IMRE = inlet_IMRE + IMRE_key(i) * IMRE_Sto(i,:) * IMRE_FC(i);  %kmol/hr
    else
        exit_IMRE = exit_IMRE + IMRE_key(i) * IMRE_Sto(i,:) * IMRE_FC(i);
    end
end
% (b) Filter Mass Balance
for i = 1:component3
    if i <= 9
        exit_IMRE_solid(i) = (1-Filter_IMRE_liq)*exit_IMRE(i);
    else
        exit_IMRE_solid(i) = exit_IMRE(i);
    end
end
for i = 1:component3
    if i == 1
        W_impurity = MW_IMRE(i) * exit_IMRE_solid(i);
    else
        W_impurity = W_impurity + MW_IMRE(i) * exit_IMRE_solid(i); %tonne/hr
    end
end
for i = 1:9
    exit_IMRE_liquid(i) = Filter_IMRE_liq * exit_IMRE(i); %electrolyte part (H+,Cl-, Ca2+, Mg2+, Al3+, Fe2+, water, Na+, OH-)
end
V_exit_IMRE_liquid = V_exit_HCl * Filter_IMRE_liq;

% [4] Carbonation Process
% (a) Inlet Flue Gas
Component_carbon = 14;
Component_flue = 4;
Wf_flue = [Wf_flue_H2O, Wf_flue_CO2, Wf_flue_N2, Wf_flue_O2]; % H2O, CO2, N2, O2/무조건 H2O 먼저
V_flue = V_flue_ref * target_product_amount; %cum/hr
Flue_mass_inlet = 0.44 * target_product_amount / Wf_flue(2);
Molar_weight_flue = [0.018, 0.0440088, 0.0280, 0.0322]; %tonne/kmol
W_CO2 = Flue_mass_inlet * Wf_flue(2);
for i = 1:Component_carbon
    if i == 1
        W_flue_inlet(i) = Flue_mass_inlet * Wf_flue(i); %tonne/hr
        C_CR_inlet(i) = exit_NaOH_Car(i) + exit_IMRE_liquid(7) + W_flue_inlet(i)/Molar_weight_flue(i); % kmol/hr
        % H2O, CO2, N2, O2, H+,Cl-, Ca2+, Mg2+, Al3+, Fe2+, Na+, OH-, Ca(OH)2, CaCO3
    elseif 2 <= i && i<= 4
        W_flue_inlet(i) = Flue_mass_inlet * Wf_flue(i);
        C_CR_inlet(i) = W_flue_inlet(i)/Molar_weight_flue(i);
    elseif Component_flue < i && i <= 10
        W_flue_inlet(i) = 0;
        C_CR_inlet(i) = exit_IMRE_liquid(i-Component_flue);
    elseif 10 < i && i <= 12
        W_flue_inlet(i) = 0;
        C_CR_inlet(i) = exit_IMRE_liquid(i+1-Component_flue) + exit_NaOH_Car(i-6);
    else
        W_flue_inlet(i) = 0;
        C_CR_inlet(i) = 0;
    end
end
% (b) Carbonation Reactor
CR_key = [C_CR_inlet(12), 0];
CR_FC = [0.5, 1.0];
CR_Sto = [0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, -2, 1, 0; 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1];
Num_CR_Rxn = size(CR_Sto,1);
for i = 1: Num_CR_Rxn
    if i == 1
       exit_CR = C_CR_inlet + CR_FC(i) * CR_key(i) * CR_Sto(i,:);
       CR_key(2) = exit_CR(2);
    else
       exit_CR = exit_CR + CR_FC(i) * CR_key(i) * CR_Sto(i,:);
    end
end
% (c) Gas-Liquid Separation (DRUM1)
for i = 1:Component_carbon
    if 2 <= i && i<=Component_flue
        exit_CR_vapor(i) = exit_CR(i);
        exit_CR_liquid(i) = 0;
    else
        exit_CR_vapor(i) = 0;
        exit_CR_liquid(i) = exit_CR(i);
    end
end
% (d) Centrifuge Mass Balance
Filter_CR_liq = 0.99999; % Liquid to Liquid
Filter_CR_sol = 0.99999; % Solid to Solid
for i = 1:Component_carbon
    if i <= Component_carbon-2
        exit_CR_recycle(i) = Filter_CR_liq * exit_CR_liquid(i);
        exit_CR_nCaCO3(i) = (1-Filter_CR_liq) * exit_CR_liquid(i);
    else
        exit_CR_recycle(i) = (1-Filter_CR_sol) * exit_CR_liquid(i);
        exit_CR_nCaCO3(i) = Filter_CR_sol * exit_CR_liquid(i);
    end
end
% (e) Product Streams
C_nCaCO3 = exit_CR_nCaCO3(Component_carbon);   % kmol/hr
W_nCaCO3 = C_nCaCO3 * Molar_weight_nCaCO3;     % tonne/hr

% [5] Recycling Process
% (a) Evaporation 1
Cf_EVP1 = [1-0.987602, 0, 0.7658680484, 0.5766565604, 0, 0,0, 0, 0, 0, 0, 0, 0, 0];   % Evaporating fraction to gas 
% H2O, CO2, N2, O2, H+,Cl-, Ca2+, Mg2+, Al3+, Fe2+, Na+, OH-, Ca(OH)2, CaCO3
exit_EVP1_gas = Cf_EVP1 .* exit_CR_recycle;
exit_EVP1_liquid = exit_CR_recycle - exit_EVP1_gas;
V_outlet_EVP1_liquid = exit_EVP1_liquid(1) * Molar_weight_water; % cum/hr
% (b) Evaporation 2
Cf_EVP2 = [0.7765026121, 0, 0.7658680484, 0.5766565604, 0, 0,0, 0, 0, 0, 0, 0, 0, 0]; % Evaporating fraction to gas 
% H2O, CO2, N2, O2, H+,Cl-, Ca2+, Mg2+, Al3+, Fe2+, Na+, OH-, Ca(OH)2, CaCO3
V_inlet_EVP2 = 220609;                                           % cum/hr
exit_EVP2_gas = Cf_EVP2 .* exit_EVP1_liquid;                     % kmol/hr
exit_EVP2_liquid = exit_EVP1_liquid - exit_EVP2_gas;
V_outlet_EVP2_liquid = exit_EVP2_liquid(1) * Molar_weight_water; % cum/hr
V_outlet_EVP2_gas = V_inlet_EVP2 - V_outlet_EVP2_liquid;
% (c) Water Split (B25)
Cf_water_split = (C_NaCl_initial(1) - exit_NaCl(1) - exit_EVP2_liquid(1))/exit_EVP2_gas(1);
exit_recycle_water = exit_EVP2_gas * Cf_water_split;
exit_water_outlet = exit_EVP2_gas - exit_recycle_water;
% (d) Recycled NaCl
recycled_NaCl_H2O = exit_EVP2_liquid(1) + exit_NaCl(1) + exit_recycle_water(1);
recycled_NaCl_NaCl = exit_EVP2_liquid(6) + exit_NaCl(2);
exit_recycle_NaCl = [recycled_NaCl_H2O, recycled_NaCl_NaCl];
W_NaCl_makeup = (C_NaCl_initial(2) - exit_recycle_NaCl(2))*Molar_weight_NaCl; % tonne/hr
% ORGANIZE PRODUCT STREAMS




%% PROCESS ENERGY BALANCE MODEL
% [1] ENERGY REQUIREMENTS: NaCl Electrolysis Process
% (a) NaCl Solution Transportation Pump (CPUMP1)
Rating_TP1 = 0.805333;                                % MJ/hr per pump, 0.2 bar
TP1_unit_frate = V_NaCl_initial_ref;                  % cum/hr capacity
Scale_factor = V_NaCl_initial/TP1_unit_frate;
num_TP1 = ceil(Scale_factor);                         % Number of TP Modules
Elec_TP1 = Rating_TP1 * num_TP1 * 0.001;              % GJ/hr electricity consumption
% (b) Heat exchanger heat duty (HEATEX1)
Heatex1_duty_ref = 2.92938;                           % Exchanging heat for 1 tonne nCaCO3 production, GJ/hr
Heatex1_duty = Heatex1_duty_ref * Scale_factor;
% (c) Heater before electrolyzer (HEATER)
HEATER1_duty_ref = 0.784579;                          % GJ/hr
LP_HEATER1_duty = HEATER1_duty_ref * Scale_factor;
% (d) Cooler after electrolyzer (COOL1)
COOL1_duty_ref = -0.59049;                            % GJ/hr
Cool_COOL1_duty = COOL1_duty_ref * Scale_factor;
% (e) NaCl electrolyzer
Elec_electrolysis_TEA = Spec_electrolysis * W_nCaCO3; % kWh/tonne nCaCO3
Elec_electrolysis = Elec_electrolysis_TEA * 0.0036;   % GJ/hr

% [2] ENERGY REQUIREMENTS: Metal Ion Extraction Process
% (a) HCl Solution Transportation Pump (CPUMP2)
Rating_TP2 = 1.44579;                    % MJ/hr per pump, 0.2 bar
TP2_unit_frate = V_exit_HCl;
num_TP2 = ceil(Scale_factor);
Elec_TP2 = Rating_TP2 * num_TP2 * 0.001;
% (b) Extraction Unit
HRT_Ex = 1;                              % residence time, hour
Ex_unit_frate = TP2_unit_frate * HRT_Ex; % cum/hr capacity
Density_Ex = 1026.05;                    % kg/cum
Volume_Ex = Ex_unit_frate * HRT_Ex;      % volume of extraction reactor, cum
l_v_Ex = 12.1737;                        % m/s, linear velocity, literature value
Diameter_Ex = ((Volume_Ex * 4)/(2.5 * pi))^(1/3);      % m
Length_Ex = L_D * Diameter_Ex;           % m
d_Ex = d_D * Diameter_Ex;                % m
b_Ex = b_D * Diameter_Ex;                % m
dr_Ex = d_Ex/2;                          % m, radius of paddle
wv_Ex = l_v_Ex /(dr_Ex * 0.10472);       % rpm, angular velocity
n_Ex = wv_Ex * 0.01666;                  % s-1, stirring rate
A_Ex = 14 + b_D * (670 * ((d_D - 0.6)^2) + 185);
B_Ex = 10^(1.3 - (4 * (b_D - 0.5)^2) - (1.14 * d_D));
P_Ex = 1.1 + 4 * b_D - 2.5 * (d_D-0.5)^2 - 7 * b_D^4;
Re_Ex = (Density_Ex * n_Ex * d_Ex * d_Ex) / Viscosity; % Reynolds number of extraction reactor
N_p_Ex = (A_Ex/Re_Ex)+(B_Ex*((1000+1.2*Re_Ex^0.66)/(1000+3.2*Re_Ex^0.66))^P_Ex)*(L_D^(0.35+b_D));
Ws_Ex = 0.001 * N_p_Ex * Density_Ex * n_Ex^3 * d_Ex^5; % shaft work of extraction reactor, kWh
Elec_Ex = Ws_Ex * 0.0036;                              % GJ/hr

% [3] ENERGY REQUIREMENTS: Impurity Removal and Separation
% (a) NaOH Transportaion Pump (CPUMP3)
Rating_TP3 = 0.665407;               % MJ/hr per pump, 0.2 bar
TP3_unit_frate = V_exit_NaOH;
num_TP3 = ceil(Scale_factor);
Elec_TP3 = Rating_TP3 * num_TP3 * 0.001;
% (b) Metal Ion Solution Transportaion Pump (CPUMP4)
Rating_TP4 = 1.43298;                % MJ/hr per pump, 0.2 bar
TP4_unit_frate = V_exit_IMRE_liquid;
num_TP4 = ceil(Scale_factor);
Elec_TP4 = Rating_TP4 * num_TP4 * 0.001;
% (c) Impurity Removal Process (IMRE)
HRT_IMRE = 1;                        % residence time, hour
IMRE_unit_frate = TP4_unit_frate + V_NaOH_IMRE * HRT_IMRE; % cum/hr capacity
Density_IMRE = 1020.05;                                    % kg/cum
Volume_IMRE = IMRE_unit_frate * HRT_IMRE;                  % volume of extraction reactor, cum
l_v_IMRE = 8.2467;                                         % m/s, linear velocity, literature value
Diameter_IMRE = ((Volume_IMRE * 4)/(2.5 * pi))^(1/3);      % m
Length_IMRE = L_D * Diameter_IMRE;                         % m
d_IMRE = d_D * Diameter_IMRE;                              % m
b_IMRE = b_D * Diameter_IMRE;                              % m
dr_IMRE = d_IMRE/2;                                        % m, radius of paddle
wv_IMRE = l_v_IMRE /(dr_IMRE * 0.10472);                   % rpm, angular velocity
n_IMRE = wv_IMRE * 0.01666;                                % s-1, stirring rate
A_IMRE = 14 + b_D * (670 * ((d_D - 0.6)^2) + 185);
B_IMRE = 10^(1.3 - (4 * (b_D - 0.5)^2) - (1.14 * d_D));
P_IMRE = 1.1 + 4 * b_D - 2.5 * (d_D-0.5)^2 - 7 * b_D^4;
Re_IMRE = (Density_IMRE * n_IMRE * d_IMRE * d_IMRE) / Viscosity; % Reynolds number of extraction reactor
N_p_IMRE = (A_IMRE/Re_IMRE)+(B_IMRE*((1000+1.2*Re_IMRE^0.66)/(1000+3.2*Re_IMRE^0.66))^P_IMRE)*(L_D^(0.35+b_D));
Ws_IMRE = 0.001 * N_p_IMRE * Density_IMRE * n_IMRE^3 * d_IMRE^5; % shaft work of extraction reactor, kWh
Elec_IMRE = Ws_IMRE * 0.0036;                                    % GJ/hr

% [4] ENERGY REQUIREMENTSL Carbonation Process
% (a) Ca Rich Ion Transportation Pump (CPUMP5)
Rating_TP5 = 1.57731;                              % MJ/hr per pump, 0.2 bar
TP5_unit_frate = V_exit_IMRE_liquid;
num_TP5 = ceil(Scale_factor);
Elec_TP5 = Rating_TP5 * num_TP5 * 0.001;           % GJ/hr
% (b) CO2 Blower (BLOWER)
Rating_Blower = 213.603;                           % MJ/hr per blower,70 kPa
Blower_unit_frate = V_flue;
num_Blower = ceil(Scale_factor);
Elec_Blower = Rating_Blower * num_Blower * 0.001;  % GJ/hr
% (c) Cooler (COOL2)
COOL2_duty_ref = -0.00040334;                      % GJ/hr
Cool_COOL2_duty = COOL2_duty_ref * Scale_factor;
V_outlet_flue = 967.019 * target_product_amount;
% (d) Carbonation Reactor (CR)
HRT_CR = 2;                                        % residence time, hour
CR_unit_frate = (TP5_unit_frate + V_NaOH_CR) * HRT_CR; % cum/hr capacity
Density_CR = 1022.13;                                  % kg/cum
Volume_CR = CR_unit_frate * HRT_CR;                    % volume of extraction reactor, cum
l_v_CR = 5.34072;                                      % m/s, linear velocity, literature value
Diameter_CR = ((Volume_CR * 4)/(2.5 * pi))^(1/3);      % m
Length_CR = L_D * Diameter_CR;                         % m
d_CR = d_D * Diameter_CR;                              % m
b_CR = b_D * Diameter_CR;                              % m
dr_CR = d_CR/2;                                        % m, radius of paddle
wv_CR = l_v_CR /(dr_CR * 0.10472);                     % rpm, angular velocity
n_CR = wv_CR * 0.01666;                                % s-1, stirring rate
A_CR = 14 + b_D * (670 * ((d_D - 0.6)^2) + 185);
B_CR = 10^(1.3 - (4 * (b_D - 0.5)^2) - (1.14 * d_D));
P_CR = 1.1 + 4 * b_D - 2.5 * (d_D-0.5)^2 - 7 * b_D^4;
Re_CR = (Density_CR * n_CR * d_CR * d_CR) / Viscosity; % Reynolds number of extraction reactor
N_p_CR = (A_CR/Re_CR)+(B_CR*((1000+1.2*Re_CR^0.66)/(1000+3.2*Re_CR^0.66))^P_CR)*(L_D^(0.35+b_D));
Ws_CR = 0.001 * N_p_CR * Density_CR * n_CR^3 * d_CR^5; % shaft work of extraction reactor, kWh
Elec_CR = Ws_CR * 0.0036;                              % GJ/hr
% (e) nCaCO3 drying (DRYER)
DRY_duty_ref = 0.0415254;                              % GJ/hr
LP_DRY_duty = DRY_duty_ref * Scale_factor;

% [5] ENERGY REQUIREMENTS: Recycling Process
% (a) Heat Exchanger Heat Duty (HEATEX2)
Heatex2_duty_ref = 13.6753;                            % Exchanging heat for 1 tonne nCaCO3 production, GJ/hr
Heatex2_duty = Heatex2_duty_ref * Scale_factor;
% (b) Heater after HEATEX (HEAT2)
HEATER2_duty_ref = 1.65881;                            % GJ/hr
LP_HEATER2_duty = HEATER2_duty_ref * Scale_factor;
% (c) Condenser after EVP1 (COND1)
COND1_duty_ref = -0.161368;                            % GJ/hr
Mole_inlet_COND1 = sum(exit_EVP1_gas, 'all');
Mole_inlet_COND1_ref = 7.64;
Cool_COND1_duty = COND1_duty_ref * (Mole_inlet_COND1/Mole_inlet_COND1_ref);
% (d) Heater after EVP1 (HEAT3)
HEATER3_duty_ref = 141.447;
Mole_inlet_HEAT3 = sum(exit_EVP1_liquid, 'all');
Mole_inlet_HEAT3_ref = 4422.9;                         % kmol/hr
LP_HEATER3_duty = HEATER3_duty_ref * (Mole_inlet_HEAT3/Mole_inlet_HEAT3_ref);
% (e) Water Transporation Pump (CPUMP6)
Rating_TP6 = 7.07222;                                  % MJ/hr per pump, 0.2 bar
TP6_unit_frate = exit_EVP2_gas(1) * Molar_weight_water;
num_TP6 = ceil(Scale_factor);
Elec_TP6 = Rating_TP6 * num_TP6 * 0.001;               % GJ/hr
% (f) Recycling NaCl Transporation Pump (CPUMP7)
Rating_TP7 = 2.74101;                                  % MJ/hr
TP7_unit_frate = V_outlet_EVP2_liquid;
num_TP7 = ceil(Scale_factor);
Elec_TP7 = Rating_TP7 * num_TP7 * 0.001;               % GJ/hr

% [6] ENERGY REQUIREMENTS: Crushing Process
Elec_crusher = 10 * Work_index * ((1/W_P80^0.5)-(1/W_F80^0.5)) * Slag_mass_inlet *3.6*0.001; %GJ/hr

% ORGANIZE ENERGY DUTIES FOR UTILITY CALCULATIONS
Elec_duty_cell = Elec_electrolysis;                 % GJ/hr
Elec_duty_other = Elec_Blower + Elec_CR + Elec_crusher + Elec_Ex + Elec_IMRE + Elec_TP1 + Elec_TP2 + Elec_TP3 + Elec_TP4 + Elec_TP5 + Elec_TP6 + Elec_TP7;
LP_duty = LP_DRY_duty + LP_HEATER1_duty + LP_HEATER2_duty;
Unit_elec_duty_cell = Elec_duty_cell / W_nCaCO3;    % GJ/ton nCaCO3
Unit_elec_duty_others = Elec_duty_other / W_nCaCO3; % GJ/ton nCaCO3
Unit_LP_duty = LP_duty / W_nCaCO3;                  % GJ/ton nCaCO3




%% TECHNOECONOMIC EVALUATION [TEA] MODEL
%============================ CAPITAL EXPENSES ===========================%
% [1] Electrolysis Process
% (a) CPUMP1
TP1_inlet_frate = TP1_unit_frate * Convert_frate;
FOB_TP1 = REFCOST_TP * (TP1_inlet_frate/REFfrate_TP)^0.59 * CEPCI/CEPCI_ref;
% (b) Electrolysis cell
FOB_cell = REFCOST_cell * Elec_electrolysis_TEA + ((LIFET/Lifetime_stack)-1) * cost_fraction * REFCOST_cell * Elec_electrolysis_TEA;
% (c) HEATEX1
Area_HEATEX1 = REFarea_HEATEX1 * Scale_factor;
FOB_HEATEX1 = REFCOST_HEATEX1 * (Area_HEATEX1/REFarea_HEATEX1)^NTH_HTX * CEPCI/CEPCI_ref2;
% (d) COOL1
Area_COOL1 = 4.54 * Scale_factor;
FOB_COOL1 = REFCOST_HEAT * (Area_COOL1/REFarea_HEAT)^NTH_HTX * CEPCI/CEPCI_ref2;
EQ_Elec = FOB_TP1 + FOB_cell + FOB_HEATEX1 +  FOB_COOL1;

% [2] Metal Ion Extraction Process
% (a) CPUMP2
TP2_inlet_frate = TP2_unit_frate * Convert_frate;
FOB_TP2 = REFCOST_TP * (TP2_inlet_frate/REFfrate_TP)^0.59 * CEPCI/CEPCI_ref * F_m_sus304_TP;
% (b) Extractor
FOB_Ex = REFCOST_reactor * (Volume_Ex/REFvolume_reactor)^NTH_reactor * CEPCI/CEPCI_ref * F_m_sus304_reactor;
% (c) Filter1
Area_filter1 = 11.32 * Scale_factor;
FOB_filter1 = REFCOST_filter * (Area_filter1/REFarea_filter) * CEPCI/CEPCI_ref2;
EQ_Ex = FOB_TP2 + FOB_Ex + FOB_filter1;

% [3] Impurity Removal and Separation
% (a) CPUMP3
TP3_inlet_frate = TP3_unit_frate * Convert_frate;
FOB_TP3 = REFCOST_TP * (TP3_inlet_frate/REFfrate_TP)^0.59 * CEPCI/CEPCI_ref * F_m_sus304_TP;
% (b) CPUMP4
TP4_inlet_frate = TP4_unit_frate * Convert_frate;
FOB_TP4 = REFCOST_TP * (TP4_inlet_frate/REFfrate_TP)^0.59 * CEPCI/CEPCI_ref * F_m_sus304_TP;
% (c) Impurity removal
FOB_IMRE = REFCOST_reactor * (Volume_IMRE/REFvolume_reactor)^NTH_reactor * CEPCI/CEPCI_ref * F_m_sus304_reactor;
% (d) Filter2
Area_filter2 = 12.07 * Scale_factor;
FOB_filter2 = REFCOST_filter * (Area_filter2/REFarea_filter) * CEPCI/CEPCI_ref2;
EQ_IMRE = FOB_TP3 + FOB_TP4 + FOB_IMRE + FOB_filter2;

% [4] Carbonation Process
% (a) CPUMP5
TP5_inlet_frate = TP5_unit_frate * Convert_frate;
FOB_TP5 = REFCOST_TP * (TP5_inlet_frate/REFfrate_TP)^0.59 * CEPCI/CEPCI_ref * F_m_sus304_TP;
% (b) Impurity removal
if Volume_CR > 100
    num_CR = ceil(Volume_CR/100);
else
    num_CR = 1;
end
FOB_CR = REFCOST_reactor * ((Volume_CR/num_CR)/REFvolume_reactor)^NTH_reactor * CEPCI/CEPCI_ref * F_m_sus304_reactor * num_CR;
% (c) Blower
Blower_inlet_frate = V_flue * Convert_frate;
FOB_Blower =  REFCOST_Blower * (Blower_inlet_frate/REFfrate_Blower) * CEPCI/CEPCI_ref2;
% (d) COOL2
Area_COOL2 = 10.50 * Scale_factor;
FOB_COOL2 = REFCOST_HEAT * (Area_COOL2/REFarea_HEAT)^NTH_HTX * CEPCI/CEPCI_ref2;
% (e) DRUM1
Scale_factor_drum = Scale_factor^(1/2.5);
Height_Drum1_ref = 5.334;       % m, from aspen plus
Diameter_Drum1_ref = 1.8288;    % m, from aspen plus
H_Drum1 = Height_Drum1_ref * Scale_factor_drum;
D_Drum1 = Diameter_Drum1_ref * Scale_factor_drum;
FOB_Drum1 = REFCOST_Drum * ((H_Drum1*D_Drum1^1.5)/REFfactor_Drum)^NTH_Drum * CEPCI/CEPCI_ref;
% (f) Centrifuge
Diameter_Centrifuge = REFdiameter_Centrifuge * Scale_factor;
FOB_Centrifuge = REFCOST_Centrifuge * (Diameter_Centrifuge/REFdiameter_Centrifuge) * REFNum_Centrifuge * CEPCI/CEPCI_ref2;
% (g) Dryer
Area_Dryer = REFArea_Dryer * Scale_factor;
FOB_Dryer = REFCOST_Dryer * (Area_Dryer/REFArea_Dryer)^NTH_Dryer * CEPCI/CEPCI_ref2;
% (h) Drum2
Height_Drum2_ref = 3.6576;     % m, from aspen plus
Diameter_Drum2_ref = 0.9144;   % m, from aspen plus
H_Drum2 = Height_Drum2_ref * Scale_factor_drum;
D_Drum2 = Diameter_Drum2_ref * Scale_factor_drum;
FOB_Drum2 = REFCOST_Drum * ((H_Drum2*D_Drum2^1.5)/REFfactor_Drum)^NTH_Drum * CEPCI/CEPCI_ref;
EQ_CR = FOB_Drum2 + FOB_Dryer + FOB_Centrifuge + FOB_Drum1 + FOB_COOL2 + FOB_TP5 + FOB_CR + FOB_Blower;

% [5] Crushing Process
FOB_Crusher = REFCOST_Crusher * (Elec_crusher /( Convert_energy*REFpower_Crusher))^NTH_Crusher * CEPCI/CEPCI_ref;
EQ_CS = FOB_Crusher;

% [6] Recycling Process
% (a) HEATEX2
Area_HEATEX2 = REFarea_HEATEX2 * Scale_factor;
FOB_HEATEX2 = REFCOST_HEATEX2 * (Area_HEATEX2/REFarea_HEATEX2)^NTH_HTX * CEPCI/CEPCI_ref2;
% (b) HEAT2
Area_HEAT2 = 4.85 * Scale_factor;
FOB_HEAT2 = REFCOST_HEAT * (Area_HEAT2/REFarea_HEAT)^NTH_HTX * CEPCI/CEPCI_ref2;
% (c) Drum3
Height_Drum3_ref = 5.4864;     % m, from aspen plus
Diameter_Drum3_ref = 1.8288;   % m, from aspen plus
H_Drum3 = Height_Drum3_ref * Scale_factor_drum;
D_Drum3 = Diameter_Drum3_ref * Scale_factor_drum;
FOB_Drum3 = REFCOST_Drum * ((H_Drum3*D_Drum3^1.5)/REFfactor_Drum)^NTH_Drum * CEPCI/CEPCI_ref;
% (d) COND1
Area_COND1 = 2.63 * Scale_factor;
FOB_COND1 = REFCOST_COND * (Area_COND1/REFarea_COND)^NTH_HTX * CEPCI/CEPCI_ref2;
% (e) HEAT3
Area_HEAT3 = 472.87 * Scale_factor;
FOB_HEAT3 = REFCOST_HEAT * (Area_HEAT3/REFarea_HEAT)^NTH_HTX * CEPCI/CEPCI_ref2;
% (f) Drum4
Height_Drum4_ref = 3.6576;     % m, from aspen plus
Diameter_Drum4_ref = 4.1148;   % m, from aspen plus
H_Drum4 = Height_Drum4_ref * Scale_factor_drum;
D_Drum4 = Diameter_Drum4_ref * Scale_factor_drum;
FOB_Drum4 = REFCOST_Drum * ((H_Drum4*D_Drum4^1.5)/REFfactor_Drum)^NTH_Drum * CEPCI/CEPCI_ref;
% (g) COOL3
Area_COOL3 = 2942.07 * Scale_factor;
FOB_COOL3 = REFCOST_HEAT * (Area_COOL3/REFarea_HEAT)^NTH_HTX * CEPCI/CEPCI_ref2;
% (h) CPUMP6
TP6_inlet_frate = TP6_unit_frate * Convert_frate;
FOB_TP6 = REFCOST_TP * (TP6_inlet_frate/REFfrate_TP)^0.59 * CEPCI/CEPCI_ref;
% (i) CPUMP7
TP7_inlet_frate = TP7_unit_frate * Convert_frate;
FOB_TP7 = REFCOST_TP * (TP7_inlet_frate/REFfrate_TP)^0.59 * CEPCI/CEPCI_ref;
EQ_EVP = FOB_TP6 + FOB_TP7 + FOB_Drum4 + FOB_COOL3 + FOB_HEAT3 + FOB_COND1 + FOB_Drum3 + FOB_HEAT2 + FOB_HEATEX2;

%%% CAPITAL COST SUMMARY
ISBL_Total = (EQ_CR + EQ_CS + EQ_Elec + EQ_EVP + EQ_Ex + EQ_IMRE) * f_IEC;
FIXED_CAPEX = ISBL_Total*(1+OSBL_OS)*(1+OSBL_DE+OSBL_CN);
CRF = (DISCO*(DISCO+1)^LIFET)/(((DISCO+1)^LIFET)-1);
ANNUALIZED_CAPEX = CRF*FIXED_CAPEX;

%=========================== OPERATING EXPENSES ==========================%
% [1] Plant Expenses
% [2] Raw Material Costs
% (a) Steel Slag Treatment Cost
COST_steel_slag = Slag_mass_inlet * (- UT_slag_treatment) * Operating_hour; % USD/yr
% (b) Unreacted Steel Slag Treatment Costs
COST_uncreated_slag = W_waste * UT_slag_treatment * Operating_hour;         % USD/yr
% (c) NaCl Make-up Costs
COST_NaCl = W_NaCl_makeup * OP_RM_NaCl * Operating_hour;                    % USD/yr
% (d) Disposal Treatment Cost
COST_Disposal = W_impurity * UT_slag_treatment * Operating_hour;            % USD/yr

% [3] Utility Costs
% (a) Electricity
COST_elec_cell = UT_ELEC * Elec_duty_cell * Operating_hour;    % USD/yr, electrolyzer
COST_elec_others = UT_ELEC * Elec_duty_other * Operating_hour; % USD/yr, other
% (b) LP Steam
COST_LP_steam = UT_LP_STEAM * LP_duty * Operating_hour;        % USD/yr
% (c) Process Water 
W_cooling_water = - (Cool_COOL1_duty + Cool_COND1_duty) / (delta_T_min * heat_capacity); % tonne/hr
COST_cooling_water = OP_RM_H2O * W_cooling_water * Operating_hour;                       % USD/yr

% [4] Fixed OPEX
%COST_FO_Labor = REFCOST_Labor * W_nCaCO3 * Operating_hour / 200000;
%COST_FO_Maintenance = ISBL_Total * f_Maintenance;
%COST_FO_admin = COST_FO_Labor * f_admin;
%COST_FO_overhead = COST_FO_Labor * f_overhead;
%COST_FO_Laboratory = ISBL_Total * f_Laboratory;
% OPERATING COST SUMMARY
%TOTAL_OPEX_FO = COST_FO_Labor + COST_FO_Maintenance + COST_FO_admin + COST_FO_overhead + COST_FO_Laboratory;

TOTAL_OPEX_Variable = COST_cooling_water + COST_LP_steam + COST_elec_others + COST_elec_cell + COST_Disposal + COST_NaCl + COST_uncreated_slag + COST_steel_slag;

%============================== COST SUMMARY =============================%
ANNUAL_PROD_COST = ANNUALIZED_CAPEX + TOTAL_OPEX_Variable; % USD/yr
COGM = ANNUAL_PROD_COST/(W_nCaCO3 * Operating_hour);       % $/ton




%% CO2 LIFE CYCLE ASSESSMENT MODEL
% [1] Direct Plant Emissions of GHG

% [2] Indirect Emissions from Energy Consumption
% (a) Indirect Emissions from Energy Consumption
EM_Elec_cell = Unit_elec_duty_cell * GWI_ELEC_GRID * 0.001;
EM_Elec_others = Unit_elec_duty_others * GWI_ELEC_GRID * 0.001;
EM_LP = Unit_LP_duty * GWI_STEAM_LP * 0.001;
% (b) Indirect Emissions from Utility Consumption
EM_Process_Water = GWI_WATER * W_cooling_water;
% (c) Indirect Emissions from Raw Material Production
Unit_disposal = W_impurity / W_nCaCO3;
Unit_slag = Slag_mass_inlet / W_nCaCO3;
Unit_waste = W_waste/ W_nCaCO3;
Unit_flue = W_CO2 / W_nCaCO3;
EM_disposal = GWI_disposal * Unit_disposal;
EM_slag = - GWI_steel_slag * Unit_slag;
EM_waste = GWI_steel_slag * Unit_waste;
EM_flue = GWI_CO2_FG * Unit_flue;

% CO2 GWI SUMMARY (in kg CO2-eq)
SPECIFIC_GWI = EM_Elec_cell + EM_Elec_others + EM_LP + EM_Process_Water + EM_disposal + EM_slag + EM_waste + EM_flue;




%% EXPORT SUSTAINABILITY CRITERIA METRICS
% Criteria 1: Technoeconomics = Cost of Goods Manufactured
% Criteria 2: Global Warming Emissions = Specific Global Warming Impact
EvalMetrics = [COGM, SPECIFIC_GWI];  
end