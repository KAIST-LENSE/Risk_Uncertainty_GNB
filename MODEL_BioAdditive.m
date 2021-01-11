%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Process + Evaluation Model for Algal Biomass-based Biofillers %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function EvalMetrics = MODEL_BioAdditive(kp)             
%% Declare Global Variables
format long g
global K_CL K_1 K_2 K_E kLa F_CUL CO2e H_ion Mu_max V_col Yc A_r I_in C1 C2
%% KEY PARAMETERS WITH UNCERTAINTY
% [1] Cultivation Process
avg_GPC = kp(1);       % Average dry weight, 0.0224 grams/10^9 cells
Mu_max = kp(2);        % Maximum Specific Growth Rate, 1.07 hr^-1
kLa = kp(3);           % Volumetric mass transfer coeff, 1.4 1/hr
P_CO2 = kp(4);         % Partial Pressure of CO2 from FG feed, 0.1358 atm
Daylight = kp(5);      % Daily fraction available for photosynthesis
I_in = kp(6);          % Average Korea incident light intensity, 750 mu*E/m^2*s

% [2] Harvest+Dewater Process
rec_CENTR = kp(7);    % Biomass recovery during centrifuge harvest, 95%
conc_CENTR = kp(8);   % Exit conc. of biomass @ Centrifuge, 25 dwt%

% [3] Blowdown + Media Recycle
% [4] Convective Flue Gas Drying

% [TEA] Techno-Economic Analysis Parameters
%%% [EQUIPMENT] Equipment Sizing Costs
C_PBR = kp(9);                  % 222,657 USD/Acre (calculated from KU Experimental Specs)
C_FG_Blower = kp(10);           % 7.5kW 3-lobe Blower, [15], 5803 USD/unit
C_Stock_Pump = kp(11);          % 5.5kW Centrifugal Pump, [15], 1826 USD/unit
C_GM_Pump = kp(12);             % 10kW Submersible Pump, [15], 20234 USD/unit
C_SEP = kp(13);                 % 7.5kW Centrifugal Separator, 62591 USD/unit
C_SS_T = kp(14);                % Cost per Tank, Carbon Steel, 54159 USD/unit
C_GM_T = kp(15);                % Cost per Tank, Carbon Steel, 478113 USD/unit
%%% [OPERATING] Plant Costs
OP_PE_MAINT = kp(16);      % $/Acre from [21]. Land maintainence + Clean-in-Place costs, 340+750
OP_PE_PBR = kp(17);        % $/Acre calculated from KU Experimental Specs
%%% [OPERATING] Raw Material Costs
OP_RM_NSOURCE = kp(18);    % $/ton of N-Nutrients (NH4Cl, [21]), 180
OP_RM_PSOURCE = kp(19);    % $/ton of P-Nutrients (KH2PO4, K2HPO4, [21]), 1400
OP_RM_SSOURCE = kp(20);    % $/ton, Industry Quote Metal Sulfate fertilizers, 150
OP_RM_CaSOURCE = kp(21);   % $/ton, Industry grade CaCl2-2H2O, 200
OP_RM_NaSOURCE = kp(22);   % $/ton, Agriculture grade EDTA-Na2, 1800
OP_RM_H2O = kp(23);                   % $/ton, 0.068 from IPON Quote
%%% [OPERATING] Utility Costs
UT_ELEC = kp(24);                     % $/GJ, 25 from IPON Quote

% [LCA] CO2 Life Cycle Assessment Parameters
%%% [MATERIAL] Global Warming Potentials for Material Consumption
GWI_NSOURCE = kp(25);                 % kg CO2-eq/ton N-Nutrients, 1860
GWI_PSOURCE = kp(26);                 % kg CO2-eq/ton P-Nutrients, 6750
GWI_SSOURCE = kp(27);                 % kg CO2-eq/ton S-Nutrients, 300
GWI_CaSOURCE = kp(28);                % kg CO2-eq/ton Ca-Nutrients, 89
GWI_NaSOURCE = kp(29);                % kg CO2-eq/ton Na-Nutrients, 53
GWI_WATER = kp(30);                   % kg CO2-eq/ton Water Consum., 0.03
%%% [UTILITY] Global Warming Potentials for Utility Consumption
GWI_ELEC_GRID = kp(31);               % 108.88 kg CO2-eq/GJ Electricity




%% NON VARYING FIXED-PARAMETERS & DESIGN SPECS
% [1] Cultivation Process
pH = 7;                % Alkaline pH control setpoint for optimal growth
V_col = 1000;          % Design Spec Vol. of Cult Bag 
K_CL = 0.0035;         % Half-TIC sat. constant, mmol/10^9 cells
K_1 = 10^-6.35;        % Dissociation constant of CO2, mol/L
K_2 = 10^-10.3;        % Dissociation constant of Bicarbonate, mol/L
K_E = 0.08;            % Light half sat. constant, mu*E/(s-10^ Cells)
A_r = 10;              % Illuminated area per unit PBR, m^2
C1 = 0.493;            % Reactor geometry constant
C2 = -0.925;           % Light wavelength constant
Yc = 1211;             % Conversion Yield of Carbon, in 10^9 cells
F_CUL = 0.0;           % Medium cycle flowrate (REF=0.12L/hr,KCRC=0.06L/hr)
VVM = 0.1;             % Flowrate of flue gas in VVM (from Korea Univ.)
RTP_FG = 0.0264385;    % RT/P of Flue Gas at 45 deg C
H = 29.41;             % Henry's Constant, atmL/mol
X_inoc = 4.858;        % Average innoculum concentration in 10^9 cells/L
TIC_0 = 0.0046;        % Initial TIC concentration in media, moles/L
HAR_T = 72;            % Semi-continuous harvest time in hours, 96 hrs 
STOCK_T = 0.5;         % Time to transfer stock soln. to GM Tank, hrs
GM_T = 1.0;            % Time to transfer growth media to Cult., hrs
REC_T = 1.0;           % Time to transfer cult. broth back to GM Tank, hrs
V_Cult_T = 1000000000; % Total cultivation volume, liters (1E6 m^3)
V_Cult_Unit = 1;       % m^3 per unit Cultivation Module
Land_Cult_Unit = 5;    % 3.3 m^2 per 1000L PBR + 1.7 m^2 Settling tank 
Har_Frac = 0.9;        % Fraction of cultivation broth harvested
FG_MM = [44.01, 18.02, 34.1, 100.09, 46.01, 32, 28.01]; 
FG_Moles = [P_CO2, 0.0818, 0.0005, 0.0154, 0.0025, 0.0354, 0.7286];

% [2] Harvest+Dewater Process
SEP_T = 10;            % Operation time/batch for culture separation, hrs

% [3] Blowdown + Media Recycle
BD_LOSS = [0.01, 0.01, 0, 0, 0, 0, 0];         % Blowdown Loss Fractions

% [4] Convective Flue Gas Dryer
Target_Moisture = 0.1;                   % 10% moisture content [A. Giostri paper]
FG_Comp = [0.0818, 0.7465, 0.1358, 0.0354, 0.0005]; % [H2O, N2, CO2, O2, H2S]
FG_MW = [18.02, 28.014, 44.01, 32, 34.1];           % [H2O, N2, CO2, O2, H2S]
FG_Temp_in = 128;                        % degC [A. Giostri paper]
FG_Temp_out = 55;                        % degC [A. Giostri paper]
FG_Pressure = 825.068;                   % mmHg
Antoine_Coeff_Water = [8.14, 1810.94, 244.485];     % Water from 99 - 374 degC
CP_FG = 1.054;                           % kJ/kgC
CP_H2O = 4.187;                          % kJ/kgC
H_lat = 2260;                            % kJ/kg
DRY_T_m = 242.7;                         % y = mx+b where y = Drying time (mins) and x = moisture content (%) [H. Hosseinizand et al., 2017]
DRY_T_b = 83.015;                        % y = mx+b where y = Drying time (mins) and x = moisture content (%) [H. Hosseinizand et al., 2017]
AREA_LOAD = 30;                          % kg/m2. Kiranoudia and Markatos (2000) suggested a maximum unit loading of 50 kg/m2 (wet)

% [TEA] Techno-Economic Analysis Parameters
%%% Chemical Engineering Plant Cost Indices (CEPCI)
CEPCI = 619.2;                        % Based on CEPCI for 2019
CEPCI_PBR = 619.2;                    % 2019 CEPCI based on [21] publish yr
CEPCI_PumpsBlowers = 541.3;           % 2016 CEPCI based on [15] publish yr
CEPCI_Separator = 541.3;              % 2016 CEPCI based on [15] publish yr
CEPCI_Upstr_Tanks = 603.1;            % 2018 CEPCI from [22]
%%% Equipment Efficiencies
NU_Pump = 0.73;                       % Efficiency for Water & Solutions
NU_Motor = 0.90;                      % Motor efficiency for 50-250kW Range
NU_HydFan = 0.70;                     % Hydraulic fan efficiency for blower, [A. Giostri et al., 2016]
NU_MechFan = 0.94;                    % Mechanical-electric fan efficiency for blower, [A. Giostri et al., 2016] 
%%% Nth-Plant Assumption Exponential Constants
NTH_PMP = 0.6;                        % From [23]
NTH_HEX = 1.2;                        % From [23], Shell and Tube
NTH_CSTR = 0.8;                       % From [23], Jacketed CSTR
NTH_DIST = 0.85;                      % From [23], Vertical Pressurized
%%% Capital Cost Lang Factors
OSBL_OS = 0.4;                        % Offsite Operations, from [23]
OSBL_DE = 0.25;                       % Design and Engineering, from [23]
OSBL_CN = 0.1;                        % Contingency, from [23]
LIFET = 30;                           % Plant Lifetime, years
DISCO = 0.07;                         % Annual Discount Rate % of TCI
%%% Equipment Costs
CAP_SS_T = 500;                       % m^3 Storage Capacity for SS Tank
CAP_GM_T = 10000;                     % m^3 Storage Capacity for GM Tank
C_Rec_Pump = C_Stock_Pump;            % 5.5kW Centrifugal Pump, [15], 1826 USD/unit
%%% [OPERATING COST] Raw Material Costs
OP_RM_FG = 0;                         % $/kg of Flue Gas

% [LCA] CO2 Life Cycle Assessment Parameters
%%% [BOUNDARY] Binary Parameters to Set LCA Boundary for Plant
B_DIR = 1;    % Include direct plant emissions? Default=1
B_ENR = 1;    % '' indirect effects from energy consumption? Default=1
B_MAT = 1;    % '' indirect effects from RawMat. consumption? Default=1
B_CaS = 0;    % '' indirect effects from plant construction & salvage?
B_PC = 0;     % '' indirect effects from product consumption?
%%% [MATERIAL] Global Warming Potentials for Material Consumption
GWI_PBR_LDPE = 2.06;                  % kg CO2-eq/kg for LDPE used for PBR
LDPE_per_Module = 1.5;                % kg LDPE per Cultivation Module (3kg per 2 years)
%%% [UTILITY] Global Warming Potentials for Utility Consumption
%%% [CONSUMPTION] Product Biocrude Consumption
GWI_Biofiller_Cons = 0;               % For biodegradable plastic additive




%% DERIVED PARAMETERS AND DESIGN SPECS
% [1] Cultivation Process
CO2e = P_CO2/H;         % CO2 conc. (moles/L) that is in eq. with gas phase
H_ion = 10^-pH;         % Concentration of Hydrogen Ions in medium
N_col = V_Cult_T/V_col; % Number of Pl104astic Columns in Cultivation Plant
% [2] Harvest+Dewater Process
% [3] Blowdown Recycle Treatment
% [4] Convective Flue Gas Drying
% [TEA] Techno-Economic Analysis Parameters
OP_Eff = HAR_T/(HAR_T+STOCK_T+SEP_T+max([GM_T,REC_T])); % Operating eff, %
Num_Har = floor((OP_Eff*8760)/HAR_T);                   % # of Harvests/yr
% [LCA] CO2 Life Cycle Assessment Parameters




%% PROCESS MASS BALANCE MODEL
% [1] Cultivation Process (Kinetic)
% (a) Solve Cultivation Mass Balance Model
NUM_INT = HAR_T*2;                          % Number of Timespan Intervals
timespan = linspace(0, HAR_T*Daylight, NUM_INT);
cult_initial = [X_inoc; TIC_0];             % Initial [TIC] is 0.0046 mol/L
[~,cult_out] = ode45(@MODEL_Kinetics_Cultivation, timespan, cult_initial);
% (b) Calculate Biomass Production
%%%% plot(timespan, cult_out(:,1))
%%%% xlabel('Cultivation Time, hours')
%%%% ylabel('PBR Cell Concentration, 10^9 cells/L')
Har_Conc = cult_out(NUM_INT,1)*avg_GPC;           % g/L conc. @ harvest
Inn_Conc = X_inoc*avg_GPC;                        % g/L conc. @ innocul.
% (c) Calculate Produced Biomass
prod_bm = (Har_Conc - Inn_Conc)*V_Cult_T/1000;    % kg biomass @ harvest
% (d) Complete Cultivation Mass Balance
%%%Initial Nutrient loading of [N, P, S, Ca, EDTA-Na], (kg/L)
nutri_conc = [0.00035 0.000432 0.00001 0.000005, 0.000005];                
%%%Specific (kg/kg) Biomass Consumption of [Water, CO2, N, P, S, Ca, EDTA-Na]
%%% Used composition data from Tokusoglu and Unal 
% CaCl2-2H2O = 110.98 + 36.03 = 147.01
% Calcium Molar mass = 40.08
% Ratio: 3.6679
% EDTA-Na2 = 292.24 + 45.98 = 338.218
% Na2 Molar Mass = 45.98
% Ratio: 7.36
cult_stoic = [0.829, 2.037, 0.238, 0.009, 0.0065, 0.00216, 0.001];   
cons_cult = cult_stoic.*prod_bm;                                  
% (f) Exit Flow Composition: [Water, Biomass, N, P, S], kg/hr 
exit_cult = [];
exit_cult(1) = (V_Cult_T - cons_cult(1))*Har_Frac/HAR_T;
exit_cult(2) = prod_bm*Har_Frac/HAR_T;
exit_cult(3:7) = ((V_Cult_T*nutri_conc)-cons_cult(3:7)).*(Har_Frac/HAR_T);
New_Innoc_Conc = cult_out(NUM_INT,1)*(1-Har_Frac);

% [2] Harvest+Dewater Process  (Black Box)
% (a) Centrifugation
bm_CENTR = rec_CENTR*exit_cult(2);             % Recovered biomass, kg/hr
water_CENTR = (bm_CENTR/conc_CENTR)-bm_CENTR;  % Exit water flowrate, kg/hr
exit_CENTR = [water_CENTR, bm_CENTR, water_CENTR*nutri_conc]; 
% (b) Recycle Culture to Blowdown
recycle_CENTR = exit_cult - exit_CENTR;        % kg/hr

% [3] Blowdown + Media Recycle (Black Box)
% (a) Exit Streams
blowdown_loss = recycle_CENTR.*BD_LOSS;                   % kg/hr
cult_recycle = recycle_CENTR - blowdown_loss;             % kg/hr
wet_bm = exit_CENTR(1:2);                                 % kg/hr 
% (b) Make-Up Stream, kg
Cult_Initial = nutri_conc.*V_Cult_T;
Cult_Consumption = cult_stoic(3:7).*prod_bm;
Cult_Final = Cult_Initial - Cult_Consumption;
Cult_Leftover = Cult_Final.*(1-Har_Frac);
Cult_Recycle = (Cult_Final.*Har_Frac) - blowdown_loss(3:7);
MU_Nutri = Cult_Initial - Cult_Leftover - Cult_Recycle;
MU_Water = cons_cult(1) + (BD_LOSS(1)*V_Cult_T*Har_Frac);

% [4] Convective Flue Gas Dryer (Black Box)
%%% Assume an adiabatic dryer where the Biomass is dried from direct
%%% convective contact with the flue gas and not via evaporation from heat
%%% conduction. The direct contact of flue gas to the Biomass does transfer
%%% heat, but it is still adiabatic because the evaporation is released with
%%% the flue gas at the outlet.
% (a) Water Evaporation Rate
Exit_Water_Flow = Target_Moisture * (wet_bm(2)/(1 - Target_Moisture));% kg/hr
Evaporation_Rate = wet_bm(1) - Exit_Water_Flow;                       % kg/hr
% (b) Humidity of Inlet
FG_Comp_Mass = (FG_MW.*FG_Comp)/sum(FG_MW.*FG_Comp);    
Mois_FG_in = FG_Comp_Mass(1);
Psat_FG_in = 10^(Antoine_Coeff_Water(1) - (Antoine_Coeff_Water(2)/(Antoine_Coeff_Water(3) + FG_Temp_in)));
Hum_in = (FG_MW(1)/mean(FG_MW.*FG_Comp)) * (Mois_FG_in*Psat_FG_in/(FG_Pressure - (Mois_FG_in*Psat_FG_in)));
% (c) Inlet FG Enthalpy (assume adiabatic process Hf_in = Hf_out)
Hf_in = ((CP_FG + CP_H2O*Hum_in)*(FG_Temp_in+273)) + (H_lat*Hum_in);  % kJ/kg, equals Hf_out
% (d) Calculate Exit Humidity from Adiabatic Assumption
Hum_out = (Hf_in - CP_FG*FG_Temp_out)/((CP_H2O*FG_Temp_out) + H_lat);
% (e) Calculate Flue Gas flowrate from Humidity Difference
FG_Flow = Evaporation_Rate/(Hum_out - Hum_in);                        % kg/hr flue gas flowrate
Prod_BM = wet_bm - [Evaporation_Rate, 0];                             % kg/hr

% ORGANIZE PRODUCT STREAMS
PROD_BIOFILLER = sum(Prod_BM);   % kg/hr
% ORGANIZE RECYCLE STREAMS
REC_WATER = cult_recycle(1);     % kg/hr
REC_NUTRI = cult_recycle(3:7);   % kg/hr




%% PROCESS ENERGY BALANCE MODEL
% [1] ENERGY REQUIREMENTS: Cultivation Process
% (a) Flue Gas Centrifugal Blower Energy (DAYTIME) (in MW) [15]
RATING_FG = 7.5;                          % kW per Blower, 0.12 bar
FG_UNIT_FWRATE = 520;                     % Nm^3/hr capacity
FWRATE_FG = VVM*V_Cult_T*0.06;            % m3/hr Flue Gas Feedrate
NUM_FG = ceil(FWRATE_FG/FG_UNIT_FWRATE);  % Number of Blower Modules
PWR_FG = NUM_FG*RATING_FG;                % kW Power Consumption
ENR_FG = PWR_FG*Daylight*HAR_T*Num_Har*0.0036;  % Daylight FG Bubbling, GJ
% (b) Air Centrifugal Blower Energy (NIGHTTIME) (in MW) [15]
PWR_Air = PWR_FG*(31/73);                 % kW Power Consumption
ENR_Air = PWR_Air*Daylight*HAR_T*Num_Har*0.0036;% Nightime FG Bubbling, GJ
% (c) Stock Solution Pump Energy (in MW) [15]
RATING_Stock = 5.5;                       % kW per Pump, 2.7 bar diffP
FWRATE_Stock = 40;                        % m^3/hr flowrate per pump
NUM_Stock = ceil((MU_Water/(1000*HAR_T))/FWRATE_Stock/STOCK_T);  
PWR_Stock = (NUM_Stock*RATING_Stock)/(NU_Pump*NU_Motor); % kW Power Cons.
ENR_Stock = PWR_Stock*STOCK_T*Num_Har*0.0036; % Stock Transf. Energy, GJ
% (d) Growth Media Pump Energy (in MW) [15]
RATING_GM = 10;                           % kW per Pump, 2.7 bar diffP
FWRATE_GM = 300;                          % m^3/hr flowrate per pump
NUM_GM = ceil((V_Cult_T*Har_Frac/1000)/FWRATE_GM/GM_T);      
PWR_GM = (NUM_GM*RATING_GM)/(NU_Pump*NU_Motor); % kW Power Cons.
ENR_GM = PWR_GM*GM_T*Num_Har*0.0036;      % GM Transf. Energy, GJ 

% [2] ENERGY REQUIREMENTS: Harvest+Dewater Process
% (a) Centrifugal Separator [15]
RATING_SEP = 7.5;                 % kW per Centrifugal Separator
CAPACITY_SEP = 4;                 % Culture processing capacity, m^3/hr
NUM_SEP = ceil((V_Cult_T*Har_Frac/(1000*HAR_T))/CAPACITY_SEP);
PWR_SEP = NUM_SEP*RATING_SEP;             % kW Separation Power Cons.
ENR_SEP = PWR_SEP*SEP_T*Num_Har*0.0036;   % Centrif. separation energy, GJ

% [3] ENERGY REQUIREMENTS: Blowdown Recycle Treatment
% (a) Recycle Growth Media Pump (w/ Filter) [15]
RATING_Rec = 5.5;                 % kW per Recycle Pump
FWRATE_Rec = 40;                  % m^3/hr flowrate processing capacity
NUM_Rec = ceil((recycle_CENTR(1)/1000)/FWRATE_Rec); 
PWR_Rec = (NUM_Rec*RATING_Rec)/(NU_Pump*NU_Motor); % kW Power Consumption
ENR_Rec = PWR_Rec*REC_T*Num_Har*0.0036;   % Broth Recycle pump energy, GJ

% [4] ENERGY REQUIREMENTS: Convective Flue Gas Dryer
% (a) Belt Dryer Conveyor Mechanical Belt
RATING_MechBelt = 90;                     % kJ/kg-H2O Evap [A. Giostri et al, 2016] 
FWRATE_MechBelt = Evaporation_Rate/3600;  % kg-H2O/s
PWR_MechBelt = RATING_MechBelt * FWRATE_MechBelt;  % kW mech belt power consumption
DRY_T = DRY_T_m*(1-conc_CENTR) + DRY_T_b; % Drying time in minutes
ENR_MechBelt = PWR_MechBelt*DRY_T*60/1000000;      % Belt mechanical energy, GJ
% (b) Flue Gas Centrifugal Blower 
P_DROP_FG = 0.5;                          % kPa FG pressure drop to drive blower
DENSITY_FG = 1.224;                       % kg/m^3 
PWR_FGBlower = ((FG_Flow*P_DROP_FG)/(DENSITY_FG*NU_HydFan*NU_MechFan))/3600;  % kW blower duty
ENR_FGBlower = PWR_FGBlower*DRY_T*60/1000000;

% ORGANIZE ENERGY DUTIES FOR UTILITY CALCULATIONS
DUTY_Elec = ENR_FG+ENR_Air+ENR_Stock+ENR_GM+ENR_SEP+ENR_Rec+ENR_MechBelt+ENR_FGBlower;




%% TECHNOECONOMIC EVALUATION [TEA] MODEL
%============================ CAPITAL EXPENSES ===========================%
% [1] Cultivation Process
% (a) Vertical Airlift Plastic Bag Cultivators [21]
%%%One unit consists of 24 panels each with volume 48m*0.7m*0.045m at fill
%%%Based on Algenol's low cost hanging bag airlift PBRS reported by [21]
%%%Each cultivation unit of 24 panels occupies 1250m^2 of land area
Num_Cult_Unit = ceil((V_Cult_T/1000)/V_Cult_Unit);    % Number of PBR Units
T_Land = Land_Cult_Unit*Num_Cult_Unit;                % Total Land Occupied
EQ_PBR_Cult = C_PBR*(T_Land/4046.86)*(CEPCI/CEPCI_PBR);
% (b) Flue Gas Centrifugal Blowers
EQ_FG_Blower = C_FG_Blower*NUM_FG*(CEPCI/CEPCI_PumpsBlowers);    
% (c) Stock Solution Centrifugal Pumps
EQ_Stock_Pump = C_Stock_Pump*NUM_Stock*(CEPCI/CEPCI_PumpsBlowers);
% (d) Growth Media Centrifugal Pumps
EQ_GM_Pump = C_GM_Pump*NUM_GM*(CEPCI/CEPCI_PumpsBlowers);       

% [2] Harvest+Dewater Process
% (a) Centrifugal Heavy Duty Separator 
EQ_SEP = C_SEP*NUM_SEP*(CEPCI/CEPCI_Separator)/SEP_T;

% [3] Blowdown + Media Recycle
% (a) Recycle Growth Media Centrifugal Pump
EQ_Rec_Pump = C_Rec_Pump*NUM_Rec*(CEPCI/CEPCI_PumpsBlowers);
% (b) Stock Solution Holding Tanks (Small Field Erected)
NUM_SS_Tank = ceil((MU_Water/1000)/CAP_SS_T);
EQ_SS_Tank = C_SS_T*NUM_SS_Tank*(CEPCI/CEPCI_Upstr_Tanks);
% (c) Growth Media Holding Tanks (Large Field Erected)
NUM_GM_Tank = ceil((V_Cult_T*Har_Frac/1000)/CAP_GM_T);
EQ_GM_Tank = C_GM_T*NUM_GM_Tank*(CEPCI/CEPCI_Upstr_Tanks);

% [4] Convective Flue Gas Dryer
% (a) Cross Sectional Area Calculation
AREA_DRY = sum(wet_bm)*(1+(wet_bm(2)/sum(wet_bm)))*(DRY_T/60)/AREA_LOAD;  %m^2 area
EQ_DRYER = 2700*AREA_DRY;          % From H. Li et al.

% CAPITAL COST SUMMARY
ISBL_Cult = EQ_PBR_Cult+EQ_FG_Blower+EQ_Stock_Pump+EQ_GM_Pump;
ISBL_HarDew = EQ_SEP;
ISBL_Media = EQ_Rec_Pump+EQ_SS_Tank+EQ_GM_Tank;
ISBL_Drying = EQ_DRYER;
ISBL_Total = ISBL_Cult + ISBL_HarDew + ISBL_Media + ISBL_Drying;
FIXED_CAPEX = ISBL_Total*(1+OSBL_OS)*(1+OSBL_DE+OSBL_CN);
CRF = (DISCO*(DISCO+1)^LIFET)/(((DISCO+1)^LIFET)-1);
ANNUALIZED_CAPEX = CRF*FIXED_CAPEX;

%=========================== OPERATING EXPENSES ==========================%
% [1] Plant Expenses
% (a) Land, Clean-In-Place & PBR Plastic Replacement Costs
COST_PE_LCP = (OP_PE_MAINT+OP_PE_PBR)*(T_Land/4046.86);
% (b) Innoculum Costs
COST_PE_INOC = sum([OP_RM_NSOURCE, OP_RM_PSOURCE, OP_RM_SSOURCE, OP_RM_CaSOURCE, OP_RM_NaSOURCE].*(nutri_conc*V_Cult_T/1000))/LIFET;

% [2] Raw Material Costs
% (a) Flue Gas Utilization Costs
MOLES_FG = (prod_bm*cult_stoic(2)*Num_Har*(1/44.01))/P_CO2;
COST_RM_FG = OP_RM_FG*MOLES_FG*sum(FG_MM.*FG_Moles);
% (b) Nutrient Stock Make-Up Costs
MU_Nutri = MU_Nutri*Num_Har;  % kg/yr Makeup Nutrients
COST_RM_NU = sum([OP_RM_NSOURCE, OP_RM_PSOURCE, OP_RM_SSOURCE, OP_RM_CaSOURCE, OP_RM_NaSOURCE].*(MU_Nutri/1000));
% (c) Water Make-up Costs 
MU_Water = MU_Water*Num_Har;  % kg/yr Makeup Water
COST_RM_H2O = OP_RM_H2O*MU_Water;      

% [3] Utility Costs
COST_UT_ELEC = UT_ELEC*DUTY_Elec;          % Yearly Electrical Costs

%======================== FIXED OPERATING EXPENSES =======================%
% Labor
LAND_ACRE = T_Land/4046.86;
INF_2014 = 1.098;
SAL_PLANT_MANAGER = 155617*INF_2014*ceil(LAND_ACRE/5000);
SAL_PLANT_ENGINEER = 82050*INF_2014*ceil(LAND_ACRE/5000);
SAL_MAINT_SUPERVISOR = 60341*INF_2014*ceil(LAND_ACRE/5000);
SAL_MODULE_OPERATOR = 38590*INF_2014*ceil(LAND_ACRE/5000);
SAL_CLERK = 38110*INF_2014*ceil(LAND_ACRE/5000);
SAL_FIELD_EMPLOYEE = 3500*LAND_ACRE*INF_2014;
LABOR = SAL_FIELD_EMPLOYEE+SAL_PLANT_MANAGER+SAL_PLANT_ENGINEER+SAL_MAINT_SUPERVISOR+SAL_MODULE_OPERATOR+SAL_CLERK;
% Maintenance
MAINT_PBR = 0.05*ISBL_Cult;
MAINT_Else = 0.03*(ISBL_HarDew + ISBL_Media + ISBL_Drying);
MAINTAINENCE = MAINT_PBR + MAINT_Else;
% Administration + Overhead
ADMIN_OVERHEAD = 0.90*LABOR;
LABORATORY = 0.01*(ISBL_Cult + ISBL_HarDew + ISBL_Media + ISBL_Drying);
TOTAL_FIXED = LABOR + MAINTAINENCE + ADMIN_OVERHEAD + LABORATORY;

% OPERATING COST SUMMARY
TOTAL_OPEX_PE = COST_PE_LCP + COST_PE_INOC;
TOTAL_OPEX_RM = COST_RM_FG + COST_RM_NU + COST_RM_H2O;
TOTAL_OPEX_UT = COST_UT_ELEC;

%============================== COST SUMMARY =============================%
ANNUAL_PROD_COST = ANNUALIZED_CAPEX + TOTAL_OPEX_PE + TOTAL_OPEX_RM + TOTAL_OPEX_UT + TOTAL_FIXED;
% Cost of Goods Manufactured, befitting TRL 3&4
COGM = ANNUAL_PROD_COST*1000/(PROD_BIOFILLER*HAR_T*Num_Har);  % USD/ton




%% CO2 LIFE CYCLE ASSESSMENT MODEL
% [1] Direct Plant Emissions of GHG
% (a) Cultivation Off Gases **ASSUME FG IS FED AS BYPASS**
CULT_CO2_In_Moles = (VVM*V_Cult_T*0.06/RTP_FG)*(P_CO2/1000)*HAR_T*Num_Har;
CULT_CO2_In = CULT_CO2_In_Moles*FG_MM(1);      % kg CO2 in Total
CULT_CO2_Consumed = cons_cult(2)*Num_Har;
% (b) Sequestration as TIC in Media
MOLES_CO2_SEQ = cult_out(NUM_INT,2)*V_Cult_T;  % CO2 sequestered as TIC
MASS_CO2_SEQ = MOLES_CO2_SEQ*(FG_MM(1)/1000)*Num_Har; % Total Sequestered
CULT_CO2_Out = CULT_CO2_In - CULT_CO2_Consumed - MASS_CO2_SEQ;
% Sum all Direct Emissions
TOT_EM_DIR = -CULT_CO2_Consumed;     % BYPASS APPROACH

% [2] Indirect Plant Emissions of GHG
% (a) Indirect Emissions from Energy Consumption
EM_ENR_ELEC = GWI_ELEC_GRID*DUTY_Elec;
TOT_EM_ENR = EM_ENR_ELEC;
% (b) Indirect Emissions from Raw Material Consumption
EM_RM_NSOURCE = ((((nutri_conc(1)*V_Cult_T)/LIFET)+MU_Nutri(1))/1000)*GWI_NSOURCE;
EM_RM_PSOURCE = ((((nutri_conc(2)*V_Cult_T)/LIFET)+MU_Nutri(2))/1000)*GWI_PSOURCE;
EM_RM_SSOURCE = ((((nutri_conc(3)*V_Cult_T)/LIFET)+MU_Nutri(3))/1000)*GWI_SSOURCE;
EM_RM_CaSOURCE = ((((nutri_conc(4)*V_Cult_T)/LIFET)+MU_Nutri(4))/1000)*GWI_CaSOURCE;
EM_RM_NaSOURCE = ((((nutri_conc(5)*V_Cult_T)/LIFET)+MU_Nutri(5))/1000)*GWI_NaSOURCE;
EM_RM_WATER = MU_Water*GWI_WATER;
EM_PBR_REPLACE = GWI_PBR_LDPE*LDPE_per_Module*((V_Cult_T/1000)/V_Cult_Unit);
TOT_EM_RM = EM_RM_NSOURCE + EM_RM_PSOURCE + EM_RM_SSOURCE + EM_RM_CaSOURCE + EM_RM_NaSOURCE + EM_RM_WATER + EM_PBR_REPLACE;
% (c) Indirect Emissions from Plant Construction and Salvage
TOT_EM_CaS = 0;             % Not Studied
% (d) Indirect Emissions from Product Consumption
TOT_EM_PC = GWI_Biofiller_Cons*(PROD_BIOFILLER*HAR_T*Num_Har);  

% CO2 GWI SUMMARY (in kg CO2-eq)
TOTAL_GWI = B_DIR*TOT_EM_DIR + B_ENR*TOT_EM_ENR + B_MAT*TOT_EM_RM + B_CaS*TOT_EM_CaS + B_PC*TOT_EM_PC;
SPECIFIC_GWI = TOTAL_GWI/(PROD_BIOFILLER*HAR_T*Num_Har);




%% ??? Breakdown
% Raw Materials-
R_NSOURCE = MU_Nutri(1)/(PROD_BIOFILLER*HAR_T*Num_Har);
R_PSOURCE = MU_Nutri(2)/(PROD_BIOFILLER*HAR_T*Num_Har);
R_SSOURCE = MU_Nutri(3)/(PROD_BIOFILLER*HAR_T*Num_Har);
R_CaSOURCE = MU_Nutri(4)/(PROD_BIOFILLER*HAR_T*Num_Har);
R_NaSOURCE = MU_Nutri(5)/(PROD_BIOFILLER*HAR_T*Num_Har);
R_WATER = MU_Water/(PROD_BIOFILLER*HAR_T*Num_Har);
R_FLUE_GAS = (CULT_CO2_In/P_CO2)/(PROD_BIOFILLER*HAR_T*Num_Har);
% Utilities
U_ELEC = DUTY_Elec/(PROD_BIOFILLER*HAR_T*Num_Har);
U_WASTE_HEAT = FG_Flow*(DRY_T/60)*HAR_T/(PROD_BIOFILLER*HAR_T*Num_Har);
% Direct Emissions




%% EXPORT SUSTAINABILITY CRITERIA METRICS
% Criteria 1: Cost of Goods Manufactured (Technoeconomic Criteria)
% Criteria 2: CO2 Global Warming Impact (Environmental Emissions Criteria)
EvalMetrics = [COGM, SPECIFIC_GWI];      
end