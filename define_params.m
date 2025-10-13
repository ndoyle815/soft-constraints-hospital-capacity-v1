% script to define model and control parameters into a mat file
clear

if not(isfolder('mats'))
    mkdir('mats')
end

% if file exists, load default values (from CHESS data)
load('./mats/Parameters.mat')
% includes:
% Dist_HM: distribution from hospital to death
% Dist_SC: distribution from symptoms to ICU
% Dist_SH: distribution from hospital to death
% Time_C:  time spent in ICU
% Time_H:  time spent in hospital (no death)
% Time_HC: time spent in hospital after ICU
% ca:      probability of critical care given symptoms
% ha:      probability of hospital (no ICU) given symptoms
% ma:      probability of death given hospital
% dt:      simulation time step

% hospital capacity
eta  = 0.468;                      % background ICU occupancy rate
Ibar = 2627/(1 - eta);         % total ICU beds (no. available to treat is (1-eta)Ibar)
etaH = 0.4597;                  % background G&A occupancy rate
Hbar = 62271/(1 - etaH);  % total G&A beds (no. available to treat is (1-eta)Hbar)
Hosp_background = etaH*Hbar;
ICU_background = eta*Ibar;   
alpha = [0.3 1];                 % relative cost of cases above capacity compared to below capacity (hard capacity)

save("./mats/Parameters.mat","Dist_HM","Dist_SC","Dist_SH","Time_C","Time_H","Time_HC","ca","ha","ma","dt",...
     "eta","Ibar","etaH","Hbar","Hosp_background","ICU_background","alpha",'-mat')

% define colours and colourmaps for plotting
% colours
vcols = [65,182,196; 29,145,192; 34,94,168; 37,52,148; 8,29,88]./255;
myred = [228 26 28]./255;
myblue = [55 126 184]./255;
mygreen = [77 175 74]./255;
mypurple = [152 78 163]./255;
myorange = [255 127 0]./255;
mywhite = [1 1 1];

% colourmaps
% purple-orange for Fig 2
dcol = 128;
colorvec = (0:1/dcol:1)';
mycolormap1 = mypurple.*(1 - colorvec) + mywhite.*colorvec;
mycolormap2 = mywhite.*(1 - colorvec) + myorange.*colorvec;
POcolormap = [mycolormap1; mycolormap2(2:end,:)];

% blue-green for Figs 3-5
mycolormap1 = myblue.*(1 - colorvec) + mywhite.*colorvec;
mycolormap2 = mywhite.*(1 - colorvec) + mygreen.*colorvec;
BGcolormap = [mycolormap1; mycolormap2(2:end,:)];

save("./mats/Cols.mat","vcols","myred","myblue","mygreen","mypurple","myorange","mywhite","POcolormap","BGcolormap",'-mat')
