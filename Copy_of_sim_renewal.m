% script to simulate an outbreak using the renewel equation transmission
% model
clear; close all

% Plotting preferences
set(0,'defaultlinelinewidth',3)
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex')
set(0,'defaultTextInterpreter','latex')
set(0,'defaultaxesfontsize',16)

% load colours
load('./mats/Cols.mat')

% load distribution parameters
para = load('./mats/Parameters.mat');

stratnames = {'Running Hot (RH)','Misguided Suppression (MS)'};

% parameters
R0 = 2;            % R_0
I0 = 10.5;           % Initial no. cases
maxtime = 1000;     % simulation time
whichR = {'Instantaneous','Case'};
Rtype = whichR{1};

% R_t timeseries
if isequal(Rtype,'Case')
    ti1 = [1 60 60+45];
    ti2 = [1 58 58+45];
else
    ts = 56;
    td = 6;
    tp = 8;
    tl = 120;
    c1 = 60;
    ti1 = [1 64+0 64+0+120+c1];
    ti2 = [1 64+7 64+7+30-5];
end
Ri1 = [R0 1 0.7];
Ri2 = [R0 0.7 0.95];

% serial interval distribution
mu = 5.4;
sig = 1.5;
wtimes = 1:maxtime;
w = gampdf(wtimes,(mu/sig)^2,1/(mu/sig^2));
w = w./sum(w);

% run models
[out1] = RENEWALmodel(Ri1,ti1,w,I0,maxtime,Rtype,para);
[out2] = RENEWALmodel(Ri2,ti2,w,I0,maxtime,Rtype,para);

%% ICU PLOTTING AND COST EVALUATION

% plotting
f7 = figure(7);
f7.Position = [100 1000 450 300];

bar(wtimes,w)
xlabel('Time $t$ (days)')
ylabel('Probability $s = t$')
title('Serial interval distribution $w_s$')
axis([0 20 0 1.2*max(w)])
grid on

% saveas(f7,'../images/F3D_renewal_wt.png')

f1 = figure(1);
f1.Position = [1000 1000 1000 225];

hold all
plot(out1.t,out1.R,'Color',mygreen)
plot(out2.t,out2.R,'Color',myblue)
legend(stratnames,'Location','north','AutoUpdate','off','FontSize',16)
yline(1,'k--','$R(t)=1$','Interpreter','latex','LineWidth',2,'FontSize',16,'LabelVerticalAlignment','top','LabelHorizontalAlignment','right','Layer','bottom')
xlabel('Time $t$ (days)')
if isequal(Rtype,'Case')
    ylab = '$R(t)^c$';
else
    ylab = '$R(t)$';
end
ylabel(ylab,'Rotation',0);
axis([0 min([500 maxtime]) 0 1.2*max([out1.R; out2.R],[],"all")])
grid on

saveas(f1,'./images/F3A_renewal_Rt.png')

f2 = figure(2);
f2.Position = [1000 400 1000 225];

hold all
plot(out1.t,out1.in_ICU,'Color',mygreen)
plot(out2.t,out2.in_ICU,'Color',myblue)
legend(stratnames,'Location','east','AutoUpdate','off','FontSize',16)
yline((1 - para.eta)*para.Ibar,'k--','Capacity','Interpreter','latex','LineWidth',2,'FontSize',16,'LabelVerticalAlignment','top','LabelHorizontalAlignment','right','Layer','bottom')
xlabel('Time $t$ (days)')
ylabel('In ICU');
axis([0 min([500 maxtime]) 0 1.3*(1-para.eta)*para.Ibar])
grid on

saveas(f2,'./images/F3B_renewal_It.png')

% changing alpha around
alphachoices = [0.1 0.3456];
whichalpha = 2;
para.alpha(1) = alphachoices(whichalpha);

% calculate costs under each constraint
vs = (1:0.01:21).*10;
Nvs = length(vs);

% hard constraint
whichconstr = 'hard';
inc_or_prev = 2;
indirect = 1;
C1_hardcosts = compute_cost(out1,para,whichconstr,vs(1),inc_or_prev,indirect);
C2_hardcosts = compute_cost(out2,para,whichconstr,vs(1),inc_or_prev,indirect);

% soft constraint
whichconstr = 'soft';
tic;
C1_softcosts = sum(compute_cost(out1,para,whichconstr,vs,inc_or_prev,indirect),2)';
C2_softcosts = sum(compute_cost(out2,para,whichconstr,vs,inc_or_prev,indirect),2)';
toc;

% plotting soft constraint costs as a function of scaling v
f3 = figure(3);
f3.Position = [1000 0 1000 225];

plot(vs,C1_softcosts,'Color',mygreen)
hold on
plot(vs,C2_softcosts,'Color',myblue)
xlim([0 max(vs)])
ax = gca;
ax.YAxis.Exponent = 0;
xlabel('Scaling of logistic curve $v$')
ylabel('Cost');
legend(stratnames,'Location','northeast','FontSize',16)
% title(['$\alpha_0 = ', num2str(para.alpha(1)), ', \alpha_1 = 1$'])
cmax = 18e3;
ylim([6e3 1.05*max([C1_softcosts(1) cmax])])
yticks(0e3:3e3:cmax)
grid on

saveas(f3,'./images/F3C_renewal_costs.png')

% save results for uncertainty analysis
save('./mats/results.mat',"out1","out2","vs","Nvs","C1_softcosts","C2_softcosts","inc_or_prev","indirect","stratnames",'-mat')
