% script to simulate and produce Figure 2 in manuscript
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

stratnames = {'Late Suppression (LS)','Early Suppression (ES)'};

% parameters
R0 = 2;            % R_0
I0 = 10;           % Initial no. cases
maxtime = 1000;     % simulation time
whichR = {'Instantaneous','Case'};
Rtype = whichR{1};

% R_t timeseries
ti1 = [1 60 90];
ti2 = [1 60 120 250];

Ri1 = [R0 1.2 0.7];
Ri2 = [R0 0.7 1.2 0.7];

% serial interval distribution
mu = 5.4;
sig = 1.5;
wtimes = 1:maxtime;
w = gampdf(wtimes,(mu/sig)^2,1/(mu/sig^2));
w = w./sum(w);

% run models
[out1] = RENEWALmodel(Ri1,ti1,w,I0,maxtime,Rtype,para);
[out2] = RENEWALmodel(Ri2,ti2,w,I0,maxtime,Rtype,para);

% plotting Fig 2A,B
f1 = figure(1);
f1.Position = [1000 1000 1100 225];

hold all
plot(out1.t,out1.R,'Color',mypurple)
plot(out2.t,out2.R,'Color',myorange)
legend(stratnames,'Location','north','AutoUpdate','off','FontSize',16)
yline(1,'k--','$R(t)=1$','Interpreter','latex','LineWidth',2,'FontSize',16,'LabelVerticalAlignment','top','LabelHorizontalAlignment','right','Layer','bottom')
xlabel('Time $t$ (days)')
ylabel('$R(t)$','Rotation',0);
axis([0 min([500 maxtime]) 0 1.2*max([out1.R; out2.R],[],"all")])
grid on

saveas(f1,'./images/F2A_renewal_Rt.png')

f2 = figure(2);
f2.Position = [1000 400 1100 225];

hold all
plot(out1.t,out1.in_ICU,'Color',mypurple)
plot(out2.t,out2.in_ICU,'Color',myorange)
legend(stratnames,'Location','east','AutoUpdate','off','FontSize',16)
yline((1 - para.eta)*para.Ibar,'k--','Capacity','Interpreter','latex','LineWidth',2,'FontSize',16,'LabelVerticalAlignment','top','LabelHorizontalAlignment','right','Layer','bottom')
xlabel('Time $t$ (days)')
ylabel('In ICU');
axis([0 min([500 maxtime]) 0 1.3*(1-para.eta)*para.Ibar])
grid on

[sum(out1.new_ICU) sum(out2.new_ICU)]

saveas(f2,'./images/F2B_renewal_It.png')

%% Cost Evaluation

% vectors for parameters v (scaling logistic curve) and k (relative cost of lockdown)
dv = 0.1;
dk = 0.05;
vs = 0:dv:100;
ks = 0:dk:50;
nv = length(vs);
nk = length(ks);

% store total costs for two scenarios with different constraints
NO_LOCKDOWN_hardcosts = zeros(nv,nk);
NO_LOCKDOWN_softcosts = zeros(nv,nk);
LOCKDOWN_hardcosts    = zeros(nv,nk);
LOCKDOWN_softcosts    = zeros(nv,nk);
inc_or_prev = 2;
indirect = 1;

% costs of control
NO_LOCKDOWN_cost_of_control = 0;
% additional control effort of Lockdown strategy is proportional to
% (RIT^2) over the intervention period
% RIT = reduction in transmission
RIT = 1 - Ri2(2)/Ri2(1);
duration = ti2(3) - ti2(2);
LOCKDOWN_cost_of_control = ks.*duration*RIT^2;

% costs of disease (hard constraint)
whichconstr = 'hard';

NO_LOCKDOWN_cost_of_inf_hard = sum(compute_cost(out1,para,whichconstr,vs(1),inc_or_prev,indirect),2);
NO_LOCKDOWN_hardcosts = repmat(NO_LOCKDOWN_cost_of_inf_hard + NO_LOCKDOWN_cost_of_control, nv, nk);

LOCKDOWN_cost_of_inf_hard = sum(compute_cost(out2,para,whichconstr,vs(1),inc_or_prev,indirect),2);
LOCKDOWN_hardcosts = repmat(LOCKDOWN_cost_of_inf_hard,1,nk) + repmat(LOCKDOWN_cost_of_control,nv,1);

% costs of disease (soft constraint)
whichconstr = 'soft';

NO_LOCKDOWN_cost_of_inf_soft = sum(compute_cost(out1,para,whichconstr,vs,inc_or_prev,indirect),2);
NO_LOCKDOWN_softcosts = repmat(NO_LOCKDOWN_cost_of_inf_soft + NO_LOCKDOWN_cost_of_control,1,nk);

LOCKDOWN_cost_of_inf_soft = sum(compute_cost(out2,para,whichconstr,vs,inc_or_prev,indirect),2);
LOCKDOWN_softcosts = repmat(LOCKDOWN_cost_of_inf_soft,1,nk) + repmat(LOCKDOWN_cost_of_control,nv,1);


% store all costs to access in loop
allcosts = zeros(nv,nk,4);
allcosts(:,:,1) = NO_LOCKDOWN_hardcosts;
allcosts(:,:,2) = LOCKDOWN_hardcosts;
allcosts(:,:,3) = NO_LOCKDOWN_softcosts;
allcosts(:,:,4) = LOCKDOWN_softcosts;

%% Plotting Figure 2C using imagesc

f4 = figure(4);
f4.Position = [200 200 550 450];
colormap(POcolormap)

imagesc(ks,vs,allcosts(:,:,4)<allcosts(:,:,3))
set(gca,'OuterPosition',[0.01 0.01 0.9 0.94])
set(gca,'YDir','Normal')
clim([0 1])
% for ticks, ticklabels, grid
vspec = vs(1:round(nv/10):nv); 
kspec = ks(1:round(nk/10):nk); 
yticks(vspec)
xticks(kspec)
yticklabels(vspec)
xticklabels(kspec)
xtickangle(0)
xlabel('$k$ scaling cost of Strategy ES')
ylabel('$v$ scaling logistic constraint')

% mark k threshold for step constraint in logistic constraint plot
kthresh = find(LOCKDOWN_hardcosts(100,:) > NO_LOCKDOWN_hardcosts(100,:), 1, 'first');
xline(ks(kthresh),'-k','LineWidth',2,'DisplayName',['$k_c = ',num2str(ks(kthresh)),'$'])

% add gridlines
yline(vspec(2:end-1), 'k--', 'Alpha', 0.3)
xline(kspec(2:end-1), 'k--', 'Alpha', 0.3)

h = axes(f4,'visible','off'); 
h.Title.Visible = 'on';
h.XLabel.Visible = 'on';
h.YLabel.Visible = 'on';

tl = title('Optimal, Soft Constraint');
tl.Position(2) = 0.97;
tl.Position(1) = 0.45;

% save figure
saveas(f4,strcat('./images/F2C_lockdowndilemma_simplemodel.png'))


% save incidence and prevalence for uncertain Hc script
save('./mats/lockdowndilemma.mat',"out1","out2","dv","dk","vs","ks","nv","nk","RIT","duration","inc_or_prev","indirect","stratnames",'-mat')
