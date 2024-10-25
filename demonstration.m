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

% parameters
R0 = 2;            % R_0
I0 = 10;           % Initial no. cases
maxtime = 300;     % simulation time
Rtype = 'Instantaneous';

% R_t timeseries (exceed example)
tvec = [1 75];
Rvec = [R0 0.8];

% R_t timeseries (don't exceed example)
tvec = [1 66];
Rvec = [R0 0.95];

% serial interval distribution
mu = 5.4;
sig = 1.5;
wtimes = 1:maxtime;
w = gampdf(wtimes,(mu/sig)^2,1/(mu/sig^2));
w = w./sum(w);

% run models
[out] = RENEWALmodel(Rvec,tvec,w,I0,maxtime,Rtype,para);

%% ICU PLOTTING AND COST EVALUATION

% plotting
f1 = figure(1);
f1.Position = [500 1000 1000 225*5];
hold all

subplot(5,1,1)
plot(out.t,out.R,'Color',myred)
yline(1,'k--','$R_t=1$','Interpreter','latex','LineWidth',2,'FontSize',16,'LabelVerticalAlignment','top','LabelHorizontalAlignment','right','Layer','bottom')
xlabel('Time $t$ (days)')
ylabel('$R_t$','Rotation',0);
axis([0 maxtime 0 1.2*max(out.R)])
grid on

subplot(5,1,2)
hold all
plot(out.t,out.in_ICU,'Color',myred)
yline((1 - para.eta)*para.Ibar,'k--','Capacity','Interpreter','latex','LineWidth',2,'FontSize',16,'LabelVerticalAlignment','top','LabelHorizontalAlignment','right','Layer','bottom')
xlabel('Time $t$ (days)')
ylabel('In ICU');
axis([0 maxtime 0 1.2*max(out.in_ICU)])
grid on

% calculate costs under each constraint
vs = [1 2 3 5 10].*10;
Nvs = length(vs);

% hard constraint
whichconstr = 'hard';
inc_or_prev = 2;
indirect = 1;
hardcosts_direct = compute_cost(out,para,whichconstr,vs(1),inc_or_prev,0);
hardcosts_total = compute_cost(out,para,whichconstr,vs(1),inc_or_prev,indirect);
hardcosts_indirect = hardcosts_total - hardcosts_direct;

% soft constraint
whichconstr = 'soft';
softcosts_direct = compute_cost(out,para,whichconstr,vs,inc_or_prev,0);
softcosts_total = compute_cost(out,para,whichconstr,vs,inc_or_prev,indirect);
softcosts_indirect = softcosts_total - softcosts_direct;

% plotting direct costs
subplot(5,1,3)
hold all

for v = 1:Nvs
    dispname = ['Logistic $\alpha$, $v = ', ' ', num2str(vs(v)), '$'];
    vcols(v,:)
    plot(out.t,softcosts_direct(v,:),'Color',vcols(v,:),'DisplayName',dispname)
end
plot(out.t,hardcosts_direct,'Color',myred,'DisplayName','Step $\alpha$')
xlim([0 maxtime])
ax = gca;
ax.YAxis.Exponent = 0;
xlabel('Time $t$ (days)')
ylabel('Direct Cost')
grid on

% plotting indirect costs
subplot(5,1,4)
hold all

for v = 1:Nvs
    dispname = ['Logistic $\alpha$, $v = ', ' ', num2str(vs(v)), '$'];
    plot(out.t,softcosts_indirect(v,:),'Color',vcols(v,:),'DisplayName',dispname)
end
plot(out.t,hardcosts_indirect,'Color',myred,'DisplayName','Step $\alpha$')
xlim([0 maxtime])
ax = gca;
ax.YAxis.Exponent = 0;
xlabel('Time $t$ (days)')
ylabel('Indirect Cost')
grid on

% plotting total costs
subplot(5,1,5)
hold all

for v = 1:Nvs
    dispname = ['Logistic $\alpha$, $v = ', ' ', num2str(vs(v)), '$'];
    plot(out.t,softcosts_total(v,:),'Color',vcols(v,:),'DisplayName',dispname)
end
plot(out.t,hardcosts_total,'Color',myred,'DisplayName','Step $\alpha$')
xlim([0 maxtime])
ax = gca;
ax.YAxis.Exponent = 0;
xlabel('Time $t$ (days)')
ylabel('Total Cost')
grid on

% 
% saveas(f3,['./renewal_images/F3C_renewal_costs_alpha0_',num2str(para.alpha(1)),'.png'])
% 
% % save results for uncertainty analysis
% save('./mats/results.mat',"out1","out2","vs","Nvs","C1_softcosts","C2_softcosts","inc_or_prev","indirect")
