% script to run model simulations
clear; close all

% Plotting preferences
set(0,'defaultlinelinewidth',3)
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex')
set(0,'defaultTextInterpreter','latex')
set(0,'defaultaxesfontsize',16)

% load colours
load('./mats/Cols.mat')

% load default parameters and coefficients
para = load('./mats/Parameters.mat');
para.eta = 0.5;
% prevalence vector
ICUmax = 5000;
ICU = (0:ICUmax/1000:ICUmax) + para.eta*para.Ibar;

% hard constraint
alpha_hard = para.alpha(1).*ones(1,length(ICU));
alpha_hard(ICU>=para.Ibar) = para.alpha(2);

% PLOT ALPHA LOGISTIC FUNCTION
vvec = [1 2 3 5 10].*1e1;

f1 = figure(1);
f1.Position = [400 1200 1100 350];
hold all
patch([0 0 para.eta para.eta], [0 2*para.alpha(2) 2*para.alpha(2) 0], 'k', 'FaceAlpha',0.4, 'EdgeColor','none', 'DisplayName','Background $\eta H_c$')
for v = 1:length(vvec)
    alpha_logistic = para.alpha(1) + (para.alpha(2)-para.alpha(1))./(1 + exp(-vvec(v).*(ICU./para.Ibar - 1)));
    dispname = ['Logistic $\alpha$, $v = ', ' ', num2str(vvec(v)), '$'];
    plot(ICU./para.Ibar,alpha_logistic,'Color',vcols(v,:),'DisplayName',dispname)
end
plot(ICU./para.Ibar,alpha_hard,'Color',myred,'DisplayName','Step $\alpha$')
axis([0 1.5 0 1.2*para.alpha(2)])
legend('Location','west','AutoUpdate','off','FontSize',16)
x1 = xline(para.Ibar./para.Ibar,'k--','Capacity $(H_c)$','Interpreter','latex','FontSize',16,'LabelOrientation','horizontal','LabelVerticalAlignment','top','LabelHorizontalAlignment','left','LineWidth',2,'Layer','bottom');
y1 = yline(para.alpha(1),'k--','$\alpha_0$','Interpreter','latex','FontSize',16,'LabelHorizontalAlignment','right','LineWidth',2,'Layer','bottom');
y2 = yline(para.alpha(2),'k--','$\alpha_1$','Interpreter','latex','FontSize',16,'LabelHorizontalAlignment','right','LineWidth',2,'Layer','bottom');
% scatter(0.925,para.alpha(1)*1.23,20,'k','filled','d')
xticks(0:0.1:1.5)
% xlabel('Total ICU Occupancy rate $(\textrm{ICU}(t) + \eta \bar{I}) / \bar{I}$')
xlabel('Total ICU Occupancy rate $(H(t) + \eta H_c) / H_c$')
title('Cost per new infections $\alpha(t)$')
grid on

saveas(f1,'./images/F1_alphat.png')


% plotting indirect cost
Cindirect = max(0, min(ICU-para.Ibar, para.eta*para.Ibar));
logfunc = 1./(1 + exp(-10.*(ICU./para.Ibar - 1)));
Cindirect1 = max(0, min((ICU-para.Ibar).*logfunc, para.eta*para.Ibar));
Cindirect2 = max(0, min(ICU-para.Ibar, para.eta*para.Ibar)).*logfunc;

f2 = figure(2);
f2.Position = [400 100 900 300];

hold all
patch([0 0 para.eta para.eta], [0 para.Ibar para.Ibar 0], 'k', 'FaceAlpha',0.4, 'EdgeColor','none', 'DisplayName','Background Occupancy $\eta H_c$')
plot(ICU./para.Ibar,Cindirect,'Color',myred)
% plot(ICU./para.Ibar,Cindirect1,'Color',myblue)
% plot(ICU./para.Ibar,Cindirect2,'Color',mygreen)
axis([0 1.5 0 1.2*para.eta*para.Ibar])
x1 = xline(para.Ibar./para.Ibar,'k--','Capacity $(H_c)$','Interpreter','latex','FontSize',16,'LabelOrientation','horizontal','LabelVerticalAlignment','top','LabelHorizontalAlignment','left','LineWidth',2,'Layer','bottom');
xticks(0:0.1:1.5)
% xlabel('Occupancy rate $(\textrm{ICU}(t) + \eta \bar{I}) / \bar{I}$')
xlabel('Total ICU Occupancy rate $(H(t) + \eta H_c) / H_c$')
title('Indirect cost $\omega(t) = \max( 0, \min(H(t) + \eta H_c - H_c, \eta H_c ) )$')
grid on

saveas(f2,'./images/indirectcost_omegat.png')
