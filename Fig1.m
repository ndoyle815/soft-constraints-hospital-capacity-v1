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
% para.alpha(1) = 0.281;

% prevalence vector
ICUmax = 5000;
ICU = (0:ICUmax/1000:ICUmax) + para.eta*para.Ibar;

% hard constraint
alpha_hard = para.alpha(1).*ones(1,length(ICU));
alpha_hard(ICU>=para.Ibar) = para.alpha(2);

% PLOT ALPHA LOGISTIC FUNCTION
vvec = [1.5 2 3 5 10].*1e1;

f1 = figure(1);
f1.Position = [400 1200 1000 300];
set(gca,'InnerPosition',[0.075 0.2 0.9 0.7])
hold all
patch([0 0 para.eta para.eta], [0 2*para.alpha(2) 2*para.alpha(2) 0], 'k', 'FaceAlpha',0.4, 'EdgeColor','none', 'DisplayName','Non-disease patients')
for v = 1:length(vvec)
    alpha_logistic = para.alpha(1) + (para.alpha(2)-para.alpha(1))./(1 + exp(-vvec(v).*(ICU./para.Ibar - 1)));
    dispname = ['Soft $\alpha$, $v = ', ' ', num2str(vvec(v)), '$'];
    plot(ICU./para.Ibar,alpha_logistic,'Color',vcols(v,:),'DisplayName',dispname)
end
plot(ICU./para.Ibar,alpha_hard,'Color',myred,'DisplayName','Hard $\alpha$')
axis([0 1.5 0 1.2*para.alpha(2)])
legend('Location','west','AutoUpdate','off','FontSize',16)
x1 = xline(para.Ibar./para.Ibar,'k--','','Interpreter','latex','FontSize',16,'LabelOrientation','horizontal','LabelVerticalAlignment','top','LabelHorizontalAlignment','left','LineWidth',2,'Layer','bottom');
y1 = yline(para.alpha(1),'k--','$\alpha_0$','Interpreter','latex','FontSize',16,'LabelHorizontalAlignment','right','LineWidth',2,'Layer','bottom');
y2 = yline(para.alpha(2),'k--','$\alpha_1$','Interpreter','latex','FontSize',16,'LabelHorizontalAlignment','right','LineWidth',2,'Layer','bottom');

% fit curve to bravata points
% myrichards = @(a,c,v,x)(para.alpha(1) + (para.alpha(2)-para.alpha(1))./(1 + a.*exp(-v.*(x - c))).^(1/a));
% xdata = para.eta + (1 - para.eta).*[0.125 0.375 0.625 0.875 1.125]';
% %xdata = para.eta + (4114/para.Ibar).*[0.125 0.375 0.625 0.875 1.125]';
% ydata = para.alpha(1).*[1 1.1 1.15 1.67 2.35]';
% myfit = fit(xdata,ydata,myrichards,'Lower',[0 0 0])
% plot(ICU./para.Ibar, myrichards(myfit.a, myfit.c, myfit.v, ICU./para.Ibar), "Color", mypurple)

% scatter(0.925,para.alpha(1)*1.23,40,'k','filled','d')
% scatter(xdata, ydata, 40, 'k', 'filled', 's')

xticks(0:0.1:1.5)
xlabel('Occupancy rate $\textrm{ICU}^o(t) / \textrm{ICU}_c$')
% xlabel('Total ICU Occupancy rate $(H(t) + \eta H_c) / H_c$')
title('Daily cost per patient $\alpha(t)$')
grid on

saveas(f1,'./images/F1_alphat.png')

% for presentation
f2 = figure(2);
f2.Position = [400 200 400 300];
set(gca,'InnerPosition',[0.1 0.2 0.8 0.65])
hold all
patch([0 0 para.eta para.eta], [0 2*para.alpha(2) 2*para.alpha(2) 0], 'k', 'FaceAlpha',0.4, 'EdgeColor','none', 'DisplayName','Background')
plot(ICU./para.Ibar,alpha_hard,'Color',myred,'DisplayName','Hard $\alpha$')
axis([0.6 1.4 para.alpha(1)-0.25 para.alpha(2)+0.25])
xline(para.Ibar./para.Ibar,'k--','','Interpreter','latex','FontSize',16,'LabelOrientation','horizontal','LabelVerticalAlignment','top','LabelHorizontalAlignment','left','LineWidth',2,'Layer','bottom');
yticks(para.alpha)
yticklabels({'$\alpha_0$', '$\alpha_1$'})
xlabel('Occupancy rate $\textrm{ICU}^o(t) / \textrm{ICU}_c$')
title('$\alpha(t)$')
grid on

saveas(f2,'./hardconstrex.png')

% for presentation
f3 = figure(3);
f3.Position = [900 200 400 300];
set(gca,'InnerPosition',[0.1 0.2 0.8 0.65])
hold all
patch([0 0 para.eta para.eta], [0 2*para.alpha(2) 2*para.alpha(2) 0], 'k', 'FaceAlpha',0.4, 'EdgeColor','none', 'DisplayName','Background')
v = 3;
alpha_logistic = para.alpha(1) + (para.alpha(2)-para.alpha(1))./(1 + exp(-vvec(v).*(ICU./para.Ibar - 1)));
%     dispname = ['Soft $\alpha$, $v = ', ' ', num2str(vvec(v)), '$'];
plot(ICU./para.Ibar,alpha_logistic,'Color',vcols(v,:),'DisplayName',dispname)
axis([0.6 1.4 para.alpha(1)-0.25 para.alpha(2)+0.25])
xline(para.Ibar./para.Ibar,'k--','','Interpreter','latex','FontSize',16,'LabelOrientation','horizontal','LabelVerticalAlignment','top','LabelHorizontalAlignment','left','LineWidth',2,'Layer','bottom');
yticks(para.alpha)
yticklabels({'$\alpha_0$', '$\alpha_1$'})
xlabel('Occupancy rate $\textrm{ICU}^o(t) / \textrm{ICU}_c$')
title('$\alpha(t)$')
grid on

saveas(f3,'./softconstrex.png')
