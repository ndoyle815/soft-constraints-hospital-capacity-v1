% script to generate schematic of alpha(t) curve
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

% prevalence vector
ICUmax = 5000;
ICU = (0:ICUmax/1000:ICUmax) + para.eta*para.Ibar;

% hard constraint
alpha_hard = para.alpha(1).*ones(1,length(ICU));
alpha_hard(ICU>=para.Ibar) = para.alpha(2);

% PLOT ALPHA LOGISTIC FUNCTION
vvec = [1.5 2 3 5 10].*1e1;

pct = 100;

f1 = figure(1);
f1.Position = [400 1200 1000 300];
set(gca,'InnerPosition',[0.075 0.2 0.9 0.7])
hold all
patch([0 0 pct.*para.eta 100.*para.eta], [0 2*para.alpha(2) 2*para.alpha(2) 0], 'k', 'FaceAlpha',0.4, 'EdgeColor','none', 'DisplayName','Non-disease patients')
for v = 1:length(vvec)
    alpha_logistic = para.alpha(1) + (para.alpha(2)-para.alpha(1))./(1 + exp(-vvec(v).*(ICU./para.Ibar - 1)));
    dispname = ['Soft $\alpha$, $v = ', ' ', num2str(vvec(v)), '$'];
    plot(pct.*ICU./para.Ibar,alpha_logistic,'Color',vcols(v,:),'DisplayName',dispname)
end
plot(pct.*ICU./para.Ibar,alpha_hard,'Color',myred,'DisplayName','Hard $\alpha$')
axis([0 pct.*1.4 0 1.2*para.alpha(2)])
legend('Location','west','AutoUpdate','off','FontSize',14)
x1 = xline(pct.*para.Ibar./para.Ibar,'k:','','Interpreter','latex','FontSize',16,'LabelOrientation','horizontal','LabelVerticalAlignment','top','LabelHorizontalAlignment','left','LineWidth',2,'Layer','bottom');
y1 = yline(para.alpha(1),'k--','$\alpha_0$','Interpreter','latex','FontSize',16,'LabelHorizontalAlignment','right','LineWidth',2,'Layer','bottom');
y2 = yline(para.alpha(2),'k--','$\alpha_1$','Interpreter','latex','FontSize',16,'LabelHorizontalAlignment','right','LineWidth',2,'Layer','bottom');
xticks(pct.*(0:0.1:1.4))
xlabel('ICU Occupancy rate $(\%)$')
ylabel('Daily cost $\alpha(t)$')
title('Daily cost per patient $\alpha(t)$')
grid on

saveas(f1,'./images/F1_alphat.png')
