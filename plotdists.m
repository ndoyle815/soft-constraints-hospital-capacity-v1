% script to plot model distributions for supplementary material
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

% serial interval distribution
mu = 5.4;
sig = 1.5;
wtimes = 1:500;
w = gampdf(wtimes,(mu/sig)^2,1/(mu/sig^2));
w = w./sum(w);

f1 = figure(1);
f1.Position = [200 800 1100 1100];

subplot(2,2,1)
bar(wtimes,w)
xlabel('Time $t$ (days)')
ylabel('Probability $s = t$')
title('Serial interval distribution $w(s)$')
axis([0 20 0 0.3])
grid on

subplot(2,2,2)
bar(0:length(para.Dist_SC)-1, para.Dist_SC)
xlabel('Time between Symptoms and ICU $t$ (days)')
ylabel('Probability $s = t$')
title('Distribution $D^{I \rightarrow H}(s)$')
xlim([-1 length(para.Dist_SC)])
grid on

subplot(2,2,3)
bar(1:length(para.Time_C), para.Time_C)
xlabel('Time in ICU $t$ (days)')
ylabel('Probability duration $s \geq t$')
title('Distribution $T^H(s)$')
axis([-1 60 0 1.1])
grid on

subplot(2,2,4)
newdist = [0 para.Time_C(1:end-1) - para.Time_C(2:end)];
newdist = para.Time_C./sum(para.Time_C);

bar(1:length(newdist), newdist)
xlabel('Time $t$ (days)')
ylabel('Probability $s = t$')
title('Weighted ICU time distribution $\bar{T}^H(s)$')
xlim([-1 60])
grid on

saveas(gcf,'./images/FS1_distributions.png')
