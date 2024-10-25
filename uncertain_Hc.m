% script to test whether the two constraints lead to a different decision
% on the "Lockdown" vs "No Lockdown" question
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

% load results
load('./mats/results.mat')
clear C1_softcosts C2_softcosts

% mean and variance on w vectors
vmeans = vs;
background = para.eta*para.Ibar;
Hmax = round(2*(1-para.eta)*para.Ibar,-3);
Nhs = Nvs;
Hs = (0:Hmax/(Nhs-1):Hmax);
Hs = (1000:3000/(Nhs-1):4000);

% store cost statistics
C1 = zeros(Nvs,Nhs);
C2 = zeros(Nvs,Nhs);
Stratboundary = zeros(1,Nhs);
inc_or_prev = 2;
whichconstr = 'soft';

% calculate soft costs
tic;
for h = 1:Nhs
    if mod(h,100) == 0
        h
    end

    para.Ibar = background + Hs(h);

    % cost for each strategy, all scaling factors v
    C1(:,h) = sum(compute_cost(out1,para,whichconstr,vs',inc_or_prev,indirect),2);
    C2(:,h) = sum(compute_cost(out2,para,whichconstr,vs',inc_or_prev,indirect),2);

    Stratboundary(h) = min([find(C1(:,h) < C2(:,h),1,'first') Nvs+1]);

end
toc;

%% Plotting

% min and max expected costs for colour range
Cmin = min((C2 - C1),[],'all');
Cmax = max((C2 - C1),[],'all');
Cabs = min(-Cmin,Cmax);

f1 = figure(1);
f1.Position = [100 1000 550 450];
colormap(BGcolormap)

imagesc((C2 - C1))
hold on
plot(1:Nhs, Stratboundary, 'r')
set(gca,'OuterPosition',[0.01 0.01 0.9 0.94])
set(gca,'YDir','normal')
clim([-Cabs Cabs])
ylabel('$v$','Rotation',0)
xlabel('Available ICU beds $(1 - \eta) H_c$')
yticks(1:round(Nvs/10):Nvs)
xticks(1:round(Nhs/5):Nhs)
yticklabels(vs(1:round(Nvs/10):Nvs))
xticklabels(Hs(1:round(Nhs/5):Nhs))
xtickangle(0)
ylim([0 find(vmeans==110)])
title('$C_{MS}(v) - C_{RH}(v)$')

h = axes(gcf,'visible','off'); 
h.Title.Visible = 'on';
h.XLabel.Visible = 'on';
h.YLabel.Visible = 'on';

dh = 0.02;
c = colorbar(h,'Position',[0.86 0.15 0.02 0.68],'FontSize',16,'TickLabelInterpreter','Latex');  % attach colorbar to h
colormap(c);
clim(h, [-Cabs Cabs]);

saveas(f1,'./images/F5_uncertain_h.png')


f2 = figure(2);
f2.Position = [900 1000 550 450];
colormap(BGcolormap)

Cmin = min((C2 - C1)./C1,[],'all');
Cmax = max((C2 - C1)./C1,[],'all');
Cabs = min(-Cmin,Cmax);

imagesc((C2 - C1)./C1)
hold on
plot(1:Nhs, Stratboundary, 'r')
set(gca,'OuterPosition',[0.01 0.01 0.9 0.94])
set(gca,'YDir','normal')
clim([-Cabs Cabs])
ylabel('$v$','Rotation',0)
xlabel('Available ICU beds $(1 - \eta) H_c$')
yticks(1:round(Nvs/10):Nvs)
xticks(1:round(Nhs/5):Nhs)
yticklabels(vs(1:round(Nvs/10):Nvs))
xticklabels(Hs(1:round(Nhs/5):Nhs))
xtickangle(0)
ylim([0 find(vmeans==110)])
title('$(C_{MS}(v) - C_{RH}(v))/C_{RH}(v)$')

h = axes(gcf,'visible','off'); 
h.Title.Visible = 'on';
h.XLabel.Visible = 'on';
h.YLabel.Visible = 'on';

c = colorbar(h,'Position',[0.86 0.15 0.02 0.68],'FontSize',16,'TickLabelInterpreter','Latex');  % attach colorbar to h
colormap(c);
clim(h, [-Cabs Cabs]);

saveas(f2,'./images/uncertain_h2.png')
