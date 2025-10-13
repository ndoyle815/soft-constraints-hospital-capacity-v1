% script to to simulate and produce Figure S5 in manuscript supplement
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
clear RH_softcosts MS_softcosts

% mean and variance on w vectors
vmeans = vs;
background = para.eta*para.Ibar;
Hmax = round(2*(1-para.eta)*para.Ibar,-3);
Nhs = Nvs;
Hs = (1000:4000);
Nhs = length(Hs);

% store cost statistics
RH = zeros(Nvs,Nhs);
MS = zeros(Nvs,Nhs);
Stratboundary = zeros(1,Nhs);
inc_or_prev = 2;
whichconstr = 'soft';
apply_discounting = 0;

% calculate soft costs
tic;
for h = 1:Nhs
    if mod(h,100) == 0
        h
    end

    para.Ibar = background + Hs(h);

    % cost for each strategy, all scaling factors v
    RH(:,h) = sum(compute_cost(outRH,para,whichconstr,vs',inc_or_prev,apply_discounting),2);
    MS(:,h) = sum(compute_cost(outMS,para,whichconstr,vs',inc_or_prev,apply_discounting),2);

    Stratboundary(h) = min([find(RH(:,h) < MS(:,h),1,'first') Nvs+1]);

end
toc;

%% Plotting

% min and max expected costs for colour range
Cmin = min((RH - MS),[],'all');
Cmax = max((RH - MS),[],'all');
Cabs = min(-Cmin,Cmax);

f1 = figure(1);
f1.Position = [100 1000 500 400];
set(gca,'InnerPosition',[0.175 0.15 0.75 0.75])

Ecolormapidxs = [1 25 49 73 97 161 185 209 223 257];
BGEcolormap = BGcolormap(Ecolormapidxs(end:-1:1),:);
colormap(BGEcolormap(5:end,:))

[M1,Cont1] = contourf(RH - MS,-4000:1000:5000,'-','LineWidth',1,'LabelColor','k','LabelSpacing',400);
clabel(M1,Cont1,'FontSize',0.9*16,'Interpreter','latex')
set(gca,'YDir','normal')
ylabel('Soft constraint scaling $v$')
xlabel('Available ICU beds')
yticks(1:round(Nvs/10):Nvs)
xticks(1:round(Nhs/6):Nhs)
yticklabels(vs(1:round(Nvs/10):Nvs))
xticklabels(Hs(1:round(Nhs/6):Nhs))
xtickangle(0)
ylim([0 find(vmeans==110)])
title('Cost difference')


saveas(f1,'./images/supplement/FS4_scen2_uncertainICU.png')
