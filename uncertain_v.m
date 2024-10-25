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

% mean and variance on w vectors
vmeans = vs;
vstds = 0:100/(Nvs-1):100;
Nvm = length(vmeans);
Nvv = length(vstds);

% store cost statistics
EC1   = zeros(Nvm,Nvv);
EC2   = zeros(Nvm,Nvv);
EsqC1 = zeros(Nvm,Nvv);
EsqC2 = zeros(Nvm,Nvv);
vv_Stratboundary = zeros(1,Nvv);
Prob_S2opt = zeros(Nvm,Nvv);
vv_StratProbboundary = zeros(1,Nvv);

% calculate soft costs
tic;
for vv = 1:Nvv
    if mod(vv,100) == 0
        vv
    end

    for vm = 1:Nvm

       % generate the normal distribution
       if vv > 1
           tmp_normdist = normpdf(vs,vmeans(vm),vstds(vv));
       else
           % if variance is zero, dirac distribution
           tmp_normdist = double(vs == vmeans(vm));
       end
       tmp_normdist = tmp_normdist./sum(tmp_normdist);  % normalised

       % expected cost of each strategy
       EC1(vm,vv) = sum(C1_softcosts.*tmp_normdist);
       EC2(vm,vv) = sum(C2_softcosts.*tmp_normdist);

       % second moment cost of each strategy
       EsqC1(vm,vv) = sum((C1_softcosts.^2).*tmp_normdist);
       EsqC2(vm,vv) = sum((C2_softcosts.^2).*tmp_normdist);

       % compute probability that Strategy 1 is optimal
       where_S2opt = C2_softcosts < C1_softcosts;
       Prob_S2opt(vm,vv) = sum(tmp_normdist(where_S2opt));
   
    end

    vv_Stratboundary(vv) = min([find(EC1(:,vv) < EC2(:,vv),1,'first') Nvm+1]);
    vv_StratProbboundary(vv) = min([find(Prob_S2opt(:,vv) < 0.5,1,'first') Nvm+1]);

end
toc;

% variance of cost for each strategy
VarC1 = EsqC1 - EC1.^2;
VarC2 = EsqC2 - EC2.^2;

%% Plotting

% min and max expected costs for colour range
Emin = min((EC2 - EC1),[],'all');
Emax = max((EC2 - EC1),[],'all');
Eabs = max(-Emin,Emax);

f1 = figure(1);
f1.Position = [100 1000 500 450];
colormap(BGcolormap)

imagesc((EC2 - EC1))
hold on
plot(1:Nvv, vv_Stratboundary, 'r')
hold on
plot(find(ismember(3.*vstds,vmeans)), find(ismember(vmeans,3.*vstds)), 'k')
text(800,1000,'$\mu_v = 3 \sigma_v$','FontSize',16,'Color','k')
set(gca,'OuterPosition',[0.01 0.01 0.9 0.94])
set(gca,'YDir','normal')
clim([-Eabs Eabs])
ylabel('$\mu_v$','Rotation',0);
xlabel('$\sigma_v$')
yticks(1:round(Nvm/10):Nvm)
xticks(1:round(Nvv/10):Nvv)
yticklabels(vmeans(1:round(Nvm/10):Nvm))
xticklabels(vstds(1:round(Nvv/10):Nvv))
xtickangle(0)
ylim([0 find(vmeans==110)])
xlim([0 find(vstds==30)])
title('$E[C_{MS}(v)] - E[C_{RH}(v)]$');

h = axes(gcf,'visible','off'); 
h.Title.Visible = 'on';
h.XLabel.Visible = 'on';
h.YLabel.Visible = 'on';

c = colorbar(h,'Position',[0.86 0.15 0.02 0.68],'FontSize',16,'TickLabelInterpreter','Latex');  % attach colorbar to h
colormap(c);
clim(h, [-Eabs Eabs]);

saveas(f1,'./images/F4C_uncertain_v_expcost.png')


f2 = figure(2);
f2.Position = [900 1000 500 450];

colormap(BGcolormap(end:-1:1,:))

imagesc(Prob_S2opt)
hold on
plot(1:Nvv, vv_StratProbboundary, 'r')
hold on
plot(find(ismember(3.*vstds,vmeans)), find(ismember(vmeans,3.*vstds)), 'k')
text(800,1000,'$\mu_v = 3 \sigma_v$','FontSize',16,'Color','k')
set(gca,'OuterPosition',[0.01 0.01 0.9 0.94])
set(gca,'YDir','normal')
clim([0 1])
ylabel('$\mu_v$','Rotation',0)
xlabel('$\sigma_v$')
yticks(1:round(Nvm/10):Nvm)
xticks(1:round(Nvv/10):Nvv)
yticklabels(vmeans(1:round(Nvm/10):Nvm))
xticklabels(vstds(1:round(Nvv/10):Nvv))
xtickangle(0)
ylim([0 find(vmeans==110)])
xlim([0 find(vstds==30)])
title('$P[C_{MS}(v) < C_{RH}(v)]$');

h = axes(gcf,'visible','off'); 
h.Title.Visible = 'on';
h.XLabel.Visible = 'on';
h.YLabel.Visible = 'on';

c = colorbar(h,'Position',[0.86 0.15 0.02 0.68],'FontSize',16,'TickLabelInterpreter','Latex');  % attach colorbar to h
colormap(c);
clim(h, [0 1]);

saveas(f2,'./images/F4D_uncertain_v_proboopt.png')
