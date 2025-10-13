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
Nvv = 121;
vmeans = vs;
vstds = 0:30/(Nvv-1):30;
Nvm = length(vmeans);
Nsamples = 5;
whichconstr = 'soft';
apply_discounting = 0;
myquantiles = [0.5 0.75 0.9 0.95 0.99];
Nq = length(myquantiles);
vstds_quantiles = 0:15:30;
Nvv_quants = length(vstds_quantiles);
vstds_quantiles_idxs = find(ismember(vstds,vstds_quantiles));

% store cost statistics
[ERH, EMS, EsqRH, EsqMS, RH90th, MS90th] = deal(zeros(Nvm,Nvv));
[RHquants, MSquants] = deal(zeros(Nq,Nvm,Nvv_quants));

% calculate soft costs
tic;
for vv = 1:Nvv
    if mod(vv,10) == 0
        vv
    end

    for vm = 1:Nvm

       % generate the normal distribution
       if vv > 1
           tmp_normdist = normpdf(vs,vmeans(vm),vstds(vv));
           tmp_normdist = tmp_normdist./sum(tmp_normdist);  % normalise
       else
           % if variance is zero, dirac distribution (sums to 1)
           tmp_normdist = double(vs == vmeans(vm));
       end   

       % expected cost of each strategy
       ERH(vm,vv) = sum(RH_softcosts.*tmp_normdist);
       EMS(vm,vv) = sum(MS_softcosts.*tmp_normdist);

       % generate the cost distribution
       [RHcostsstorted, RHcostssortedidx] = sort(RH_softcosts);
       RHvcostdist = cumsum(tmp_normdist(RHcostssortedidx));
       v90idx = find(RHvcostdist<0.9,1,'last');
       if isempty(v90idx)
           RH90th(vm,vv) = RH_softcosts(1);
       else
           tmpv = vs(RHcostssortedidx);
           RH90th(vm,vv) = RH_softcosts(find(vs==tmpv(v90idx)));
       end

       [~, MScostssortedidx] = sort(MS_softcosts);
       MSvcostdist = cumsum(tmp_normdist(MScostssortedidx));
       v90idx = find(MSvcostdist<0.9,1,'last');
       if isempty(v90idx)
           MS90th(vm,vv) = MS_softcosts(1);
       else
           tmpv = vs(MScostssortedidx);
           MS90th(vm,vv) = MS_softcosts(find(vs==tmpv(v90idx)));
       end

       if ismember(vv,vstds_quantiles_idxs)
           vqidx = find(ismember(vstds_quantiles_idxs,vv));
           vsamples = randsample(vs,Nsamples,'true',tmp_normdist)';
           RHs = sum(compute_cost(outRH,para,whichconstr,vsamples,inc_or_prev,apply_discounting),2);
           MSs = sum(compute_cost(outMS,para,whichconstr,vsamples,inc_or_prev,apply_discounting),2);
    
           RHquants(:,vm,vqidx) = quantile(RHs,myquantiles);
           MSquants(:,vm,vqidx) = quantile(MSs,myquantiles);
       end

       % second moment cost of each strategy
       EsqRH(vm,vv) = sum((RH_softcosts.^2).*tmp_normdist);
       EsqMS(vm,vv) = sum((MS_softcosts.^2).*tmp_normdist);
   
    end

end
toc;

% variance of cost for each strategy
VarRH = EsqRH - ERH.^2;
VarMS = EsqMS - EMS.^2;

%% Plotting

% min and max expected costs for colour range
Emin = min((ERH - EMS),[],'all');
Emax = max((ERH - EMS),[],'all');
Eabs = max(-Emin,Emax);

f1 = figure(1);
f1.Position = [200 1000 500 400];
set(gca,'InnerPosition',[0.175 0.15 0.75 0.75])
hold all

Ecolormapidxs = [1 25 49 73 97 161 185 209 223 257];
BGEcolormap = BGcolormap(Ecolormapidxs(end:-1:1),:);
colormap(BGEcolormap(5:end,:))

[M1,RH] = contourf(ERH - EMS,-4000:1000:5000,'-','LineWidth',1,'LabelColor','k','LabelSpacing',300);
clabel(M1,RH,'FontSize',0.9*16,'Interpreter','latex')
% plot(find(ismember(3.*vstds,vmeans)), find(ismember(vmeans,3.*vstds)), 'Color',myorange, 'LineWidth',2, 'LineStyle',':')
set(gca,'YDir','normal')
ylabel('$\mu_v$','Rotation',0);
xlabel('$\sigma_v$')
yticks(1:round(Nvm/10):Nvm)
xticks(1:round(Nvv/6):Nvv)
yticklabels(vmeans(1:round(Nvm/10):Nvm))
xticklabels(vstds(1:round(Nvv/6):Nvv))
xtickangle(0)
ylim([1 find(vmeans>=110,1,'first')])
xlim([1 find(vstds>=30,1,'first')])
% title('$E[C_{MS}(v)] - E[C_{RH}(v)]$');
title('Expected MS benefit')

saveas(f1,'./images/F4C_uncertain_v_expcost.png')


f2 = figure(2);
f2.Position = [800 1000 500 400];
set(gca,'InnerPosition',[0.175 0.15 0.75 0.75])
hold all

colormap(BGEcolormap(5:end,:))

[M1,RH] = contourf(RH90th - MS90th,-4000:1000:5000,'-','LineWidth',1,'LabelColor','k','LabelSpacing',300);
clabel(M1,RH,'FontSize',0.9*16,'Interpreter','latex')
% plot(find(ismember(3.*vstds,vmeans)), find(ismember(vmeans,3.*vstds)), 'Color',myorange, 'LineWidth',2, 'LineStyle',':')
set(gca,'YDir','normal')
ylabel('$\mu_v$','Rotation',0);
xlabel('$\sigma_v$')
yticks(1:round(Nvm/10):Nvm)
xticks(1:round(Nvv/6):Nvv)
yticklabels(vmeans(1:round(Nvm/10):Nvm))
xticklabels(vstds(1:round(Nvv/6):Nvv))
xtickangle(0)
ylim([1 find(vmeans>=110,1,'first')])
xlim([1 find(vstds>=30,1,'first')])
% title('$E[C_{MS}(v)] - E[C_{RH}(v)]$');
title('MS benefit at 90th percentile')

saveas(f2,'./images/F4D_decisionbounds.png')


%% For each quantile, find the decision boundary
dec_boundary = zeros(Nq,Nvv_quants);

for q = 1:Nq

    Cdiff_quantile = reshape(MSquants(q,:,:) - RHquants(q,:,:),Nvm,Nvv_quants);

    for vv = 1:Nvv_quants
    
        % check whether the boundary exists
        bound_idx = find(Cdiff_quantile(:,vv)>=0,1,'first');
    
        if isempty(bound_idx)
            bound_idx = Nvm;
        end
    
        dec_boundary(q,vv) = vmeans(bound_idx);
    end
end

% plotting
f3 = figure(3);
f3.Position = [1400 1000 500 400];
set(gca,'InnerPosition',[0.15 0.125 0.75 0.75])
hold all

markers = {'o','s','D','^','pentagram'};
quantcols = [252,146,114; 251,106,74; 239,59,44; 203,24,29; 153,0,13]./255;

for q = 1:Nq
    plot(vstds_quantiles, dec_boundary(q,:), 'k-', 'LineWidth', 2, 'Marker',markers{q}, ...
         'MarkerFaceColor',quantcols(q,:), 'MarkerSize',10, 'DisplayName',['$q = ',' ',num2str(myquantiles(q)),'$'])
end

legend('FontSize',13,'Location','southeast','AutoUpdate','off','NumColumns',2)
patch([0 0 11 11],[25 35 35 25],myblue,'FaceAlpha',0.8,'Edgecolor','none')
patch([0 0 11 11],[110 120 120 110],mygreen,'FaceAlpha',0.8,'Edgecolor','none')
text(1,30,'MS optimal','FontSize',16)
text(1,115,'RH optimal','FontSize',16)
set(gca,'InnerPosition',[0.175 0.15 0.75 0.75])
set(gca,'YDir','normal')
ylabel('$\mu_v$','Rotation',0)
xlabel('$\sigma_v$')
yticks(vmeans(1:round(Nvm/10):Nvm))
xticks(vstds(1:round(Nvv/10):Nvv))
xtickangle(0)
ylim([min(vmeans) 130])
xlim([min(vstds) 30])
title('Decision boundaries');
grid on

% saveas(f2,'./images/F4D_decisionbounds.png')
