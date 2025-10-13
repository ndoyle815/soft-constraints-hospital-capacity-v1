% script to to simulate and produce Figure 4C,D in manuscript
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
title('Expected MS benefit')

saveas(f1,'./images/F4C_uncertain_v_expcost.png')


f2 = figure(2);
f2.Position = [800 1000 500 400];
set(gca,'InnerPosition',[0.175 0.15 0.75 0.75])
hold all

colormap(BGEcolormap(5:end,:))

[M1,RH] = contourf(RH90th - MS90th,-4000:1000:5000,'-','LineWidth',1,'LabelColor','k','LabelSpacing',300);
clabel(M1,RH,'FontSize',0.9*16,'Interpreter','latex')
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
title('MS benefit at 90th percentile')

saveas(f2,'./images/F4D_decisionbounds.png')
