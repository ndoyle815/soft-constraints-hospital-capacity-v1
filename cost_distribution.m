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

% generate distribution of v
[~,minidx] = min(abs(C1_softcosts - C2_softcosts));  % pick mu_v so means are as close as possible
mu_v = vs(minidx);
mu_v = 50;
sigma_v = 15;
vdist = normpdf(vs,mu_v,sigma_v);
vdist = vdist./sum(vdist);

% generate samples
Nsamples = 5000;
vsamples = randsample(vs,Nsamples,'true',vdist)';

% calculate soft costs by sampling from distribution
whichconstr = 'soft';

C1_mu_v = sum(compute_cost(out1,para,whichconstr,mu_v,inc_or_prev,indirect),2);
C1s     = sum(compute_cost(out1,para,whichconstr,vsamples,inc_or_prev,indirect),2);

C2_mu_v = sum(compute_cost(out2,para,whichconstr,mu_v,inc_or_prev,indirect),2);
C2s     = sum(compute_cost(out2,para,whichconstr,vsamples,inc_or_prev,indirect),2);

% plotting
f1 = figure(1);
f1.Position = [200 1000 600 200];

plot(vs,vdist,'r')
xlabel('$v$');
% xl.Position(2) = xl.Position(2) + 0.1;
ylabel('Probability')
title(['$v \sim \mathcal{N}(\mu_v = ',num2str(mu_v),', \sigma_v = ',num2str(sigma_v),')$'])
xlim([0 max(vs)])
grid on

saveas(f1,'./images/vdist.png')

%% 

f2 = figure(2);
f2.Position = [800 600 900 600];

subplot(2,1,1)
histogram(C1s,'FaceColor',mygreen)
xline(C1_mu_v,'r','LineWidth',3)
xlabel('Cost $C_{\mathcal{S}_1}(v)$')
ylabel('Frequency')
xlim([min([C1s; C2s]) max([C1s; C2s])])
legend(stratnames{1},'Location','north')
title(['$v \sim \mathcal{N}(\mu_v = ',num2str(mu_v),', \sigma_v = ',num2str(sigma_v),')$'])
grid on

subplot(2,1,2)
histogram(C2s,'FaceColor',myblue)
xline(C2_mu_v,'r','LineWidth',3)
xlabel('Cost $C_{\mathcal{S}_2}(v)$')
ylabel('Frequency')
%xlim([min([C1s; C2s]) max([C1s; C2s])])
legend(stratnames{2},'Location','north')
grid on

saveas(f2,'./images/cost_distribution.png')

%% Sample and plot cost distributions in a loop for different sigma

% vector of standard deviations
% sigs = 5:5:ceil(mu_v/3);
sigs = [2 4 8 16];
nsigs = length(sigs);

% generate normally distributed v for each sigma
vdists = normpdf(vs, mu_v, sigs');

% plot v distributions
f3 = figure(3);
f3.Position = [200 650 nsigs*250 225];
hold all

for s = 1:nsigs
    subplot(1,nsigs,s)
    plot(vs,vdists(s,:),'r')
    axis([0 max(vs) 0 0.21])
    xl = xlabel('$v$');
    % xl.Position(2) = xl.Position(2) + 0.1;
    if s == 1
        yl = ylabel('PDF');
        yl.Position(1) = yl.Position(1) - 0.15;
    else
        yticklabels([])
    end
    title(['$\sigma_v = $', ' ', num2str(sigs(s))])
end

saveas(f3,'./images/F4A_v_distributions.png')

% plot cost distributions
f4 = figure(4);
f4.Position = [200 300 nsigs*250 225];
hold all

for s = 1:nsigs
    
    % sample and compute cost distributions
    sigma_v = sigs(s);

    % generate samples
    vsamples = randsample(vs,Nsamples,'true',vdists(s,:))';

    C1s = sum(compute_cost(out1,para,whichconstr,vsamples,inc_or_prev,indirect),2);
    C2s = sum(compute_cost(out2,para,whichconstr,vsamples,inc_or_prev,indirect),2);

    % plotting
    cd Violinplot-Matlab-master/

    subplot(1,nsigs,s)
    violinplot([C1s C2s], {'RH','MS'}, 'ViolinColor', [mygreen; myblue], 'ShowMean', true)
    axis([0.6 2.4 10e3 20e3])
    xticks([1 2])
    ax = gca;
    ax.YAxis.Exponent = 0;
    if s == 1
        yl = ylabel('Cost $C_{\mathcal{S}}(v)$');
        yl.Position(1) = yl.Position(1) - 0.15;
    else
        yticklabels([])
    end
    title(['$\sigma_v = $', ' ', num2str(sigma_v)])
    cd ..

end

saveas(f4,'./images/F4B_cost_distributions.png')
