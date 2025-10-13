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
apply_discounting = 0;

% load results
load('./mats/results.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TO PRODUCE SENSITIVITY ANALYSIS FIGURES, MODIFY THE BELOW VARIABLE
% KEY:
% 0: Fig 4a,b
% 1: Fig S6a,b  (mu_v = 75)
% 2: Fig S6c,d  (mu_v = 100)
sensitivity_analysis = 0;

mus = [50 75 100];
mu_v = mus(1);

if sensitivity_analysis == 1
    mu_v = mus(2);

elseif sensitivity_analysis == 2
    mu_v = mus(3);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate distribution of v
[~,minidx] = min(abs(RH_softcosts - MS_softcosts));  % pick mu_v so means are as close as possible

sigma_v = 15;
vdist = normpdf(vs,mu_v,sigma_v);
vdist = vdist./sum(vdist);

% generate samples
Nsamples = 5000;
vsamples = randsample(vs,Nsamples,'true',vdist)';

% calculate soft costs by sampling from distribution
whichconstr = 'soft';

RH_mu_v = sum(compute_cost(outRH,para,whichconstr,mu_v,inc_or_prev,apply_discounting),2);
RHs     = sum(compute_cost(outRH,para,whichconstr,vsamples,inc_or_prev,apply_discounting),2);

MS_mu_v = sum(compute_cost(outMS,para,whichconstr,mu_v,inc_or_prev,apply_discounting),2);
MSs     = sum(compute_cost(outMS,para,whichconstr,vsamples,inc_or_prev,apply_discounting),2);

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

%% 

f2 = figure(2);
f2.Position = [800 600 900 600];

subplot(2,1,1)
histogram(RHs,'FaceColor',mygreen)
xline(RH_mu_v,'r','LineWidth',3)
xlabel('Cost $C_{\mathcal{S}_1}(v)$')
ylabel('Frequency')
xlim([min([RHs; MSs]) max([RHs; MSs])])
legend(stratnames{1},'Location','north')
title(['$v \sim \mathcal{N}(\mu_v = ',num2str(mu_v),', \sigma_v = ',num2str(sigma_v),')$'])
grid on

subplot(2,1,2)
histogram(MSs,'FaceColor',myblue)
xline(MS_mu_v,'r','LineWidth',3)
xlabel('Cost $C_{\mathcal{S}_2}(v)$')
ylabel('Frequency')
%xlim([min([RHs; MSs]) max([RHs; MSs])])
legend(stratnames{2},'Location','north')
grid on

% if mu_v == 50
%     saveas(f1,'./images/vdist.png')
%     saveas(f2,'./images/cost_distribution.png')
% end

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
    % xl.Position(2) = xl.Position(2) + 0.05;
    if s == 1
        yl = ylabel('PDF');
        yl.Position(1) = yl.Position(1) - 0.15;
    else
        yticklabels([])
    end
    title(['$\sigma_v = $', ' ', num2str(sigs(s))])

    ax = gca;
    ax.InnerPosition(2) = 0.24;
    ax.InnerPosition(4) = 0.64;
end


% plot cost distributions
f4 = figure(4);
f4.Position = [200 300 nsigs*250 225];
hold all

tic;
for s = 1:nsigs
    
    % sample and compute cost distributions
    sigma_v = sigs(s);

    % generate samples
    vsamples = randsample(vs,Nsamples,'true',vdists(s,:))';

    RHs = sum(compute_cost(outRH,para,whichconstr,vsamples,inc_or_prev,apply_discounting),2);
    MSs = sum(compute_cost(outMS,para,whichconstr,vsamples,inc_or_prev,apply_discounting),2);

    % plotting
    cd Violinplot-Matlab-master/

    subplot(1,nsigs,s)
    violinplot([RHs MSs], {'RH','MS'}, 'ViolinColor', [mygreen; myblue], 'ShowMean', true);
    axis([0.6 2.4 8e3 14e3])
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
toc;

if mu_v == 50
    saveas(f3,'./images/F4A_v_distributions.png')
    saveas(f4,'./images/F4B_cost_distributions.png')
else
    saveas(f3,['./images/supplement/FS5_v_dists_mu', num2str(mu_v), '.png'])
    saveas(f4,['./images/supplement/FS5_cost_dists_mu', num2str(mu_v), '.png'])
end

%% Plot big one for presentation

% set(0,'defaultaxesfontsize',24)
% 
% f6 = figure(6);
% f6.Position = [200 200 500 500];
% 
% cd Violinplot-Matlab-master/
% violinplot([RHs MSs], {'RH','MS'}, 'ViolinColor', [mygreen; myblue], 'ShowMean', true);
% axis([0.6 2.4 8e3 14e3])
% xticks([1 2])
% ax = gca;
% ax.YAxis.Exponent = 0;
% % yl = ylabel('Cost $C_{\mathcal{S}}(v)$');
% % yl.Position(1) = yl.Position(1) - 0.15;
% xlabel('Control Strategy')
% title(['Cost distribution'])
% cd ..
% 
% saveas(f6,'~/Pictures/test.png')