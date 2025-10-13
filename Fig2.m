% script to simulate and produce Figure 2 in manuscript
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

new_background = para.eta*para.Ibar;   % doesnt change
bedincrease = 1;
new_available = bedincrease*(1 - para.eta)*para.Ibar;
para.Ibar = new_background + new_available;
para.eta = new_background/para.Ibar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TO PRODUCE SENSITIVITY ANALYSIS FIGURES, MODIFY THE BELOW VARIABLE
% KEY:
% 0: Fig 2
% 1: Fig S7a  (alpha_0 = 0.2)
% 2: Fig S7b  (alpha_0 = 0.5)
% 3: Fig S10a (alpha_1 = 5 alpha_0)
% 4: Fig S10b (alpha_1 = 2 alpha_0)
sensitivity_analysis = 4;

% TO PRODUCE FS3 DISTRIBUTIONS ON BEDS AVAILABLE, SET THE BELOW VARIABLE
% TO 1
plot_distributions = 0;

alpha0s = [0.2 0.5];
alpha1s = [para.alpha(1)*5 para.alpha(1)*2];

if sensitivity_analysis == 1
    para.alpha(1) = alpha0s(1);

elseif sensitivity_analysis == 2
    para.alpha(1) = alpha0s(2);

elseif sensitivity_analysis == 3
    para.alpha(2) = alpha1s(1);

elseif sensitivity_analysis == 4
    para.alpha(2) = alpha1s(2);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stratnames = {'Late Suppression (LS)','Early Suppression (ES)'};

% parameters
R0 = 2;            % R_0
I0 = 10;           % Initial no. cases
maxtime = 750;     % simulation time
whichR = {'Instantaneous','Case'};
Rtype = whichR{1};

% R_t timeseries
ti1 = [1 70 140];
Ri1 = [R0 0.7 0.7];
ti2 = [1 63 140];
Ri2 = [R0 0.91 0.91];

% serial interval distribution
mu = 5.4;
sig = 1.5;
w = discretise_distribution(1:maxtime,'Gamma',mu,sig);

% run models
[outLS] = RENEWALmodel(Ri1,ti1,w,I0,maxtime,Rtype,para);
[outES] = RENEWALmodel(Ri2,ti2,w,I0,maxtime,Rtype,para);

[sum(outLS.new_ICU) sum(outES.new_ICU)]


% Cost Evaluation

% vectors for parameters v (scaling logistic curve) and k (relative cost of lockdown)
dv = 0.2;
dk = 0.05;
vs = 10:dv:210;
ks = 0:dk:50;
nv = length(vs);
nk = length(ks);

if bedincrease < 1
    mult = 3;
else
    mult = 1;
end
dk = mult*dk;
ks = mult.*ks;

% store total costs for two scenarios with different constraints
inc_or_prev = 2;
apply_discounting = 1;

% costs of control
LS_cost_of_control = 0;
% additional control effort of Lockdown strategy is proportional to
% (RIT^2) over the intervention period
% RIT = reduction in transmission
RIT = 1 - Ri2(2)/Ri2(1);
duration = ti2(3) - ti2(2);
ES_cost_of_control = ks.*duration*RIT^2;

% costs of disease (hard constraint)
whichconstr = 'hard';

LS_cost_of_inf_hard = sum(compute_cost(outLS,para,whichconstr,vs(1),inc_or_prev,apply_discounting),2);
LS_hardcosts = repmat(LS_cost_of_inf_hard + LS_cost_of_control, nv, nk);

ES_cost_of_inf_hard = sum(compute_cost(outES,para,whichconstr,vs(1),inc_or_prev,apply_discounting),2);
ES_hardcosts = repmat(ES_cost_of_inf_hard,1,nk) + repmat(ES_cost_of_control,nv,1);

% costs of disease (soft constraint)
whichconstr = 'soft';

LS_cost_of_inf_soft = sum(compute_cost(outLS,para,whichconstr,vs,inc_or_prev,apply_discounting),2);
LS_softcosts = repmat(LS_cost_of_inf_soft + LS_cost_of_control,1,nk);

ES_cost_of_inf_soft = sum(compute_cost(outES,para,whichconstr,vs,inc_or_prev,apply_discounting),2);
ES_softcosts = repmat(ES_cost_of_inf_soft,1,nk) + repmat(ES_cost_of_control,nv,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting

% R(t) timeseries
f1 = figure(1);
f1.Position = [1300 1200 600 300];
set(gca,'InnerPosition',[0.16 0.21 0.8 0.75])
% f1.Position = [1300 1200 1000 250];
% set(gca,'InnerPosition',[0.09 0.25 0.89 0.7])
hold all
plot(outLS.t,outLS.R,'Color',mypurple)
plot(outES.t,outES.R,'Color',myorange)
legend(stratnames,'Location','northeast','AutoUpdate','off','FontSize',16)
yline(1,'k--','$R(t)=1$','Interpreter','latex','LineWidth',2,'FontSize',16,'LabelVerticalAlignment','top','LabelHorizontalAlignment','right','Layer','bottom')
xlabel('Time $t$ (days)')
ylab = ylabel('$R(t)$','Rotation',0);
ylab.Position(1) = ylab.Position(1) + 15;
ylab.Position(2) = ylab.Position(2) - 0.2;
axis([0 min([350 maxtime]) 0 1.2*max([outLS.R; outES.R],[],"all")])
grid on

% ICU prevalence
f2 = figure(2);
f2.Position = [1300 1200 600 300];
set(gca,'InnerPosition',[0.16 0.21 0.8 0.75])
% f2.Position = [1300 1200 1000 250];
% set(gca,'InnerPosition',[0.09 0.25 0.89 0.7])
hold all
patch([0 0 maxtime maxtime], [0 para.eta para.eta 0].*para.Ibar, 'k', 'FaceAlpha',0.4, 'EdgeColor','none', 'DisplayName','Background occupancy')
plot(outLS.t,outLS.in_ICU,'Color',mypurple)
plot(outES.t,outES.in_ICU,'Color',myorange)
yline(para.Ibar,'k--','Total beds ICU$_c$','Interpreter','latex','LineWidth',2,'FontSize',16,'LabelVerticalAlignment','top','LabelHorizontalAlignment','right','Layer','bottom')
xlabel('Time $t$ (days)')
ylabel('In ICU');
% title(['Available bed change: $',num2str(100*(bedincrease-1)),'\%$']);
axis([0 min([350 maxtime]) 0 max([1.3*para.Ibar outLS.in_ICU outES.in_ICU])])
grid on

% costs
f3 = figure(3);
f3.Position = [1300 1200 600 300];
set(gca,'InnerPosition',[0.16 0.21 0.8 0.75])
% f3.Position = [1300 1200 1000 250];
% set(gca,'InnerPosition',[0.09 0.25 0.89 0.7])
hold all
plot(vs,LS_cost_of_inf_soft,'Color',mypurple)
plot(vs,ES_cost_of_inf_soft,'Color',myorange)
yline(LS_cost_of_inf_hard,'Color',mypurple,'LineStyle','-.','LineWidth',3,'Layer','bottom')
yline(ES_cost_of_inf_hard,'Color',myorange,'LineStyle','--','LineWidth',3,'Layer','bottom')
xticks(0:50:200)
axis([0 max(vs) 0.95*min(LS_cost_of_inf_hard, ES_cost_of_inf_hard) 1.05*max([LS_cost_of_inf_soft; ES_cost_of_inf_soft])])
% axis([0 max(vs) 2750 7750])
ax = gca;
ax.YAxis.Exponent = 0;
% xticks([vs(1:round(nv/5):end)])
xlabel('Soft constraint scaling $v$')
ylabel('Cost')
grid on


% save figures
if sensitivity_analysis
    labs = {['alpha0_',num2str(para.alpha(1))],['alpha0_',num2str(para.alpha(1))],['alpha1_',num2str(para.alpha(2))],['alpha1_',num2str(para.alpha(2))]};
    SAfigID = labs{sensitivity_analysis};
    saveas(f3,['./images/supplement/FS6_scen1_Costs_',SAfigID,'.png'])

else
    saveas(f1,'./images/F2A_scen1_Rt.png')
    saveas(f2,'./images/F2B_scen1_It.png')
    saveas(f3,'./images/F2C_scen1_Costs.png')
    
    % save incidence and prevalence for uncertain Hc script
    save('./mats/lockdowndilemma.mat',"outLS","outES","dv","dk","vs","ks","nv","nk","RIT","duration","inc_or_prev","stratnames","apply_discounting",'-mat')
end


% saveas(f,['./testing_images/F2_bedchange_',num2str(round(100*(bedincrease-1))),'.png'])
% saveas(f4,['./images/F2BC_bedchange_',num2str(round(100*(bedincrease-1))),'.png'])
% saveas(f5,['./images/F2BC_nok_bedchange_',num2str(round(100*(bedincrease-1))),'.png'])

%% Produce Fig 4D: Uncertainty in available beds

% we will assume the background number of beds occupied (not %) known
background = para.eta*para.Ibar;

% quantifying distribution on hospital capacity
Ibars = 1000:4000;

% define means and standard deviations
mu_lo = 0.8*para.Ibar*(1-para.eta);
mu_md = 1.0*para.Ibar*(1-para.eta);
mu_hi = 1.2*para.Ibar*(1-para.eta);

muidx = find(Ibars == round(mu_md));

sigs = 0:100:400;
mus = mu_md.*ones(size(sigs));
ndists = length(sigs);

% matrix of different normal pdfs
Ibar_dists = zeros(ndists,length(Ibars));
Ibar_dists(1,muidx) = 1;
CI95 = zeros(ndists,2);
CI95(1,:) = [mu_md mu_md];

for s = 2:ndists
    Ibar_dists(s,:) = discretise_distribution(Ibars,'Normal',mus(s),sigs(s));
    CI95(s,:) = Ibars([find(cumsum(Ibar_dists(s,:))>0.025,1,'first') find(cumsum(Ibar_dists(s,:))>0.975,1,'first')]);
end

% 95% CI for dists
CI95

% begin with matrices of the probability of Lockdown being preferable
% decision (lower cost) for each v and k under both constraints, append
% to it iteratively via the Ibar distribution
clear ES_softcosts ES_hardcosts LS_softcosts LS_hardcosts
[EC_LS_hard, EC_ES_hard, Prob_lockdown_hard] = deal(zeros(1,ndists));
[EC_LS_soft, EC_ES_soft, Prob_lockdown_soft] = deal(zeros(nv,ndists));
Costdiff = zeros(nv,length(Ibars));

tic
for H = 1:length(Ibars)
    if mod(H,100) == 0
        H
    end
    para.Ibar = Ibars(H) + background;

    % costs of control
    % LS_cost_of_control = 0;
    % ES_cost_of_control = ks.*duration*RIT^2.*0;
    
    % costs of disease (hard constraint)
    whichconstr = 'hard';
    
    LS_hardcosts = sum(compute_cost(outLS,para,whichconstr,vs(1),inc_or_prev,apply_discounting),2);
    ES_hardcosts = sum(compute_cost(outES,para,whichconstr,vs(1),inc_or_prev,apply_discounting),2);
    
    % costs of disease (soft constraint)
    whichconstr = 'soft';
    
    LS_softcosts = sum(compute_cost(outLS,para,whichconstr,vs,inc_or_prev,apply_discounting),2);
    ES_softcosts = sum(compute_cost(outES,para,whichconstr,vs,inc_or_prev,apply_discounting),2);


    % find v indices where Lockdown is optimal
    ES_optimal_hard = ES_hardcosts < LS_hardcosts;
    ES_optimal_soft = ES_softcosts < LS_softcosts;

    % cost difference (for supplementary figure)
    Costdiff(:,H) = LS_softcosts - ES_softcosts;

    % trivial calculation for first distribution
    if H == muidx
        EC_LS_hard(1)           = LS_hardcosts;
        EC_ES_hard(1)           = ES_hardcosts;
        EC_LS_soft(:,1)         = LS_softcosts;
        EC_ES_soft(:,1)         = ES_softcosts;
        Prob_lockdown_hard(1)   = ES_optimal_hard;
        Prob_lockdown_soft(:,1) = ES_optimal_soft;
    end

    % accumulate metrics using the pdf
    for s = 2:ndists

        % expected cost of each strategy
        EC_LS_hard(s) = EC_LS_hard(s) + LS_hardcosts.*Ibar_dists(s,H);
        EC_ES_hard(s)    = EC_ES_hard(s)    + ES_hardcosts.*Ibar_dists(s,H);
        EC_LS_soft(:,s) = EC_LS_soft(:,s) + LS_softcosts.*Ibar_dists(s,H);
        EC_ES_soft(:,s)    = EC_ES_soft(:,s)    + ES_softcosts.*Ibar_dists(s,H);

        % probability lockdown is optimal under each constraint
        Prob_lockdown_hard(s) = Prob_lockdown_hard(s) + ES_optimal_hard.*Ibar_dists(s,H);
        Prob_lockdown_soft(:,s) = Prob_lockdown_soft(:,s) + ES_optimal_soft.*Ibar_dists(s,H);
    end

end
toc

%% Plotting Fig 4D

f4 = figure(4);
f4.Position = [300 300 600 300];
set(gca,'InnerPosition',[0.16 0.21 0.8 0.75])
hold all

mysigmacols = [252,146,114; 251,106,74; 239,59,44; 203,24,29; 153,0,13]./255;
mysigmacols = [107,174,214; 66,146,198; 33,113,181; 8,81,156; 8,48,107]./255;
% mysigmacols = [27,158,119; 217,95,2; 102,166,30; 166,118,29; 102,102,102]./255;


for s = 1:ndists
    whatweplotting_soft = (EC_LS_soft(:,s) - EC_ES_soft(:,s));%./EC_ES_soft(:,s);
    whatweplotting_hard = (EC_LS_hard(s) - EC_ES_hard(s));%./EC_ES_hard(s);
    if s == 1
        dispname = '$\sigma_c \rightarrow 0$';
        myaxmin = 1.05*min(whatweplotting_hard,[],"all");
    else
        dispname = ['$\sigma_c = ',num2str(sigs(s)),'$'];
    end
    h(s) = plot(vs,whatweplotting_soft,'Color',mysigmacols(length(mysigmacols)-ndists+s,:),'DisplayName',dispname);
    yline(whatweplotting_hard,'Color',mysigmacols(length(mysigmacols)-ndists+s,:),'LineStyle',':','LineWidth',2,'Layer','bottom')
end
yline(0,'Color','k','LineStyle','-','LineWidth',2,'Layer','bottom')
% patch([0 0 max(vs) max(vs)], [0 1e4 1e4 0],myorange,'EdgeColor','none','FaceAlpha',0.1)
% patch([0 0 max(vs) max(vs)], [0 -1e4 -1e4 0],mypurple,'EdgeColor','none','FaceAlpha',0.1)
legend(h,'Location','northeast','AutoUpdate','off','FontSize',14.5,'NumColumns',3)
xticks(0:50:200)
% text(50,600,sprintf("\\alpha \\in %c", 120124),'Interpreter','tex')
axis([0 max(vs) myaxmin 1.5*max(whatweplotting_soft,[],"all")])
% axis([0 max(vs) -500 1000])
ax = gca;
ax.YAxis.Exponent = 0;
xlabel('Soft constraint scaling $v$')
ylab = ylabel('Expected cost difference');
% ylab = ylabel('\\tex[B][B]\{$E[C_{LS}(v)] - E[C_{ES}(v)]$\}','Interpreter','none');
% ylab = ylabel('\tex[B][B]{$E[C_{LS}(v)]$}','FontSize',15);
% ylab = ylabel('p2','FontSize',15,'Interpreter','none');
ylab.Position(2) = ylab.Position(2) - 100;
% title('$E[C_{LS}(v)] - E[C_{ES}(v)]$')
grid on


f5 = figure(5);
f5.Position = [100 1200 500 400];
set(gca,'InnerPosition',[0.175 0.15 0.75 0.75])
colormap(POcolormap([1:96 162:257],:))
hold all

[cont1,cont2] = contourf(Costdiff,-1500:500:1500,'-','LineWidth',1,'LabelColor','k','LabelSpacing',300);
clabel(cont1,cont2,'FontSize',0.9*16,'Interpreter','latex')
set(gca,'YDir','Normal')
% clim([0 100])
clim([-1500 1500])
% xline(Ibars(find(LS_overwhelm>0,1,'last')), 'Color',mypurple, 'LineStyle','--', 'LineWidth',4)
% xline(Ibars(find(ES_overwhelm>0,1,'last')), 'Color',myorange, 'LineStyle','--', 'LineWidth',4)
% for ticks, ticklabels, grid
vspec = (1:round(nv/20):nv); 
Hspec = (1:round(length(Ibars)/12):length(Ibars)); 
yticks(vspec(1:2:end))
xticks(Hspec(1:2:end))
yticklabels(vs(vspec(1:2:end)))
xticklabels(Ibars(Hspec(1:2:end)))
xtickangle(0)
% axis([Ibars(1) Ibars(end) vs(1) vs(end)])
ylim([1 501])
xlabel('Available ICU beds')
ylabel('Soft constraint scaling $v$')
title('Cost difference')


if plot_distributions

    f7 = figure(7);
    f7.Position = [1000 300 600 300];
    set(gca,'InnerPosition',[0.16 0.21 0.8 0.75])
    hold all

    for s = 2:ndists
        plot(Ibars, Ibar_dists(s,:), 'Color',mysigmacols(length(mysigmacols)-ndists+s,:), 'DisplayName', ['$\sigma_c = ',num2str(sigs(s)),'$'])
    end
    legend('Location','northeast','FontSize',16)
    axis([Ibars(1) Ibars(end) 0 0.004])
    xlabel('Available beds')
    ylabel('PMF')
    ax = gca;
    ax.YAxis.Exponent = 0;
    % title(['$\sigma_c = ',num2str(sigs(s)),'$']);
    grid on

    saveas(f7,'./images/supplement/FS3_ICUbeddists.png')

end    

if sensitivity_analysis
    labs = {['alpha0_',num2str(para.alpha(1))],['alpha0_',num2str(para.alpha(1))],['alpha1_',num2str(para.alpha(2))],['alpha1_',num2str(para.alpha(2))]};
    SAfigID = labs{sensitivity_analysis};
    saveas(f4,['./images/supplement/FS6_scen1_uncertainty_',SAfigID,'.png'])

else
    saveas(f4,'./images/F2D_scen1_uncertainICU.png')
    saveas(f5,'./images/supplement/FS4_scen1_uncertainICU.png')
    
end
