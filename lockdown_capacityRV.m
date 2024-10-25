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

% load default parameters
if not(isfolder('mats'))
    disp('No parameters saved: Running define_params.m')
    define_params
end
para = load('./mats/Parameters.mat');

% obtain data from lockdown dilemma, capacity known
load('./mats/lockdowndilemma.mat');

% we will assume the background number of beds occupied (not %) known
background = para.eta*para.Ibar;

% quantifying distribution on hospital capacity
dH = para.Ibar/750;
Ibars = (0.4*para.Ibar:dH:1.6*para.Ibar).*(1-para.eta);

sigmax = round(para.Ibar*(1-para.eta)/6,-2);
sigs = 0:sigmax/20:sigmax;
sigs(1) = sigs(2)/2;
nsig = length(sigs);

% matrix of different normal pdfs
Ibar_dists = zeros(length(sigs),length(Ibars));
for s = 1:nsig
    Ibar_dists(s,:) = normpdf(Ibars,para.Ibar*(1-para.eta),sigs(s));
end

% ... plus new distributions on hospital capacity
sigma_Ibar = 250;
mu_lo = 0.8*para.Ibar*(1-para.eta);
mu_md = 1.0*para.Ibar*(1-para.eta);
mu_hi = 1.2*para.Ibar*(1-para.eta);

Ibars_dist_lo = normpdf(Ibars,mu_lo,sigma_Ibar);  % underestimate
Ibars_dist_md = normpdf(Ibars,mu_md,sigma_Ibar);  % correct mean
Ibars_dist_hi = normpdf(Ibars,mu_hi,sigma_Ibar);  % overestimate

mus = [mu_md.*ones(1,length(sigs)) mu_lo mu_md mu_hi];
sigs = [sigs sigma_Ibar sigma_Ibar sigma_Ibar];
Ibar_dists = [Ibar_dists; Ibars_dist_lo; Ibars_dist_md; Ibars_dist_hi];

Ibar_dists = Ibar_dists./sum(Ibar_dists,2);  % normalisation
ndists = size(Ibar_dists,1);

% 95% CI for dist 23 (to be included as Fig 2d)
Ibars([find(cumsum(Ibar_dists(end-1,:))>0.025,1,'first') find(cumsum(Ibar_dists(end-1,:))>0.975,1,'first')])

% begin with matrices of the probability of Lockdown being preferable
% decision (lower cost) for each v and k under both constraints, append
% to it iteratively via the Ibar distribution
[EC_NO_LOCKDOWN_hard, EC_NO_LOCKDOWN_soft, EC_LOCKDOWN_hard, EC_LOCKDOWN_soft, Prob_lockdown_hard, Prob_lockdown_soft] = deal(zeros(nv,nk,ndists));

tic
for H = 1:length(Ibars)
    if mod(H,100) == 0
        H
    end
    para.Ibar = Ibars(H) + background;

    % costs of control
    NO_LOCKDOWN_cost_of_control = 0;
    LOCKDOWN_cost_of_control = ks.*duration*RIT^2;
    
    % costs of disease (hard constraint)
    whichconstr = 'hard';
    
    NO_LOCKDOWN_cost_of_inf_hard = sum(compute_cost(out1,para,whichconstr,vs(1),inc_or_prev,indirect),2);
    NO_LOCKDOWN_hardcosts = repmat(NO_LOCKDOWN_cost_of_inf_hard + NO_LOCKDOWN_cost_of_control, nv, nk);
    
    LOCKDOWN_cost_of_inf_hard = sum(compute_cost(out2,para,whichconstr,vs(1),inc_or_prev,indirect),2);
    LOCKDOWN_hardcosts = repmat(LOCKDOWN_cost_of_inf_hard,1,nk) + repmat(LOCKDOWN_cost_of_control,nv,1);
    
    % costs of disease (soft constraint)
    whichconstr = 'soft';
    
    NO_LOCKDOWN_cost_of_inf_soft = sum(compute_cost(out1,para,whichconstr,vs,inc_or_prev,indirect),2);
    NO_LOCKDOWN_softcosts = repmat(NO_LOCKDOWN_cost_of_inf_soft + NO_LOCKDOWN_cost_of_control,1,nk);
    
    LOCKDOWN_cost_of_inf_soft = sum(compute_cost(out2,para,whichconstr,vs,inc_or_prev,indirect),2);
    LOCKDOWN_softcosts = repmat(LOCKDOWN_cost_of_inf_soft,1,nk) + repmat(LOCKDOWN_cost_of_control,nv,1);


    % find v and k indices where Lockdown is optimal
    LOCKDOWN_optimal_hard = LOCKDOWN_hardcosts < NO_LOCKDOWN_hardcosts;
    LOCKDOWN_optimal_soft = LOCKDOWN_softcosts < NO_LOCKDOWN_softcosts;

    % accumulate metrics using the pdf
    for s = 1:ndists

        % expected cost of each strategy
        EC_NO_LOCKDOWN_hard(:,:,s) = EC_NO_LOCKDOWN_hard(:,:,s) + NO_LOCKDOWN_hardcosts.*Ibar_dists(s,H);
        EC_LOCKDOWN_hard(:,:,s)    = EC_LOCKDOWN_hard(:,:,s)    + LOCKDOWN_hardcosts.*Ibar_dists(s,H);
        EC_NO_LOCKDOWN_soft(:,:,s) = EC_NO_LOCKDOWN_soft(:,:,s) + NO_LOCKDOWN_softcosts.*Ibar_dists(s,H);
        EC_LOCKDOWN_soft(:,:,s)    = EC_LOCKDOWN_soft(:,:,s)    + LOCKDOWN_softcosts.*Ibar_dists(s,H);

        % probability lockdown is optimal under each constraint
        Prob_lockdown_hard(:,:,s) = Prob_lockdown_hard(:,:,s) + LOCKDOWN_optimal_hard.*Ibar_dists(s,H);
        Prob_lockdown_soft(:,:,s) = Prob_lockdown_soft(:,:,s) + LOCKDOWN_optimal_soft.*Ibar_dists(s,H);
    end

end
toc

% save tensor for analysis
save('./mats/lockdowndilemma.mat',"EC_NO_LOCKDOWN_hard","EC_LOCKDOWN_hard","EC_NO_LOCKDOWN_soft","EC_LOCKDOWN_soft",...
     "Prob_lockdown_hard","Prob_lockdown_soft","Ibar_dists","Ibars","mus","sigs","background",'-mat','-append')

%% PLOTTING COST SURFACES

% 1,6,21: correct mean, different sigma
% 22,23,24: under/correct/over estimate mean, specific sigma
distidxs = [1 6 21 22 23 24];

figno = 1;
for nd = distidxs

    % plotting distribution PDF
    figure(figno);
    set(gcf,'Position',[125*(figno-1) 1200 550 150])
    
    plot(Ibars, Ibar_dists(nd,:), 'r', 'DisplayName', ['$\mu = ',num2str(round(mus(nd))),', \sigma^2 = ',num2str(sigs(nd)),'$'])
    xlim([0 Ibars(end)])
    xlabel('Capacity $H_c$')
    ylabel('PDF')
    title('Distribution $H_c \sim \mathcal{N}(\mu, \sigma)$')
    legend('Location','west')
    grid on
    
    saveas(gcf, ['./images/LKdilemma_mu',num2str(round(mus(nd))),'_sigma',num2str(sigs(nd)),'_Ibardist.png'])
    

    % plotting Prob[Lockdown optimal]
    figure(figno+1);
    set(gcf,'Position',[125*(figno-1) 200 550 450]);
    colormap(POcolormap)
    
    imagesc(ks,vs,Prob_lockdown_soft(:,:,nd))
    set(gca,'OuterPosition',[0.01 0.01 0.9 0.94])
    set(gca,'YDir','Normal')
    clim([0 1])
    yticks(vs(1:round(nv/10):nv))
    xticks(ks(1:round(nk/10):nk))
    yticklabels(vs(1:round(nv/10):nv))
    xticklabels(ks(1:round(nk/10):nk))
    xtickangle(0)
    xlabel('$k$ scaling cost of Strategy ES')
    ylabel('$v$ scaling logistic constraint')

    h = axes(gcf,'visible','off'); 
    h.Title.Visible = 'on';
    h.XLabel.Visible = 'on';
    h.YLabel.Visible = 'on';
    
    c = colorbar(h,'Position',[0.86 0.15 0.02 0.68],'FontSize',16,'TickLabelInterpreter','Latex');  % attach colorbar to h
    colormap(c);
    clim(h,[0,1]);
    c.Ticks = 0:0.2:1;
    c.Label.String = "Probability Strategy ES Optimal";
    c.Label.Rotation = 270;
    c.Label.VerticalAlignment = "bottom";
    c.Label.Interpreter = "latex";
    
    if nd == distidxs(5)
        tl = title('Optimal, Soft Constraint');
        tl.Position(2) = 0.97;
        tl.Position(1) = 0.45;
        saveas(gcf,strcat('./images/F2D_lockdowndilemma_Ibaruncertainty_sigma',num2str(sigs(nd)),'.png'))
    else
        title(['$\mu = ',num2str(round(mus(nd))),', \sigma^2 = ',num2str(sigs(nd)),'$']);
        saveas(gcf,strcat('./images/LKdilemma_mu',num2str(round(mus(nd))),'_sigma',num2str(sigs(nd)),'_ProbLockOpt.png'))
    end


    figno = figno + 2;

end

%% Plotting difference in expected costs of each strategy
close all

figno = 1;
for nd = distidxs

    for hs = 0:1

        % plotting difference in expected costs (hard constraint)
        if hs == 0
            ECDIFF = EC_LOCKDOWN_hard(:,:,nd) - EC_NO_LOCKDOWN_hard(:,:,nd);
            plotname = ['./images/LKdilemma_mu',num2str(round(mus(nd))),'_sigma',num2str(sigs(nd)),'_Ecosthard.png'];
        else
            ECDIFF = EC_LOCKDOWN_soft(:,:,nd) - EC_NO_LOCKDOWN_soft(:,:,nd);
            plotname = ['./images/LKdilemma_mu',num2str(round(mus(nd))),'_sigma',num2str(sigs(nd)),'_Ecostsoft.png'];
        end
        Emax = max(ECDIFF,[],"all");
        Emin = min(ECDIFF,[],"all");
        Eabs = max(-Emin,Emax);
    
        figure(2*length(distidxs)+figno+hs);
        set(gcf,'Position',[125*(figno-1) 1200-1000*hs 500 450]);
        colormap(POcolormap(end:-1:1,:))
        
        imagesc(ks,vs,ECDIFF)
        set(gca,'OuterPosition',[0.01 0.01 0.9 0.94])
        set(gca,'YDir','Normal')
        clim([-Eabs Eabs])
        yticks(vs(1:round(nv/10):nv))
        xticks(ks(1:round(nk/10):nk))
        yticklabels(vs(1:round(nv/10):nv))
        xticklabels(ks(1:round(nk/10):nk))
        xtickangle(0)
        xlabel('$k$ scaling cost of Lockdown')
        ylabel('$v$ scaling logistic constraint')
        title(['$\mu = ',num2str(round(mus(nd))),', \sigma^2 = ',num2str(sigs(nd)),'$']);
        
        h = axes(gcf,'visible','off'); 
        h.Title.Visible = 'on';
        h.XLabel.Visible = 'on';
        h.YLabel.Visible = 'on';
        
        c = colorbar(h,'Position',[0.86 0.15 0.02 0.68],'FontSize',16,'TickLabelInterpreter','Latex');  % attach colorbar to h
        colormap(c);
        clim(h,[-Eabs Eabs]);
        % c.Ticks = 0:0.2:1;
        c.Label.Interpreter = "latex";
        c.Label.String = 'Difference in $E[C_{\mathcal{S}}(v)]$';
        c.Label.Rotation = 270;
        c.Label.VerticalAlignment = "bottom";

        saveas(gcf,plotname)

    end

    figno = figno + 2;

end
