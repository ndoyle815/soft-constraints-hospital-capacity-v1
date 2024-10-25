% script to test whether the two constraints lead to a different decision
% on the "Lockdown" vs "No Lockdown" question
clear; close all

% Plotting preferences
set(0,'defaultlinelinewidth',3)
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex')
set(0,'defaultTextInterpreter','latex')
set(0,'defaultaxesfontsize',16)
% colours
wcols = [65,182,196; 29,145,192; 34,94,168; 37,52,148; 8,29,88]./255;
myred = [228 26 28]./255;
myblue = [55 126 184]./255;
mygreen = [77 175 74]./255;
mypurple = [152 78 163]./255;
myorange = [255 127 0]./255;

% load distribution parameters
para = load('./mats/Parameters.mat');

% load results
load('./mats/results.mat')

% random variable is vhat = 1/10*v
vhats = 10./vs;

% mean and variance on v vectors
alphas = 25/(Nvs-1):25/(Nvs-1):10;
betas  = (25/(Nvs-1):25/(Nvs-1):10) + 1;
Na = length(alphas);
Nb = length(betas);

% store cost statistics
[EC1, EC2, EsqC1, EsqC2, Prob_S2opt] = deal(zeros(Na,Nb));
vv_Stratboundary = zeros(1,Nb);
vv_StratProbboundary = zeros(1,Nb);

% generate beta distributions
% betadists = zeros(Na,Nb,Nvs);
% tic;
% for v = 1:Nvs
%     betadists(:,:,v) = betapdf(vhats(v), repmat(alphas,Na,1), repmat(betas,Nb,1));
% end
% toc;

% calculate soft costs
tic;
for b = 1:Nb
    if mod(b,100) == 0
        b
    end

    for a = 1:Na

       % generate the beta distribution
       tmp_betadist = betapdf(vhats,alphas(a),betas(b));
       tmp_betadist = tmp_betadist./sum(tmp_betadist);  % normalised

       % expected cost of each strategy
       EC1(a,b) = sum(C1_softcosts.*tmp_betadist);
       EC2(a,b) = sum(C2_softcosts.*tmp_betadist);

       % second moment cost of each strategy
       EsqC1(a,b) = sum((C1_softcosts.^2).*tmp_betadist);
       EsqC2(a,b) = sum((C2_softcosts.^2).*tmp_betadist);

       % compute probability that Strategy 1 is optimal
       where_S2opt = C2_softcosts < C1_softcosts;
       Prob_S2opt(a,b) = sum(tmp_betadist(where_S2opt));
   
    end

    vv_Stratboundary(b) = min([find(EC2(:,b) < EC1(:,b),1,'first') Na+1]);
    vv_StratProbboundary(b) = min([find(Prob_S2opt(:,b) > 0.5,1,'first') Na+1]);

end
toc;

% variance of cost for each strategy
VarC1 = EsqC1 - EC1.^2;
VarC2 = EsqC2 - EC2.^2;

%% Plotting
mu_vhat = 10/68.4;

% min and max expected costs for colour range
Emin = min([EC1; EC2],[],'all');
Emax = max([EC1; EC2],[],'all');

f1 = figure(1);
f1.Position = [100 1000 450 450];

Emin = min((EC2 - EC1)./EC1,[],'all');
Emax = max((EC2 - EC1)./EC1,[],'all');
Ediff = Emax - Emin;

mywhite = [1 1 1];

dcol1 = round(256*(-Emin/(Emax - Emin)));
dcol1 = 128;
colorvec1 = (0:1/dcol1:1)';
mycolormap1 = myblue.*(1 - colorvec1) + mywhite.*colorvec1;

dcol2 = 256 - dcol1;
dcol2 = 128;
colorvec2 = (0:1/dcol2:1)';
mycolormap2 = mywhite.*(1 - colorvec2) + mygreen.*colorvec2;
mycolormap = [mycolormap1; mycolormap2(2:end,:)];
colormap(mycolormap)

imagesc((EC2 - EC1)./EC1)
hold on
plot(1:Nb, vv_Stratboundary, 'r')
hold on
% plot(find(ismember((mu_vhat.*betas)./(1 - mu_vhat),alphas)), find(ismember(alphas,(mu_vhat.*betas)./(1 - mu_vhat))), 'k')
% text(800,1000,'$\mu_v = 3 \sigma_v$','FontSize',16,'Color','k')
subax = gca;
% subax.Position(1) = subax.Position(1) - 0.04;
set(subax,'YDir','normal')
clim([-Emax Emax])
ylab = ylabel('$\alpha_v$','Rotation',0);
xlabel('$\beta_v$')
yticks(round(Na/10):round(Na/10):Na)
xticks(round(Nb/10):round(Nb/10):Nb)
yticklabels(alphas(round(Na/10):round(Na/10):Na))
xticklabels(betas(round(Nb/10):round(Nb/10):Nb))
xtickangle(0)
tl = title('$(E[C_{\mathcal{S}_2}(v)] - E[C_{\mathcal{S}_1}(v)])/E[C_{\mathcal{S}_1}(v)]$');
% tl.Position(2) = tl.Position(2) + 75;

h = axes(gcf,'visible','off'); 
h.Title.Visible = 'on';
h.XLabel.Visible = 'on';
h.YLabel.Visible = 'on';

c = colorbar(h,'Position',[0.86 0.15 0.02 0.68],'FontSize',16,'TickLabelInterpreter','Latex');  % attach colorbar to h
colormap(c);
clim(h, [-Emax Emax]);

set(subax,'InnerPosition',[0.16 0.13 0.66 0.73])

% saveas(f1,'./images/F4C_uncertain_v_expcost.png')


f2 = figure(2);
f2.Position = [900 1000 450 450];

% Emin = min((EC1 - EC2)./EC2,[],'all');
% Emax = max((EC1 - EC2)./EC2,[],'all');
% Ediff = Emax - Emin;
% 
% mywhite = [1 1 1];
% 
% dcol1 = round(256*(-Emin/(Emax - Emin)));
% dcol1 = 128;
% colorvec1 = (0:1/dcol1:1)';
% mycolormap1 = myblue.*(1 - colorvec1) + mywhite.*colorvec1;
% 
% dcol2 = 256 - dcol1;
% dcol2 = 128;
% colorvec2 = (0:1/dcol2:1)';
% mycolormap2 = mywhite.*(1 - colorvec2) + mygreen.*colorvec2;
% mycolormap = [mycolormap1; mycolormap2(2:end,:)];
colormap(mycolormap(end:-1:1,:))

imagesc(Prob_S2opt)
hold on
plot(1:Nb, vv_StratProbboundary, 'r')
hold on
% plot(find(ismember(3.*vstds,vmeans)), find(ismember(vmeans,3.*vstds)), 'k')
% text(800,1000,'$\mu_v = 3 \sigma_v$','FontSize',16,'Color','k')
subax = gca;
% subax.Position(1) = subax.Position(1) - 0.04;
set(subax,'YDir','normal')
clim([0 1])
ylab = ylabel('$\alpha_v$','Rotation',0);
xlabel('$\beta_v$')
yticks(round(Na/10):round(Na/10):Na)
xticks(round(Nb/10):round(Nb/10):Nb)
yticklabels(alphas(round(Na/10):round(Na/10):Na))
xticklabels(betas(round(Nb/10):round(Nb/10):Nb))
xtickangle(0)
tl = title('$P[C_{\mathcal{S}_2}(v) < C_{\mathcal{S}_1}(v)]$');
% tl.Position(2) = tl.Position(2) + 75;

h = axes(gcf,'visible','off'); 
h.Title.Visible = 'on';
h.XLabel.Visible = 'on';
h.YLabel.Visible = 'on';

c = colorbar(h,'Position',[0.86 0.15 0.02 0.68],'FontSize',16,'TickLabelInterpreter','Latex');  % attach colorbar to h
colormap(c);
clim(h, [0 1]);

set(subax,'InnerPosition',[0.16 0.13 0.66 0.73])

% saveas(f2,'./images/F4D_uncertain_v_proboopt.png')
