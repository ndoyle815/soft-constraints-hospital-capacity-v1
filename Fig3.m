% script to simulate an outbreak using the renewel equation transmission
% model
clear; close all

% Plotting preferences
set(0,'defaultlinelinewidth',3)
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex')
set(0,'defaultTextInterpreter','latex')
set(0,'defaultaxesfontsize',16)
figwidth = 750;  % width of Fig 3

% load colours
load('./mats/Cols.mat')

% load distribution parameters
para = load('./mats/Parameters.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TO PRODUCE SENSITIVITY ANALYSIS FIGURES, MODIFY THE BELOW VARIABLE
% KEY:
% 0: Fig 3
% 1: Fig S8  (alpha_0 = 0.2)
% 2: Fig S9  (alpha_0 = 0.5)
% 3: Fig S11 (alpha_1 = 5 alpha_0)
% 4: Fig S12 (alpha_1 = 2 alpha_0)
% 5: Fig S14 (with vs without discounting)
sensitivity_analysis = 5;

alpha0s = [0.2 0.5];
alpha1s = [para.alpha(1)*5 para.alpha(1)*2];
tSA = 0;  % (mis)timing of Strategy MS varies slightly with cost bounds 

if sensitivity_analysis == 1
    para.alpha(1) = alpha0s(1);
    tSA = -1;

elseif sensitivity_analysis == 2
    para.alpha(1) = alpha0s(2);
    tSA = 1;

elseif sensitivity_analysis == 3
    para.alpha(2) = alpha1s(1);

elseif sensitivity_analysis == 4
    para.alpha(2) = alpha1s(2);
    tSA = 1;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stratnames = {'Running Hot (RH)','Mistimed Suppression (MS)'};

% parameters
R0 = 2;            % R_0
I0 = 10;           % Initial no. cases
maxtime = 750;     % simulation time
whichR = {'Instantaneous','Case'};
Rtype = whichR{1};

% R_t timeseries
if isequal(Rtype,'Case')
    ti1 = [1 60 60+45];
    ti2 = [1 58 58+45];
else
    ts = 64;
    td = 8 + tSA;  % varies with sensitivity analysis
    tl = 150;
    ti1 = [1 ts ts+tl];
    ti2 = [1 ts+td ts+td+25];
end
Ri1 = [R0 1 0.7];
Ri2 = [R0 0.7 0.95];

% Discretise a gamma distributed serial interval distribution using the method of
% Cori et al.;
mu = 5.4;
sig = 1.5;
w = discretise_distribution(1:maxtime,'Gamma',mu,sig);

% run models
[outRH] = RENEWALmodel(Ri1,ti1,w,I0,maxtime,Rtype,para);
[outMS] = RENEWALmodel(Ri2,ti2,w,I0,maxtime,Rtype,para);


% schematic for presentation
Ri3 = [3.02 0.5:0.035:0.85 1 0.7];
ti3 = [1 40:50 51 230];
[out3] = RENEWALmodel(Ri3,ti3,w,I0,maxtime,Rtype,para);

f6 = figure(6);
f6.Position = [100 100 600 300];
set(gca,'InnerPosition',[0.1 0.2 0.85 0.69])
hold all
% patch([0 0 maxtime maxtime], [0 para.eta para.eta 0].*para.Ibar, 'k', 'FaceAlpha',0.4, 'EdgeColor','none', 'DisplayName','Background occupancy')
plot(out3.t, out3.in_ICU - para.eta*para.Ibar,'Color',myred)
yline(1.01*max(out3.in_ICU - para.eta*para.Ibar),'k--','Total beds','Interpreter','latex','LineWidth',2,'FontSize',16,'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','right','Layer','bottom')
xlabel('Time $t$')
title('In ICU');
axis([0 min([350 maxtime]) 0 1.1*max(out3.in_ICU - para.eta*para.Ibar)])
grid on
saveas(f6,'./runninghotex.png')

%% ICU PLOTTING AND COST EVALUATION

% plotting
f1 = figure(1);
f1.Position = [700 400 1000 250];
set(gca,'InnerPosition',[0.09 0.25 0.89 0.7])

hold all
plot(outRH.t,outRH.R,'Color',mygreen)
plot(outMS.t,outMS.R,'Color',myblue)
if figwidth > 700
    legend(stratnames,'Location','north','AutoUpdate','off','FontSize',16)
end
yline(1,'k--','$R(t)=1$','Interpreter','latex','LineWidth',2,'FontSize',16,'LabelVerticalAlignment','top','LabelHorizontalAlignment','right','Layer','bottom')
xlabel('Time $t$ (days)')
if isequal(Rtype,'Case')
    ylab = '$R(t)^c$';
else
    ylab = '$R(t)$';
end
ylabel(ylab,'Rotation',0);
axis([0 min([350 maxtime]) 0 1.2*max([outRH.R; outMS.R],[],"all")])
grid on


f2 = figure(2);
f2.Position = [700 400 1000 250];
set(gca,'InnerPosition',[0.09 0.25 0.89 0.7])
% f2.Position = [1300 1200 600 300];
% set(gca,'InnerPosition',[0.16 0.21 0.8 0.75])

hold all
patch([0 0 maxtime maxtime], [0 para.eta para.eta 0].*para.Ibar, 'k', 'FaceAlpha',0.4, 'EdgeColor','none', 'DisplayName','Background occupancy')
plot(outRH.t,outRH.in_ICU,'Color',mygreen)
plot(outMS.t,outMS.in_ICU,'Color',myblue)
if figwidth > 1250
    legend(stratnames,'Location','east','AutoUpdate','off','FontSize',16)
end
yline(para.Ibar,'k--','Total beds $\textrm{ICU}_c$','Interpreter','latex','LineWidth',2,'FontSize',16,'LabelVerticalAlignment','top','LabelHorizontalAlignment','right','Layer','bottom')
xlabel('Time $t$ (days)')
ylabel('In ICU');
axis([0 min([350 maxtime]) 0 1.3*para.Ibar])
grid on


% calculate costs under each constraint
vs = 10:0.2:210;
Nvs = length(vs);

% hard constraint
whichconstr = 'hard';
inc_or_prev = 2;
apply_discounting = 0;
RH_hardcosts = sum(compute_cost(outRH,para,whichconstr,vs(1),inc_or_prev,apply_discounting),2)';
MS_hardcosts = sum(compute_cost(outMS,para,whichconstr,vs(1),inc_or_prev,apply_discounting),2)';

% soft constraint
whichconstr = 'soft';
RH_softcosts = sum(compute_cost(outRH,para,whichconstr,vs,inc_or_prev,apply_discounting),2)';
MS_softcosts = sum(compute_cost(outMS,para,whichconstr,vs,inc_or_prev,apply_discounting),2)';

% compute threshold v
[~,minidx] = min(abs(RH_softcosts - MS_softcosts));
vs(minidx)
[sum(outRH.new_ICU) sum(outMS.new_ICU)]./1e3

% plotting soft constraint costs as a function of scaling v
f3 = figure(3);
f3.Position = [700 400 1000 250];
set(gca,'InnerPosition',[0.09 0.25 0.89 0.7])
% f3.Position = [1300 1200 600 300];
% set(gca,'InnerPosition',[0.16 0.21 0.8 0.75])
hold all

plot(vs,RH_softcosts,'Color',mygreen)
plot(vs,MS_softcosts,'Color',myblue)
yline(RH_hardcosts,'Color',mygreen,'LineStyle','-.','LineWidth',3,'Layer','bottom')
yline(MS_hardcosts,'Color',myblue,'LineStyle','--','LineWidth',3,'Layer','bottom')
xlim([0 max(vs)])
ax = gca;
ax.YAxis.Exponent = 0;
xlabel('Soft constraint scaling $v$')
ylabel('Cost')
xticks(0:50:200)
if figwidth > 1250
    legend(stratnames,'Location','northeast','AutoUpdate','off','FontSize',16)
end
% title(['$\alpha_0 = ', num2str(para.alpha(1)), ', \alpha_1 = 1$'])
if sensitivity_analysis == 0
    ylim([8e3 15e3])
    yticks(0:3e3:15e3)
end
grid on


% meanwhile, what's going on with hospitals?
f4 = figure(4);
f4.Position = [700 400 1000 250];
set(gca,'InnerPosition',[0.09 0.25 0.89 0.6])

hold all
ax = gca;
ax.YAxis.Exponent = 0;
patch([0 0 maxtime maxtime], [0 para.etaH para.etaH 0].*para.Hbar, 'k', 'FaceAlpha',0.4, 'EdgeColor','none', 'DisplayName','Background occupancy')
plot(outRH.t,outRH.in_Hosp,'Color',mygreen)
plot(outMS.t,outMS.in_Hosp,'Color',myblue)
% legend(stratnames,'Location','east','AutoUpdate','off','FontSize',16)
yline(para.Hbar,'k--','Total beds $\textrm{H}_c$','Interpreter','latex','LineWidth',2,'FontSize',16,'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','right','Layer','bottom')
xlabel('Time $t$ (days)')
title('In hospital, non-ICU');
axis([0 min([350 maxtime]) round(para.etaH*para.Hbar,-4) round(para.Hbar,-4)])
grid on

% discounting
if sensitivity_analysis == 5
    apply_discounting = 1;
    RH_discounting = sum(compute_cost(outRH,para,whichconstr,vs,inc_or_prev,apply_discounting),2)';
    MS_discounting = sum(compute_cost(outMS,para,whichconstr,vs,inc_or_prev,apply_discounting),2)';
    
    
    % plotting soft constraint costs as a function of scaling v
    f14 = figure(14);
    f14.Position = [700 400 1000 250];
    set(gca,'InnerPosition',[0.09 0.25 0.89 0.7])
    hold all
    
    plot(vs,RH_softcosts,'Color',mygreen,'LineWidth',1,'LineStyle','-','DisplayName','RH, no discounting')
    plot(vs,MS_softcosts,'Color',myblue,'LineWidth',1,'LineStyle','-','DisplayName','MS, no discounting')
    plot(vs,RH_discounting,'Color',mygreen,'LineStyle','--','DisplayName','RH, discounting')
    plot(vs,MS_discounting,'Color',myblue,'LineStyle','--','DisplayName','MS, discounting')
    xlim([0 max(vs)])
    ax = gca;
    ax.YAxis.Exponent = 0;
    xlabel('Soft constraint scaling $v$')
    ylabel('Cost')
    legend('Location','northeast','FontSize',16,'NumColumns',2)
    ylim([8e3 15e3])
    yticks(0:3e3:15e3)
    grid on
    
    saveas(f14,'./images/supplement/FS14_discounting.png')

end

% save figures
if sensitivity_analysis > 0 && sensitivity_analysis < 5
    labs = {['alpha0_',num2str(para.alpha(1))],['alpha0_',num2str(para.alpha(1))],['alpha1_',num2str(para.alpha(2))],['alpha1_',num2str(para.alpha(2))]};
    SAfigID = labs{sensitivity_analysis};
    saveas(f1,['./images/supplement/FS8_scen2_Rt_',SAfigID,'.png'])
    saveas(f2,['./images/supplement/FS8_scen2_It_',SAfigID,'.png'])
    saveas(f3,['./images/supplement/FS8_scen2_Costs_',SAfigID,'.png'])

elseif sensitivity_analysis == 0
    saveas(f1,'./images/F3A_scen2_Rt.png')
    saveas(f2,'./images/F3B_scen2_It.png')
    saveas(f3,'./images/F3C_scen2_Costs.png')
    saveas(f4,'./images/supplement/FS2_hospitaloccupancy.png')

    % save results for uncertainty analysis
    save('./mats/results.mat',"outRH","outMS","vs","Nvs","RH_softcosts","MS_softcosts","inc_or_prev","stratnames",'-mat')

end



%% Plot regret

% Regretcost = zeros(Nvs);
% 
% % RH_Inferredcost = repmat(RH_softcosts,Nvs,1);
% % RH_Truecost = repmat(RH_softcosts',1,Nvs);
% % MS_Inferredcost = repmat(MS_softcosts,Nvs,1);
% % MS_Truecost = repmat(MS_softcosts,Nvs,1);
% 
% % [inferred, actual] = deal(zeros(Nvs));
% 
% for i = 1:Nvs
%     for j = 1:Nvs
%         % what do we infer?
%         inferred = RH_softcosts(i) - MS_softcosts(i);
%         actual = RH_softcosts(j) - MS_softcosts(j);
%         % if inferred > 0 and actual > 0: we made the right decision (MS optimal). No regret
%         % if inferred > 0 and actual < 0: we made the wrong decision (RH
%         % actually optimal). Regret is difference (add absolute values)
%         % if inferred < 0 and actual > 0: we made the wrong decision (MS
%         % actually optimal). Regret is difference (add absolute values)
%         % if inferred < 0 and actual < 0: we made the right decision (RH optimal). No regret
%         if ~isequal(sign(inferred),sign(actual))
%             Regretcost(i,j) = actual;
%         end
%     end
% end
% 
% vmarks = 1:round(Nvs/5):Nvs;
% contvec = horzcat(-6000:1000:-2000, -1500:500:1500, 2000:1000:6000);
% % contvec = horzcat(0:500:1000, 2000:1000:6000);
% decision_boundary = find(RH_softcosts>MS_softcosts,1,"last");
% arr = 1:13:110;
% 
% mycolormap = BGcolormap(horzcat(arr, fliplr(repmat(258,1,length(arr)) - arr)),:);
% 
% f5 = figure(5);
% f5.Position = [800 1000 500 400];
% set(gca,'InnerPosition',[0.15 0.125 0.75 0.75])
% hold all
% 
% colormap(BGcolormap)
% [cont1,cont2] = contourf(Regretcost,contvec,'-','LineWidth',1,'LabelColor','k','LabelSpacing',50);
% clabel(cont1,cont2,'FontSize',0.9*16,'Interpreter','latex')
% % imagesc(Regretcost)
% patch([0 0 decision_boundary decision_boundary], [0 decision_boundary decision_boundary 0], 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3)
% patch([decision_boundary decision_boundary Nvs Nvs], [decision_boundary  Nvs Nvs decision_boundary], 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3)
% xline(decision_boundary, 'r-', '$v_c$','LineWidth',2,'Interpreter','latex','FontSize',16)
% yline(decision_boundary, 'r-', '$v_c$','LineWidth',2,'Interpreter','latex','FontSize',16)
% set(gca,'YDir','Normal')
% clim([-6000 6000])
% xlabel('Inferred $v$')
% ylabel('True $v$')
% xticks(vmarks)
% yticks(vmarks)
% xticklabels(vs(vmarks))
% yticklabels(vs(vmarks))
% title('Regret')
% 
% colorbar(gca)