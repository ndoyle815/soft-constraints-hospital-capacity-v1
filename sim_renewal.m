% script to simulate an outbreak using the renewel equation transmission
% model
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

% parameters
R0 = 2;            % R_0
I0 = 10;           % Initial no. cases
maxtime = 1500;     % simulation time
whichR = {'Instantaneous','Case'};
Rtype = whichR{1};

% R_t timeseries
if isequal(Rtype,'Case')
    ti1 = [1 60 60+45];
    ti2 = [1 58 58+45];
else
    ts = 56;
    td = 6;
    tp = 8;
    tl = 120;
    ti1 = [1 ts+tp    ts+tp+tl+60];
    ti2 = [1 ts+td+tp ts+td+tp+tl/4];
end
Ri1 = [R0 1 0.8];
Ri2 = [R0 0.8 0.95];

% serial interval distribution
mu = 5.4;
sig = 1.5;
wtimes = 1:maxtime;
w = gampdf(wtimes,(mu/sig)^2,1/(mu/sig^2));
w = w./sum(w);

% run models
[out1] = RENEWALmodel(Ri1,ti1,w,I0,maxtime,Rtype,para);
[out2] = RENEWALmodel(Ri2,ti2,w,I0,maxtime,Rtype,para);

sum(out1.new_ICU)
sum(out2.new_ICU)

%% ICU PLOTTING AND COST EVALUATION

% plotting
f7 = figure(7);
f7.Position = [100 100 450 300];

bar(wtimes,w)
xlabel('Time $t$ (days)')
ylabel('Probability $s = t$')
title('Serial interval distribution $w_s$')
axis([0 20 0 1.2*max(w)])
grid on

saveas(f7,'../images/F3D_renewal_wt.png')

f1 = figure(1);
f1.Position = [100 100 900 300];

hold all
plot(out1.t,out1.R,'Color',mygreen)
plot(out2.t,out2.R,'Color',myblue)
legend({'Strategy $\mathcal{S}_1$','Strategy $\mathcal{S}_2$'},'Location','northeast')
xlabel('Time $t$ (days)')
if isequal(Rtype,'Case')
    ylab = '$R_t^c$';
else
    ylab = '$R_t$';
end
yl = ylabel(ylab,'Rotation',0);
%yl.Position(1) = yl.Position(1) - 0.0;
title('Renewal simulations')
axis([0 min([540 maxtime]) 0 1.2*max(max([out1.R; out2.R]))])
grid on

saveas(f1,'../images/F3A_renewal_Rt.png')

f2 = figure(2);
f2.Position = [100 100 900 300];

hold all
% plot(alltimes,ExpIts(2,:),'Color',mygreen)
% plot(alltimes,ExpIts(1,:),'Color',myblue)
plot(out1.t,out1.in_ICU,'Color',mygreen)
plot(out2.t,out2.in_ICU,'Color',myblue)
legend({'Strategy $\mathcal{S}_1$','Strategy $\mathcal{S}_2$'},'Location','east','AutoUpdate','off')
%bar(alltimes,Its(1,:),'FaceColor',myblue,'FaceAlpha',0.4,'BarWidth',1)
%bar(alltimes,Its(2,:),'FaceColor',mygreen,'FaceAlpha',0.4,'BarWidth',1)
yline((1 - para.eta)*para.Ibar,'k--','Capacity','Interpreter','latex','LineWidth',2,'FontSize',16,'LabelVerticalAlignment','top','LabelHorizontalAlignment','right')
xlabel('Time $t$ (days)')
% yl = ylabel({'In','ICU'},'Rotation',0);
yl = ylabel('In ICU');
%yl.Position(1) = yl.Position(1) - 0.0;
xlim([0 min([540 maxtime])])
grid on

saveas(f2,'../images/F3B_renewal_It.png')

% changing alpha around
alphachoices = [0.1 0.3456];
whichalpha = 2;
para.alpha(1) = alphachoices(whichalpha);

% calculate costs under each constraint
vs = (1:0.01:21)'.*10;
Nvs = length(vs);

% hard constraint
whichconstr = 'hard';
inc_or_prev = 0;
indirect = 0;
C1_hardcosts = compute_cost(out1,para,whichconstr,vs(1),inc_or_prev,indirect);
C2_hardcosts = compute_cost(out2,para,whichconstr,vs(1),inc_or_prev,indirect);

% soft constraint
whichconstr = 'soft';
tic;
C1_softcosts = sum(compute_cost(out1,para,whichconstr,vs,2,indirect),2)';
C2_softcosts = sum(compute_cost(out2,para,whichconstr,vs,2,indirect),2)';
toc;
tic;
C1_softcosts1 = compute_cost(out1,para,whichconstr,vs,1,indirect);
C2_softcosts1 = compute_cost(out2,para,whichconstr,vs,1,indirect);
toc;
tic;
C1_softcosts2 = compute_cost(out1,para,whichconstr,vs,2,indirect);
C2_softcosts2 = compute_cost(out2,para,whichconstr,vs,2,indirect);
toc;

% plotting (cost of both curves as a function of logistic scaling)
vidx = 2000;
figure(8)
subplot(1,2,1)
plot(1:maxtime,C1_softcosts1(vidx,:),'k:', 1:maxtime,C1_softcosts2(vidx,:),'r')

subplot(1,2,2)
plot(1:maxtime,C2_softcosts1(vidx,:),'k:', 1:maxtime,C2_softcosts2(vidx,:),'r')

f3 = figure(3);
f3.Position = [600 600 900 300];

plot(vs,C1_softcosts,'Color',mygreen)
hold on
plot(vs,C2_softcosts,'Color',myblue)
xlim([0 max(vs)])
ax = gca;
ax.YAxis.Exponent = 0;
xlabel('Scaling of logistic curve $w$')
yl = ylabel('Cost');
% yl.Position(1) = yl.Position(1) - 0.15;
legend({'Strategy $\mathcal{S}_1$','Strategy $\mathcal{S}_2$'},'Location','northeast')
title(['$\alpha_0 = ', num2str(para.alpha(1)), ', \alpha_1 = 1$'])
cmax = 16e3;
ylim([0 max([C1_softcosts(1) cmax])])
yticks(0:3e3:cmax)
grid on

saveas(f3,['../images/F3C_renewal_costs_alpha0_',num2str(para.alpha(1)),'.png'])

% save results for uncertainty analysis
vs = vs';
save('./mats/results.mat',"out1","out2","vs","Nvs","C1_softcosts","C2_softcosts")

%% HOSPITAL PLOTTING AND COST EVALUATION
Hbar1 = 9000;

% plotting
f5 = figure(5);
f5.Position = [100 100 900 300];

hold all
% plot(alltimes,ExpIts(2,:),'Color',mygreen)
% plot(alltimes,ExpIts(1,:),'Color',myblue)
plot(out1.t,out1.in_Hosp,'Color',mygreen)
plot(out2.t,out2.in_Hosp,'Color',myblue)
legend({'R = (2,1,0.7)','R = (2,0.7,0.9)'},'Location','east','AutoUpdate','off')
%bar(alltimes,Its(1,:),'FaceColor',myblue,'FaceAlpha',0.4,'BarWidth',1)
%bar(alltimes,Its(2,:),'FaceColor',mygreen,'FaceAlpha',0.4,'BarWidth',1)
yline(Hbar1,'k--','Capacity','Interpreter','latex','LineWidth',2,'FontSize',16,'LabelVerticalAlignment','top','LabelHorizontalAlignment','right')
xlabel('Time $t$ (days)')
yl = ylabel({'In','ICU'},'Rotation',0);
yl.Position(1) = yl.Position(1) - 00;
%title(['Cases $I_t$, using ', Rtype, ' Reproduction Number'])
xlim([0 min([540 maxtime])])
grid on

% % define hospital capacity
para.eta = 0.4597;
Hbar = Hbar1/(1 - para.eta);

% calculate costs under each constraint
alpha = [para.ma 1];
%alpha = [1 10];

C1_softcosts = zeros(size(vs));
C2_softcosts = zeros(size(vs));

% hard constraint
chard1 = (alpha(1).*ones(size(out1.in_Hosp)) + (alpha(2) - alpha(1)).*(out1.in_Hosp+para.eta*Hbar>Hbar)).*out1.new_Hosp;
chard2 = (alpha(1).*ones(size(out2.in_Hosp)) + (alpha(2) - alpha(1)).*(out2.in_Hosp+para.eta*Hbar>Hbar)).*out2.new_Hosp;

% soft constraint
for v = 1:Nvs
    C1_softcosts(v) = sum((alpha(1) + (alpha(2) - alpha(1))./(1 + exp(-vs(v)*((out1.in_Hosp+para.eta*Hbar)./Hbar - 1)))).*out1.new_Hosp);
    C2_softcosts(v) = sum((alpha(1) + (alpha(2) - alpha(1))./(1 + exp(-vs(v)*((out2.in_Hosp+para.eta*Hbar)./Hbar - 1)))).*out2.new_Hosp);
end


% plotting (cost of both curves as a function of logistic scaling)
f6 = figure(6);
f6.Position = [600 600 900 300];

plot(vs,C1_softcosts,'Color',mygreen)
hold on
plot(vs,C2_softcosts,'Color',myblue)
xlim([0 max(vs)])
xlabel('Scaling of logistic curve $w$')
ylabel('Cost')
legend({'Scenario 1','Scenario 2'},'Location','north')
grid on

% %saveas(f6,'./images/F3D_renewal_costs.png')
