%Make KDE plots of parameter densities as an alternative to histograms, as
%suggested by Ariel.

%Version 2 also shows the priors.
%The prior is the same for the forward and restoration methods, so I'm just
%plotting the forward version.

folder = 'DMCBootInflSigma06_Results';
% folder = 'DMCInfl_Results';
N = 200;
fontsize = 16;

param_names = {'Fault1 Strike','Fault1 Dip','HorizonA Depth','HorizonA to HorizonB Thickness',...
    'Fault1 Maximum Displacement','Fault1 Asymmetry','Fault1 Range','Fault1 Displacement Horizontal Axis',...
    'Fault1 Displacement Vertical Axis'};
param_labels = {'Strike (\circ)','Dip (\circ)','Horizon A Restored Depth (m)','Thickness Horizon A to B (m)',...
    'Maximum Displacement (m)','Displacement Asymmetry','Reverse Drag Radius(m)','Displacement Ellipse Along-Strike Semiaxis (m)',...
    'Displacement Ellipse Along-Dip Semiaxis (m)'}; %These are the names for the plot.
true_vals = [70,55,1100,300,250,0.6,700,1000,500];
% mins = [60,40,1075,280,50,0.0,100,250,250];
% maxs = [80,65,1125,320,500,1.0,1500,1500,1500];
mins = [60,40,1000,100,50,0.0,100,250,250];
maxs = [80,65,1200,400,500,1.0,1500,1500,1500];

load(['./',folder,'/N',num2str(N),'_Forward.mat'],'params_initial_raw','params_final_raw','info')
% load('N200_Forward_scattered.mat');
gen_params_forward = zeros(info.nparams.general,N);
gen_params_forward_prior = zeros(info.nparams.general,N);
parfor (n = 1:N,4)
    gen_params_forward(:,n) = unbounded2bounded(params_final_raw(1:info.nparams.general,n),...
        info.gen_param_mins',info.gen_param_maxs');
    gen_params_forward_prior(:,n) = unbounded2bounded(params_initial_raw(1:info.nparams.general,n),...
        info.gen_param_mins',info.gen_param_maxs');
end
info_forward = info;
clear params_final_raw info

load(['./',folder,'/N',num2str(N),'_Restoration.mat'],'params_initial_raw','params_final_raw','info')
% load('N200_Restoration_scattered.mat');
gen_params_restoration = zeros(info.nparams.general,N);
gen_params_restoration_prior = zeros(info.nparams.general,N);
parfor (n = 1:N,4)
    gen_params_restoration(:,n) = unbounded2bounded(params_final_raw(1:info.nparams.general,n),...
        info.gen_param_mins',info.gen_param_maxs');
    gen_params_restoration_prior(:,n) = unbounded2bounded(params_initial_raw(1:info.nparams.general,n),...
        info.gen_param_mins',info.gen_param_maxs');
end
info_restoration = info;
clear params_final_raw info

letters = {'A','B','C','D','E','F','G','H','I'};
for iP = 1:length(param_names)
    subplot(3,3,iP)
    ind1 =find(strcmp(info_forward.gen_param_names,param_names{iP}));
    ind2 =find(strcmp(info_restoration.gen_param_names,param_names{iP}));
    pts = mins(iP):(maxs(iP)-mins(iP))/499:maxs(iP);
    [f1,xi1] = ksdensity(gen_params_forward(ind1,:),pts);
    [f1_prior,xi1_prior] = ksdensity(gen_params_forward_prior(ind1,:),pts);
    plot(xi1_prior,f1_prior,'k:','LineWidth',2)
    hold on
    plot(xi1,f1,'b-','LineWidth',2)
    [f2,xi2] = ksdensity(gen_params_restoration(ind2,:),pts);
    plot(xi2,f2,'r-','LineWidth',2)
    yl = ylim;
    plot(mean(gen_params_forward(ind1,:))*ones(1,2),yl,'b--','LineWidth',2)
    plot(mean(gen_params_restoration(ind2,:))*ones(1,2),yl,'r--','LineWidth',2)
    plot([true_vals(iP),true_vals(iP)],yl,'k-','LineWidth',2);
    hold off
    ylim(yl);
    xlabel(param_labels{iP});
    if mod(iP,3)==1
        ylabel('Probability Density')
    end
    xlim([mins(iP),maxs(iP)]);
    xl = xlim;
    text(xl(1)+(xl(2)-xl(1))*0.05,yl(2)-(yl(2)-yl(1))*0.1,letters{iP},'FontSize',fontsize);
    xlim(xl);
    if iP==9
        legend('Prior','Forward Posterior','Restoration Posterior','Forward Mean','Restoration Mean','True Value')
    end
    set(gca,'FontSize',fontsize);
    disp([param_labels{iP},' Forward: ',num2str(mean(gen_params_forward(ind1,:))),' +/- ',...
        num2str(std(gen_params_forward(ind1,:))),' Restoration: ',num2str(mean(gen_params_restoration(ind2,:))),...
        ' +/- ',num2str(std(gen_params_restoration(ind2,:)))]);
end