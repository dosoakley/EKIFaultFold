%Make KDE plots of parameter densities as an alternative to histograms, as
%suggested by Ariel.

%Version 2 also shows the priors.
%The prior is the same for the forward and restoration methods, so I'm just
%plotting the forward version.

% folder = 'DMCBootInflSigma06_Results';
% folder = 'DMCInfl_Results';
N = 200;
fontsize = 16;

param_names = {'Fault1 Maximum Displacement','Fault1 Asymmetry'};
param_labels = {'Maximum Displacement (m)','Displacement Asymmetry'}; %These are the names for the plot.
true_vals = [250,0.6];
mins = [50,0.0];
maxs = [500,1.0];

load('N200_Forward_scattered.mat');
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

load('N200_Restoration_scattered.mat');
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

letters = {'A','B'};
for iP = 1:length(param_names)
    subplot(2,1,iP)
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
%     yl(2) = yl(2)*1.2; %Exend this slightly.
    plot(mean(gen_params_forward(ind1,:))*ones(1,2),yl,'b--','LineWidth',2)
    plot(mean(gen_params_restoration(ind2,:))*ones(1,2),yl,'r--','LineWidth',2)
    plot([true_vals(iP),true_vals(iP)],yl,'k-','LineWidth',2);
    hold off
    ylim(yl);
    xlabel(param_labels{iP});
    ylabel('Probability Density')
    xlim([mins(iP),maxs(iP)]);
    xl = xlim;
%     text(xl(1)+(xl(2)-xl(1))*0.05,yl(2)-(yl(2)-yl(1))*0.1,letters{iP},'FontSize',fontsize);
    xlim(xl);
    if iP==2
        legend('Prior','Forward Posterior','Restoration Posterior','Forward Mean','Restoration Mean','True Value')
    end
    set(gca,'FontSize',fontsize);
    disp([param_labels{iP},' Forward: ',num2str(mean(gen_params_forward(ind1,:))),' +/- ',...
        num2str(std(gen_params_forward(ind1,:))),' Restoration: ',num2str(mean(gen_params_restoration(ind2,:))),...
        ' +/- ',num2str(std(gen_params_restoration(ind2,:)))]);
end