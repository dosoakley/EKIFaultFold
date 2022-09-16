%Make a figure to compare the localization / inflation methods and
%different ensemble sizes.
%Plot mean and standard deviation of dummy parameter.

methods = {'Restoration','Forward'};
folders = {'DMC_Results','DMCBootSigma06_Results','DMCInfl_Results','DMCBootInflSigma06_Results'};
titles = {'EKI-DMC','With Localization','With Inflation','With Localization and Inflation'};
fontsize = 17; %20;

% Nvals = [50,100,200,500,1000,2000];
Nvals = [100,200,500,1000,2000];
colors = {'k','b'};
dummy_mean_final = zeros(2,4,length(Nvals));
dummy_std_final = zeros(2,4,length(Nvals));
asym_mean_final = zeros(2,4,length(Nvals));
asym_std_final = zeros(2,4,length(Nvals));

p = [];
for iM = 1:length(methods)
    for iF = 1:length(folders)
        for jN = 1:length(Nvals)
            %Load the data.
            load(['.\',folders{iF},'\N',num2str(Nvals(jN)),'_',methods{iM},'.mat'],'params_final_raw','info')
            dummy_param_final = params_final_raw(end,:);
            dummy_mean_final(iM,iF,jN) = mean(dummy_param_final);
            dummy_std_final(iM,iF,jN) = std(dummy_param_final);
        end
        %Plot the dummy parameter standard deviation.
        subplot(2,4,iF)
        pl = plot(Nvals,squeeze(dummy_std_final(iM,iF,:)),[colors{iM},'-x']);
        if iF==4
            p(iM) = pl;
        end
        hold on
        if iM==2
            p_i = plot(Nvals,ones(size(Nvals)),'r-');
%             xlabel('Ensemble Size')
            if iF==4
%                 legend('Final St. Dev. (Restoration)','Final St. Dev. (Forward)','Initial St. Dev.')
                legend([p_i,p(1),p(2)],{'Initial','Final (Restoration)','Final (Forward)'})
            end
        end
        if iF==1
            ylabel('Dummy Parameter St. Dev.')
        end
        ylim([0,2])
        title(titles{iF})
%         xticks(Nvals(2:end))
        xticks(Nvals)
        xtickangle(-60)
        set(gca,'FontSize',fontsize)
        %Plot the dummy parameter mean.
        subplot(2,4,4+iF)
        plot(Nvals,squeeze(dummy_mean_final(iM,iF,:)),[colors{iM},'-x']);
        hold on
        if iM==2
            plot(Nvals,zeros(size(Nvals)),'r-');
            xlabel('Ensemble Size')
%             if iF==4
%                 legend('Final Mean (Restoration)','Final Mean (Forward)','Initial Mean')
%             end
        end
        if iF==1
            ylabel('Dummy Parameter Mean')
        end
        ylim([-1,1])
%         title(titles{iF})
%         xticks(Nvals(2:end))
        xticks(Nvals)
        xtickangle(-60)
        set(gca,'FontSize',fontsize) 
    end
end
