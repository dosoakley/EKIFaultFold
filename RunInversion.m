%Run the Test5 Ensemble Kalman Inversion.
%This version is meant to be run with EnKF_Combined_Code8 or EnKF_Combined_Code9.
%For the Dunlop versions (which didn't work well), you need to use 8.
%For the restoration version, you need to use 9.
%Other versions should run the same with either version of the code.

% options_file_name = 'Options_EmeraldField_Small_forward.m';
% options_file_name = 'Options_EmeraldField_Small_restoration_RMSPrior.m';
% options_file_name = 'Options_EmeraldField_Small_restoration_RMSPrior_centered.m';
% options_file_name = 'Options_EmeraldField_Small_forward_v2.m';
options_file_name = 'Options_EmeraldField_Small_restoration_v2.m';

run(options_file_name);
[params_initial_raw,params_final_raw,info] = RunEKI(opt); %Run the inversion.

% save('ResultsN200.mat','params_initial_raw','params_final_raw','info','opt');
% save('ResultsN200_restoration_RMSPrior.mat');
% save('ResultsN200_restoration_RMSPrior_centered.mat');
% save('ResultsN200_forward_v2.mat','params_initial_raw','params_final_raw','info','opt');
save('ResultsN200_restoration_v2.mat','params_initial_raw','params_final_raw','info','opt');

%%
%Display the general parameters.
gen_params = unbounded2bounded(params_final_raw(1:info.nparams.general,:),info.gen_param_mins',info.gen_param_maxs');
for i = 1:info.nparams.general
    disp([info.gen_param_names{i},': ',num2str(mean(gen_params(i,:))),' +/- ',num2str(std(gen_params(i,:))),...
        '     Limits: ',num2str(info.gen_param_mins(i)),' to ',num2str(info.gen_param_maxs(i))])
end