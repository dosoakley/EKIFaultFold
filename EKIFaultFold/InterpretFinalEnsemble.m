% load('ResultsN200.mat');
load('ResultsN200_restoration.mat');

%Transform the final parameters back to the bounded domain.
gen_params_final = zeros(info.nparams.general,opt.N);
faults_final = repmat(fault_class.empty(),opt.N,opt.nfaults);
horiz_params_final = zeros(info.nparams.horiz,opt.N);
parfor i = 1:opt.N
    [gen_params_final(:,i),faults_final(i,:),horiz_params_final(:,i)] = InterpretModel(params_final_raw(:,i),opt,info);
end

% save('ResultsN200_Interp.mat','info','opt','gen_params_final','faults_final','horiz_params_final');
save('ResultsN200_restoration_Interp.mat','info','opt','gen_params_final','faults_final','horiz_params_final');