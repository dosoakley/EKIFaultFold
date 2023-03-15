%Run the Test6 Ensemble Kalman Inversion.

%I am running this using EnKF_Combined_Code15.
%The data being fit are horizon points, fault points, and restored-state
%elevation misfits.

%Loop through all the difference inversions we want to do.
rng default; %Set the random number seed here so that results are repeatable.

method = 'Restoration';

options_file_name = ['Options_Test6_Dense_',method,'.m'];
run(options_file_name);

%Make lists of the different settings to try.
infl = [false,true,false,true];
boot = [false,false,true,true];
folder_names = {'DMC_Results','DMCInfl_Results','DMCBootSigma06_Results','DMCBootInflSigma06_Results'};
N = [50,100,200,500,1000,2000];

for i = 1:4
    for j = 1:length(N)
        opt.infl = infl(i);
        opt.boot = boot(i);
        opt.N = N(j);
        opt.inflation_z_file = ['./PriorEnsembles/N',num2str(opt.N),'_',method];
        opt.boot_inds_file = ['./PriorEnsembles/N',num2str(opt.N),'_',method];
        opt.prior_path = ['./PriorEnsembles/N',num2str(opt.N),'_',method];
        opt.rlzts_path = ['./PriorEnsembles/N',num2str(opt.N),'_',method];
        disp(['Starting',folder_names{i},'/N',num2str(N(j))])
        [params_initial_raw,params_final_raw,info] = RunEKI(opt); %Run the inversion.
        filename = [folder_names{i},'/N',num2str(N(j)),'_',method,'.mat'];
        save(filename,'params_initial_raw','params_final_raw','info','opt');
    end
end