%Create prior ensembles, inflation z, and bootstrap index values to use 
%with the Run_All_Inversions_forward and restoration scripts.
%This should make results more reproduceable.

rng default

nparams_restoration = 2621; %Including dummy
nparams_forward = 1737; %Including dummy
nparams_gen_restoration = 24;
ndata_restored = 882;
ndata_fault = 85;
ndata_horiz = 701;
sigma_data_restored = 10;
sigma_data_fault = 10;
sigma_data_horiz = 10;
Nz_infl = 100;
Nvals = [50,100,200,500,1000,2000];
Nboot = 50;
maxit = 100;

for N = Nvals

    %Inflation z
    z = LHCube_UncorNorm(zeros(1,Nz_infl),ones(1,Nz_infl),N,'Iterations',Nz_infl);
    z = z-mean(z,2); %Shift the mean of each row to 0.
    z = z.*1./std(z,[],2); %Scale so the std of each row is exactly the expected value. (This requires mean is 0.)
    save('inflationz.mat','z');
    
    %Bootstrapping indices.
    boot_inds = zeros(Nboot,N,maxit);
    for i = 1:maxit
       boot_inds(:,:,i) = randi([1,N],Nboot,N); 
    end
    
    %Data:
    sigma_data = [sigma_data_restored*ones(1,ndata_restored),sigma_data_fault*ones(1,ndata_fault),...
        sigma_data_horiz*ones(1,ndata_horiz)]';
    ndata_total = ndata_restored+ndata_fault+ndata_horiz;
    rlzts = LHCube_UncorNorm(zeros(1,ndata_total),sigma_data,N,'Iterations',100);
    rlzts = rlzts-mean(rlzts,2); %Shift the mean of each row to 0.
    rlzts = rlzts.*sigma_data./std(rlzts,[],2); %Scale so the std of each row is exactly the expected value. (This requires mean is 0.)
    perturbations = true; %Tells that these are perturbations, not realizations of the data themselves.
    
    %Restoration parameters:
    sigma_params = [1.8*ones(1,nparams_gen_restoration),ones(1,nparams_restoration-nparams_gen_restoration)]'; %+1 to include dummy param.
    params_initial = LHCube_UncorNorm(zeros(1,nparams_restoration),sigma_params,N,'Iterations',100);
    %Shift to the desired mean and standard deviation.
    params_initial = params_initial-mean(params_initial,2); %Shift the mean of each row to 0.
    params_initial = params_initial.*sigma_params./std(params_initial,[],2); %Scale so the std of each row is exactly the expected value. (This requires mean is 0.)
    save(['PriorEnsembles/N',num2str(N),'_Restoration.mat'],'params_initial','rlzts','perturbations','z','boot_inds')
    
    %Forward parameters:
    params_initial = params_initial([1:18,21:1297,1739:2179,end],:); %Skip the block offsets that are specific to the restoration method. Skip the second fault block horizon params.
    save(['PriorEnsembles/N',num2str(N),'_Forward.mat'],'params_initial','rlzts','perturbations','z','boot_inds')
end
