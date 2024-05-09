function [out_params,EKI_info] = EKI_DMC(func,N,inparams,data,xi,Gamma,maxit,options)
%This function implements regularizing iterative Ensemble Kalman inversion 
%with the data misfit controller, as described by Iglesias and Yang (2021).
%It optionally allows the use of bootstrapping-based screening (Zhang and 
%Oliver (2010), covariance inflation with the method described in Evensen 
%(2009) and Lacerda et al. (2019), subspace pseudo-inversion (Evensen,
%2004; 2009), and subspace pseudo-inversion with a low-rank approximation
%of the error covariance matrix (Evensen, 2004; 2009).
%Arguments:

%I'm not sure if I should update the model in the last step (as I'm doing 
%right now) or quit as soon as the convergence criterion is reached .

arguments
    func function_handle %The function to be evaluated. It should return a 1D vector of misfits of length ndata.
    N (1,1) {mustBeInteger} %Number of ensemble members.
    inparams (:,:) {mustBeReal} %An array of initial parameter values of size (nparams,N).
    data (:,1) {mustBeReal} %An array of data of size (ndata,1).
    xi (:,:) {mustBeReal} %xi = An array of perturbations of the data of size (ndata,N).
    Gamma (:,:) {mustBeReal} %Gamma = Covariance matrix for data perturbations.
    maxit (1,1) {mustBeInteger} %maxit = Maximum number of iterations allowed, regardless of convergence.
    options.boot_screen (1,1) logical = false %Whether to use bootstrapping-based screening of K.
    options.Nboot (1,1) {mustBeInteger} = 100 %Number of bootstrap samples for the Kalman gain matrix, only used if boot_screen is true.
    options.sigma_alpha (1,1) {mustBeReal} = 0.6 %Parameter in Eqn. 4 of Zhang and Oliver (2010).
    options.load_boot_inds = false %Whether to load bootstrap sample indices from a file. (Can be useful for reproducibility.)
    options.boot_inds_file = '' %Path to a .mat file from which to load the bootstrap sample indices if load_boot_inds is true.
    options.inflate (1,1) logical = false %Whether to use covariance inflation.
    options.Nz (1,1) {mustBeInteger} = 100 %Number of white noise model parameters to use for calculating the covariance inflation factor.
    options.persistent_z = false %Whether to keep using the same dummy parameters for inflation or regenerate them at each iteration.
    options.load_z = false %Whether to load initial z values from a file. (Can be useful for reproducibility).
    options.z_file = '' %Path to a .mat file from which to load z values if load_z is true.
    options.algorithm {mustBeMember(options.algorithm,['exact','subspace','lowrank'])} = 'exact' %Exact inversion, subspace pseudoinverse, or subspace with low rank error matrix approximation.
    options.Nthreads = 4 %Number of threads to use for parfor loops.
    options.dataset_start_inds (:,1) {mustBeReal} = [] %An array indicating start indices of different datasets for which different RMS errors will be calculated.
end

nparams = size(inparams,1);
M = size(data,1); %Number of data points.
if options.inflate
    params = [inparams;zeros(options.Nz,N)]; %Initialize the parameter values.
else
    params = inparams;
end
clear inparams %Free up memory

if any(Gamma==0)
    Gamma = sparse(Gamma); %Gamma is in many cases sparse, so do this to save space.
end
Gamma_inv = Gamma^-1; %This only has to be calculated once.

n=0;
tic;
t1_cpu = cputime;
RMS = zeros(1,maxit);
alpha_all = zeros(1,maxit);
gamma_all = zeros(1,maxit);
t = 0;
while (n<maxit)
    n = n+1;
    disp(['Iteration ',num2str(n)])
    %Run through the ensemble of models.
    G = zeros(M,N); %G will hold predicted data values.
    parfor (j = 1:N,options.Nthreads)
%     for j = 1:N
%         disp(j)
        G(:,j) = func(params(:,j)); %Need to figure out how to deal with inflation params. Maybe two cases as for out_params. Maybe separate z from params.
    end
    %Calculate alpha.
    Phi = zeros(1,N);
    parfor (j = 1:N,options.Nthreads)
        Phi(j) = 0.5*(data-G(:,j))' * Gamma_inv * (data-G(:,j))
        %The equation used to calculate the least-squares misfit term, Phi, is Eqn. 13 of Ma et al. (2017).
        %This is (I am pretty sure) the same Eqn. 3 of Iglesias and Yang (2021), but is an easier form to understand and implement.
    end
    alpha = 1/(min(max(M/(2*mean(Phi)),sqrt(M/(2*var(Phi)))),1-t)); %Eqn. 14
    %If using inflation, make the vector of white noise to be used for determining the inflation factor.
    if options.inflate && (n==1 || ~options.persistent_z)
        if options.load_z
            load(options.z_file,'z')
            z = z(1:options.Nz,1:N,n);
        else
            z = lhsnorm(zeros(1,options.Nz),eye(options.Nz),N)'; %Generate z from a standard normal distribution using Latin hypercube sampling.
        end
        z = z-mean(z,2); %Shift z so the mean of each row is exactly zero.
        z = z./std(z,[],2); %Scale z so the std of each row is exactly one.
        params(nparams+1:end,:) = z;
    end
    %Calculate the Kalman gain matrix.
    switch options.algorithm
        case('exact')
            K = KalmanGain(G,params,Gamma,N,alpha);
        case('subspace')
            K = KalmanGain_sub(G,params,Gamma,N,M,alpha);
        case('lowrank')
            K = KalmanGain_sub_low(G,params,xi,N,M,alpha);
    end
    %If using bootstrap-based screening, calculate bootstrapped mean and variance of K and use it to screen K.
    if options.boot_screen
        Kvar = zeros(size(K));
        if options.load_boot_inds
            load(options.boot_inds_file,'boot_inds')
            boot_inds = boot_inds(1:options.Nboot,1:N,n);
        end
%         parfor (j = 1:options.Nboot,options.Nthreads) %Bootstrap K.
        for j = 1:options.Nboot %Bootstrap K.
            if options.load_boot_inds
                inds = boot_inds(j,:);
            else
                inds = randi([1,N],1,N);
            end
            switch options.algorithm
                case('exact')
                    Kboot = KalmanGain(G(:,inds),params(:,inds),Gamma,N,alpha);
                case('subspace')
                    Kboot = KalmanGain_sub(G(:,inds),params(:,inds),Gamma,N,M,alpha);
                case('lowrank')
                    Kboot = KalmanGain_sub_low(G(:,inds),params(:,inds),xi(:,inds),N,M,alpha);
                case default
                    error('Unrecognized Algorithm');
            end
            Kvar = Kvar+(Kboot-K).^2; %As in Eqn. 1 of Zhang and Oliver (2010), use true K for calculating the variance.
        end
        Kvar = Kvar./options.Nboot; %Eqn. 1 of Zhang and Oliver (2010). Note that they use Nboot, not Nboot-1.
        R = Kvar./K.^2; %From Zhang and Oliver, 2010.
        clear Kvar Kboot
        alpha_boot = 1./(1+R.^2*(1+1/options.sigma_alpha^2)); %Eqn. 4 of Zhang and Oliver (2010).
        clear R
        alpha_boot(isnan(alpha_boot)) = 0; %Fixed value parameters (which I shouldn't really have anyway) give NaN, so fix that.
        K = alpha_boot.*K; %Screened Kalman gain matrix. (Unnumbered Eqns. at end of section 2 of Zhang and Oliver (2010).)
        clear alpha_boot
    end
    %Update the parameters in all ensemble members.
%     parfor (j = 1:N,options.Nthreads)
    for j = 1:N
        params(:,j) = params(:,j)+K*(data+sqrt(alpha)*xi(:,j)-G(:,j));  %Update params (Eqn. 4)
    end
    t = t+alpha^-1; %Update t. (Eqn. 15)
    %If using covariance inflation, apply it now.
    if options.inflate
        z = params(nparams+1:end,:);
        sigma_z = sqrt((1/(N-1))*sum((z-mean(z,2)).^2,2)); %Eqn. 37 in Lacerda et al. (2019).
        gamma_k = 1/((1/options.Nz)*sum(sigma_z)); %Eqn. 36 in Lacerda et al. (2019).
        params = gamma_k*(params-mean(params,2))+mean(params,2);
    end
    %Calculate RMS and display progress so far.
    RMS(n) = sqrt(mean((G-data).^2,'all')); %RMS error.
    alpha_all(n) = alpha;
    disp(['RMS Misfit: ',num2str(RMS(1:n))])
    disp(['alpha: ',num2str(alpha_all(1:n))])
    if options.inflate
        gamma_all(n) = gamma_k;
        disp(['gamma: ',num2str(gamma_all(1:n))])
    end
    if ~isempty(options.dataset_start_inds)
        ndatasets = length(options.dataset_start_inds);
        dataset_end_inds = [options.dataset_start_inds(2:end)-1;M];
        RMS_datasets = zeros(1,ndatasets);
        RMS_datasets2 = zeros(1,ndatasets);
        for i = 1:ndatasets
            inds = options.dataset_start_inds(i):dataset_end_inds(i);
            RMS_datasets(i) = sqrt(mean((G(inds,:)-data(inds)).^2,'all')); %RMS for all points
            RMS_datasets2(i) = sqrt(mean((mean(G(inds,:),2)-data(inds)).^2)); %RMS for the dataset relative to mean model.
        end
        disp(['RMS over all models by dataset: ',num2str(RMS_datasets)]);
        disp(['RMS for the average model by dataset: ',num2str(RMS_datasets2)]);
    end
    %Check for convergence, using the principle given in section 1.3 of Iglesias and Yang (2021).
    %From looking at Iglesias and Yang (2021), I have determined that a
    %final update should be done before stopping. Thus the sum of all
    %alphas will equal 1.
    if t==1
        if options.inflate
            out_params = params(1:end-options.Nz,:);
        else
            out_params = params;
        end
        clear params
        break
    end

end
time = toc;
t2_cpu = cputime;
cpu_time = t2_cpu-t1_cpu;

%Summarize the results.
if n<maxit
    disp(['Converged after ',num2str(n),' iterations'])
else
    disp('Failed to converge.')
    if options.inflate
        out_params = params(1:end-options.Nz,:);
    else
        out_params = params;
    end
    clear params
end
disp(['Wall Time = ',num2str(time),' seconds']);
disp(['CPU Time = ',num2str(cpu_time),' seconds']);
disp(['Final RMS Error = ',num2str(RMS(n))])
EKI_info.RMS = RMS(1:n);
EKI_info.alpha = alpha_all(1:n);
EKI_info.gamma = gamma_all(1:n);
EKI_info.wall_time = time;
EKI_info.cpu_time = cpu_time;
end

function K = KalmanGain(G,params,Gamma,N,alpha)
%Kalman gain type matrix of Iglesias and Yang (2021).
CnGG = cov(G');   %Eqn. 16 of Iglesias and Yang (2021).
K = (1/(N-1))*(params-mean(params,2))*((G-mean(G,2))'/(CnGG+alpha*Gamma)); %This if equivalent to CnuG*inv(CnGG+Gamma) but faster.
end

function K = KalmanGain_sub(G,params,Gamma,N,M,alpha)
%Subspace pseudo inversion.
S = sqrt(1/(N-1))*(G-mean(G,2)); %S*S'=CnGG.
[U0e,Sigma0,~] = svd(S,'econ');
U0 = sparse(M,M);
if N<M
    U0(:,1:N) = U0e;
else
   U0 = U0e; 
end
sigma = diag(Sigma0);
tol = max(size(S))*eps(max(sigma)); %Same as the default tolerance if I called rank(S).
clear S U0e Sigma0
p = sum(sigma>tol); %This is the rank of S. This is faster than calling rank(S) because that would do the svd again.
Sigma0Plus = sparse(1:p,1:p,1./sigma(1:p),N,M);
X0 = full(Sigma0Plus*U0'*(alpha*Gamma)*U0*Sigma0Plus');
[Z1,Lambda1] = eig(X0);
X1 = U0*Sigma0Plus'*Z1;
Cplus = X1/(eye(N)+Lambda1)*X1';
K = (1/(N-1))*(params-mean(params,2))*((G-mean(G,2))'*Cplus); %This is equivalent to CnuG*Cplus but faster.
end

function K = KalmanGain_sub_low(G,params,E,N,M,alpha) %E will be xi from above.
%Subspace pseudo inversion with low-rank approximation of Gamma.
S = sqrt(1/(N-1))*(G-mean(G,2)); %S*S'=CnGG.
[U0e,Sigma0,~] = svd(S,'econ');
U0 = sparse(M,M);
if N<M
    U0(:,1:N) = U0e;
else
    U0 = U0e;
end
sigma = diag(Sigma0);
tol = max(size(S))*eps(max(sigma)); %Same as the default tolerance if I called rank(S).
p = sum(sigma>tol); %This is the rank of S. This is faster than calling rank(S) because that would do the svd again.
Sigma0Plus = sparse(1:p,1:p,1./sigma(1:p),N,M);
X0 = Sigma0Plus*U0'*sqrt(alpha)*E;
[U1,Sigma1,~] = svd(X0);
X1 = U0*Sigma0Plus'*U1;
Cplus = (X1/(eye(N)+Sigma1.^2)*X1')/(N-1); %Since Sigma1 is a diagonal matrix, Sigma1^2 or Sigma1.^2 should be the same. .^2 is probably faster.
K = (1/(N-1))*(params-mean(params,2))*((G-mean(G,2))'*Cplus); %This is equivalent to CnuG*Cplus but faster.
end