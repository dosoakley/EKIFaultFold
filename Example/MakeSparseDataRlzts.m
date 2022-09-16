%Make prior data realizations for the sparse data so that too can be
%reproduceable.

rng default

ndata_restored = 882;
ndata_fault = 14;
ndata_horiz = 97;
sigma_data_restored = 10;
sigma_data_fault = 10;
sigma_data_horiz = 10;
N = 200;

sigma_data = [sigma_data_restored*ones(1,ndata_restored),sigma_data_fault*ones(1,ndata_fault),...
    sigma_data_horiz*ones(1,ndata_horiz)]';
ndata_total = ndata_restored+ndata_fault+ndata_horiz;
rlzts = LHCube_UncorNorm(zeros(1,ndata_total),sigma_data,N,'Iterations',100);
rlzts = rlzts-mean(rlzts,2); %Shift the mean of each row to 0.
rlzts = rlzts.*sigma_data./std(rlzts,[],2); %Scale so the std of each row is exactly the expected value. (This requires mean is 0.)
perturbations = true; %Tells that these are perturbations, not realizations of the data themselves.

save('SparseDataRlzts.mat','rlzts','perturbations')