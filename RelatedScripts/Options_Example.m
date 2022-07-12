%Options file for Test6.

%Options for the Ensemble Kalman inversion.
opt.N = 100; %Number of models in the ensemble.
opt.max_iterations = 30; %Maximum number of iterations for which to run the ensemble Kalman filter.
opt.boot = true; %If true, use EKI_DMC_boot_Infl. (Requires regularize and DMC to be true to.)
opt.nboot = 50; %Number of bootstrap samples. Use only if boot is true.
opt.sigma_alpha = 0.6; %sigma_alpha for bootstrap screening. Use only if boot is true.
opt.infl = true; %If true, use covariance inflation.
opt.nz_infl = 100; %Number of dummy parameters for determining inflation factor. Use only if infl is true.
opt.algorithm = 'subspace'; %Options are 'exact', 'subspace', or 'lowrank'.
opt.Nthreads = 8; %Number of threads to use for parallel processing.

%Choose the structural modelling method: restoration or forward modelling.
opt.method = 'restoration'; %Options are 'forward' or 'restoration'.

%Define some basic information about the model.
opt.horizon_names = {'HorizonA','HorizonB'}; %Horizons should be in order, either shallowest to deepest or vice versa.
opt.nhorizons = length(opt.horizon_names);
opt.n_fault_blocks = 2; %Number of fault blocks. (Only used for restoration.
opt.nfaults = 1; %Number of faults in the model.

%Define paths tothe reference model and the data.
opt.data_path = '..\ForwardModel\DataDense';
opt.ref_model_path = 'ReferenceModel_Deformed_Test6.mat';

%Choose what properties to fit for
opt.fit_disp = true; %Fit for displacement on the fault.
opt.fit_surf = true; %Fit for the fault surface geometry.
opt.fit_horizons = true; %Fit for the horizon geometry.

%Choose whether to fit for the residual terms or just the trend for each
%gridded property. The trends alone are an elliptical displacement field, a
%planar fault, and a planar horizon. Note that these will only be used if
%the corresponding entry in the block above is also true.
opt.fit_disp_residuals = true;
opt.fit_surf_residuals = true;
opt.fit_horizons_residuals = true;

%Choose which types of data to use.
opt.restored_elevation_data = true; %Difference between restored point elevations and the horizontal trend for each horizon.
opt.fault_points_data = true; %Points on the fault.
opt.horizon_points_data = true; %Points on the horizons.

%Specify the uncertainty in the data of each data type.
opt.sigma_misfit_restored = 10; %Uncertainty in restored-state misfits.
opt.sigma_fault = 10; %Uncertainty in the fault surface points.
opt.sigma_horizon = 10; %Uncertainty in the horizon points.
opt.correlated_data = false; %If true, use the range to create spatially correlated data realizations.
opt.range_misfit_restored_data = 1e3; %Correlation lengths to use if correlated_data is true.
opt.range_fault_data = 1e3;
opt.range_horizon_data = 2e3;

%Options for the fault surface realizations.
%Each of these should be of size nfaults x 2.
opt.surf_range_lims = repmat([1000,5000],opt.nfaults,1); %Range of allowed values for surface variogram range.
opt.surf_std_lims = repmat([20,150],opt.nfaults,1); %Range of allowed standard deviation for surface Gaussian random fields.
opt.zlims = [800,1600]; %z coordinates (depth) of each fault center.
opt.strike_lims = [60,80]; %strike of each fault.
opt.dip_lims = [40,65]; %dip of each fault.

%Options for the fault displacement realizations.
%Each of these should be of size nfaults x 2.
opt.disp_range_lims = [100,1000]; %Range of allowed values for displacement variogram range.
opt.disp_std_lims = [10,60]; %Range of allowed standard deviation for displacement Gaussian random fields.
opt.disp_max_lims = [100,350]; %Range of allowed maximum displacement values.
opt.disp_axis_lims = [250,1500]; %Range of allowed axis lengths for ellipse.
opt.u0_lims = [-200,200]; %Displacement ellipse center u coordinate.
opt.v0_lims = [-200,200]; %Displacement ellipse center v c00rdinate.

%Options for the horizon realizations:
opt.horiz_range_lims = repmat([1e3,5e3],opt.nhorizons,1); %Range of allowed values for surface correlation length.
opt.horiz_std_lims = repmat([1,100],opt.nhorizons,1); %Range of allowed standard deviation for surface Gaussian random fields.
opt.horiz1_depth_lims = [800,1400]; %Depth of the first horizon (either shallowest or deepest).
opt.horiz_thickness_lims = [100,600]; %Thickness of each horizon (except the last). Positive if first horizon is shallowest, negative if its deepest.
opt.slopex_lims = [tand(-5),tand(5)]; %Limits of the restored-state slope in the x direction.
opt.slopey_lims = [tand(-5),tand(5)]; %Limits of the restored-state slope in the y direction.
opt.block_offset_lims = repmat([-100,100],opt.n_fault_blocks,1); %Vertical offsets of each fault block relative to the restored state horizon trends.
%Note: opt.block_offset_lims is only used if opt.method = 'restoration'. 
%If this is the case but opt.fit_disp = true, then at least one set of
%block_offset_lims should be set to 0 (min and max).

%Define limits for asymmetry and range, which are the two kinematic models
%that are fit for. Each of these should be of size nfaults x 2.
opt.asym_lims = [0,1];
opt.range_lims = [100,1.5e3];

%Fault names.
%These must be in the order in which the faults were added (opposite of the order in which they will be removed).
opt.fault_names = {'Fault1'};

%Fault block matrix.
%This tells relationship of each fault block to each fault.
%0 = no relationship, -1 = FW side, 1 = HW side.
%Row specifies the fault block and column the fault.
%This information comes from the FaultBlocks section in fault_model.rms in
%the reference model.
opt.fault_blocks_relationships = [1;-1];

%Fault truncations matrix.
%This tells how faults truncate each other.
%Row specified the fault being truncated and column the fault that
%truncates it.
%0 = no truncatione rule, -1 = truncated in FW, 1 = truncated in HW.
opt.fault_truncations = [0];

%Set up the horizontal grid.
%xstep should divide evenly into xmax-xmin, and ystep should divide evenly
%into ymax-ymin.
[opt.xmin,opt.xmax] = deal(-1e3,1e3);
[opt.ymin,opt.ymax] = deal(-1e3,1e3);
[opt.xstep,opt.ystep] = deal(100,100); %x and y step sizes for grid.

%Optionally, trim the horizon grids to within some distance of the 
%reference model points in that fault block to reduce the total number of
%model parameters.
opt.trim_horiz_grids = false;
opt.trim_horiz_dist = 500;