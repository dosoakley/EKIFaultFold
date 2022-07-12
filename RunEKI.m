function [params_initial,params_final,info] = RunEKI(opt)
%This runs the Ensemble Kalman inversion using the opt structure, with values 
%for various user-specified option. This must be in the form created by the
%example file.

%Check that the modelling method is either forward or restoration.
if ~(strcmp(opt.method,'forward') || strcmp(opt.method,'restoration'))
    error('Unrecognized Modelling Method')
end

%Load the reference model.
load(opt.ref_model_path,'faults_ref','horiz_ref');
info.faults_ref = copy(faults_ref);
clear faults_ref

%Create the horizontal grid.
[xgrid, ygrid] = meshgrid(opt.xmin:opt.xstep:opt.xmax, opt.ymin:opt.ystep:opt.ymax);
info.nx = size(xgrid,2);
info.ny = size(xgrid,1);
info.xgrid = xgrid(:);
info.ygrid = ygrid(:);
clear xgrid ygrid
info.n_horiz_pts = length(info.xgrid); %The number of data points per horizon.

%If necessary, create smaller horizontal grids for each fault block for restoration.
if strcmp(opt.method,'restoration')
    info.grid_pts_mask = cell(opt.n_fault_blocks,1);
    info.n_horiz_pts_block = zeros(opt.n_fault_blocks,1);
    if opt.trim_horiz_grids
        [x,y] = deal([],[]);
        [FaultBlock,active] = deal([],logical([]));
        for i = 1:opt.nhorizons
            x = [x;horiz_ref.x{i}];
            y = [y;horiz_ref.y{i}];
            FaultBlock = [FaultBlock;horiz_ref.FaultBlock{i}];
            active = [active;horiz_ref.active{i}];
        end
        for i = 1:opt.n_fault_blocks
            %This double for loop is slow.
            min_dist = zeros(size(info.xgrid));
            mask = active & FaultBlock==i;
            parfor j = 1:length(min_dist)
                min_dist(j) = min(sqrt((info.xgrid(j)-x(mask)).^2+(info.ygrid(j)-y(mask)).^2));
            end
            mask2 = min_dist <= opt.trim_horiz_dist;
            info.grid_pts_mask{i} = mask2;
            info.n_horiz_pts_block(i) = sum(mask2);
        end
        clear active FaultBlock x y z mask2 min_dist
    else
        for i = 1:opt.n_fault_blocks
            info.grid_pts_mask{i} = true(size(info.xgrid));
            info.n_horiz_pts_block(i) = info.n_horiz_pts;
        end
    end
end

%Interpolate the reference horizon model onto the grid.
info.horiz_pts_ref = zeros(opt.nhorizons*info.n_horiz_pts,1);
ind_start = 1;
for i = 1:opt.nhorizons
    if strcmp(opt.method,'restoration')
        for j = 1:opt.n_fault_blocks
            npts = info.n_horiz_pts_block(j);
            mask = horiz_ref.FaultBlock{i}==j;
            interpolant = scatteredInterpolant(horiz_ref.x{i}(mask),horiz_ref.y{i}(mask),horiz_ref.z{i}(mask));
            info.horiz_pts_ref(ind_start:ind_start+npts-1) = interpolant(info.xgrid(info.grid_pts_mask{j}),info.ygrid(info.grid_pts_mask{j}));
            ind_start = ind_start+npts;
        end
    else
        npts = info.n_horiz_pts;
        interpolant = scatteredInterpolant(horiz_ref.x{i},horiz_ref.y{i},horiz_ref.z{i});
        info.horiz_pts_ref(ind_start:ind_start+npts-1) = interpolant(info.xgrid,info.ygrid);
        ind_start = ind_start+npts;
    end
end
clear horiz_ref

%Calculate the matrix of distances between grid points.
info.dist_mat = zeros(info.n_horiz_pts,info.n_horiz_pts); %Inter-point distances.
for i = 1:info.n_horiz_pts-1 %Compute distances. For i==j, dist_mat will be 0, so we don't change it.
    dists = sqrt((info.xgrid(i+1:info.n_horiz_pts)-info.xgrid(i)).^2+(info.ygrid(i+1:info.n_horiz_pts)-info.ygrid(i)).^2);
    info.dist_mat(i,i+1:info.n_horiz_pts) = dists;
    info.dist_mat(i+1:info.n_horiz_pts,i) = dists; %Since the distances are the same in either direction.
end

%Create the data covariance matrix and data realizations:
if opt.restored_elevation_data
    %Restored-state data:
    %The "data" are the misfits of the restored points from a horizontal surface.
    info.ndata.restored_elev = opt.nhorizons*info.n_horiz_pts; %The total number of restored-state data points.
    if opt.correlated_data
        Gamma_d_restored_elev = [];
        for i = 1:nhorizons
            Gamma_d_misfit_part = (opt.sigma_misfit_restored^2)*CorrSpher(info.xgrid,info.ygrid,opt.range_misfit_restored_data);
            if i==1
                Gamma_d_restored_elev = Gamma_d_misfit_part;
            else
                Gamma_d_restored_elev = [Gamma_d_restored_elev,zeros(length(Gamma_d_restored_elev),npts);zeros(npts,length(Gamma_d_restored_elev)),Gamma_d_misfit_part];
            end
        end
        clear Gamma_d_misfit_part
    else
        Gamma_d_restored_elev = (opt.sigma_misfit_restored^2)*eye(info.ndata.restored_elev);
    end
    restored_elev_data_pts = zeros(info.ndata.restored_elev,1); %The expected value is 0, since it's measured as the misfit between model and data.
else
    info.ndata.restored_elev = 0;
    restored_elev_data_pts = [];
    Gamma_d_restored_elev = [];
end
if opt.fault_points_data
    %Fault surface data.
    %Since the fault reference plane is not known ahead of time, the "data"
    %are actually the misfits between the fault plane and the actual data, 
    %with expected value of 0.
    fault_data = [];
    Gamma_d_fault = [];
    info.ndata.fault = 0;
    for i = 1:opt.nfaults
        file = fopen([opt.data_path,'\',info.faults_ref(i).name,'_Data.txt']);
        fault_data = [fault_data;textscan(file,'%f%f%f')];
        npts_fault = length(fault_data{i,1});
        if opt.correlated_data
            [x,y,z] = deal(fault_data{i,1},fault_data{i,2},fault_data{i,3});
            fault_uv = info.faults_ref(i).M\[x,y,z,ones(size(x))]';
            [u,v] = deal(fault_uv(1,:)',fault_uv(2,:)');
            Gamma_d_fault_part = (opt.sigma_fault^2)*CorrSpher(u,v,opt.range_fault_data);
        else
            Gamma_d_fault_part = (opt.sigma_fault^2)*eye(npts_fault);
        end
        if i == 1
            Gamma_d_fault = Gamma_d_fault_part;
        else
            Gamma_d_fault = [Gamma_d_fault,zeros(length(Gamma_d_fault),npts_fault);...
                zeros(npts_fault,length(Gamma_d_fault)),Gamma_d_fault_part];
        end
        fclose(file);
        info.ndata.fault = info.ndata.fault+npts_fault;
    end
    fault_data_pts = zeros(info.ndata.fault,1); %Since the fault data are measured as misfit from the surface, the expected value of the misfit is 0.
    clear Gamma_d_fault_part x y z u v file
else
    fault_data_pts = [];
    info.ndata.fault = 0;
    fault_data = [];
    Gamma_d_fault = [];
end
if opt.horizon_points_data
    %Deformed-state horizons data.
    x_horiz_data = [];
    y_horiz_data = [];
    id_horiz_data = [];
    horiz_data_pts = [];
    Gamma_d_horiz = [];
    for i = 1:opt.nhorizons
        file = fopen([opt.data_path,'\',opt.horizon_names{i},'_Data.txt']);
        horizon_data = textscan(file,'%f%f%f');
        [x,y,z] = deal(horizon_data{1},horizon_data{2},horizon_data{3});
        mask = (x>opt.xmin & x<opt.xmax & y>opt.ymin & y<opt.ymax);
        x_horiz_data = [x_horiz_data;x(mask)];
        y_horiz_data = [y_horiz_data;y(mask)];
        id_horiz_data = [id_horiz_data;i*ones(sum(mask),1)];
        horiz_data_pts = [horiz_data_pts;z(mask)];
        if opt.correlated_data
            Gamma_d_horiz_part = (opt.sigma_horizon^2)*CorrSpher(x(mask),y(mask),opt.range_horizon_data);
        else
            Gamma_d_horiz_part = (opt.sigma_horizon^2)*eye(sum(mask));
        end
        if i == 1
            Gamma_d_horiz = Gamma_d_horiz_part;
        else
            Gamma_d_horiz = [Gamma_d_horiz,zeros(length(Gamma_d_horiz),sum(mask));zeros(sum(mask),length(Gamma_d_horiz)),Gamma_d_horiz_part];
        end
        fclose(file);
    end
    clear Gamma_d_horiz_part x y z mask horizon_data file
    info.ndata.horiz = length(x_horiz_data);
else
    info.ndata.horiz = 0;
    x_horiz_data = [];
    y_horiz_data = [];
    id_horiz_data = [];
    horiz_data_pts = [];
    Gamma_d_horiz = [];
end
%Combine the three data types:
ndata_cumsum = cumsum([info.ndata.restored_elev,info.ndata.fault,info.ndata.horiz]);
data = zeros(ndata_cumsum(3),1);
data(1:ndata_cumsum(1)) = restored_elev_data_pts;
data(ndata_cumsum(1)+1:ndata_cumsum(2)) = fault_data_pts;
data(ndata_cumsum(2)+1:ndata_cumsum(3)) = horiz_data_pts;
Gamma_d = zeros(ndata_cumsum(3));
Gamma_d(1:ndata_cumsum(1),1:ndata_cumsum(1)) = Gamma_d_restored_elev;
Gamma_d(ndata_cumsum(1)+1:ndata_cumsum(2),ndata_cumsum(1)+1:ndata_cumsum(2)) = Gamma_d_fault;
Gamma_d(ndata_cumsum(2)+1:ndata_cumsum(3),ndata_cumsum(2)+1:ndata_cumsum(3)) = Gamma_d_horiz;
%Create realizations of the data:
if opt.rlzts_from_file
    load(opt.rlzts_path,'rlzts','perturbations')
    if size(rlzts,2)>opt.N
        rlzts = rlzts(:,1:opt.N);
    elseif size(rlzts,2)<opt.N
        error(['Loaded data realizations array must contain at least ',num2str(opt.N),' entries.'])
    end 
    if perturbations
        %perturbations should be true or false to tell if these are perturbations around data (true) or the data + perturbations (false).
        rlzts = data+rlzts;
    end
else
    if opt.correlated_data
        %rlzts = mvnrnd(data,Gamma_d,opt.N)';
        rlzts = lhsnorm(data,Gamma_d,opt.N)'; %Latin hypercube. This can get too memory intensive.
    else
        sigma = sqrt(diag(Gamma_d));
        rlzts = LHCube_UncorNorm(data,sigma,opt.N);
    end
end
%Clear unneeded values to save memory.
clear restored_data_pts fault_data_pts Gamma_d_restored Gamma_d_fault Gamma_d_horiz
if any(Gamma_d==0)
    Gamma_d = sparse(Gamma_d);
end

%Determine the numbers of different kinds of parameters.
info.nparams.general = opt.nhorizons+2; %nhorizons bed depths + 2 slope parameters.
info.nparams.fault = 0;
if opt.fit_surf
    if opt.fit_surf_residuals
        if strcmp(opt.variogram,'matern')
            info.nparams.general = info.nparams.general+6*opt.nfaults;
        else
            info.nparams.general = info.nparams.general+5*opt.nfaults;
        end
        for i = 1:opt.nfaults
            info.nparams.fault = info.nparams.fault+info.faults_ref(i).N;
        end
    else
        info.nparams.general = info.nparams.general+3*opt.nfaults;
    end
end
info.nparams.disp = 0;
if opt.fit_disp
    if opt.fit_disp_residuals
        if strcmp(opt.variogram,'matern')
            info.nparams.general = info.nparams.general+10*opt.nfaults;
        else    
            info.nparams.general = info.nparams.general+9*opt.nfaults;
        end
        for i = 1:opt.nfaults
            info.nparams.disp = info.nparams.disp+info.faults_ref(i).N;
        end
    else
        info.nparams.general = info.nparams.general+7*opt.nfaults;
    end
end
if opt.fit_horizons
    if opt.fit_horizons_residuals
        if strcmp(opt.variogram,'matern')
            info.nparams.general = info.nparams.general+3*opt.nhorizons;
        else
            info.nparams.general = info.nparams.general+2*opt.nhorizons;
        end
        if strcmp(opt.method,'forward')
            info.nparams.horiz = opt.nhorizons*length(info.xgrid);
        else
            info.nparams.horiz = opt.nhorizons*sum(info.n_horiz_pts_block);
        end
    else
        info.nparams.horiz = 0;
    end
    if strcmp(opt.method,'restoration')
        info.nparams.general = info.nparams.general+opt.n_fault_blocks*opt.fit_horizons;
    end
else
    info.nparams.horiz = 0;
end

%Create arrays of mins and maxs for the general parameters (kinematic
%parameters, horizon restored-state depth, and hierarchical parameters).
[info.gen_param_mins,info.gen_param_maxs] = deal([],[]);
info.gen_param_names = {};
info.gen_param_mins = [info.gen_param_mins,opt.horiz1_depth_lims(1),opt.horiz_thickness_lims(:,1)',opt.slopex_lims(1),opt.slopey_lims(1)];
info.gen_param_maxs = [info.gen_param_maxs,opt.horiz1_depth_lims(2),opt.horiz_thickness_lims(:,2)',opt.slopex_lims(2),opt.slopey_lims(2)];
info.gen_param_names = [info.gen_param_names,{[opt.horizon_names{1},' Depth']}];
for i = 1:opt.nhorizons-1
    info.gen_param_names = [info.gen_param_names,{[opt.horizon_names{i},' to ',opt.horizon_names{i+1},' Thickness']}];
end
info.gen_param_names = [info.gen_param_names,{'Slope x','Slope y'}];
if opt.fit_surf
    for i = 1:opt.nfaults
        info.gen_param_mins = [info.gen_param_mins,opt.surf_offset_lims(i,1),opt.strike_lims(i,1),opt.dip_lims(i,1)];
        info.gen_param_maxs = [info.gen_param_maxs,opt.surf_offset_lims(i,2),opt.strike_lims(i,2),opt.dip_lims(i,2)];
        info.gen_param_names = [info.gen_param_names,{[opt.fault_names{i},' Offset'],[opt.fault_names{i},' Strike'],...
            [opt.fault_names{i},' Dip']}];
        if opt.fit_surf_residuals
            info.gen_param_mins = [info.gen_param_mins,opt.surf_range_lims(i,1),opt.surf_std_lims(i,1)];
            info.gen_param_maxs = [info.gen_param_maxs,opt.surf_range_lims(i,2),opt.surf_std_lims(i,2)];
            info.gen_param_names = [info.gen_param_names,{[opt.fault_names{i},' Surface Range'],[opt.fault_names{i},' Surface St Dev']}];
            if strcmp(opt.variogram,'matern') || strcmp(opt.variogram,'dunlop')
                info.gen_param_mins = [info.gen_param_mins,opt.surf_nu_lims(i,1)];
                info.gen_param_maxs = [info.gen_param_maxs,opt.surf_nu_lims(i,2)];
                info.gen_param_names = [info.gen_param_names,{[opt.fault_names{i},' Surface Smoothness']}];
            end
        end
    end
end
if opt.fit_disp
    for i = 1:opt.nfaults
        info.gen_param_mins = [info.gen_param_mins,opt.disp_axis_lims_horiz(i,1),opt.disp_axis_lims_vert(i,1),opt.disp_max_lims(i,1),...
            opt.u0_lims(i,1),opt.v0_lims(i,1),opt.asym_lims(i,1),opt.range_lims(i,1)];
        info.gen_param_maxs = [info.gen_param_maxs,opt.disp_axis_lims_horiz(i,2),opt.disp_axis_lims_vert(i,2),opt.disp_max_lims(i,2),...
            opt.u0_lims(i,2),opt.v0_lims(i,2),opt.asym_lims(i,2),opt.range_lims(i,2)];
        info.gen_param_names = [info.gen_param_names,{[opt.fault_names{i},' Displacement Horizontal Axis'],...
            [opt.fault_names{i},' Displacement Vertical Axis'],[opt.fault_names{i},' Maximum Displacement'],...
            [opt.fault_names{i},' Displacement Center u'],[opt.fault_names{i},' Displacement Center v'],...
            [opt.fault_names{i},' Asymmetry'],[opt.fault_names{i},' Range']}];
        if opt.fit_disp_residuals
            info.gen_param_mins = [info.gen_param_mins,opt.disp_range_lims(i,1),opt.disp_std_lims(i,1)];
            info.gen_param_maxs = [info.gen_param_maxs,opt.disp_range_lims(i,2),opt.disp_std_lims(i,2)];
            info.gen_param_names = [info.gen_param_names,{[opt.fault_names{i},' Displacement Range'],[opt.fault_names{i},' Displacement St Dev']}];
            if strcmp(opt.variogram,'matern') || strcmp(opt.variogram,'dunlop')
                info.gen_param_mins = [info.gen_param_mins,opt.disp_nu_lims(i,1)];
                info.gen_param_maxs = [info.gen_param_maxs,opt.disp_nu_lims(i,2)];
                info.gen_param_names = [info.gen_param_names,{[opt.fault_names{i},' Displacement Smoothness']}];
            end
        end
    end
end
if opt.fit_horizons
    if strcmp(opt.method,'restoration')
        for i = 1:opt.n_fault_blocks
            info.gen_param_mins = [info.gen_param_mins,opt.block_offset_lims(i,1)];
            info.gen_param_maxs = [info.gen_param_maxs,opt.block_offset_lims(i,2)];
            info.gen_param_names = [info.gen_param_names,{['Fault Block ',num2str(i),' Offset']}];
        end
    end
    if opt.fit_horizons_residuals
        for i = 1:opt.nhorizons
            info.gen_param_mins = [info.gen_param_mins,opt.horiz_range_lims(i,1),opt.horiz_std_lims(i,1)];
            info.gen_param_maxs = [info.gen_param_maxs,opt.horiz_range_lims(i,2),opt.horiz_std_lims(i,2)];
            info.gen_param_names = [info.gen_param_names,{[opt.horizon_names{i},' Range'],[opt.horizon_names{i},' StDev']}];
            if strcmp(opt.variogram,'matern') || strcmp(opt.variogram,'dunlop')
                info.gen_param_mins = [info.gen_param_mins,opt.horiz_nu_lims(i,1)];
                info.gen_param_maxs = [info.gen_param_maxs,opt.horiz_nu_lims(i,2)];
                info.gen_param_names = [info.gen_param_names,{[opt.horizon_names{i},' Smoothness']}];
            end
        end
    end
end

if opt.priors_from_file == true %Load the prior parameters from a Matlab file.
    load(opt.prior_path,'params_initial');
    if size(params_initial,2)>opt.N
        params_initial = params_initial(:,1:opt.N);
    elseif size(params_initial,2)<opt.N
        error('Loaded parameters array must contain at least opt.N entries.')
    end 
else %Generate a random prior.
    %Draw sample of all parameters using Latin hypercube sampling.
    n_field_params = info.nparams.fault+info.nparams.disp+info.nparams.horiz;
    n_params_total = info.nparams.general+n_field_params+1; %+1 to include dummy param. If not wanting to do that, use above line.
    sigma = [1.8*ones(1,info.nparams.general),ones(1,n_field_params+1)]'; %+1 to include dummy param.
    params_initial = LHCube_UncorNorm(zeros(1,n_params_total),sigma,opt.N);
    %If using a centered hierarchical parameterization, deal with that.
    if opt.centered_hierarchical
        parfor i = 1:opt.N
            params_initial(:,i) = InterpretCentered(params_initial(:,i),opt,info);
            %This converts the non-centered hierarchical parameterization to a
            %centered one.
        end
    end
end

%Define the function that will used to run the model with the EKI.
func = @(x)RunModel(x,opt,info,fault_data,x_horiz_data,y_horiz_data,horiz_data_pts,id_horiz_data);

%Run the Inversion.
dataset_start_inds = [1,1+ndata_cumsum(1:end-1)];
dataset_start_inds = dataset_start_inds([true,diff(dataset_start_inds)~=0]); %Eliminate any datasets that we're not using.
[params_final,EKI_info] = EKI_DMC(func,opt.N,params_initial,data,rlzts-data,Gamma_d,opt.max_iterations,...
    'boot_screen',opt.boot,'Nboot',opt.nboot,'sigma_alpha',opt.sigma_alpha,'load_boot_inds',opt.load_boot_inds,...
    'boot_inds_file',opt.boot_inds_file,'inflate',opt.infl,'Nz',opt.nz_infl,'persistent_z',opt.persistent_z,...
    'load_z',opt.load_inflation_z,'z_file',opt.inflation_z_file,'algorithm',opt.algorithm,'Nthreads',opt.Nthreads,...
    'dataset_start_inds',dataset_start_inds);

%Add output from EKI_DMC to the info structure.
info.RMS = EKI_info.RMS;
info.alpha = EKI_info.alpha;
info.gamma = EKI_info.gamma;
info.wall_time = EKI_info.wall_time;
info.cpu_time = EKI_info.cpu_time;
end