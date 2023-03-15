function params = InterpretCentered(hierarchical_params,opt,info)
%Make an initial realization for a centered hierarchical parameterization.
%This is based on InterpretModel.

params = hierarchical_params; %The general params part will be returned the same, while the rest will be changed.
gen_params = unbounded2bounded(hierarchical_params(1:info.nparams.general),...
    info.gen_param_mins',info.gen_param_maxs');
n = opt.nhorizons+2; %This keeps track of the last parameter dealt with to know where to start next.
nparams_part = info.nparams.general; %This keeps track of the last parameter for adding new field parameters to params.
if opt.fit_surf
    ind_start = nparams_part+1;
    for i = 1:opt.nfaults
        n = n+3;
        if opt.fit_surf_residuals
            [surf_range,stdev] = deal(gen_params(n+1),gen_params(n+2));
            n = n+2;
            if strcmp(opt.variogram,'matern')
                surf_nu = gen_params(n+1);
                n = n+1;
                vparams = [surf_range,surf_nu];
            else %Circular variogram.
                vparams = surf_range;
            end
            params(ind_start:ind_start+info.faults_ref(i).N-1) = stdev*MakeGRF(info.faults_ref(i).u(:),info.faults_ref(i).v(:),vparams,...
                hierarchical_params(ind_start:ind_start+info.faults_ref(i).N-1),opt.variogram,opt.rlzt_method);
            ind_start = ind_start+info.faults_ref(i).N;
        end
    end
    nparams_part = nparams_part+info.nparams.fault;
end
if opt.fit_disp
    ind_start = nparams_part+1;
    for i = 1:opt.nfaults
        n = n+7;
        if opt.fit_disp_residuals
            [disp_range,stdev] = deal(gen_params(n+1),gen_params(n+2));
            n = n+2;
            if strcmp(opt.variogram,'matern')
                disp_nu = gen_params(n+1);
                n = n+1;
                vparams = [disp_range,disp_nu];
            else %Circular variogram.
                vparams = disp_range;
            end
            params(ind_start:ind_start+info.faults_ref(i).N-1) = stdev*MakeGRF(info.faults_ref(i).u(:),info.faults_ref(i).v(:),vparams,...
                hierarchical_params(ind_start:ind_start+info.faults_ref(i).N-1),opt.variogram,opt.rlzt_method);
            ind_start = ind_start+info.faults_ref(i).N;
        end
    end
    nparams_part = nparams_part+info.nparams.disp;
end
if opt.fit_horizons>0
    if strcmp(opt.method,'restoration') %Fault block offset is the same for all horizons, so we only have to do this once.
        n = n+opt.n_fault_blocks;
    end
    for i = 1:opt.nhorizons
        if opt.fit_horizons_residuals
            [horiz_range,stdev] = deal(gen_params(n+1),gen_params(n+2));
            n = n+2;
            if strcmp(opt.variogram,'matern')
                horiz_nu = gen_params(n+1);
                n = n+1;
                vparams = [horiz_range,horiz_nu];
            else %Circular variogram.
                vparams = horiz_range;
            end
        end
        if strcmp(opt.method,'forward')
            if opt.fit_horizons_residuals
                ind_start = nparams_part+(i-1)*info.n_horiz_pts+1;
                params(ind_start:ind_start+info.n_horiz_pts-1) = stdev*MakeGRF(info.xgrid,info.ygrid,vparams,...
                    hierarchical_params(ind_start:ind_start+info.n_horiz_pts-1),opt.variogram,opt.rlzt_method);
            end
        else
            for j = 1:opt.n_fault_blocks
                if opt.fit_horizons_residuals
                    ind_start = nparams_part+(i-1)*sum(info.n_horiz_pts_block)+sum(info.n_horiz_pts_block(1:j-1))+1;
                    params(ind_start:ind_start+info.n_horiz_pts_block(j)-1) = stdev*MakeGRF(info.xgrid(info.grid_pts_mask{j}),info.ygrid(info.grid_pts_mask{j}),vparams,...
                        hierarchical_params(ind_start:ind_start+info.n_horiz_pts_block(j)-1),opt.variogram,opt.rlzt_method);
                end
            end
        end
    end
end
end