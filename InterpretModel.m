function [gen_params,faults,horiz_params] = InterpretModel(hierarchical_params,opt,info)
%Interpret the model parameters from the form that the EnKF is fitting to
%the actual meaningful form.
faults = copy(info.faults_ref);
gen_params = unbounded2bounded(hierarchical_params(1:info.nparams.general),...
    info.gen_param_mins',info.gen_param_maxs');
restored_depths = cumsum(gen_params(1:opt.nhorizons));
slopes = gen_params(opt.nhorizons+(1:2));
thicknesses = diff(restored_depths);
n = opt.nhorizons+2; %This keeps track of the last parameter dealt with to know where to start next.
nparams_part = info.nparams.general; %This keeps track of the last parameter for adding new field parameters to params.
if opt.fit_surf
    ind_start = nparams_part+1;
    for i = 1:opt.nfaults
        [offset,strike,dip] = deal(gen_params(n+1),gen_params(n+2)*pi/180,...
            gen_params(n+3)*pi/180);
        n = n+3;
        [cx,cy,cz] = deal(faults(i).M(1,4),faults(i).M(2,4),faults(i).M(3,4));
        MakeM(faults(i),cx,cy,cz,strike,dip);
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
            SurfFromAux(faults(i),offset,vparams,stdev,hierarchical_params(ind_start:ind_start+faults(i).N-1),opt)
            ind_start = ind_start+faults(i).N;
        else
            faults(i).f(:) = offset;
        end
    end
    nparams_part = nparams_part+info.nparams.fault;
end
if opt.fit_disp
    ind_start = nparams_part+1;
    for i = 1:opt.nfaults
        [lu,lv,max_val,u0,v0,asymmetry,range] = deal(gen_params(n+1),gen_params(n+2),...
            gen_params(n+3),gen_params(n+4),gen_params(n+5),gen_params(n+6),gen_params(n+7));
        n = n+7;
        [faults(i).asymmetry,faults(i).Range] = deal(asymmetry,range);
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
            DispFromAux(faults(i),lu,lv,max_val,u0,v0,vparams,stdev,hierarchical_params(ind_start:ind_start+faults(i).N-1),opt)
            ind_start = ind_start+faults(i).N;
        else
            [faults(i).lu,faults(i).lv] = deal(lu,lv);
            faults(i).max_val = max_val;
            [faults(i).u0,faults(i).v0] = deal(u0,v0);
            faults(i).displacement = EllipticalTrend(faults(i),faults(i).u,faults(i).v);
            CalculateThicknessAttribute(faults(i))
        end
    end
    nparams_part = nparams_part+info.nparams.disp;
end
if opt.nfaults>1 && (opt.fit_surf || opt.fit_disp)
    faults = ApplyTruncations(faults,opt.fault_truncations);
end
if opt.fit_horizons>0
    if strcmp(opt.method,'restoration') %Fault block offset is the same for all horizons, so we only have to do this once.
        fault_block_offsets = gen_params(n+1:n+opt.n_fault_blocks);
        n = n+opt.n_fault_blocks;
        horiz_params = zeros(opt.nhorizons*sum(info.n_horiz_pts_block),1); %This will equal info.nparams.horiz if fit_horizons_residuals is true.
    else
        horiz_params = zeros(opt.nhorizons*length(info.xgrid),1); %This will equal info.nparams.horiz if fit_horizons_residuals is true.
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
            if i == 1
                x = info.xgrid-(max(info.xgrid)+min(info.xgrid))/2; %Center x and y.
                y = info.ygrid-(max(info.ygrid)+min(info.ygrid))/2;
                trend = restored_depths(i)+slopes(1)*x+slopes(2)*y;
            else
                start_end_inds = (i-2)*info.n_horiz_pts+[1,info.n_horiz_pts];
                trend = horiz_params(start_end_inds(1):start_end_inds(2)) + thicknesses(i-1); %Previous horizon plus expected thickness.
            end
            if opt.fit_horizons_residuals
                ind_start = nparams_part+(i-1)*info.n_horiz_pts+1;
                if opt.centered_hierarchical
                    horiz_params((i-1)*info.n_horiz_pts+1:i*info.n_horiz_pts) = ...
                        trend+hierarchical_params(ind_start:ind_start+info.n_horiz_pts-1);
                else
                    u = MakeGRF(info.xgrid,info.ygrid,vparams,...
                        hierarchical_params(ind_start:ind_start+info.n_horiz_pts-1),opt.variogram,opt.rlzt_method);
                    horiz_params((i-1)*info.n_horiz_pts+1:i*info.n_horiz_pts) = trend+stdev*u(:);
                end
            else
                horiz_params((i-1)*info.n_horiz_pts+1:i*info.n_horiz_pts) = trend;
            end
        else
            for j = 1:opt.n_fault_blocks
                if i == 1
                    x = info.xgrid(info.grid_pts_mask{j})-(max(info.xgrid(info.grid_pts_mask{j}))+min(info.xgrid(info.grid_pts_mask{j})))/2; %Center x and y.
                    y = info.ygrid(info.grid_pts_mask{j})-(max(info.ygrid(info.grid_pts_mask{j}))+min(info.ygrid(info.grid_pts_mask{j})))/2;
                    trend = restored_depths(i)+slopes(1)*x+slopes(2)*y+fault_block_offsets(j);
                else
                    ind_start_lasti = (i-2)*sum(info.n_horiz_pts_block)+sum(info.n_horiz_pts_block(1:j-1))+1;
                    trend = horiz_params(ind_start_lasti:ind_start_lasti+info.n_horiz_pts_block(j)-1)+thicknesses(i-1)+fault_block_offsets(j); %Previous horizon plus expected thickness.
                end
                if opt.fit_horizons_residuals
                    ind_start = nparams_part+(i-1)*sum(info.n_horiz_pts_block)+sum(info.n_horiz_pts_block(1:j-1))+1;
                    ind_start_horiz = ind_start-nparams_part;
                    if opt.centered_hierarchical
                        horiz_params(ind_start_horiz:ind_start_horiz+info.n_horiz_pts_block(j)-1) = ...
                            trend+hierarchical_params(ind_start:ind_start+info.n_horiz_pts_block(j)-1);
                    else
                        u = MakeGRF(info.xgrid(info.grid_pts_mask{j}),info.ygrid(info.grid_pts_mask{j}),vparams,...
                            hierarchical_params(ind_start:ind_start+info.n_horiz_pts_block(j)-1),opt.variogram,opt.rlzt_method);
                        horiz_params(ind_start_horiz:ind_start_horiz+info.n_horiz_pts_block(j)-1) = trend+stdev*u(:);
                    end
                else
                    horiz_params(ind_start_horiz:ind_start_horiz+info.n_horiz_pts_block(j)-1) = trend;
                end
            end
        end
    end
else
    horiz_params = [];
end
end