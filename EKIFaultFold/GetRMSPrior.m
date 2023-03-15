%Read the model realizations made in RMS and convert them to a prior
%ensemble for EKI.

%A lot at the beginning of this is copied from RunEKI.

%Note about what the options file should look like:
%1) Variogram should be spherical, not Matern. If non-centered,
%hierarchical parameters should be converted to Gaussian random fields by
%Cholesky decomposition of the covariance matrix.
%2) Everything including residuals should be fit for.
%3) Method should be restoration.

%Load the realizations.
% load('D:\EmeraldFieldRealizations\Realizations.mat');
load('..\..\EmeraldFieldRealizations\Realizations.mat');

%Run the options file we want to use.
% run('Options_EmeraldField_Small_restoration_RMSPrior.m')
run('Options_EmeraldField_Small_restoration_RMSPrior_v2.m')

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

%Figure out the numbers of parameters.
info.nparams.general = opt.nhorizons+2; %nhorizons bed depths + 2 slope parameters.
info.nparams.general = info.nparams.general+5*opt.nfaults; %For fault surface
info.nparams.general = info.nparams.general+9*opt.nfaults; %For displacement
info.nparams.general = info.nparams.general+2*opt.nhorizons+opt.n_fault_blocks*opt.fit_horizons; %For horizons
info.nparams.fault = 0;
for i = 1:opt.nfaults
    info.nparams.fault = info.nparams.fault+info.faults_ref(i).N;
end
info.nparams.disp = 0;
for i = 1:opt.nfaults
    info.nparams.disp = info.nparams.disp+info.faults_ref(i).N;
end
info.nparams.horiz = opt.nhorizons*sum(info.n_horiz_pts_block);

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
for i = 1:opt.nfaults
    info.gen_param_mins = [info.gen_param_mins,opt.surf_offset_lims(i,1),opt.strike_lims(i,1),opt.dip_lims(i,1)];
    info.gen_param_maxs = [info.gen_param_maxs,opt.surf_offset_lims(i,2),opt.strike_lims(i,2),opt.dip_lims(i,2)];
    info.gen_param_names = [info.gen_param_names,{[opt.fault_names{i},' Offset'],[opt.fault_names{i},' Strike'],...
        [opt.fault_names{i},' Dip']}];
    info.gen_param_mins = [info.gen_param_mins,opt.surf_range_lims(i,1),opt.surf_std_lims(i,1)];
    info.gen_param_maxs = [info.gen_param_maxs,opt.surf_range_lims(i,2),opt.surf_std_lims(i,2)];
    info.gen_param_names = [info.gen_param_names,{[opt.fault_names{i},' Surface Range'],[opt.fault_names{i},' Surface St Dev']}];
end
for i = 1:opt.nfaults
    info.gen_param_mins = [info.gen_param_mins,opt.disp_axis_lims_horiz(i,1),opt.disp_axis_lims_vert(i,1),opt.disp_max_lims(i,1),...
        opt.u0_lims(i,1),opt.v0_lims(i,1),opt.asym_lims(i,1),opt.range_lims(i,1)];
    info.gen_param_maxs = [info.gen_param_maxs,opt.disp_axis_lims_horiz(i,2),opt.disp_axis_lims_vert(i,2),opt.disp_max_lims(i,2),...
        opt.u0_lims(i,2),opt.v0_lims(i,2),opt.asym_lims(i,2),opt.range_lims(i,2)];
    info.gen_param_names = [info.gen_param_names,{[opt.fault_names{i},' Displacement Horizontal Axis'],...
        [opt.fault_names{i},' Displacement Vertical Axis'],[opt.fault_names{i},' Maximum Displacement'],...
        [opt.fault_names{i},' Displacement Center u'],[opt.fault_names{i},' Displacement Center v'],...
        [opt.fault_names{i},' Asymmetry'],[opt.fault_names{i},' Range']}];
    
    info.gen_param_mins = [info.gen_param_mins,opt.disp_range_lims(i,1),opt.disp_std_lims(i,1)];
    info.gen_param_maxs = [info.gen_param_maxs,opt.disp_range_lims(i,2),opt.disp_std_lims(i,2)];
    info.gen_param_names = [info.gen_param_names,{[opt.fault_names{i},' Displacement Range'],[opt.fault_names{i},' Displacement St Dev']}];
end
for i = 1:opt.n_fault_blocks
    info.gen_param_mins = [info.gen_param_mins,opt.block_offset_lims(i,1)];
    info.gen_param_maxs = [info.gen_param_maxs,opt.block_offset_lims(i,2)];
    info.gen_param_names = [info.gen_param_names,{['Fault Block ',num2str(i),' Offset']}];
end
if opt.fit_horizons_residuals
    for i = 1:opt.nhorizons
        info.gen_param_mins = [info.gen_param_mins,opt.horiz_range_lims(i,1),opt.horiz_std_lims(i,1)];
        info.gen_param_maxs = [info.gen_param_maxs,opt.horiz_range_lims(i,2),opt.horiz_std_lims(i,2)];
        info.gen_param_names = [info.gen_param_names,{[opt.horizon_names{i},' Range'],[opt.horizon_names{i},' StDev']}];
    end
end

gen_params = zeros(info.nparams.general,opt.N);
fault_params = zeros(info.nparams.fault,opt.N);
disp_params = zeros(info.nparams.disp,opt.N);
horiz_params = zeros(info.nparams.horiz,opt.N);
for n = 1:opt.N
    disp(n)
    
    %Get the fault surface parameters.
    ng = opt.nhorizons+2; %This keeps track of the last general parameter dealt with to know where to start next.
    ind_start = 1;
    faults_new = [];
    for i = 1:opt.nfaults
        fault = copy(faults{n}(i));
        
%         %There was an error in ReadRMS_Realizations that switched asymmetry
%         %and range, so I'm fixing it here instead of rerunning
%         %ReadRMS_Realizations, which I can't do on my home computer.
%         asym = fault.Range;
%         fault.Range = fault.asymmetry;
%         fault.asymmetry = asym;
%         %This has been fixed now and ReadRMS_Realizations has been rerun.
        
        fault_new = copy(info.faults_ref(i));
        %Transform the fault to the reference fault center.
        fault_xyz = fault.M*[fault.u(:),fault.v(:),fault.f(:),ones(size(fault.u(:)))]';
        fault_new.M = fault.M;
        fault_new.M(1:3,4) = info.faults_ref(i).M(1:3,4);
        fault_uv = fault_new.M\fault_xyz;
        %Interpolate onto the new fault grid.
        fault_new.f =  InterpScatteredData(fault_uv(1,:)',fault_uv(2,:)',fault_uv(3,:)',fault_new.u,fault_new.v);
        fault_new.displacement =  InterpScatteredData(fault_uv(1,:)',fault_uv(2,:)',fault.displacement(:),fault_new.u,fault_new.v);
        fault_new.ThicknessAttribute = InterpScatteredData(fault_uv(1,:)',fault_uv(2,:)',fault.ThicknessAttribute(:),fault_new.u,fault_new.v);
        fault_new.asymmetry = fault.asymmetry;
        fault_new.Range = fault.Range;
        %Find the strike and dip.
        strike = atan2d(fault_new.M(1,1),fault_new.M(2,1));
        dip = atan2d(-fault_new.M(3,2),-fault_new.M(3,3));
        %Find the surface trend (offset) and make the rest be the residuals.
        offset = mean(fault_new.f(:));
        gen_params(ng+1:ng+3,n) = [offset,strike,dip];
        ng = ng+3;
        residuals = fault_new.f(:)-offset;
        %Fit a variogram to the surface.
        pts = [fault_new.u(:),fault_new.v(:)];
        data = fault_new.f(:);
        v = variogram(pts,data,'nrbins',50);
        a0 = 1e3; %Initial range value.
        c0 = 1e2; %Initial sill value.
        [range,sill,~,~] = variogramfit(v.distance,v.val,a0,c0,v.num,'model','spherical','plotit',false);
        gen_params(ng+1:ng+2,n) = [range,sqrt(sill)];
        ng = ng+2;
        if opt.centered_hierarchical
            fault_params(ind_start:ind_start+faults{n}(i).N-1,n) = residuals;
        else
            c = sill*CorrSpher(fault_new.u(:),fault_new.v(:),range);
            L = chol(c,'lower');
            aux = L\residuals;
            fault_params(ind_start:ind_start+fault_new.N-1,n) = aux;
        end
        faults_new = [faults_new,fault_new];
        ind_start = ind_start+faults{n}(i).N;
    end
        
    %Get the displacement parameters.
    ind_start = 1;
    for i = 1:opt.nfaults
        fault_new = faults_new(i);
%         %Find the displacement trend.
% %         options = optimoptions('fsolve','Display','off','Algorithm','levenberg-marquardt');
%         fun = @(x) EllipseMisfit(fault_new.u(:),fault_new.v(:),fault_new.displacement(:),x);
%         x0 = [faults{n}(i).lu,faults{n}(i).lv,faults{n}(i).max_val,faults{n}(i).u0,faults{n}(i).v0];
% %         ellipse_params = fsolve(fun,x0,options);
%         options = optimoptions('fmincon','Display','off');
%         lb = [0,0,0,-1e5,-1e5];
%         ub = [1e5,1e5,1e5,1e5,1e5];
%         if faults{n}(i).max_val<0
%             x0(3) = max(faults{n}(i).displacement);
%         end
%         if any((x0<lb | x0>ub))
%             disp('Out of Bounds')
%         end
%         ellipse_params = fmincon(fun,x0,[],[],[],[],lb,ub,[],options);
        %Use Havana's values. Don't fit.
        ellipse_params = [faults{n}(i).lu,faults{n}(i).lv,faults{n}(i).max_val,faults{n}(i).u0,faults{n}(i).v0];
%         if faults{n}(i).max_val<0
%             disp('<0')
%         end
        gen_params(ng+1:ng+7,n) = [ellipse_params,faults{n}(i).asymmetry,faults{n}(i).Range];
        ng = ng+7;
        trend = EllipseTrend(fault_new.u(:),fault_new.v(:),ellipse_params);
        residuals = fault_new.displacement(:)-trend;
        [fault_new.lu,fault_new.lv] = deal(ellipse_params(1),ellipse_params(2));
        fault_new.max_val = ellipse_params(3);
        [fault_new.u0,fault_new.v0] = deal(ellipse_params(1),ellipse_params(2));
        %Fit a variogram to the displacement, excluding anything outside the tip line.
        mask = fault_new.ThicknessAttribute>0;
        pts = [fault_new.u(mask),fault_new.v(mask)];
        data = residuals(mask(:));
        v = variogram(pts,data,'nrbins',50);
        a0 = 1e3; %Initial range value.
        c0 = 1e2; %Initial sill value.
        [range,sill,~,~] = variogramfit(v.distance,v.val,a0,c0,v.num,'model','spherical','plotit',false);
        if sill == 0
            sill = 0.1; %0 causes problems, so just chage it.
        end
        gen_params(ng+1:ng+2,n) = [range,sqrt(sill)];
        ng = ng+2;
        if opt.centered_hierarchical
            disp_params(ind_start:ind_start+faults{n}(i).N-1,n) = residuals;
        else
            c = sill*CorrSpher(fault_new.u(:),fault_new.v(:),range);
            L = chol(c,'lower');
            aux = L\residuals;
            disp_params(ind_start:ind_start+fault_new.N-1,n) = aux;
        end
        ind_start = ind_start+faults{n}(i).N;
        faults_new(i) = fault_new;
    end
    
    %Get the horizon parameters.
    ind_start = 1;
    for i = 1:opt.nhorizons
        [x,y,z] = deal(horizons{n}.x{i},horizons{n}.y{i},horizons{n}.z{i});
        FaultBlock = horizons{n}.FaultBlock{i}+1;
        [x_interp,y_interp,z_interp,block_interp] = deal([],[],[],[]);
        for j = 1:opt.n_fault_blocks
            mask = FaultBlock==j;
            interpolant = scatteredInterpolant(x(mask),y(mask),z(mask));
            z_interp_block = interpolant(info.xgrid(info.grid_pts_mask{j}),info.ygrid(info.grid_pts_mask{j}));
            [x_interp,y_interp,z_interp,block_interp] = deal([x_interp;info.xgrid(info.grid_pts_mask{j})],...
                [y_interp;info.ygrid(info.grid_pts_mask{j})],[z_interp;z_interp_block],[block_interp;j*ones(size(z_interp_block))]);
        end
        %Estimate the restored state depths / thicknesses:
        block = AssignFaultBlock(x_interp,y_interp,z_interp,faults_new,opt.fault_blocks_relationships);
        active = block==block_interp;
        [x0,y0,z0,~] = MoveModel(x_interp(active)',y_interp(active)',z_interp(active)',faults_new,opt.method,opt.fault_deformation);
        xcenter = (max(info.xgrid)+min(info.xgrid))/2;
        ycenter = (max(info.ygrid)+min(info.ygrid))/2;
        xg = x0-xcenter; %Center x and y.
        yg = y0-ycenter;
        if i==1
            A = [xg',yg',ones(size(xg'))];
            b = A\z0';
            gen_params(1,n) = b(3);
            gen_params(opt.nhorizons+(1:2),n) = [b(1:2)];
        else
            gen_params(i,n) = mean(z0-(b(1)*xg+b(2)*yg+b(3))); %Distance of this restored state plane above the last one.
        end
        %Estimate the block offsets:
        %(We're just using the first horizon to do this.)
        if i == 1
            offset = zeros(1,opt.n_fault_blocks);
            for j = 1:opt.n_fault_blocks
                mask = block_interp==j;
                offset(j) = mean(z_interp(mask)-(b(1)*(x_interp(mask)-xcenter)+b(2)*(y_interp(mask)-ycenter)+b(3)));
                gen_params(ng+1,n) = offset(j);
                ng = ng+1;
            end
        end
        %Calculate the residuals.
%         interpolant = scatteredInterpolant(x_interp(active),y_interp(active),z_interp(active));
%         zgrid = interpolant(info.xgrid,info.ygrid);
%         block_grid = AssignFaultBlock(info.xgrid,info.ygrid,zgrid,faults_new,opt.fault_blocks_relationships);
        residuals = zeros(size(z_interp));
        for j = 1:opt.n_fault_blocks
            if i == 1
                residuals(block_interp==j) = z_interp(block_interp==j)-(b(1)*(x_interp(block_interp==j)-xcenter)+...
                    b(2)*(y_interp(block_interp==j)-ycenter)+b(3)+offset(j));
            else
                residuals(block_interp==j) = z_interp(block_interp==j)-z_interp_last(block_interp==j);
            end
        end
        z_interp_last = z_interp;
%         if i == 1
%             residuals = zeros(size(zgrid));
%             for j = 1:opt.n_fault_blocks
%                 residuals(block_grid==j) = zgrid(block_grid==j)-(b(1)*(info.xgrid-xcenter)+b(2)*(info.ygrid-ycenter)+b(3)+offset(j));
%             end
%         else
%             residuals = zgrid-zgrid_last;
%         end
%         zgrid_last = zgrid;
        %Estimate the variogram:
        %This isn't exactly right because of the faults, but hopefully it's
        %close enough. The offsets should help to account for that to some
        %extent.
        pts = [x_interp,y_interp];
        data = residuals;
        v = variogram(pts,data,'nrbins',50);
        a0 = 1e3; %Initial range value.
        c0 = 1e2; %Initial sill value.        
        [range,sill,~,~] = variogramfit(v.distance,v.val,a0,c0,v.num,'model','spherical','plotit',false);
        if sill == 0
            sill = 0.1; %0 causes problems, so just chage it.
        end
        gen_params(ng+1:ng+2,n) = [range,sqrt(sill)];
        ng = ng+2;
        npts = length(residuals);
        if opt.centered_hierarchical
            horiz_params(ind_start:ind_start+npts-1,n) = residuals;
        else
            for j = 1:opt.n_fault_blocks
                c = sill*CorrSpher(x_interp(block_interp==j),y_interp(block_interp==j),range);
                L = chol(c,'lower');
                aux = L\residuals(block_interp==j);
                ind_start_j = (i-1)*sum(info.n_horiz_pts_block)+sum(info.n_horiz_pts_block(1:j-1))+1;
                horiz_params(ind_start_j:ind_start_j+info.n_horiz_pts_block(j)-1,n) = aux;
            end
        end
        ind_start = ind_start+npts;
    end
end

for i = 1:info.nparams.general
    disp([info.gen_param_names{i},': ',num2str(mean(gen_params(i,:))),' +/- ',num2str(std(gen_params(i,:))),...
        '     Range: ',num2str(min(gen_params(i,:))),' to ',num2str(max(gen_params(i,:)))])
end

%Convert to unbounded domain.
gen_params = bounded2unbounded(gen_params,info.gen_param_mins',info.gen_param_maxs');
params_initial = [gen_params;fault_params;disp_params;horiz_params];
% save('RMSPrior.mat','params_initial');
save('RMSPrior_v2.mat','params_initial');

function misfit = EllipseMisfit(u,v,d,params)
trend = EllipseTrend(u,v,params);
misfit = sum((trend-d).^2);
end

function trend = EllipseTrend(u,v,params)
[lu,lv,max_val,u0,v0] = deal(params(1),params(2),params(3),params(4),params(5));
r = sqrt(((u-u0)./lu).^2+((v-v0)./lv).^2);
mu0 = zeros(size(r));
mu0(r<=1) = 2*(1-r(r<=1)).*sqrt(((1+r(r<=1)).^2)./4-r(r<=1).^2); %Georgsen (2012) Eqn. 3.
trend =  max_val*mu0; %Georgsen (2012) Eqn. 4
end
