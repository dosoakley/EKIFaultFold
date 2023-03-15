function predicted_data = RunModel(modelparams,opt,info,fault_data,x_horiz_data,...
    y_horiz_data,z_horiz_data,id_horiz_data)
%This is a function to create the horizon and fault model and predict the 
%values at the data points.
%Disable the duplicate point interpolation warning.
warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');
warning('off','MATLAB:scatteredInterpolant:InterpEmptyTri2DWarnId');
%Read the model parameters and assign them to a new fault structure:
[gen_params,faults,horiz_params] = InterpretModel(modelparams,opt,info);
restored_depths = cumsum(gen_params(1:opt.nhorizons));
slopes = gen_params(opt.nhorizons+(1:2));
if opt.fit_horizons
    model_horizons = horiz_params;
    clear horiz_params
else
    model_horizons = info.horiz_pts_ref;
end

%Calculate the predicted data.
predicted_data = []; %Initialize.
if info.ndata.restored_elev>0 %Fitting for restored-state horizon geometry.
    %Evaluate how well the restored-state horizons fit a flat plane.
    misfitj = zeros(opt.nhorizons*info.n_horiz_pts,1);
    for k = 1:opt.nhorizons
        if strcmp(opt.method,'forward')
            %Horizon parameters define the restored state directly.
            z_restored_grid = model_horizons(1+(k-1)*info.n_horiz_pts:k*info.n_horiz_pts); %z for this horizon.
        else
            %Restore the points.
            [x,y,z] = deal([],[],[]);
            for j = 1:opt.n_fault_blocks
                xj = info.xgrid(info.grid_pts_mask{j})';
                yj = info.ygrid(info.grid_pts_mask{j})';
                ind_start = (k-1)*sum(info.n_horiz_pts_block)+sum(info.n_horiz_pts_block(1:j-1))+1;
                zj = model_horizons(ind_start:ind_start+info.n_horiz_pts_block(j)-1)';
                [xj,yj,zj] = deal(xj',yj',zj');
                model_blocks = AssignFaultBlock(xj,yj,zj,faults,opt.fault_blocks_relationships);
                active = model_blocks==j;
                [x,y,z] = deal([x;xj(active)],[y;yj(active)],[z;zj(active)]);
            end
            [x,y,z,~,success] = MoveModel(x',y',z',faults,opt.method,opt.fault_deformation,[opt.sigma_misfit_restored,opt.sigma_fault]);
            z_restored_grid = InterpScatteredData(x(success)',y(success)',z(success)',info.xgrid,info.ygrid);
        end
        xg = info.xgrid-(max(info.xgrid)+min(info.xgrid))/2; %Center x and y.
        yg = info.ygrid-(max(info.ygrid)+min(info.ygrid))/2;
        misfitj(1+(k-1)*info.n_horiz_pts:k*info.n_horiz_pts) = z_restored_grid-(restored_depths(k)+slopes(1)*xg+slopes(2)*yg);
    end
    predicted_data = [predicted_data;misfitj];
end
if info.ndata.fault>0
    if strcmp(opt.method,'forward') %faults needs to be changed from restored to deformed state.
        [~,~,~,faults_deformed,~] = MoveModel([],[],[],faults,opt.method,opt.fault_deformation,[]);
    end
    for i = 1:opt.nfaults
        fault_data_uv = faults(i).M\[fault_data{i,1},fault_data{i,2},fault_data{i,3},ones(size(fault_data{i,1}))]';
        fault_data_pts = fault_data_uv(3,:)';
        fault_data_u = fault_data_uv(1,:)';
        fault_data_v = fault_data_uv(2,:)';
        if strcmp(opt.method,'forward')
            model_surf_data = InterpolateSurface(faults_deformed(i),fault_data_u,fault_data_v); %Calculate the fault surface at the fault data points.
        else %faults is already in deformed state
            model_surf_data = InterpolateSurface(faults(i),fault_data_u,fault_data_v); %Calculate the fault surface at the fault data points.
        end
        model_surf_misfit = model_surf_data-fault_data_pts;
        predicted_data = [predicted_data;model_surf_misfit]; %Start the predicted data vector.
    end
end
if  info.ndata.horiz>0
    %Calculate the horizon surfaces at the horizon data points.
    for i = 1:opt.nhorizons
        if strcmp(opt.method,'forward')
            %Forward model the points.
            x = info.xgrid';
            y = info.ygrid';
            z = model_horizons(1+(i-1)*info.n_horiz_pts:i*info.n_horiz_pts)';
            [x,y,z,faults_deformed,~] = MoveModel(x,y,z,faults,opt.method,opt.fault_deformation,[]);
            [x,y,z] = deal(x',y',z');
            %Assign fault blocks to model and data points.
            %Note: For the data points, this uses the true values, not the
            %realization, since that's simpler.
            block_model = AssignFaultBlock(x,y,z,faults_deformed,opt.fault_blocks_relationships);
            mask1 = id_horiz_data==i;
            block_data = AssignFaultBlock(x_horiz_data(mask1),y_horiz_data(mask1),z_horiz_data(mask1),faults_deformed,opt.fault_blocks_relationships);
            %Loop through the fault blocks, doing the interpolation
            %separately in each fault block.
            model_horiz_data = zeros(sum(id_horiz_data==i),1);
            for j = 1:opt.n_fault_blocks
                mask2 = block_model==j;
                mask3 = mask1;
                mask3(mask1) = block_data==j;
                if sum(mask2)>0
                    model_horiz_data(block_data==j) = InterpScatteredData(x(mask2),y(mask2),z(mask2),x_horiz_data(mask3),y_horiz_data(mask3));
                else
                    %Occasionally, we can get no model points in the fault
                    %block. We need some value, so in this case, we will
                    %just interpolate from all model points without regard
                    %to fault block.
                    model_horiz_data(block_data==j) = InterpScatteredData(x,y,z,x_horiz_data(mask3),y_horiz_data(mask3)); 
                end
            end
%             model_horiz_data = InterpScatteredData(x,y,z,x_horiz_data(mask1),y_horiz_data(mask1)); 
        else
            %Horizon parameters define the deformed state directly.
            mask1 = id_horiz_data==i;
            block_data = AssignFaultBlock(x_horiz_data(mask1),y_horiz_data(mask1),z_horiz_data(mask1),faults,opt.fault_blocks_relationships);
            model_horiz_data = zeros(sum(id_horiz_data==i),1);
            for j = 1:opt.n_fault_blocks
                x = info.xgrid(info.grid_pts_mask{j});
                y = info.ygrid(info.grid_pts_mask{j});
                ind_start = (i-1)*sum(info.n_horiz_pts_block)+sum(info.n_horiz_pts_block(1:j-1))+1;
                z = model_horizons(ind_start:ind_start+info.n_horiz_pts_block(j)-1);
                mask3 = mask1;
                mask3(mask1) = block_data==j;
                model_horiz_data(block_data==j) = InterpScatteredData(x,y,z,x_horiz_data(mask3),y_horiz_data(mask3));
            end
        end
        predicted_data = [predicted_data;model_horiz_data];
    end
end
if any(isnan(predicted_data))
    disp('Error: NaNs')
end
end