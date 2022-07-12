%Run the Test6 Ensemble Kalman Inversion.

%I am running this using EnKF_Combined_Code14.
%The data being fit are horizon points, fault points, and restored-state.

% options_file_name = 'Options_Test6_Scattered_Forward.m';
% options_file_name = 'Options_Test6_Dense_Forward.m';
options_file_name = 'Options_Test6_Dense_Restoration.m';
% options_file_name = 'Options_Test6_Dense_Forward_old_lims.m';

run(options_file_name);
[params_initial_raw,params_final_raw,info] = RunEKI(opt); %Run the inversion.

%%
% run(options_file_name); %Do this to get the values that are contained in the options file.

%Disable the duplicate point interpolation warning that will occur a lot
%during the plotting otherwise.
warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');

%Transform the final parameters back to the bounded domain.
N = size(params_final_raw,2); %Ensemble size
gen_params_initial = zeros(info.nparams.general,N);
gen_params_final = zeros(info.nparams.general,N);
faults_init = fault_class.empty();
faults_final = fault_class.empty();
horiz_params_init = zeros(info.nparams.horiz,N);
horiz_params_final = zeros(info.nparams.horiz,N);
parfor i = 1:N
% for i = 1:N
    [gen_params_initial(:,i),faults_init(i),horiz_params_init(:,i)] = InterpretModel(params_initial_raw(:,i),opt,info);
    [gen_params_final(:,i),faults_final(i),horiz_params_final(:,i)] = InterpretModel(params_final_raw(:,i),opt,info);
end

%Display the general parameters.
for i = 1:info.nparams.general
    disp([info.gen_param_names{i},': ',num2str(mean(gen_params_final(i,:))),' +/- ',num2str(std(gen_params_final(i,:)))])
end

%Get the grid dimensions.
[nxgrid,nygrid] = deal(length(opt.xmin:opt.xstep:opt.xmax),length(opt.ymin:opt.ystep:opt.ymax));

%Plot means and their changes from the original mean:
if opt.fit_horizons
    %Plot the horizons.
    xgrid_lims = [min(info.xgrid,[],'all'),max(info.xgrid,[],'all')];
    ygrid_lims = [min(info.ygrid,[],'all'),max(info.ygrid,[],'all')];
    for i = 1:opt.nhorizons
        %Final:
        zinterp_final = zeros(size(info.xgrid,1),N);
        zinterp_final_moved = zeros(size(info.xgrid,1),N);
        parfor j = 1:N
%         for j = 1:N
            if strcmp(opt.method,'forward')
                x = info.xgrid;
                y = info.ygrid;
                z = horiz_params_final(1+(i-1)*info.n_horiz_pts:i*info.n_horiz_pts,j);
                zinterp_final(:,j) = z;
            else
                x = zeros(sum(info.n_horiz_pts_block),1);
                y = zeros(sum(info.n_horiz_pts_block),1);
                FaultBlock = zeros(sum(info.n_horiz_pts_block),1);
                for k = 1:opt.n_fault_blocks
                    x(sum(info.n_horiz_pts_block(1:k-1))+1:sum(info.n_horiz_pts_block(1:k))) = info.xgrid(info.grid_pts_mask{k});
                    y(sum(info.n_horiz_pts_block(1:k-1))+1:sum(info.n_horiz_pts_block(1:k))) = info.ygrid(info.grid_pts_mask{k});
                    FaultBlock(sum(info.n_horiz_pts_block(1:k-1))+1:sum(info.n_horiz_pts_block(1:k))) = k-1;
                end
                z = horiz_params_final(1+(i-1)*sum(info.n_horiz_pts_block):i*sum(info.n_horiz_pts_block),j); %z for all fault block grids.
                active = true(size(FaultBlock));
                [x,y,z] = FilterHWFW(x,y,z,active,FaultBlock,faults_final(j),opt.fault_blocks_relationships);
                interpolant = scatteredInterpolant(x,y,z);
                zinterp_final(:,j) = interpolant(info.xgrid,info.ygrid);
            end
            %Move the points
            if strcmp(opt.method,'forward')
                [x,y,z] = MovePts_forward(faults_final(j),x',y',z');
            else
                [x,y,z,fval] = MovePts_restore(faults_final(j),x',y',z');
                [x,y,z] = deal(x(fval<opt.sigma_misfit_restored),y(fval<opt.sigma_misfit_restored),z(fval<opt.sigma_misfit_restored));
            end
            interpolant2 = scatteredInterpolant(x',y',z');
            zinterp_final_moved(:,j) = interpolant2(info.xgrid,info.ygrid);
        end
        mean_z_final = mean(zinterp_final,2);
        mean_z_final = reshape(mean_z_final,[nxgrid,nygrid]);
        mean_z_final_moved = mean(zinterp_final_moved,2);
        mean_z_final_moved = reshape(mean_z_final_moved,[nxgrid,nygrid]);
        %Initial:
        zinterp_init = zeros(size(info.xgrid,1),N);
        zinterp_init_moved = zeros(size(info.xgrid,1),N);
        for j = 1:N
            if strcmp(opt.method,'forward')
                x = info.xgrid;
                y = info.ygrid;
                z = horiz_params_init(1+(i-1)*info.n_horiz_pts:i*info.n_horiz_pts,j);
            else
                x = zeros(sum(info.n_horiz_pts_block),1);
                y = zeros(sum(info.n_horiz_pts_block),1);
                FaultBlock = zeros(sum(info.n_horiz_pts_block),1);
                for k = 1:opt.n_fault_blocks
                    x(sum(info.n_horiz_pts_block(1:k-1))+1:sum(info.n_horiz_pts_block(1:k))) = info.xgrid(info.grid_pts_mask{k});
                    y(sum(info.n_horiz_pts_block(1:k-1))+1:sum(info.n_horiz_pts_block(1:k))) = info.ygrid(info.grid_pts_mask{k});
                    FaultBlock(sum(info.n_horiz_pts_block(1:k-1))+1:sum(info.n_horiz_pts_block(1:k))) = k-1;
                end
                z = horiz_params_init(1+(i-1)*sum(info.n_horiz_pts_block):i*sum(info.n_horiz_pts_block),j); %z for all fault block grids.
                active = true(size(FaultBlock));
                [x,y,z] = FilterHWFW(x,y,z,active,FaultBlock,faults_init(j),opt.fault_blocks_relationships);
            end
            interpolant = scatteredInterpolant(x,y,z);
            zinterp_init(:,j) = interpolant(info.xgrid,info.ygrid);
            %Move the points
            if strcmp(opt.method,'forward')
                [x,y,z] = MovePts_forward(faults_init(j),x',y',z');
            else
                [x,y,z,fval] = MovePts_restore(faults_init(j),x',y',z');
                [x,y,z] = deal(x(fval<opt.sigma_misfit_restored),y(fval<opt.sigma_misfit_restored),z(fval<opt.sigma_misfit_restored));
            end
            interpolant2 = scatteredInterpolant(x',y',z');
            zinterp_init_moved(:,j) = interpolant2(info.xgrid,info.ygrid);
        end
        mean_z_init = mean(zinterp_init,2);
        mean_z_init = reshape(mean_z_init,[nxgrid,nygrid]);
        mean_z_init_moved = mean(zinterp_init_moved,2);
        mean_z_init_moved = reshape(mean_z_init_moved,[nxgrid,nygrid]);
        %True
        if strcmp(opt.method,'restoration')
            x = zeros(sum(info.n_horiz_pts_block),1);
            y = zeros(sum(info.n_horiz_pts_block),1);
            FaultBlock = zeros(sum(info.n_horiz_pts_block),1);
            for k = 1:opt.n_fault_blocks
                x(sum(info.n_horiz_pts_block(1:k-1))+1:sum(info.n_horiz_pts_block(1:k))) = info.xgrid(info.grid_pts_mask{k});
                y(sum(info.n_horiz_pts_block(1:k-1))+1:sum(info.n_horiz_pts_block(1:k))) = info.ygrid(info.grid_pts_mask{k});
                FaultBlock(sum(info.n_horiz_pts_block(1:k-1))+1:sum(info.n_horiz_pts_block(1:k))) = k-1;
            end
            z =  info.horiz_pts_ref(1+(i-1)*sum(info.n_horiz_pts_block):i*sum(info.n_horiz_pts_block)); %z for all fault block grids.
            active = true(size(FaultBlock));
            [x,y,z] = FilterHWFW(x,y,z,active,FaultBlock,info.faults_ref,opt.fault_blocks_relationships);
            interpolant = scatteredInterpolant(x,y,z);
            true_z = interpolant(info.xgrid,info.ygrid);
            true_z = reshape(true_z,[nxgrid,nygrid]);
        end
        figure(1)
        subplot(4,opt.nhorizons,i)
        image(xgrid_lims,ygrid_lims,mean_z_init,'CDataMapping','scaled')
        set(gca,'YDir','normal')
        colorbar
        title(['Initial Mean Depth ',opt.horizon_names{i}])
        %     caxis([min(mean_z_init,[],'all'),max(mean_z_init,[],'all')])
        subplot(4,opt.nhorizons,i+opt.nhorizons)
        image(xgrid_lims,ygrid_lims,mean_z_final,'CDataMapping','scaled')
        set(gca,'YDir','normal')
        colorbar
        title(['Final Mean Depth ',opt.horizon_names{i}])
        caxis([min(mean_z_final,[],'all'),max(mean_z_final,[],'all')])
        subplot(4,opt.nhorizons,i+2*opt.nhorizons)
        image(xgrid_lims,ygrid_lims,mean_z_final-mean_z_init,'CDataMapping','scaled')
        set(gca,'YDir','normal')
        colorbar
        title(['Change in Mean Depth',opt.horizon_names{i}])
        if strcmp(opt.method,'restoration')
            subplot(4,opt.nhorizons,i+3*opt.nhorizons)
            image(xgrid_lims,ygrid_lims,true_z,'CDataMapping','scaled')
            set(gca,'YDir','normal')
            colorbar
            title(['True Depth',opt.horizon_names{i}])
            caxis([min(mean_z_final,[],'all'),max(mean_z_final,[],'all')])
        end
        figure(2)
        subplot(4,opt.nhorizons,i)
        image(xgrid_lims,ygrid_lims,mean_z_init,'CDataMapping','scaled')
        set(gca,'YDir','normal')
        colorbar
        title(['Initial Mean Depth Moved',opt.horizon_names{i}])
        %     caxis([min(mean_z_init,[],'all'),max(mean_z_init,[],'all')])
        subplot(4,opt.nhorizons,i+opt.nhorizons)
        image(xgrid_lims,ygrid_lims,mean_z_final_moved,'CDataMapping','scaled')
        set(gca,'YDir','normal')
        colorbar
        title(['Final Mean Depth Moved ',opt.horizon_names{i}])
        caxis([min(mean_z_final_moved,[],'all'),max(mean_z_final_moved,[],'all')])
        subplot(4,opt.nhorizons,i+2*opt.nhorizons)
        image(xgrid_lims,ygrid_lims,mean_z_final_moved-mean_z_init_moved,'CDataMapping','scaled')
        set(gca,'YDir','normal')
        colorbar
        title(['Change in Mean Depth Moved',opt.horizon_names{i}])
    end
    
    %Horizons Standard deviation:
    xgrid_lims = [min(info.xgrid,[],'all'),max(info.xgrid,[],'all')];
    ygrid_lims = [min(info.ygrid,[],'all'),max(info.ygrid,[],'all')];
    for i = 1:opt.nhorizons
        %Final:
        std_z_final = std(zinterp_final,[],2);
        std_z_final = reshape(std_z_final,[nxgrid,nygrid]);
        %Initial:
        std_z_init = std(zinterp_init,[],2);
        std_z_init = reshape(std_z_init,[nxgrid,nygrid]);
        %Here it looks like I don't need to make mean_z be mean_z'. Must be
        %something different about how these grids are organized vs. the ones
        %for disp and fdata.
        figure(3)
        subplot(3,opt.nhorizons,i)
        image(xgrid_lims,ygrid_lims,std_z_init,'CDataMapping','scaled')
        set(gca,'YDir','normal')
        colorbar
        title(['Initial St.Dev. Depth ',opt.horizon_names{i}])
        subplot(3,opt.nhorizons,i+opt.nhorizons)
        image(xgrid_lims,ygrid_lims,std_z_final,'CDataMapping','scaled')
        set(gca,'YDir','normal')
        colorbar
        title(['Final St.Dev. Depth ',opt.horizon_names{i}])
        subplot(3,opt.nhorizons,i+2*opt.nhorizons)
        image(xgrid_lims,ygrid_lims,std_z_final-std_z_init,'CDataMapping','scaled')
        set(gca,'YDir','normal')
        colorbar
        title(['Change in St.Dev. Depth',opt.horizon_names{i}])
    end
end

%Contour the fault tip line.
figure(4)
c = contour(info.faults_ref(1).u,info.faults_ref(1).v,info.faults_ref(1).ThicknessAttribute,[0,0]);
trueTipContour_u = c(1,2:end);
trueTipContour_v = c(2,2:end);

if opt.fit_surf
    %Fault Surface:
    mean_fdata_final = mean(cat(3,faults_final.f),3);
    mean_fdata_init = mean(cat(3,faults_init.f),3);
    u_fdata_lims = [min(info.faults_ref(1).u,[],'all'),max(info.faults_ref(1).u,[],'all')];
    v_fdata_lims = [min(info.faults_ref(1).v,[],'all'),max(info.faults_ref(1).v,[],'all')];
    figure(4)
    subplot(1,4,1)
    image(u_fdata_lims,v_fdata_lims,mean_fdata_init','CDataMapping','scaled')
    hold on
    plot(trueTipContour_u,trueTipContour_v,'k')
    set(gca,'YDir','normal')
    colorbar
    title('Initial Mean Surface')
    caxis([-150,150])
    subplot(1,4,2)
    image(u_fdata_lims,v_fdata_lims,mean_fdata_final','CDataMapping','scaled')
    hold on
    plot(trueTipContour_u,trueTipContour_v,'k')
    set(gca,'YDir','normal')
    colorbar
    title('Final Mean Surface')
    caxis([-150,150])
    subplot(1,4,3)
    image(u_fdata_lims,v_fdata_lims,(mean_fdata_final-mean_fdata_init)','CDataMapping','scaled')
    hold on
    plot(trueTipContour_u,trueTipContour_v,'k')
    set(gca,'YDir','normal')
    colorbar
    title('Change in Mean Surface ')
    subplot(1,4,4)
    true_fdata = info.faults_ref(1).f; %Not really the true one but doesn't have the interpolation errors in fault_true.
    image(u_fdata_lims,v_fdata_lims,true_fdata','CDataMapping','scaled')
    hold on
    plot(trueTipContour_u,trueTipContour_v,'k')
    set(gca,'YDir','normal')
    colorbar
    title('True Model Surface')
%     title('Reference Model Surface')
    caxis([-150,150])

    %Fault Surface Standard Deviation:
    std_fdata_final = std(cat(3,faults_final.f),[],3);
    std_fdata_init = std(cat(3,faults_init.f),[],3);
    u_fdata_lims = [min(info.faults_ref(1).u,[],'all'),max(info.faults_ref(1).u,[],'all')];
    v_fdata_lims = [min(info.faults_ref(1).v,[],'all'),max(info.faults_ref(1).v,[],'all')];
    figure(5)
    subplot(1,3,1)
    image(u_fdata_lims,v_fdata_lims,std_fdata_init','CDataMapping','scaled')
    hold on
    plot(trueTipContour_u,trueTipContour_v,'k')
    set(gca,'YDir','normal')
    colorbar
    title('Initial StDev Surface')
    caxis([0,50])
    subplot(1,3,2)
    image(u_fdata_lims,v_fdata_lims,std_fdata_final','CDataMapping','scaled')
    hold on
    plot(trueTipContour_u,trueTipContour_v,'k')
    set(gca,'YDir','normal')
    colorbar
    title('Final StDev Surface')
    caxis([0,50])
    subplot(1,3,3)
    image(u_fdata_lims,v_fdata_lims,(std_fdata_final-std_fdata_init)','CDataMapping','scaled')
    hold on
    plot(trueTipContour_u,trueTipContour_v,'k')
    set(gca,'YDir','normal')
    colorbar
    title('Change in StDev Surface ')
end

%Displacement mean:
if opt.fit_disp
    mean_disp_final = mean(cat(3,faults_final.displacement),3);
    mean_disp_init = mean(cat(3,faults_init.displacement),3);
    u_disp_lims = [min(info.faults_ref(1).u,[],'all'),max(info.faults_ref(1).u,[],'all')];
    v_disp_lims = [min(info.faults_ref(1).v,[],'all'),max(info.faults_ref(1).v,[],'all')];
    figure(6)
    subplot(1,4,1)
    image(u_disp_lims,v_disp_lims,mean_disp_init','CDataMapping','scaled')
    hold on
    plot(trueTipContour_u,trueTipContour_v,'k')
    set(gca,'YDir','normal')
    colorbar
    title('Initial Mean Displacement')
    caxis([0,300])
%     caxis([0,500])
    subplot(1,4,2)
    image(u_disp_lims,v_disp_lims,mean_disp_final','CDataMapping','scaled')
    hold on
    plot(trueTipContour_u,trueTipContour_v,'k')
    set(gca,'YDir','normal')
    colorbar
    title('Final Mean Displacement')
    caxis([0,300])
%     caxis([0,500])
    subplot(1,4,3)
    image(u_disp_lims,v_disp_lims,(mean_disp_final-mean_disp_init)','CDataMapping','scaled')
    hold on
    plot(trueTipContour_u,trueTipContour_v,'k')
    set(gca,'YDir','normal')
    colorbar
    title('Change in Mean Displacement')
    subplot(1,4,4)
    true_disp = info.faults_ref(1).displacement;
    image(u_disp_lims,v_disp_lims,true_disp','CDataMapping','scaled')
    hold on
    plot(trueTipContour_u,trueTipContour_v,'k')
    set(gca,'YDir','normal')
    colorbar
    title('True Model Displacement')
    caxis([0,300])
%     caxis([0,500])
end

%Displacement Standard Deviation:
if opt.fit_disp
    std_disp_final = std(cat(3,faults_final.displacement),[],3);
    std_disp_init = std(cat(3,faults_init.displacement),[],3);
    u_disp_lims = [min(info.faults_ref(1).u,[],'all'),max(info.faults_ref(1).u,[],'all')];
    v_disp_lims = [min(info.faults_ref(1).v,[],'all'),max(info.faults_ref(1).v,[],'all')];
    figure(7)
    subplot(1,3,1)
    image(u_disp_lims,v_disp_lims,std_disp_init','CDataMapping','scaled')
    hold on
    plot(trueTipContour_u,trueTipContour_v,'k')
    set(gca,'YDir','normal')
    colorbar
    title('Initial StDev Displacement')
    caxis([0,100])
    subplot(1,3,2)
    image(u_disp_lims,v_disp_lims,std_disp_final','CDataMapping','scaled')
    hold on
    plot(trueTipContour_u,trueTipContour_v,'k')
    set(gca,'YDir','normal')
    colorbar
    title('Final StDev Displacement')
    caxis([0,100])
    subplot(1,3,3)
    image(u_disp_lims,v_disp_lims,(std_disp_final-std_disp_init)','CDataMapping','scaled')
    hold on
    plot(trueTipContour_u,trueTipContour_v,'k')
    set(gca,'YDir','normal')
    colorbar
    title('Change in StDev Displacement')
end

%Check how the dummy parameter has changed.
dummy_param_final = params_final_raw(end,:); %The dummy parameter isn't transformed, so raw is fine.
dummy_param_init = params_initial_raw(end,:);
figure(8)
histogram(dummy_param_init,-3:0.5:3);
hold on
histogram(dummy_param_final,-3:0.5:3);
hold off
legend('Initial','Final')
title('Dummy Parameter')

%Plot histograms of asymmetry and range and see how they changed.
if opt.fit_disp
    figure(9)
    subplot(1,2,1)
    ind_asym = find(strcmp(info.gen_param_names,'Fault1 Asymmetry'));
    ind_range = find(strcmp(info.gen_param_names,'Fault1 Range'));
    histogram(gen_params_initial(ind_asym,:),info.gen_param_mins(ind_asym):0.05:info.gen_param_maxs(ind_asym));
    hold on
    histogram(gen_params_final(ind_asym,:),info.gen_param_mins(ind_asym):0.05:info.gen_param_maxs(ind_asym));
    hold off
    legend('Initial','Final')
    title('Asymmetry')
    subplot(1,2,2)
    histogram(gen_params_initial(ind_range,:),info.gen_param_mins(ind_range):50:info.gen_param_maxs(ind_range));
    hold on
    histogram(gen_params_final(ind_range,:),info.gen_param_mins(ind_range):50:info.gen_param_maxs(ind_range));
    hold off
    legend('Initial','Final')
    title('Range')
end

%Plot histograms of the fault maximum displacement.
figure(10)
ind_max_disp = find(strcmp(info.gen_param_names,'Fault1 Maximum Displacement'));
histogram(gen_params_initial(ind_max_disp,:),info.gen_param_mins(ind_max_disp):25:info.gen_param_maxs(ind_max_disp))
hold on
histogram(gen_params_final(ind_max_disp,:),info.gen_param_mins(ind_max_disp):25:info.gen_param_maxs(ind_max_disp))
hold off
legend('Initial','Final')
title('Fault1 Maximum Displacement')

%Get strike and dip from fault M matrices. (Only relevant for
%_fixed_strike_dip.
dip = zeros(1,opt.N);
strike = zeros(1,opt.N);
for i = 1:opt.N
dip(i) = asind(-faults_final(i).M(3,2));
strike(i) = atan2d(faults_final(i).M(1,1),faults_final(i).M(2,1));
end
disp(['Strike = ',num2str(mean(strike)),' +/- ',num2str(std(strike))]);
disp(['Dip = ',num2str(mean(dip)),' +/- ',num2str(std(dip))]);