%Plot the mean parameter fields.

%Note: This is set up on the assumption of the forward modelling method.
%Some changes would need to be made for restoration.

%For now I'm just plotting one horizon, but that can be easily changed.

fontsize = 16;
cmap = 'jet';

% folder = 'DMCBootInflSigma06_Results';
% N = 200;
% load(['./',folder,'/N',num2str(N),'_forward.mat'],'params_final_raw','info','opt')
load('N200_Forward_scattered.mat');

%Load the reference model horizons.
load('ReferenceModel_Deformed_Test6.mat')
horiz_ref_deformed = horiz_ref;
load('ReferenceModel_Flat_Test6.mat')
horiz_ref_restored = horiz_ref;

%Get the grid dimensions.
[nxgrid,nygrid] = deal(length(opt.xmin:opt.xstep:opt.xmax),length(opt.ymin:opt.ystep:opt.ymax));

%Transform the final parameters back to the bounded domain.
N = size(params_final_raw,2); %Ensemble size
gen_params_final = zeros(info.nparams.general,N);
faults_final = fault_class.empty();
horiz_params_final = zeros(info.nparams.horiz,N);
parfor i = 1:N
% for i = 1:N
    [gen_params_final(:,i),faults_final(i),horiz_params_final(:,i)] = InterpretModel(params_final_raw(:,i),opt,info);
end

figure(1)
%Top horizon restored and deformed:
z_restored = zeros(size(info.xgrid,1),N);
zinterp_deformed = zeros(size(info.xgrid,1),N);
parfor j = 1:opt.N
    x = info.xgrid;
    y = info.ygrid;
    z = horiz_params_final(1:info.n_horiz_pts,j);
    z_restored(:,j) = z;
    [x,y,z] = MovePts_forward(faults_final(j),x',y',z');
    interpolant2 = scatteredInterpolant(x',y',z');
    zinterp_deformed(:,j) = interpolant2(info.xgrid,info.ygrid);
end
mean_z_restored = mean(z_restored,2);
mean_z_restored = reshape(mean_z_restored,[nxgrid,nygrid]);
std_z_restored = std(z_restored,[],2);
std_z_restored = reshape(std_z_restored,[nxgrid,nygrid]);
mean_z_deformed = mean(zinterp_deformed,2);
mean_z_deformed = reshape(mean_z_deformed,[nxgrid,nygrid]);
std_z_deformed = std(zinterp_deformed,[],2);
std_z_deformed = reshape(std_z_deformed,[nxgrid,nygrid]);
% min_depth = min(min(mean_z_restored,[],'all'),min(mean_z_deformed,[],'all'));
% max_depth = max(max(mean_z_restored,[],'all'),max(mean_z_deformed,[],'all'));
% min_depth_std = min(min(std_z_restored,[],'all'),min(std_z_deformed,[],'all'));
% max_depth_std = max(max(std_z_restored,[],'all'),max(std_z_deformed,[],'all'));
min_depth = 1070;
max_depth = 1190;
min_depth_std = 0;
max_depth_std = 70;
%True
true_z_restored = horiz_ref_restored.z{1}';
true_z_restored = reshape(true_z_restored,[nxgrid,nygrid]);
x = horiz_ref_deformed.x{1};
y = horiz_ref_deformed.y{1};
z = horiz_ref_deformed.z{1};
FaultBlock = horiz_ref_deformed.FaultBlock{1};
block = AssignFaultBlock(x,y,z,info.faults_ref,opt.fault_blocks_relationships);
active = FaultBlock==block;
[x,y,z] = deal(x(active),y(active),z(active));
interpolant = scatteredInterpolant(x,y,z);
true_z_deformed = interpolant(info.xgrid,info.ygrid);
true_z_deformed = reshape(true_z_deformed,[nxgrid,nygrid]);
xgrid_lims = [min(info.xgrid,[],'all'),max(info.xgrid,[],'all')];
ygrid_lims = [min(info.ygrid,[],'all'),max(info.ygrid,[],'all')];
subplot(2,3,1)
image(xgrid_lims,ygrid_lims,true_z_deformed,'CDataMapping','scaled')
set(gca,'YDir','normal')
colormap(gca,cmap)
colorbar
% title(['True Deformed Depth',opt.horizon_names{1}])
axis equal
caxis([min_depth,max_depth])
xlabel('Distance East (m)')
ylabel('Distance North (m)')
set(gca,'FontSize',fontsize);
subplot(2,3,2)
image(xgrid_lims,ygrid_lims,mean_z_deformed,'CDataMapping','scaled')
set(gca,'YDir','normal')
colormap(gca,cmap)
colorbar
% title(['Ensemble Mean Deformed Depth ',opt.horizon_names{1}])
axis equal
caxis([min_depth,max_depth])
xlabel('Distance East (m)')
ylabel('Distance North (m)')
set(gca,'FontSize',fontsize);
subplot(2,3,3)
image(xgrid_lims,ygrid_lims,std_z_deformed,'CDataMapping','scaled')
set(gca,'YDir','normal')
colormap(gca,cmap)
colorbar
% title(['St.Dev. Deformed Depth ',opt.horizon_names{1}])
axis equal
caxis([min_depth_std,max_depth_std])
xlabel('Distance East (m)')
ylabel('Distance North (m)')
set(gca,'FontSize',fontsize);

%Plot the fault surface
%Contour the fault tip line.
figure(3)
c = contour(info.faults_ref(1).u,info.faults_ref(1).v,info.faults_ref(1).ThicknessAttribute,[0,0]);
trueTipContour_u = c(1,2:end);
trueTipContour_v = c(2,2:end);
figure(1)

%Fault Surface:
mean_fdata_final = mean(cat(3,faults_final.f),3);
std_fdata_final = std(cat(3,faults_final.f),[],3);
u_fdata_lims = [min(info.faults_ref(1).u,[],'all'),max(info.faults_ref(1).u,[],'all')];
v_fdata_lims = [min(info.faults_ref(1).v,[],'all'),max(info.faults_ref(1).v,[],'all')];
subplot(2,3,4)
true_fdata = info.faults_ref(1).f; %Not really the true one but doesn't have the interpolation errors in fault_true.
% true_fdata(info.faults_ref(1).ThicknessAttribute<0) = nan;
image(u_fdata_lims,v_fdata_lims,true_fdata','CDataMapping','scaled')
hold on
plot(trueTipContour_u,trueTipContour_v,'k')
set(gca,'YDir','normal')
colormap(gca,cmap)
colorbar
% title('True Fault Surface')
caxis([-150,150])
axis equal
xlabel('Along Strike Distance (m)')
ylabel('Along Dip Distance (m)')
set(gca,'FontSize',fontsize);
hold off
subplot(2,3,5)
image(u_fdata_lims,v_fdata_lims,mean_fdata_final','CDataMapping','scaled')
hold on
plot(trueTipContour_u,trueTipContour_v,'k')
set(gca,'YDir','normal')
colormap(gca,cmap)
colorbar
% title('Ensemble Mean Fault Surface')
caxis([-150,150])
axis equal
xlabel('Along Strike Distance (m)')
ylabel('Along Dip Distance (m)')
set(gca,'FontSize',fontsize);
hold off
subplot(2,3,6)
image(u_fdata_lims,v_fdata_lims,std_fdata_final','CDataMapping','scaled')
hold on
plot(trueTipContour_u,trueTipContour_v,'k')
set(gca,'YDir','normal')
colormap(gca,cmap)
colorbar
% title('St. Dev. Fault Surface')
axis equal
caxis([0,260])
xlabel('Along Strike Distance (m)')
ylabel('Along Dip Distance (m)')
set(gca,'FontSize',fontsize);
hold off

