%Make a reference model.
%Load the reference fault model made from the ForwardModel folder.
%Create a grid of horizon points at constant restored-state elevation for 
%the horizons.

%Name of the file to save the reference model to:
save_file_name = 'ReferenceModel_Deformed_Test6.mat';

%Load the reference fault.
load('./ForwardModel/fault.mat','fault');
faults_ref = copy(fault);

%Create the horizontal grid.
[xmin,xmax] = deal(-1e3,1e3);
[ymin,ymax] = deal(-1e3,1e3);
[xstep,ystep] = deal(100,100); %x and y step sizes for grid.
[xgrid, ygrid] = meshgrid(xmin:xstep:xmax, ymin:ystep:ymax);
[nxgrid,nygrid] = size(xgrid);
xgrid = xgrid(:);
ygrid = ygrid(:);
n_horiz_pts = length(xgrid); %The number of data points per horizon.

%Names of the horizons:
horizon_names = {'HorizonA','HorizonB'};

%Make the reference horizons model.
nhorizons = length(horizon_names);
horiz_ref.x = cell(1,nhorizons);
horiz_ref.y = cell(1,nhorizons);
horiz_ref.z = cell(1,nhorizons);
horiz_ref.FaultBlock = cell(1,nhorizons);
horiz_ref.active = cell(1,nhorizons);
filename = './ForwardModel/horizons.mat';
load(filename);
for i = 1:nhorizons
    horiz_ref.x{i} = horizons{i}(1,:)';
    horiz_ref.y{i} = horizons{i}(2,:)';
    horiz_ref.z{i} = horizons{i}(3,:)';
    data_xyz = [horiz_ref.x{i}';horiz_ref.y{i}';horiz_ref.z{i}';ones(size(horiz_ref.x{i}'))];
    data_uvw = fault.M\data_xyz;
    data_f = InterpolateSurface(fault,data_uvw(1,:),data_uvw(2,:));
    FaultBlocks = zeros(size(horiz_ref.x{i}));
    FaultBlocks(data_uvw(3,:)>=data_f) = 1; %HW
    FaultBlocks(data_uvw(3,:)<data_f) = 2; %FW
    horiz_ref.FaultBlock{i} = FaultBlocks;
    horiz_ref.active{i} = true(size(horiz_ref.x{i})); %Since we didn't create any inactive points, all points are active.
end

save(save_file_name,'faults_ref','horiz_ref')