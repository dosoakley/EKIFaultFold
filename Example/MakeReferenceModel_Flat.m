%Make a reference model.
%Load the reference fault model made from the ForwardModel folder.
%Create a grid of horizon points at constant restored-state elevation for 
%the horizons.

%Name of the file to save the reference model to:
save_file_name = 'ReferenceModel_Flat_Test6.mat';

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
flat_depths = [1100,1400];
nhorizons = length(horizon_names);
horiz_ref.x = cell(1,nhorizons);
horiz_ref.y = cell(1,nhorizons);
horiz_ref.z = cell(1,nhorizons);
horiz_ref.FaultBlock = cell(1,nhorizons);
horiz_ref.active = cell(1,nhorizons);
for i = 1:nhorizons
    horiz_ref.x{i} = xgrid;
    horiz_ref.y{i} = ygrid;
    horiz_ref.z{i} = flat_depths(i)*ones(n_horiz_pts,1);
    horiz_ref.FaultBlock{i} = zeros(n_horiz_pts,1);
    horiz_ref.active{i} = true(n_horiz_pts,1);
end

save(save_file_name,'faults_ref','horiz_ref')