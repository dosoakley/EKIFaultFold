%Make the forward model and save it.

%The EnKF_Combined_Code8 folder should be added to the Matlab path to run
%this correctly.

rng default; %Set the random number seed here so that results are repeatable.

%Create the undeformed flat horizons.
depths = [1100,1400]; %The depths of horizons A and B.
x1 = -1000:100:1000;
y1 = -1000:100:1000;
[xx1, yy1] = meshgrid(x1, y1);
zz1 = depths(1)*ones(size(xx1));
x2 = x1;
y2 = y1;
[xx2, yy2] = meshgrid(x2, y2);
zz2 = depths(2)*ones(size(xx2));
horiz_pts1 = [xx1(:),yy1(:),zz1(:)]';
horiz_pts2 = [xx2(:),yy2(:),zz2(:)]';
horizons_flat = {horiz_pts1,horiz_pts2};

%Define the elliptical fault. (These are the same parameters that are in
%the ell_fault.txt file for Havana in Test4.)
name = 'Fault1';
[cx,cy,cz] = deal(0.0,0.0,1200);
[strike,dip] = deal(70.0,55.0);
max_displacement = 250;
asymmetry = 0.6;
fault_length = 2000.0; %These are full length and height of the fault, so lu and lv will be half this.
fault_height = 1000.0;
range = 700.0; %Radius over which slip dies out away from the fault.
slip_sense = 1; %1 for normal; 0 for reverse.

%Create the elliptical fault.
fault = fault_class();
fault.name = name;
MakeM(fault,cx,cy,cz,strike*pi/180,dip*pi/180);
fault.NU = 26;
fault.NV = 16;
fault.dU = 100;
fault.dV = 100;
fault.LU = (fault.NU-1)*fault.dU;
fault.LV = (fault.NV-1)*fault.dV;
fault.U0 = -fault.LU/2;
fault.V0 = -fault.LV/2;
fault.N = fault.NU*fault.NV;
[fault.u,fault.v] = ndgrid(fault.U0:fault.dU:(fault.U0+fault.LU),fault.V0:fault.dV:(fault.V0+fault.LV));
fault.f = zeros(fault.NU,fault.NV);
fault.lu = fault_length/2;
fault.lv = fault_height/2;
fault.max_val = max_displacement;
fault.u0 = 0; %Center is at (u=0,v=0).
fault.v0 = 0;
fault.slip_type = slip_sense;
fault.displacement = EllipticalTrend(fault,fault.u,fault.v);
CalculateThicknessAttribute(fault);
fault.asymmetry = asymmetry;
fault.Range = range;

% % %Modify the fault geometry into a kind of s shape, or alternatively, a
% c shape. (Currenlty using c shape.)
% % fault.f = fault.f+200*sign(fault.u).*(fault.u./max(fault.u)).^2;
fault.f = fault.f+250*(fault.u./max(fault.u)).^2;
fault.f = fault.f-mean(fault.f); %Recenter it around 0.

% %Modify the fault geometry with a spherical variogram.
% surf_std = 200; %Variogram standard deviation
% surf_Lc = 3000; %Variogram range
% r = CorrSpher(fault.u(:),fault.v(:),surf_Lc);
% [L,p] = chol(r,'lower');
% if p > 0
%     error ('Cf not positive definite');
% end
% perturbations = surf_std*L*randn(size(fault.u(:)));
% mask = fault.ThicknessAttribute>0;
% fault.f(mask) = fault.f(mask)+perturbations(mask);
% %This doesn't actually work as well as I'd like.

% %Modify the displacement field by perturbation with a spherical variogram.
disp_std = 20; %Variogram standard deviation
disp_Lc = 600; %Variogram range
r = CorrSpher(fault.u(:),fault.v(:),disp_Lc);
[L,p] = chol(r,'lower');
if p > 0
    error ('Cf not positive definite');
end
perturbations = disp_std*L*randn(size(fault.u(:)));
mask = fault.ThicknessAttribute>0;
fault.displacement(mask) = fault.displacement(mask)+perturbations(mask);

%Forward model the deformation.
horizons = horizons_flat;
for i = 1:2
    [horizons{i}(1,:),horizons{i}(2,:),horizons{i}(3,:)] = ...
        MovePts_forward(fault,horizons{i}(1,:),...
        horizons{i}(2,:),horizons{i}(3,:));
end

%Find the fault lines for each horizon. These will be used to identify
%regions where the horizon is absent and should not be included in well
%data.
fault_blocks_relationships = [1;-1];
n_fault_blocks = 2;
fault_lines = {[],[],[],[]}; %This holds hw and fw separately.
fault_lines_combined = {[],[]}; %This combines hw and fw fault lines into a single polygon in 2D space.
for i = 1:length(horizons)
    bed_uv = fault.M\[horizons{i}(1,:);horizons{i}(2,:);horizons{i}(3,:);ones(size(horizons{i}(1,:)))];
    blocks_pts = AssignFaultBlock(horizons{i}(1,:)',horizons{i}(2,:)',horizons{i}(3,:)',fault,fault_blocks_relationships);
%     f = InterpolateSurface(fault,bed_uv(1,:),bed_uv(2,:),false);
%     dist = bed_uv(3,:)-f; %Distance from the fault.
    %t = InterpolateThicknessAttribute(fault,bed_uv(1,:),bed_uv(2,:));
    for k = 1:n_fault_blocks
        mask = blocks_pts==k; %& t'>=0;
        zI = scatteredInterpolant(horizons{i}(1,mask)',horizons{i}(2,mask)',horizons{i}(3,mask)','linear','nearest');
        [xx1, yy1] = meshgrid(x1, y1);
        zz1 = zI(xx1,yy1);
        grid_uv = fault.M\[xx1(:),yy1(:),zz1(:),ones(size(zz1(:)))]';
        ff1 = InterpolateSurface(fault,grid_uv(1,:),grid_uv(2,:)); %I'm using true here b/c AssignFaultBlock above uses true.
        dist = grid_uv(3,:)-ff1; %Distance from the fault.
        dd1 = reshape(dist,size(xx1));
        C = contour(xx1,yy1,dd1,[0,0]);
        contour_line = C(:,2:end);
        if C(2,1)<size(contour_line,2) %If there is more than one 0 contour line, this can happen.
            contour_line = contour_line(:,1:C(2,1));
            warning('Multiple 0 contour lines for fault line. Using the first one only.')
        end
        [x_line,y_line] = deal(contour_line(1,:),contour_line(2,:));
%         zI = scatteredInterpolant(horizons{i}(1,mask)',horizons{i}(2,mask)',horizons{i}(3,mask)');
        z_line = zI(x_line,y_line);
        line_uv = fault.M\[x_line;y_line;z_line;ones(size(z_line))];
        t_line = InterpolateThicknessAttribute(fault,line_uv(1,:),line_uv(2,:));
        ind = 2*(i-1)+k;
        if k==1
            fault_lines{ind} = [x_line(t_line>=0);y_line(t_line>=0);z_line(t_line>=0)];
        else
            d1 = sqrt((fault_lines{ind-1}(1,1)-x_line(1))^2+(fault_lines{ind-1}(2,1)-y_line(1))^2);
            d2 = sqrt((fault_lines{ind-1}(1,end)-x_line(1))^2+(fault_lines{ind-1}(2,end)-y_line(1))^2);
            if d1>d2 %New start point is closer to last end point.
                fault_lines{ind} = [x_line(t_line>=0);y_line(t_line>=0);z_line(t_line>=0)];
            else %New start point is closer to last start point.
                fault_lines{ind} = fliplr([x_line(t_line>=0);y_line(t_line>=0);z_line(t_line>=0)]);
            end
        end
    end
    fault_lines_combined{i} = [fault_lines{2*(i-1)+1},fault_lines{2*i},fault_lines{2*(i-1)+1}(:,1)]; %Add the first point to the end to make it a complete polygon.
end
save('fault_lines.mat','fault_lines','fault_lines_combined');

%Filter horizon points near the fault.
filter_distance = 25; %50;
for i = 1:2
    horiz_pts_uv = fault.M\[horizons{i};ones(1,size(horizons{i},2))];
    f = InterpolateSurface(fault,horiz_pts_uv(1,:),horiz_pts_uv(2,:));
    t = InterpolateThicknessAttribute(fault,horiz_pts_uv(1,:),horiz_pts_uv(2,:));
    mask = t>0 & abs(horiz_pts_uv(3,:)-f)<=filter_distance;
    horizons{i} = horizons{i}(:,~mask);
%     disp(sum(mask))
end

%Save the fault.
save('fault.mat','fault');

%Save the horizons.
save('horizons.mat','horizons');

%Save the horizon data as text files, with added noise.
sigma_horiz = 10; %Noise level.
horiz_letters = {'A','B'};
for i = 1:length(horiz_letters)
    outfile = fopen(['./DataDense/Horizon',horiz_letters{i},'_Data.txt'],'w');
    x = horizons{i}(1,:);
    y = horizons{i}(2,:);
    z = horizons{i}(3,:);
    z = z+sigma_horiz*randn(size(z)); %Add error.
    npts = length(x);
    for j = 1:npts
        fprintf(outfile,'%f %f %f\n',x(j),y(j),z(j));
    end
    fclose(outfile);
end

%Interpolate the fault data onto a grid and save those points as the fault
%data points.
%This is the same grid used for the beds data points with added noise.
sigma_fault = 10; %Noise level.
mask = fault.ThicknessAttribute>0;
fault_xyz = fault.M*[fault.u(mask),fault.v(mask),fault.f(mask),ones(sum(mask,'all'),1)]';
FInterpolant = scatteredInterpolant(fault_xyz(1,:)',fault_xyz(2,:)',fault_xyz(3,:)','linear','none');
x1 = -1000:100:1000;
y1 = -1000:100:1000;
[xx1, yy1] = meshgrid(x1, y1);
[xx1,yy1] = deal(xx1(:),yy1(:));
Fdata = FInterpolant(xx1,yy1);
%Add noise and remove NaN values:
mask2 = ~isnan(Fdata);
pts_uvw = fault.M\[xx1(mask2),yy1(mask2),Fdata(mask2),ones(sum(mask2),1)]';
pts_uvw(3,:) = pts_uvw(3,:)+sigma_fault*randn(size(pts_uvw(3,:))); %Add error in the fault-normal direction.
pts_xyz = fault.M*pts_uvw;
%Save
file = fopen('./DataDense/Fault1_Data.txt','w');
% for i = 1:numel(Fdata)
%     if ~isnan(Fdata(i))
%         fprintf(file,'%f %f %f\n',xx1(i),yy1(i),Fdata(i));
%     end
% end
for i = 1:size(pts_xyz,2)
    fprintf(file,'%f %f %f\n',pts_xyz(1,i),pts_xyz(2,i),pts_xyz(3,i));
end
fclose(file);

%Extract scattered data at a series of randomly chosen well points.
num_wells = 50;
% Ewells = -1000+2000*rand(num_wells,1); %Uniform distribution.
% Nwells = -1000+2000*rand(num_wells,1);
sigma_well_dist = 1000; %600; %Well locations will be drawn from a normal distribution with this standard deviation, centered around 0.
Ewells = sigma_well_dist*randn(num_wells,1); %Normal distribution.
Nwells = sigma_well_dist*randn(num_wells,1);
while any(abs(Ewells)>1e3 | abs(Nwells)>1e3)
    mask = abs(Ewells)>1e3 | abs(Nwells)>1e3;
    Ewells(mask) = sigma_well_dist*randn(sum(mask),1);
    Nwells(mask) = sigma_well_dist*randn(sum(mask),1);
end
horizA = horizons{1};
horizB = horizons{2};
AInterpolant = scatteredInterpolant(horizA(1,:)',horizA(2,:)',horizA(3,:)','linear','none');
BInterpolant = scatteredInterpolant(horizB(1,:)',horizB(2,:)',horizB(3,:)','linear','none');
picksA = zeros(1,num_wells);
picksB = zeros(1,num_wells);
picksF = zeros(1,num_wells);
for i = 1:num_wells
    picksA(i) = AInterpolant(Ewells(i),Nwells(i));
    picksB(i) = BInterpolant(Ewells(i),Nwells(i));
    picksF(i) = FInterpolant(Ewells(i),Nwells(i));
end

%Check if the wells are in the area where a horizon is not present.
in = {[],[]};
for i = 1:length(horizons)
   in{i} = inpolygon(Ewells,Nwells,fault_lines_combined{i}(1,:),fault_lines_combined{i}(2,:));
end

%Save the well locations.
save('wells.mat','Ewells','Nwells');

%Save the well picks.
horiz_picks = [picksA;picksB];
for i = 1:length(horiz_letters)
    outfile = fopen(['./DataSparse/Horizon',horiz_letters{i},'_Data.txt'],'w');
    for j = 1:num_wells
        if ~in{i}(j)
            fprintf(outfile,'%f %f %f\n',Ewells(j),Nwells(j),horiz_picks(i,j));
        end
    end
    fclose(outfile);
end
outfile = fopen('./DataSparse/Fault1_Data.txt','w');
for j = 1:num_wells
    if ~isnan(picksF(j))
        fprintf(outfile,'%f %f %f\n',Ewells(j),Nwells(j),picksF(j));
    end
end
fclose(outfile);

%Plot the locations of the wells relative to the fault.
figure(1)
scatter(fault_xyz(1,:),fault_xyz(2,:),[],fault_xyz(3,:))
colorbar
hold on
mask = ~isnan(picksF);
scatter(Ewells(mask),Nwells(mask),[],picksF(mask));
plot(Ewells,Nwells,'k.')
hold off
axis([-1e3,1e3,-1e3,1e3])

%Plot the model.
figure(2)
mask = fault.ThicknessAttribute>0;
fault_xyz = fault.M*[fault.u(mask),fault.v(mask),fault.f(mask),ones(sum(mask,'all'),1)]';
scatter3(fault_xyz(1,:),fault_xyz(2,:),-fault_xyz(3,:))
hold on
scatter3(horizA(1,:),horizA(2,:),-horizA(3,:))
scatter3(horizB(1,:),horizB(2,:),-horizB(3,:))
hold off