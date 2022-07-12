function [fault] = ReadHavanaFault(filename,fault_name,asymmetry,Range)
%Read the .grid file produced by Havana for a fault and save it to a
%fault_class object.
%The displacement is reinterpolated onto the same grid as the surface and 
%tip line.
%Since asymmetry and Range are not included in the .grid file, they need to
%be assigned manually in the call to this function.

fault = fault_class();
fault.name = fault_name;
file = fopen(filename);
while ~feof(file)
   line = fgetl(file);
   if strcmp(line,'TransformMatrix4x4')
       fault.M = cell2mat(textscan(file,'%f%f%f%f',4));
   elseif strcmp(line,'# Fault Surface: ')
       fault_surface = true;
       displacement_field = false;
   elseif strcmp(line,'# Displacement Field: ')
       fault_surface = false;
       displacement_field = true;
   elseif length(line)>2 && strcmp(line(1:2),'NU')
       if fault_surface
           surf.NU = str2double(line(4:end));
       elseif displacement_field
           disp.NU = str2double(line(4:end));
       end
   elseif length(line)>2 && strcmp(line(1:2),'NV')
       if fault_surface
           surf.NV = str2double(line(4:end));
       elseif displacement_field
           disp.NV = str2double(line(4:end));
       end
   elseif length(line)>4 && strcmp(line(1:4),'Umin')
       surf.U0 = str2double(line(6:end));
   elseif length(line)>4 && strcmp(line(1:4),'Vmin')
       surf.V0 = str2double(line(6:end));
   elseif length(line)>2 && strcmp(line(1:2),'U0')
       disp.U0 = str2double(line(4:end));
   elseif length(line)>2 && strcmp(line(1:2),'V0')
       disp.V0 = str2double(line(4:end));
   elseif length(line)>2 && strcmp(line(1:2),'LU')
       surf.LU = str2double(line(4:end));
   elseif length(line)>2 && strcmp(line(1:2),'LV')
       surf.LV = str2double(line(4:end));
   elseif length(line)>8 && strcmp(line(1:8),'U-length')
       disp.LU = str2double(line(10:end));
   elseif length(line)>8 && strcmp(line(1:8),'V-length')
       disp.LV = str2double(line(10:end));
   elseif strcmp(line,'# Center point for elliptic trend ')
       fgetl(file); %Skip one line.
       temp = textscan(file,'%s%f',3);
       disp.ctr = temp{2};
   elseif strcmp(line,'#Trend Length ')
       disp.trend_length = str2double(fgetl(file));
   elseif strcmp(line,'#Trend height ')
       disp.trend_height = str2double(fgetl(file));
   elseif strcmp(line,'#Max trend ')
       disp.max_trend = str2double(fgetl(file));
   elseif strcmp(line,'Data')
       data = cell2mat(textscan(file,'%f',surf.NU*surf.NV));
   elseif strcmp(line,'AttributeThickness')
       ThicknessAttribute = cell2mat(textscan(file,'%f',surf.NU*surf.NV));
   elseif strcmp(line,'AttributeDisplacement')
       Displacement = cell2mat(textscan(file,'%f',disp.NU*disp.NV));
   end
end
fclose(file);

%Get the surface U and V values.
surf.dU = surf.LU/(surf.NU-1); %U increment
surf.dV = surf.LV/(surf.NV-1);
[surf.u,surf.v] = ndgrid(surf.U0:surf.dU:(surf.U0+surf.LU),...
    surf.V0:surf.dV:(surf.V0+surf.LV));
surf.data = reshape(data,[surf.NU,surf.NV]); %Put this in the same format as u and v.
surf.ThicknessAttribute = reshape(ThicknessAttribute,[surf.NU,surf.NV]); %Put this in the same format as u and v.

%Get the displacement U and V values.
disp.dU = disp.LU/(disp.NU-1); %U increment
disp.dV = disp.LV/(disp.NV-1);
[disp.u,disp.v] = ndgrid(disp.U0:disp.dU:(disp.U0+disp.LU),...
    disp.V0:disp.dV:(disp.V0+disp.LV));
disp.displacement = reshape(Displacement,[disp.NU,disp.NV]); %Put this in the same format as u and v.

%Reinterpolate displacement onto the same grid as surface.
disp_reduced_grid = interp2(disp.u',disp.v',disp.displacement',surf.u,surf.v,'linear',0);
[disp.NU,disp.NV] = deal(surf.NU,surf.NV);
[disp.U0,disp.V0] = deal(surf.U0,surf.V0);
[disp.LU,disp.LV] = deal(surf.LU,surf.LV);
[disp.dU,disp.dV] = deal(surf.dU,surf.dV);
[disp.u,disp.v] = deal(surf.u,surf.v);
disp.displacement = disp_reduced_grid;
    
%Add fields to the fault object.
[fault.NU,fault.NV] = deal(surf.NU,surf.NV);
[fault.U0,fault.V0] = deal(surf.U0,surf.V0);
[fault.LU,fault.LV] = deal(surf.LU,surf.LV);
[fault.dU,fault.dV] = deal(surf.dU,surf.dV);
[fault.u,fault.v] = deal(surf.u,surf.v);
fault.f = surf.data;
fault.ThicknessAttribute = surf.ThicknessAttribute;
fault.displacement = disp.displacement;
fault.displacement(surf.ThicknessAttribute<0) = 0;

%Record the total number of U,V points too.
fault.N = fault.NU*fault.NV;

%Determine the sense of slip.
%This follows the convention of the FaultSet.txt file that "Slip type"is 0
%for normal and 1 for reverse.
if disp.max_trend>=0
    fault.slip_type = 1; %Normal fault
else
    fault.slip_type = 0; %Reverse fault
end

%Set the parameters of the displacement ellipse.
fault.lu = disp.trend_length/2;
fault.lv = disp.trend_height/2;
fault.max_val = disp.max_trend;
fault.u0 = disp.ctr(1);
fault.v0 = disp.ctr(2);
%disp.ctr has 3 values labelled x, y, and z in the .grid file. But from
%looking at some examples, it looks like z is always 0 and the values of x
%and y are such that they only make sense if interpreted as (u,v)
%coordinates in the fault-aligned coordinate system rather than as world
%coordinates.

%Set the asymmetry and range.
fault.asymmetry = asymmetry;
fault. Range = Range;
end