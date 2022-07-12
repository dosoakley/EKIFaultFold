function [x_lines,y_lines,z_lines,fault_num,hw] = GetFaultLines(faults,x,y,z)
%Find where a surface intersects one or more faults.
%Create lines offset very slightly into the HW and FW from the
%intersection.

%Currently, I only use 1 contour line, but there's no real reason not to
%use them all if there is more than one.

%Interpolate on fault grid.
%x, y, and z just need to be row vectors.
%This is simpler, but the downside is that fault lines won't extend beyond the edges of the fault grid.
[x_lines,y_lines,z_lines,fault_num] = deal([],[],[],[]);
hw = true(0,0);
for i = 1:length(faults)
    bed_uv = faults(i).M\[x;y;z;ones(size(x))];
    f = InterpolateSurface(faults(i),bed_uv(1,:),bed_uv(2,:)); %I'm using true here b/c AssignFaultBlock above uses true.
    dist = bed_uv(3,:)-f; %Distance from the fault.
%     dist_f = InterpScatteredData(bed_uv(1,:)',bed_uv(2,:)',dist',faults(i).u,faults(i).v);
%     I = scatteredInterpolant(bed_uv(1,:)',bed_uv(2,:)',dist','linear','none');
    I = scatteredInterpolant(bed_uv(1,:)',bed_uv(2,:)',dist','linear','linear');
    dist_f = I(faults(i).u,faults(i).v);
    dist_f(faults(i).ThicknessAttribute<0) = NaN; %These won't get included in the contouring this way.
    uvec = faults(i).U0+(0:faults(i).dU:faults(i).LU);
    vvec = faults(i).V0+(0:faults(i).dV:faults(i).LV);
    if ~isempty(dist_f)
        c = contourc(uvec,vvec,dist_f',[0,0]);
        if ~isempty(c)
            n = 1;
            contour_line = [];
            while n<size(c,2) %Ideally, there should only be 1 contour, but in case there is more than 1, we do this.
                contour_line = [contour_line,c(:,n+1:n+c(2,n))]; %Add the points in this contour.
                n = n+c(2,n)+1;
            end
            [u_line,v_line] = deal(contour_line(1,:),contour_line(2,:));
            w_line = InterpolateSurface(faults(i),u_line,v_line);
            uvw_line2 = [[u_line,u_line];[v_line,v_line];[w_line+1,w_line-1];ones(1,2*length(u_line))]; %2 lines, in hw and fw.
            xyz_line2 = faults(i).M*uvw_line2;
            [x_lines,y_lines,z_lines] = deal([x_lines,xyz_line2(1,:)],[y_lines,xyz_line2(2,:)],[z_lines,xyz_line2(3,:)]);
            fault_num = [fault_num,i*ones(1,size(xyz_line2,2))];
            hw = [hw,true(1,length(u_line)),false(1,length(u_line))]; %Tells if points are in the hw (true) or fw (false).
        end
    end
end

% %Interpolate on horizon grid.
% %x, y, and z should be able to be reshaped into a grid in meshgrid format.
% [x_lines,y_lines,z_lines] = deal([],[],[]);
% for i = 1:length(faults)
%     [xvec,yvec] = deal(unique(x,'sorted'),unique(y,'sorted'));
%     bed_uv = faults(i).M\[x;y;z;ones(size(x))];
%     f = InterpolateSurface(faults(i),bed_uv(1,:),bed_uv(2,:)); %I'm using true here b/c AssignFaultBlock above uses true.
%     dist = bed_uv(3,:)-f; %Distance from the fault.
%     dd1 = reshape(dist,[length(yvec),length(xvec)]);
%     c = contourc(xvec,yvec,dd1,[0,0]);
%     if ~isempty(c)
%         n = 1;
%         contour_line = [];
%         while n<size(c,2) %Ideally, there should only be 1 contour, but in case there is more than 1, we do this.
%             contour_line = [contour_line,c(:,n+1:n+c(2,n))]; %Add the points in this contour.
%             n = n+c(2,n)+1;
%         end
%         [x_line,y_line] = deal(contour_line(1,:),contour_line(2,:));
%         z_line = InterpScatteredData(x',y',z',x_line,y_line);
%         uvw_line = faults(i).M\[x_line;y_line;z_line;ones(size(x_line))];
%         uvw_line2 = [[uvw_line(1:2,:),uvw_line(1:2,:)];[uvw_line(3,:)+1,uvw_line(3,:)-1];[uvw_line(4,:),uvw_line(4,:)]]; %2 lines, in hw and fw.
%         xyz_line2 = faults(i).M*uvw_line2;
%         [x_lines,y_lines,z_lines] = deal([x_lines,xyz_line2(1,:)],[y_lines,xyz_line2(2,:)],[z_lines,xyz_line2(3,:)]);
%     end
% end

%The second method seems to be slower.
end