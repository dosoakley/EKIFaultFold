function faults = ApplyTruncations(faults,truncations)
%Apply truncation relations to a group of faults.
%faults is an array of fault objects of the class fault_class.
%truncations is a square array in which the row specifies the fault being 
%truncated and column the fault that truncates it, with the meaning of the
%values being: 0 = no truncation rule, -1 = truncated in FW, 1 = truncated in HW.

nfaults = length(faults);
for i = 1:nfaults %Fault being truncated.
    for j = 1:nfaults %Fault doing the truncating.
        if truncations(i,j) ~= 0
            %Create a gridded surface for the displacement (u,v) grid, setting it to 0
            xyz_data_surf = faults(i).M*[faults(i).u(:),faults(i).v(:),faults(i).f(:),ones(faults(i).N,1)]';
            uv_data_surf = (faults(j).M\xyz_data_surf); %Truncated fault (i) surface points converted to truncating fault (j) coordinate system.
            f_surf = InterpolateSurface(faults(j),uv_data_surf(1,:),uv_data_surf(2,:));
            hw_surf = uv_data_surf(3,:)>=f_surf; %Points in the truncated fault (i) that are in the hanging wall of the truncating fault (j).
            a_surf = InterpolateThicknessAttribute(faults(j),uv_data_surf(1,:),uv_data_surf(2,:))>0; %Determine if the nearest part of the truncating surface is active.
            if truncations(i,j)==1
                faults(i).displacement(hw_surf & a_surf) = 0; %In hw of truncating fault and within its tip line.
                faults(i).ThicknessAttribute(hw_surf & a_surf) = -1;
            elseif truncations(i,j)==-1
                faults(i).displacement(~hw_surf & a_surf) = 0; %In fw of truncating fault and within its tip line.
                faults(i).ThicknessAttribute(~hw_surf & a_surf) = -1;
            else
                disp('Error: Unrecognized value in fault truncations matrix.')
            end
        end
    end
end
end