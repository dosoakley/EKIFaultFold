function [x,y,z,faults_out,success] = MoveModel(x,y,z,faults_in,method,fault_deformation,restoration_max_fval)
%Restore or forward model a model involving multiple faults.
%The faults in faults_in should be in order from oldest to youngest.
%(x,y,z) are points (such horizon points) to be moved.
%method should be either 'forward' or 'restoration'.
%fault_deformation should be either true or false.
%restoration_max_fval should have two entries: the first for horizon points
%and the second for fault points. In restoration, points with fval greater
%than these values are labeled as false in the success array, so that they
%can be rejected later.


faults_out = copy(faults_in);
success = true(size(x)); %Initialize. (For forward method, it's always true; for restoration it may change.)
if strcmp(method,'forward')
    for i = 1:length(faults_in)
        if ~isempty(z)
            [x,y,z] = MovePts_forward(faults_out(i),x,y,z);
        end
        if fault_deformation
            for j = 1:(i-1)
                fault_xyz = faults_out(j).M*[faults_out(j).u(:),faults_out(j).v(:),faults_out(j).f(:),ones(faults_out(j).N,1)]';
                [fault_xyz(1,:),fault_xyz(2,:),fault_xyz(3,:)] = MovePts_forward(faults_out(i),fault_xyz(1,:),fault_xyz(2,:),fault_xyz(3,:));
                fault_uvw = faults_out(j).M\fault_xyz;
                faults_out(j).f = InterpScatteredData(fault_uvw(1,:)',fault_uvw(2,:)',fault_uvw(3,:)',faults_out(j).u,faults_out(j).v);
                faults_out(j).displacement = InterpScatteredData(fault_uvw(1,:)',fault_uvw(2,:)',faults_out(j).displacement(:),faults_out(j).u,faults_out(j).v);
                faults_out(j).ThicknessAttribute = InterpScatteredData(fault_uvw(1,:)',fault_uvw(2,:)',faults_out(j).ThicknessAttribute(:),faults_out(j).u,faults_out(j).v);
                u0v0_xyz = faults_out(j).M*[faults_out(j).u0,faults_out(j).v0,InterpolateSurface(faults_out(j),faults_out(j).u0,faults_out(j).v0),1]';
                [u0v0_xyz(1),u0v0_xyz(2),u0v0_xyz(3)] = MovePts_forward(faults_out(i),u0v0_xyz(1),u0v0_xyz(2),u0v0_xyz(3));
                u0v0_uvw = faults_out(j).M\u0v0_xyz;
                [faults_out(j).u0,faults_out(j).v0] = deal(u0v0_uvw(1),u0v0_uvw(2));
            end
        end
    end
elseif strcmp(method,'restoration')
    for i = length(faults_in):-1:1
        if ~isempty(z)
            [x,y,z,fval] = MovePts_restore(faults_out(i),x,y,z);
            success = success & abs(fval)<=restoration_max_fval(1);
        end
        if fault_deformation
            for j = (i-1):-1:1
                fault_xyz = faults_out(j).M*[faults_out(j).u(:),faults_out(j).v(:),faults_out(j).f(:),ones(faults_out(j).N,1)]';
                [fault_xyz(1,:),fault_xyz(2,:),fault_xyz(3,:),fval] = MovePts_restore(faults_out(i),fault_xyz(1,:),fault_xyz(2,:),fault_xyz(3,:));
                success_fault = abs(fval)<=restoration_max_fval(2);
                fault_uvw = faults_out(j).M\fault_xyz;
                faults_out(j).f = InterpScatteredData(fault_uvw(1,success_fault)',fault_uvw(2,success_fault)',...
                    fault_uvw(3,success_fault)',faults_out(j).u,faults_out(j).v);
                faults_out(j).displacement = InterpScatteredData(fault_uvw(1,success_fault)',fault_uvw(2,success_fault)',faults_out(j).displacement(success_fault)',faults_out(j).u,faults_out(j).v);
                faults_out(j).ThicknessAttribute = InterpScatteredData(fault_uvw(1,success_fault)',fault_uvw(2,success_fault)',faults_out(j).ThicknessAttribute(success_fault)',faults_out(j).u,faults_out(j).v);
                u0v0_xyz = faults_out(j).M*[faults_out(j).u0,faults_out(j).v0,InterpolateSurface(faults_out(j),faults_out(j).u0,faults_out(j).v0),1]';
                [u0v0_xyz(1),u0v0_xyz(2),u0v0_xyz(3),fval] = MovePts_restore(faults_out(i),u0v0_xyz(1),u0v0_xyz(2),u0v0_xyz(3));
                if abs(fval)<=restoration_max_fval(2)
                    u0v0_uvw = faults_out(j).M\u0v0_xyz;
                    [faults_out(j).u0,faults_out(j).v0] = deal(u0v0_uvw(1),u0v0_uvw(2));
                else
                    %I'm interpolating the change in u or v, rather than
                    %the absolute value of u or v so that if changes are 
                    %small, I don't get a large change for a (u0,v0)
                    %located outside the grid.
                    u0 = faults_out(j).u0+InterpScatteredData(faults_out(j).u(success_fault)',faults_out(j).v(success_fault)',fault_uvw(1,success_fault)'-faults_out(j).u(success_fault)',faults_out(j).u0,faults_out(j).v0);
                    v0 = faults_out(j).v0+InterpScatteredData(faults_out(j).u(success_fault)',faults_out(j).v(success_fault)',fault_uvw(2,success_fault)'-faults_out(j).v(success_fault)',faults_out(j).u0,faults_out(j).v0);
                    [faults_out(j).u0,faults_out(j).v0] = deal(u0,v0);
                end
            end
        end
    end
else
    error('Unrecognized method.')
end

end
