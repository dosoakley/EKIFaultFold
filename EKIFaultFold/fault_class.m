classdef fault_class < matlab.mixin.Copyable
    %This is a class for a fault, based on the way that faults are
    %represented in Havana and RMS and using the elliptical displacement
    %model of Georgsen et al. (2012).
    properties
        %Basic Set Up:
        name %Fault name (string)
        M %4x4 transformation array
        NU %Number of elements in U direction
        NV %Number of elements in V direction
        U0 %u coordinate of left edge of fault.
        V0 %v coordinate of bottom edge of fault.
        LU %Length of fault in U direction.
        LV %Length of fault in V direction.
        dU %Grid step size in u direction
        dV %Grid step size in v direction
        N %Total number of grid points.
        %Fault surface:
        u %Grid point u coordinates (NUxNV)
        v %Grid point v coordinates (NUxNV)
        f %Distances of fault above reference plane at grid points (NUxNV)
        %Fault displacement:
        displacement %Fault displacement values at grid points. (NUxNV)
        ThicknessAttribute %Level set function. Tip line is at 0 contour. (NUxNV)
        lu %Displacement ellipse length in u coordinates
        lv %Displacement ellipse length in v coordinates
        max_val %Displacement ellipse maximum value
        u0 %Displacement ellipse center u coordinate
        v0 %Displacement ellipse center v coordinate
        slip_type %1 for normal faults, 0 for reverse faults
        asymmetry %Displacement asymmetry. In the range from 0 (all FW) to 1 (all HW).
        Range %Distance normal to fault surface at which displacement dies out (reverse drag radius)
    end
    methods
        function MakeM(F,cx,cy,cz,strike,dip)
            %Make the transformation array for the fault (F)
            %(cx,cy,cz) = Center coordinates of reference plane
            %(strike,dip) = Orientation of reference plane
            F.M = [sin(strike),-cos(strike)*cos(dip),cos(strike)*sin(dip),cx;...
            cos(strike),sin(strike)*cos(dip),-sin(strike)*sin(dip),cy;...
            0,-sin(dip),-cos(dip),cz;...
            0,0,0,1];
        end
        function SurfFromAux(F,trend,vparams,stdev,aux,opt)
            %Create the fault surface from auxiliary variables
            if opt.centered_hierarchical
                GRFrlzt = aux;
            else
                GRFrlzt = stdev*MakeGRF(F.u(:),F.v(:),vparams,aux,opt.variogram,opt.rlzt_method); %Realization of the Gaussian random field.
            end
            fdata = trend+GRFrlzt;
            F.f = reshape(fdata,[F.NU,F.NV]);
        end
        function trend = EllipticalTrend(F,u,v)
            %Calculate the value of the elliptical displacement trend function at a set of points.
            r = sqrt(((u-F.u0)./F.lu).^2+((v-F.v0)./F.lv).^2);
            mu0 = zeros(size(r));
            mu0(r<=1) = 2*(1-r(r<=1)).*sqrt(((1+r(r<=1)).^2)./4-r(r<=1).^2); %Georgsen (2012) Eqn. 3.
            trend =  F.max_val*mu0; %Georgsen (2012) Eqn. 4
            if F.slip_type == 0 %Reverse fault
                trend = - trend;
            end
        end
        function DispFromAux(F,lu,lv,max_val,u0,v0,vparams,stdev,aux,opt)
            %Create the fault displacement from auxiliary variables
            [F.lu,F.lv] = deal(lu,lv);
            F.max_val = max_val;
            [F.u0,F.v0] = deal(u0,v0);
            trend = EllipticalTrend(F,F.u(:),F.v(:));
            CalculateThicknessAttribute(F)
            if opt.centered_hierarchical
                GRFrlzt = aux;
            else
                GRFrlzt = stdev*MakeGRF(F.u(:),F.v(:),vparams,aux,opt.variogram,opt.rlzt_method); %Realization of the Gaussian random field.
            end
            displ = trend+GRFrlzt;
            displ(F.ThicknessAttribute<0) = 0; %Set displacement outside the tip line to 0.
            if F.slip_type == 1 %Normal faults
                displ = abs(displ); %Slip <0 isn't allowed for normal faults.
            else %Reverse faults
                displ = -abs(displ); %Slip >0 isn't allowed for reverse faults.
            end
            F.displacement = reshape(displ,[F.NU,F.NV]);
        end
        function f = InterpolateSurface(F,uq,vq)
            I = griddedInterpolant(F.u,F.v,F.f,'linear','linear'); %Interpolant for the fault surface.
            f = I(uq,vq);
        end
        function d = InterpolateDisplacement(F,uq,vq)
            %Interpolate the displacement at the points (uq,vq).
            d = zeros(size(uq));
            D = griddedInterpolant(F.u,F.v,F.displacement,'linear','nearest');
            mask1 = uq>=F.U0 & uq<=F.U0+F.LU & vq>=F.V0 & vq<=F.V0+F.LV; %Points in the fault grid area.
            d(mask1) = D(uq(mask1),vq(mask1)); 
            mask2 = ~mask1;
            if any(mask2) %If some points are not inside the grid.
                t = InterpolateThicknessAttribute(F,uq(mask2),vq(mask2));
                mask2(mask2) = t>0; %Only use active points (inside the tip line).
                d(mask2) = EllipticalTrend(F,uq(mask2),vq(mask2)); %Outside of the grid area, we just assume the elliptical trend holds without perturbation.
            end
        end
        function t = InterpolateThicknessAttribute(F,uq,vq)
            T = griddedInterpolant(F.u,F.v,F.ThicknessAttribute,'linear','nearest');
            t = T(uq,vq); 
            mask = t>0 & ~(uq>=F.U0 & uq<=F.U0+F.LU & vq>=F.V0 & vq<=F.V0+F.LV); %Active points outside the fault grid area
            if any(mask)
                %If some points are not inside the grid and the nearest
                %neighbor interpolation made them active, we want to check
                %that they are actually inside the elliptical fault tip
                %line.
                r = sqrt(((uq(mask)-F.u0)./F.lu).^2+((vq(mask)-F.v0)./F.lv).^2); %Same as in EllipticalTrend. Goes to 1 at the tip line.
                mask(mask) = r>1; %Only true for points outside the tip line.
                t(mask) = -1; %These are points outside the tip line that were incorrectly set to >0, so we are fixing them.
            end
        end
        function [x,y,z] = MovePts_forward(F,x,y,z)
            %Move the points
            %fault = The fault object on which to move the points.
            %R = The range of the kinematic model
            %gamma = The asymmetry of the kinematic model.
            %(x,y,z) = The points to move. They should be row vectors. z is positive down (depth)
            %A note on sign conventions: This script expects the fault transformation
            %matrix fault.M in the form used by Havana. In this way, z in world
            %coordinates points down. However, v in fault coordinates points updip, not
            %downdip as in Georgsen et al. (2012). The displacement in the fault object
            %is expected to be positive for normal sense displacement and negative for
            %reverse sense displacement. A result of all this is that the displacement
            %in the v direction needs to be opposite in sign from the Dy given by
            %Georgsen et al. (2012)'s Eqns. 7 and 8.
            
            %Convert the points to the fault-aligned coordinate system.
            npts = length(x); %x, y, z should all be the same length.
            pts_xyz = [x;y;z;ones(1,npts)];
            pts_uvw = F.M\pts_xyz;
            [up,vp,wp] = deal(pts_uvw(1,:),pts_uvw(2,:),pts_uvw(3,:));
            
            [~,vp,wp] = ForwardModel(F,up,vp,wp);
            
            %Convert back to the world coordinate system.
            % pts_uv(1,:) = u; %No need to update this because u doesn't change.
            pts_uvw(2,:) = vp;
            pts_uvw(3,:) = wp;
            pts_xyz = F.M*pts_uvw;
            [x,y,z] = deal(pts_xyz(1,:),pts_xyz(2,:),pts_xyz(3,:));
        end
        function [x,y,z,fval] = MovePts_restore(F,x,y,z)
            %Restore the points.
            %fault = The fault object on which to move the points.
            %R = The range of the kinematic model
            %gamma = The asymmetry of the kinematic model.
            %(x,y,z) = The points to move. They should be row vectors. z is positive down (depth)
            %As discussed in Georgsen et al. (2012), there isn't an exact inverse
            %for the displacement function. Therefore, similar to Georgsen et al.
            %(2012), we first move grids of points on the fault in the hanging wall and
            %footwall to make an approximate interpolated guess of the restored state
            %position of each point, and then we further refine that guess.
            %A note on sign conventions: This script expects the fault transformation
            %matrix fault.M in the form used by Havana. In this way, z in world
            %coordinates points down. However, v in fault coordinates points updip, not
            %downdip as in Georgsen et al. (2012). The displacement in the fault object
            %is expected to be positive for normal sense displacement and negative for
            %reverse sense displacement. A result of all this is that the displacement
            %in the v direction needs to be opposite in sign from the Dy given by
            %Georgsen et al. (2012)'s Eqns. 7 and 8.
            
            %Convert the points to the fault-aligned coordinate system.
            npts = length(x); %x, y, z should all be the same length.
            pts_xyz = [x;y;z;ones(1,npts)];
            pts_uvw = F.M\pts_xyz;
            [up,vp,wp] = deal(pts_uvw(1,:),pts_uvw(2,:),pts_uvw(3,:));
            
            %Create hw and fw grids and move them. Since u doesn't change, we don't
            %need to care about it. For w, since this is a 2D grid, we calculate change
            %in w rather than its value so that it can be used to make an estimate for
            %other w values than those of the grid.
            [~,hw_vgrid,hw_wgrid] = ForwardModel(F,F.u(:),F.v(:),F.f(:)+0.1);
            [~,fw_vgrid,fw_wgrid] = ForwardModel(F,F.u(:),F.v(:),F.f(:)-0.1);
            hw_v0_interp = scatteredInterpolant(F.u(:),hw_vgrid,F.v(:));
            fw_v0_interp = scatteredInterpolant(F.u(:),fw_vgrid,F.v(:));
            hw_dw_interp = scatteredInterpolant(F.u(:),hw_vgrid,(F.f(:)+0.1)-hw_wgrid);
            fw_dw_interp = scatteredInterpolant(F.u(:),fw_vgrid,(F.f(:)-0.1)-fw_wgrid);
            
            %Create initial estimates of the starting positions of each point.
            [v0,w0] = deal(zeros(1,npts),zeros(1,npts)); %u = u0 because there is no strike-slip motion.
            fp = InterpolateSurface(F,up,vp);
            hw = wp>=fp;
            v0(hw) = hw_v0_interp(up(hw),vp(hw));
            v0(~hw) = fw_v0_interp(up(~hw),vp(~hw));
            w0(hw) = wp(hw)+hw_dw_interp(up(hw),vp(hw));
            w0(~hw) = wp(~hw)+fw_dw_interp(up(~hw),vp(~hw));

            %Now narrow down each point from this starting point.
%             options = optimoptions('fminunc','Display','off','MaxFunctionEvaluations',300);
            options = optimoptions('fsolve','Display','off','Algorithm','levenberg-marquardt','FunctionTolerance',1.0);
            fval = zeros(size(up));
            for i = 1:npts
                if objfun(F,up(i),vp(i),wp(i),vp(i),wp(i))==0 %In this case the point doesn't need to move at all to reach its current position.
                    [v0(i),w0(i)] = deal(vp(i),wp(i));
                else
                    fun = @(a) objfun(F,up(i),a(1),a(2),vp(i),wp(i));
                    a0 = [v0(i),w0(i)];
                    [a,fval(i),~,~] = fsolve(fun,a0,options); 
                    [v0(i),w0(i)] = deal(a(1),a(2));
                end
            end
            
            %Find points that didn't converge to the 1 m tolerance and try
            %to improve them.
            mask = abs(fval)>1;
            if any(mask) && sum(mask)>3
                Iv = scatteredInterpolant(up(~mask)',vp(~mask)',wp(~mask)',v0(~mask)');
                Iw = scatteredInterpolant(up(~mask)',vp(~mask)',wp(~mask)',w0(~mask)');
                for i = find(mask)
                    a0 = [Iv(up(i),vp(i),wp(i)),Iw(up(i),vp(i),wp(i))];
                    niter = 0;
                    while abs(fval(i))>1 && niter<5
                        fun = @(a) objfun(F,up(i),a(1),a(2),vp(i),wp(i));
                        if ~isempty(a0)
                            [a,fval_new] = fsolve(fun,a0,options);
                            a0 = a; %Use this as the starting point for the next iteration if it hasn't converged yet.
                        else
                            disp('a0 empty') %Why is this occasionally happening?
                            disp(sum(mask))
                            fval_new = fval;
                        end
                        if fval_new<fval(i)
                            [v0(i),w0(i)] = deal(a(1),a(2));
                            fval(i) = fval_new;
                        end
                        niter = niter+1;
                    end
                end
            end
            
            %Convert back to the world coordinate system.
            % pts_uv(1,:) = u0; %No need to update this because u doesn't change.
            pts_uvw(2,:) = v0;
            pts_uvw(3,:) = w0;
            pts_xyz = F.M*pts_uvw;
            [x,y,z] = deal(pts_xyz(1,:),pts_xyz(2,:),pts_xyz(3,:));
        end
        function dist = objfun(F,u0,v0,w0,v,w)
            %Objective function to minimize with fsolve. Distance in (v,w) space
            %between forward modelled point and observed point.
            [~,vm,wm] = ForwardModel(F,u0,v0,w0); %Calculate modelled v and w values.
            dist = sqrt((vm-v).^2+(wm-w).^2);
        end
        function [u,v,w] = ForwardModel(F,u,v,w)
            %This function implements the forward deformation model of Georgsen et al.
            %(2012)
            mask = InterpolateThicknessAttribute(F,u,v)>0; %Only use points within the tip line / active region.
            fp = InterpolateSurface(F,u(mask),v(mask));
            mask1 = abs(w(mask)-fp)<=F.Range; %Only use points within the reverse drag radius.
            mask(mask) = mask1;
            hw = w(mask)>=fp(mask1);
            d = -InterpolateDisplacement(F,u(mask),v(mask));
            %d is negative for the reasons discussed in the "note on sign conventions" above.
            alpha = (1-abs(w(mask)-fp(mask1))/F.Range).^2;
            gamma = F.asymmetry*ones(size(hw));
            gamma(~hw) = F.asymmetry-1;
            v(mask) = v(mask)+gamma.*d.*alpha; %Eqns. 7 and 8 of Georgsen (2012).
            f2 = InterpolateSurface(F,u(mask),v(mask)); %Recalculate with updated vp.
            w(mask) = w(mask)+f2-fp(mask1); %Eqn. 9 of Georgsen (2012).
            %No change in u because no strike-slip deformation.
        end
        function CalculateThicknessAttribute(F)
            %Calculate the ThicknessAttribute for a Havana fault based on
            %the tip line of the elliptic trend. The distance of
            %change_length is based on what I saw in some examples from
            %Havana. Having values other than -1 and 1 is not especially
            %useful, but may help if one wants to contour
            %ThicknessAttribute.
            r = sqrt(((F.u-F.u0)./F.lu).^2+((F.v-F.v0)./F.lv).^2); %Same as in EllipticalTrend. Goes to 1 at the tip line.
            change_length = sqrt((5*F.dU./F.lu).^2+(5*F.dV./F.lv).^2); %Distance over which ThicknessAttribute change from 0 to +/- 1.
            TA = (1-r)/change_length;
            TA(TA>1) = 1;
            TA(TA<-1) = -1;
            F.ThicknessAttribute = TA;
        end
    end
end