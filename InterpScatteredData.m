function zq = InterpScatteredData(x,y,z,xq,yq)
%This function calls scatteredInterpolant but if necessary deals with
%problems of too few data or NaN values that I may otherwise encounter.
% Interpolant = scatteredInterpolant(x,y,z,'linear','nearest'); %Use nearest extrapolation to avoid huge values.););
Interpolant = scatteredInterpolant(x,y,z,'natural','nearest'); %Use nearest extrapolation to avoid huge values.););
zq = Interpolant(xq,yq);
if any(isnan(zq))
    %I'm not sure why this happens, but it does occasionally,
    %and this seems to fix it.
    mask = isnan(zq);
%         Interpolant2 = scatteredInterpolant(xq(~mask),yq(~mask),zq(~mask),'linear','nearest'); %Use nearest extrapolation to avoid huge values.
    Interpolant2 = scatteredInterpolant(xq(~mask),yq(~mask),zq(~mask),'natural','nearest'); %Use nearest extrapolation to avoid huge values.
    zq(mask) = Interpolant2(xq(mask),yq(mask));
end
if isempty(zq)
    %This can happen if there are only 2 points in x, y, and z or if
    %the points are all colinear, and maybe for other reasons.
    zq = mean(z)*ones(size(xq));
end
end