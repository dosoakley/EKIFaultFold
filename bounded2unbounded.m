function z = bounded2unbounded(x,mins,maxs)
%Transform from a bounded domain to an infinite one.
y = (x-mins)./(maxs-mins);
z = log(y./(1-y));
z(maxs==mins,:) = 0; %maxs==mins cause z = NaN, but setting it to any number will give the correct value in unbounded2bounded.
z(isinf(z))=sign(z(isinf(z)))*1e100; %Get rid of any infinities. I was using realmax() instead of 1e100, but then I still got Infs when I multiplied that by other things.
% z(abs(z)>1e100)=sign(z(abs(z)>1e100))*1e100;
end