function x = unbounded2bounded(z,mins,maxs)
%Reverse the transformation done by bounded2unbounded.
y = 1./(1+exp(-z));
x = y.*(maxs-mins)+mins;
end