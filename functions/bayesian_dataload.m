function [y,X,m,d] = bayesian_dataload(name)
data = load(['.\bayesian_data\',name,'.txt']);
y = data(:,end);
if ismember(2,y)
    y = y - 1;
end
X0 = data(:,1:end-1);
zx = normalize(X0);
m = size(zx,1);
X = [ones(m,1),zx];
d = size(X,2);
end