function rvs = multi_t_sample(u,d,mean,cholsigma,df)
% u: d(+1) \times N
if df == inf
    rvs = mean+cholsigma*norminv(u(1:d,:)); % d \times N
else
    z = norminv(u);
    uu = normcdf(z);
    X = cholsigma*norminv(uu(1:d,:));
    V = gaminv(u(end,:),df/2,2/df);
    rvs = mean+X./sqrt(V); % d \times N
end

