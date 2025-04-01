function f = multi_t_pdf(x,d,mean,invsigma,df,flog)
% x d \times N
% f N \times 1
if df == inf
    % f = log(mvnpdf(x',mean',sigma));
    % f = -1/2*(x-mean)'*invsigma*(x-mean);
    temp = (x-mean)';
    f = -1/2*sum((temp*invsigma).*temp, 2); 
else
    % logC = log(gamma((df+d)/2)) - (log(gamma(df/2))+1/2*log(det(sigma))+d/2*(log(df)+log(pi)));
    % f = logC -(df+d)/2*log(1+(x-mean)'*inv(sigma)*(x-mean)/df);
    temp = (x-mean)';
    f = -(df+d)/2*log(1+sum((temp*invsigma).*temp, 2)/df);
end
if flog == 0
    f = exp(f);
end