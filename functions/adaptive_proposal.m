function [R,y,R0] = adaptive_proposal(model,N,d,u,mu,invsigma,C,x,proposal,df)
y = zeros(d,N+1);
y(:,1) = x(:);
if proposal == "adaptive_t"
    cholsigma = C*sqrt(sd*(df-2)/df);
else
    cholsigma = C;
end
y(:,2:N+1) = multi_t_sample(u',d,mu,cholsigma,df);
logkernely = multi_t_pdf(y,d,mu,invsigma,df,1);
log_piy = model.logtargetpdf(y);
y = y'; % y N+1 \times d
R0 = log_piy-logkernely;
R0 = reshape(R0-max(R0),1,N+1);
logR = R0-log(sum(exp(R0)));
R = exp(logR);
end