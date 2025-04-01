function rsample = tMPMHnIS(L,d,sequence,model,proposal,N,f,burniniter,df,flagsequence)
x = model.initsample();
if contains(proposal,"adaptive_t")
    dd = d+1;
else
    dd = d;
end

if flagsequence == "independent"
    u = [];
    for k = 1:L
        uuk = seqfunMH(sequence,N,dd,log2(N));
        u = [u;uuk];
    end
else
    u = seqfunMH(sequence,L*(N),dd,log2(L*N));
end

u = [rand(burniniter*N,dd);u];
% prepare
xlt = zeros((L+burniniter)+1,d);
xlt(1,:) = x;
fd = size(f(x),2);
es = zeros(L+burniniter,fd);
mu_0 = model.pi0_mu;
mu_0 = mu_0(:);
d = size(mu_0,1);
sigma_0 = model.pi0_sigma;
weight_mean = zeros(d,L+burniniter+1);
weight_cov = zeros(d,d,L+burniniter+1);
weight_mean(:,1) = mu_0;
weight_cov(:,:,1) = sigma_0;
weight = zeros(L+burniniter,N+1);

% iteration
for k = 1:L+burniniter
    uk = u((k-1)*N+1:k*N,:);
    C_0 =chol(sigma_0,'lower');
    invsigma_0 = inv(sigma_0);
    [W,y,~] = adaptive_proposal(model,N,d,uk,mu_0,invsigma_0,C_0,x,proposal,df);
    mu_k = W*y; % W 1 \times N+1; y N+1 \times d
    mu_0 = mu_0+(mu_k'-mu_0)/(k+1);
    weight_mean(:,k+1) = mu_0;
    temp1 = y'-mu_0; % d\times N+1
    sigma_k = temp1*(W'.*temp1');
    sigma_0 = sigma_0+(sigma_k-sigma_0)/(k+1);
    weight_cov(:,:,k+1) = sigma_0;
    I = disc_sample(N+1,W,rand);
    weight(k,:) = W;
    es(k,:) = W*f(y);
    x = y(I,:);
    xlt(k+1,:) = x;
end
rsample.es = es;
rsample.weight = weight;
rsample.mean = mean(es(burniniter+1:end,:),1);
rsample.mu = weight_mean;
rsample.sigma = weight_cov;
rsample.sample = xlt;
end