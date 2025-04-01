function result = main_fun_nIS(d,mdl,R,L,N,sequence,...
    proposal,burniniter,f,starv,df,flagsequence)
ml = length(L);
m = length(N);
Ns = length(sequence);
Np = length(proposal);
Rs = repmat(repmat(reshape(repmat(sequence,R,1),[],1),Np,1),m,1);
Rp = repmat(reshape(repmat(proposal,Ns*R,1),[],1),m,1);
Rr = R*Np*Ns;
Nn = reshape(repmat(N,Rr,1),[],1);
bin = reshape(repmat(burniniter,Rr,1),[],1);
Ln = reshape(repmat(L,Rr,1),[],1);
Rn = Rr*m;
dd = size(f(zeros(1,d)),2);
rmean = zeros(Rn,dd);
% rsample = cell(1,Rn);
tic
parfor i = 1:Rn
    % rng(i)
    n1 = Nn(i);
    seq = Rs(i);
    pro = Rp(i);
    bi = bin(i);
    l = Ln(i);
    rsample = tMPMHnIS(l,d,seq,mdl,pro,n1,f,bi,df,flagsequence);
    rmean(i,:) = rsample.mean;
    if mod(i,50) == 0
        disp(['finish ',num2str(floor(i/Rn*100)),'%'])
    end
end
toc
%
format short g
result = cal_mean(rmean,Ns,m,ml,R,starv);