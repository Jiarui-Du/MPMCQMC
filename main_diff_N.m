clc
clear
addpath("./functions")
% choose model
modelnumber = 1; % different number for different model
% 1: Linear_boston
% 2: Probit_australian
% 3: logistic_pima
modelname = choosemodel(modelnumber);

% load data and model
filename = ['./data/',modelname,'.mat'];
load(filename)

% pilot run
f = @(x) x;
rsample1 = tMPMHnIS(2^8,d,"iid",mdl,"adaptive",2^8,f,0,inf,'dependent');
mdl.pi0_mu = rsample1.mu(:,end);
mdl.pi0_sigma = rsample1.sigma(:,:,end);

% true value
if modelnumber == 1
    starv = reshape(mdl.apprmean,1,d);
else
    nn = 14;%
    N = 2.^nn;
    m = length(N);
    n = 9;
    L = ones(1,m)*2^n;
    burniniter = floor(burnin./N);
    sequence = "FELFSR";
    resultiid = main_fun_nIS(d,mdl,R,L,N,sequence,...
        proposal,burniniter,f,[],df,flagsequence);
    starv = reshape(resultiid(1,:,:),1,d);
    filename = ['./result/',modelname,'-resultiid.mat'];
    save(filename,'resultiid','starv','L','N','R','burnin')
    % filename = ['./result/',modelname,'-resultiid.mat'];
    % load(filename)
end

% ARB-MP-MCQMC
R = 50;
nn = 5:11;%
N = 2.^nn;
m = length(N);
n = 9;
L = ones(1,m)*2^n;
burnin = 2^12;
flagsequence = 'dependent';
proposal = "adaptive";
df = inf;
burniniter = floor(burnin./N);
sequence = ["iid","FELFSR","sobol-owen"];
resultp = main_fun_nIS(d,mdl,R,L,N,sequence,...
    proposal,burniniter,f,starv,df,flagsequence);
resultpmse = resultp(3*m+1:4*m,:,:);
filename = ['./result/',modelname,'-diffN-pilot.mat'];
save(filename,'resultp','resultpmse','L','N','R','burnin')

% ARB-MP-MCQMC
sequence = "sobol-owen";
flagsequence = 'independent';
resultind = main_fun_nIS(d,mdl,R,L,N,sequence,...
    proposal,burniniter,f,starv,df,flagsequence);
resultindmse = resultind(3*m+1:4*m,:,:);
filename = ['./result/',modelname,'-diffN-pilot-ind.mat'];
save(filename,'resultind','resultindmse','L','N','R','burnin')

disp([modelname,' finish!'])
format short g