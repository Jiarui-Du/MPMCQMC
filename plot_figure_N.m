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

filename = ['./result/',modelname,'-diffN-pilot-ind.mat'];
load(filename)

filename = ['./result/',modelname,'-diffN-pilot.mat'];
load(filename)

% figure
m = size(resultpmse,1);
result_mse = [resultpmse,resultindmse];
result_rmse = sqrt(result_mse);
result_tmse = sum(result_mse,3);
mm = size(result_tmse,2);
result_tmrf = zeros(m,mm);
for i = 1:mm
    result_tmrf(:,i) = result_tmse(:,1)./result_tmse(:,i);
end
result_trmse = sqrt(result_tmse);
result_trrf = zeros(m,mm);
for i = 1:mm
    result_trrf(:,i) = result_trmse(:,1)./result_trmse(:,i);
end
index = [1,ceil(m/2),m];
% table for presentation
table = [result_trmse(index,1),result_trrf(index,2:4)]';
smv = log2(result_rmse);
tsmv = log2(result_trmse);

% convergence rate
nn = log2(N);
kk = 1;
for i = 1:mm
k(i,:) = polyfit(nn(kk:end),tsmv(kk:end,i),1);
end

% convergence plot of total rmse
figure
style = ["o","p","x","s","d","*"];
colors = ["c","b","k","r","g","m"];
for i = 1:mm
plot(nn,tsmv(:,i),'color',colors(i+1),'marker',style(1),'linestyle','-','LineWidth',2)
hold on
end
xx = nn-nn(1);
plot(nn,-1/2*xx+tsmv(1,1),'color','k','marker','none','linestyle','--','LineWidth',1.5)
hold on
plot(nn,-1*xx+tsmv(1,3),'color','k','marker','none','linestyle','-.','LineWidth',1.5)
hold on
plot(nn,-1/2*xx+tsmv(1,2),'color','k','marker','none','linestyle','--','LineWidth',1.5)
xlabel('log_2(N)');
ylabel('log_2(RMSE)');
% ylim([tsmv(end,3)-4 tsmv(1,1)+1.5])
legend(["IID","FELFSR","Sobol'-de","Sobol'-inde","Slope:-1/2","Slope:-1"],...
    'Fontsize',16,'Location','southwest','NumColumns',2)
set(gca,'FontSize',16)
set(gca,'LooseInset',get(gca,'TightInset'))
saveas(gcf,['./figs/',modelname,'_tvar_diffN.eps'],'epsc')
% end