function modelname = choosemodel(modelnumber)
allmodel = ["Linear_boston","Probit_australian","logistic_pima"];

modelname = char(allmodel(modelnumber));

switch modelname
    case 'Linear_boston'
        % (506*14) boston house data from R package MASS
        data = load('./data/bostondata.txt');
        Y = log(data(:,end));
        X = data(:,1:end-1);
        X(:,5) = (data(:,5)).^2;
        X(:,6) = (data(:,6)).^2;
        X(:,8) = log(data(:,8));
        X(:,9) = log(data(:,9));
        X(:,13) = log(data(:,13));
        m = length(X);
        X = [ones(m,1),normalize(X)];
        d = size(X,2);
        sigma = 0.5;
        scale = 1;
        mdl = model_bayeslinear(modelname,X,Y,d,sigma,scale);
        %
        filename = ['./data/',modelname,'.mat'];
        save(filename,"Y","X","d","sigma","mdl");
    case 'Probit_australian'
        % (392*9) Machine Learning Repository. Pima indians diabetes data set, 2012d.
        [Y,zX,M,d] = bayesian_dataload('australian');
        sigma = [20,5*ones(1,d-1)];
        % sigma = 100*ones(1,d);
        scale = 0.05;
        mdl = model_bayesprobit(modelname,Y,zX,M,d,sigma,scale);
        %
        filename = ['./data/',modelname,'.mat'];
        save(filename,"Y","zX","M","d","sigma","mdl");
    case 'logistic_pima'
        % (392*9)Machine Learning Repository. Pima indians diabetes data set, 2012d.
        load('./data/pima.mat')
        xx = [];
        for i = 2:8
            xx(:,i-1) = X(:,i)==0;
        end
        Ix = ~any(xx,2);
        DataX = X(Ix,:); % Eliminating outliers
        M = size(DataX,1);
        dataX = [ones(M,1),DataX];
        Y = y(Ix,:);
        d = size(dataX,2);
        zX = [dataX(:,1),normalize(dataX(:,2:end))];
        alpha = 1;
        scale = 0.5;
        mdl = model_bayeslogistic(modelname,alpha,Y,zX,M,d,scale);
        %
        filename = ['./data/',modelname,'.mat'];
        save(filename,"Y","zX","M","d","mdl");
end
