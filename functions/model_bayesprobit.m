classdef model_bayesprobit
    %MODELCLASS Summary of this class goes here
    %   Detailed explanation goes here

    properties
        ModelName      % Model name
        dimension
        numsamples
        respones
        designmatrix
        pi0_mu
        pi0_sigma
        init_sample
        pri_sigma
        proposal_sigma
        Lap_sigma
        cholC
    end

    methods
        function obj = model_bayesprobit(Modelname,y,X,m,d,sigma,scale,varargin)
            obj.ModelName  = Modelname;
            obj.dimension = d;
            obj.proposal_sigma = scale*eye(d);
            obj.cholC = chol(obj.proposal_sigma,'lower');
            obj.numsamples = m;
            obj.designmatrix = X;
            obj.respones = y;
            obj.pri_sigma = sigma;
            obj.pi0_mu = fsolve(@obj.gradlogtargetpdf,zeros(d,1),optimset('Display','off'));
            obj.Lap_sigma = inv(obj.Lapsigma(obj.pi0_mu));
            obj.pi0_sigma = inv(obj.Fisherinfo(obj.pi0_mu));
            if nargin > 7
                obj.init_sample = varargin{1};
            else
                obj.init_sample = obj.initsample();
            end
        end

        %% Initialize parameters
        function theta0 = initsample(obj)
            d = obj.dimension;
            % theta0 = obj.pi0_mu;
            theta0 = mvnrnd(obj.pi0_mu,obj.pi0_sigma);
            % theta0 = zeros(d,1);
            theta0 = reshape(theta0,1,d);
        end

        %% log_prior
        function f = logprior(obj,theta)
            d = obj.dimension;
            dd = size(theta,1);
            sigma = obj.pri_sigma;
            if dd ~= d
                theta = theta';
            end
            f = multi_t_pdf(theta,d,zeros(d,1),diag(1./sigma),inf,1);
        end

        %% log_likelihood
        function f = loglike(obj,theta)
            X = obj.designmatrix;
            Y = obj.respones;
            d = obj.dimension;
            dd = size(theta,1);
            if dd ~= d
                theta = theta';
            end
            xb = X*theta;
            pxb = normcdf(xb);
            Y1 = Y == 1;
            Y0 = Y == 0;
            f(Y1,:) = log(pxb(Y1,:));
            f(Y0,:) = log(1-pxb(Y0,:));
            f = sum(f,1);
            % f = Y'*log(pxb)+(1-Y')*log(1-pxb);
            % When pxb_i near 1 and Y_i is 1,
            % this equation will also calculate log(1-pxb_i)
            f = f';
        end

        %% target distribution
        function fy = logtargetpdf(obj,theta)
            logp = obj.logprior(theta); % theta d \times N
            logl = obj.loglike(theta); % theta d \times N
            fy = logp+logl;
        end

        %% gradtarget distribution
        function f = gradlogtargetpdf(obj,theta)
            X = obj.designmatrix;
            Y = obj.respones;
            sigma = obj.pri_sigma;
            theta = theta(:);
            xb = X*theta;
            cxb = normcdf(xb);
            f = X'*((Y-cxb)./(cxb.*(1-cxb)).*normpdf(xb))-theta./sigma';
        end

        %% grad2target distribution
        function G = Fisherinfo(obj,theta)
            X = obj.designmatrix;
            theta = theta(:);
            sigma = obj.pri_sigma;
            xb = X*theta;
            cxb = normcdf(xb);
            G = X'*((normpdf(xb)).^2./(cxb.*(1-cxb)).*X)+diag(1./sigma);
        end

        function [G,GY0,GY1] = Lapsigma(obj,theta)
            X = obj.designmatrix;
            Y = obj.respones;
            theta = theta(:);
            sigma = obj.pri_sigma;
            xb = X*theta;
            cxb = normcdf(xb);
            pxb = normpdf(xb);
            Y1 = Y == 1;
            Y0 = Y == 0;
            zcxb = cxb.*xb;
            GY1 = (X(Y1,:))'*(pxb(Y1)./(cxb(Y1)).^2.*(pxb(Y1)+zcxb(Y1)).*X(Y1,:));
            GY0 = (X(Y0,:))'*(pxb(Y0)./(1-cxb(Y0)).^2.*(pxb(Y0)+zcxb(Y0)-xb(Y0)).*X(Y0,:));
            G = GY1+GY0+diag(1./sigma);
        end

        function [GIS,G] = Lapsigma2(obj,theta)
            X = obj.designmatrix;
            Y = obj.respones;
            N = length(Y);
            theta = theta(:);
            sigma = obj.pri_sigma;
            xb = X*theta;
            cxb = normcdf(xb);
            pxb = normpdf(xb);
            d = obj.dimension;
            G = zeros(d,d,N);
            for i = 1:N
                XTXi = (X(i,:))'*X(i,:);
                pxbi = pxb(i);
                cxbi = cxb(i);
                xbi = xb(i);
                if Y(i) == 1
                    G(:,:,i) = pxbi/cxbi^2*(pxbi^2+xbi*cxbi)*XTXi;
                else
                    G(:,:,i) = pxbi/(1-cxbi)^2*(pxbi^2+xbi*cxbi-xbi)*XTXi;
                end
            end
            GIS = sum(G,3)+diag(1./sigma);
        end
    end

end

