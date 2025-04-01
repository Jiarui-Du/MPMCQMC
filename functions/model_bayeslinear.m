classdef model_bayeslinear
    %MODELCLASS Summary of this class goes here
    %   Detailed explanation goes here

    properties
        ModelName      % Model name
        dimension
        numsamples
        respones
        noise
        designmatrix
        proposal_sigma
        cholC
        sigma
        Sigma_0
        apprmean
        apprcov
        pi0_mu
        pi0_sigma
        Lap_sigma
        Fisher
        init_sample
    end

    methods
        function obj = model_bayeslinear(Modelname,X,obs,d,sigma,scale,varargin)
            obj.ModelName  = Modelname;
            obj.dimension = d;
            obj.proposal_sigma = scale*eye(d);
            obj.cholC = chol(obj.proposal_sigma,'lower');
            obj.sigma = sigma;
            obj.numsamples = size(X,1);
            obj.designmatrix = X;
            obj.respones = obs;
            obj.Sigma_0 = sigma*obj.numsamples*inv(X'*X);
            obj.apprmean = obj.target_mu();
            obj.apprcov = obj.target_sigma();
            obj.pi0_mu = fsolve(@obj.gradlogtargetpdf,zeros(d,1),optimset('Display','off'));
            obj.Fisher = obj.Fisherinfo(obj.pi0_mu);
            obj.pi0_sigma = obj.target_sigma();
            obj.Lap_sigma = obj.target_sigma();
            if nargin > 6
                obj.init_sample = varargin{1};
            else
                obj.init_sample = obj.initsample();
            end
        end

        %% Initialize parameters
        function theta0 = initsample(obj)
            d = obj.dimension;
            mu0 = obj.pi0_mu;
            sigma0 = obj.pi0_sigma;
            % theta0 = zeros(d,1);
            theta0 = mvnrnd(mu0,sigma0);
            theta0 = reshape(theta0,1,d);
        end

        %% log_prior
        function f = logprior(obj,theta)
            d = obj.dimension;
            dd = size(theta,1);
            if dd ~= d
                theta = theta';
            end
            f = multi_t_pdf(theta,d,zeros(d,1),inv(obj.Sigma_0),inf,1); 
        end

        %% log_likelihood
        function f = loglike(obj,theta)
            X = obj.designmatrix;
            y = obj.respones;
            d = obj.dimension;
            dd = size(theta,1);
            if dd ~= d
                theta = theta';
            end
            temp = (X*theta-y)';
            f = -1/(2*obj.sigma)*sum(temp.*temp, 2); 
        end

        %% target distribution
        function fy = logtargetpdf(obj,theta)
            logp = obj.logprior(theta);
            logl = obj.loglike(theta);
            fy = logp+logl;
        end

        %% gradtarget distribution
        function f = gradlogtargetpdf(obj,theta)
            X = obj.designmatrix;
            y = obj.respones;
            theta = theta(:);
            G = obj.Fisherinfo(theta);
            f = -G*theta+X'*y/obj.sigma;
        end

        %% grad2target distribution
        function G = Fisherinfo(obj,~)
            X = obj.designmatrix;
            G = X'*X/obj.sigma+inv(obj.Sigma_0);
        end
        
        function u = target_mu(obj)
            X = obj.designmatrix;
            y = obj.respones;
            XSigma = X'*X+obj.sigma*inv(obj.Sigma_0);
            u = inv(XSigma)*X'*y;
        end

        function u = target_sigma(obj)
            X = obj.designmatrix;
            XSigma = X'*X+obj.sigma*inv(obj.Sigma_0);
            u = obj.sigma*inv(XSigma);
        end
            

    end

end

