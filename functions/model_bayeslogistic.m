classdef model_bayeslogistic
    %MODELCLASS Summary of this class goes here
    %   Detailed explanation goes here

    properties
        ModelName      % Model name
        dimension
        numsamples
        respones
        designmatrix
        alpha
        pi0_mu
        pi0_sigma
        Lap_sigma
        Fisher
        init_sample
        proposal_sigma
        cholC
    end

    methods
        function obj = model_bayeslogistic(Modelname,alpha,y,X,m,d,scale,varargin)
            obj.ModelName  = Modelname;
            obj.dimension = d;
            obj.proposal_sigma = scale*eye(d);
            obj.cholC = chol(obj.proposal_sigma,'lower');
            obj.numsamples = m;
            obj.designmatrix = X;
            obj.respones = y;
            obj.alpha = alpha;
            obj.pi0_mu = fsolve(@obj.gradlogtargetpdf,zeros(d,1),optimset('Display','off'));
            obj.Fisher = obj.Fisherinfo(obj.pi0_mu);
            obj.pi0_sigma = inv(obj.Fisher);
            obj.Lap_sigma = obj.pi0_sigma;
            if nargin > 7
                obj.init_sample = varargin{1};
            else
                obj.init_sample = obj.initsample();
            end
        end

        %% Initialize parameters
        function theta0 = initsample(obj)
            % a = obj.alpha;
            d = obj.dimension;
            mu = obj.pi0_mu;
            sigma = obj.pi0_sigma;
            theta0 = mvnrnd(mu,sigma);
            % theta0 = mvnrnd(zeros(d,1),a*eye(d));
            % theta0 = reshape(theta0,1,d);
            % d = obj.dimension;
            % theta0 = obj.pi0_mu;
            % theta0 = zeros(d,1);
            theta0 = reshape(theta0,1,d);
        end

        %% log_prior
        function f = logprior(obj,theta)
            a = obj.alpha;
            d = obj.dimension;
            dd = size(theta,1);
            if dd ~= d
                theta = theta';
            end
            % f = -1/2*(theta'*theta)/a;
            % f = log(mvnpdf(theta',zeros(1,lt),a*eye(lt)));
            % theta d \times N; f N \times 1
            f = multi_t_pdf(theta,d,zeros(d,1),eye(d)/a,inf,1); 
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
            % theta d \times N; f N \times 1
            % N = size(y,2);
            % f = 0;
            % for i = 1:N
            %     f = f + y(i)*X(i,:)*theta-log(1+exp(X(i,:)*theta));
            % end
            f = theta'*X'*y-sum(log(1+exp(theta'*X')),2);
        end

        %% target distribution
        function fy = logtargetpdf(obj,theta)
            logp = obj.logprior(theta); % theta d \times N
            logl = obj.loglike(theta); % theta d \times N
            fy = logp+logl;
            % fy = logl;
        end

        %% gradtarget distribution
        function f = gradlogtargetpdf(obj,theta)
            X = obj.designmatrix;
            y = obj.respones;
            a = obj.alpha;
            theta = theta(:);
            xt = X*theta;
            f = X'*(y-1./(1+exp(-xt)))-theta/a;
            % f = X'*(y-1./(1+exp(-xt)));
        end

        %% grad2target distribution
        function G = Fisherinfo(obj,theta)
            X = obj.designmatrix;
            a = obj.alpha;
            d = obj.dimension;
            theta = theta(:);
            xt = X*theta;
            % 若theta在变，则在使用auxiliary不能使用z的几何信息，故在此固定theta
            % theta = obj.apprmean;
            p = 1./(1+exp(-xt));
            v = p.*(1-p);
            Dig = diag(v);
            G = X'*Dig*X+1/a*eye(d);
        end

    end

end

