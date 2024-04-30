function [latents] = LCM_sim(alpha,beta,data,opts)
    
    % Simple simulation of data under the latent cause model.
    %
    % USAGE: [lik, latents] = LCM_sim(alpha,beta,data,[opts])
    %
    % INPUTS:
    %   alpha - concentration parameter
    %   beta - scaling parameter parameter mapping model CR to output CR
    %   data - single-subject data
    %   opts - options structure
    %
    % OUTPUTS:
    %   lik - log-likelihood
    %   latents - structure containing latent variables:
    %               .b  - simulated beta value
    %               .alpha - simulated alpha value
    %               .CR - predicted CR based on input data sequence and
    %               simulated parameter values
    %
    % Agnes Norbury, Jan 2020: based on function LCM_lik by Sam Gershman, Jan 2019
    
    % set concentration parameter
    opts.alpha = alpha;
    
    % run particle filter
    results = LCM_infer([data.US data.CS], opts);
    
    % use linear regression to fit model output to CR
    N = length(results.V);
    X = results.V;
    b=beta;
    noise=normrnd(0,1,[N,1]);                           % as per Gershman and Hartley (2015)
    CRpred = X*b + noise;                    % predicted CR (with injected noise)
    
    % return latent variables
    if nargout > 0
        latents = results;
        latents.b = b;
        latents.alpha=alpha;
        latents.CR = CRpred;
    end
