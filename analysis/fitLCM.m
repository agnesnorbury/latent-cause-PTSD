% script to run latent cause modeling (LCM) described in Gershman, Blei and
% Niv (2010), https://github.com/sjgershm/LCM. See Gershman &
% Hartley, 2015, for example implementation in human Pavlovian fear
% conditioning (GSR) data.

clear; close all;

%-initialise file paths
%==========================================================================
homedir='/Users/agnesnorbury/Desktop/causal_gen/iCBT_gen/code_share';
%where to find raw data:
datadir=[homedir '/data'];
%where to find LCM codebase (see github link in header for source):
LCM_path=[homedir '/analysis/LCM']; addpath(LCM_path)
%where to save LCM results:
resdir=[homedir '/results']; if ~exist(resdir); mkdir(resdir); end

%-specify sessions and relevant task info
%==========================================================================
sessions={'pilot'};    
%task vars:
n_t=90;                 %number of trials
n_t_exclRecall=60;      %conditioning and extinction trials only 
n_t_block=10;           %number of trials per block
n_c=2;                  %number of CSs

%-prepare data for LCM fitting 
% [not needed here, dead code provided for reference on data generation]
%==========================================================================
% USAGE: results = LCM_fit(data,[opts])
%
% INPUTS:
%   data - [nSubjects x 1] structure containing the following fields:
%           .CR - [nTrials x 1] conditioned response
%           .CS - [nTrials x nCues] conditioned stimului **AN: nb this should be binary vectors
%           for presence/absence of each CS on each trial - not CS identity on each trial**
%           .US - [nTrials x 1] unconditioned response
%
% Sam Gershman, Jan 2019

% for sess=1:numel(sessions)
%     %get files for subjects who completed that session
%     cd([datadir '/' sessions{sess}])
%     data_f=dir('*.csv');
%     
%     %for each subject
%     for s=1:numel(data_f)
%         
%         %first, get participant ID from output file:
%         tmp=readcell([data_f(s).name]);
%         sID=cellstr(tmp{5,11});                          %cell where subject ID is stored in csv file        
%         sub_ID=regexpi(sID, '[A-Z][\w]{6}','match');     %extract 7 digit alphanumeric subject ID code (case insensitive)
%         
%         %then, extract task data:
%         tmp=readmatrix(data_f(s).name);
%         tmp=tmp(:,[5,6,13,14,15,16]);        %get cols: 1.randomization group, 2.RT, 3.rating, 4.trial_indx, 5.trial_type, 6.CS (1=CS+, 2=CS-)
%         tmp=tmp(~isnan(tmp(:,4)),:);         %remove all rows without a task trial_indx number (col 4) [nb trial_indx starts at 0]
%         
%         %and add column representing trial outcome (only loss trials are CS+-US i.e. trial_type 2)
%         for i=1:length(tmp), if tmp(i,5)==2, tmp(i,7)=1; else, tmp(i,7)=0; end; end
%         
%         %NB we will fit LCM model only to initial conditioning and *first
%         %two* extinction blocks, so model is unbiased by later trials:
%         LCMdata=tmp(1:(n_t_exclRecall-n_t_block),:);
%         data(s).subID=sub_ID{1,1};
%         data(s).randgroup=LCMdata(1,1);
%         %add binary vectors for presence/absence of each CS on each trial
%         %(row 1=CS+, row2=CS-):
%         for i=1:length(LCMdata)
%             if LCMdata(i,6)==1          
%                 data(s).CS(i,1)=1;
%                 data(s).CS(i,2)=0;
%             elseif LCMdata(i,6)==2
%                 data(s).CS(i,1)=0;
%                 data(s).CS(i,2)=1; 
%             end
%         end
%         data(s).CR=zscore(LCMdata(:,3)); %z score ratings within-ppt
%         data(s).US=LCMdata(:,7);         %0=no loss; 1=loss
%         data(s).RT=LCMdata(:,2);         %trial RT in ms
%         
%         %while we're here let's also create long format data for other analysis (e.g. in R)
%         data_long_s(1:length(LCMdata),1)=s;
%         data_long_s(:,2)=1:1:length(LCMdata);
%         data_long_s(:,3)=LCMdata(:,6);
%         data_long_s(:,4)=LCMdata(:,3);
%         data_long_s(:,5)=zscore(LCMdata(:,3));
%         data_long_s(:,6)=LCMdata(:,7);
%         data_long_s(:,7)=LCMdata(:,2);
%         data_long = vertcat(data_long, data_long_s);
%     
%     end
%
%     %save .mat data
%     save([resdir '/' sessions{sess} '_raw_dataZ'],'data')
%     
%     %save .csv (long format) data
%     headers = {'subID','trial','CS','CR','CRz','US','RT'}; headers = strjoin(headers, ',');
%     fid = fopen([resdir '/' sessions{sess} '_raw_data_long.csv'], 'w'); 
%     fprintf(fid,'%s\n', headers);
%     fclose(fid);
%     dlmwrite([resdir '/' sessions{sess} '_raw_data_long.csv'], data_long, '-append');
%     
% end
 

%-fit latent cause model to data
%==========================================================================
% USAGE: results = LCM_fit(data,[opts])
%
% opts (optional) - structure defining LCM options (see LCM_opts.m)
% DEFAULTS:
%           opts.M = 100        (number of particles)
%           opts.a = 1          (hyperparameter of beta prior: pseudo-count for feature presence)
%           opts.b = 1          (hyperparameter of beta prior: pseudo-count for feature absence)
%           opts.alpha = 0      (concentration parameter for Chinese restaurant process prior)
%           opts.stickiness = 0 (stickiness parameer for Chinese restaurant process prior)
%           opts.K = 10         (maximum number of latent causes)
%
% OUTPUTS:
%   results - [nSubjects x 1] structure containing the following fields:
%               .alpha - concentration parameter values
%               .P - posterior probability distribution over alpha
%               .lik - log-likelihood for each alpha value
%               .latents - latent variables for each alpha value (see LCM_infer)
%               .logBF - log Bayes factor for the alpha>=0 model
%                       relative to the alpha=0 model
%
% Sam Gershman, Jan 2019

for sess=1:numel(sessions)
    
    %load session data
    load([datadir '/' sessions{sess} '_raw_dataZ'])
   
    %fit LCM to data
    opts.a=1; opts.b=1; 
    %opts.M=1;       %for MAP (can run this instead of particle filter for quick first pass)
    opts.M = 500;    %for higher accuracy particle filtering model estimation
    [results] = LCM_fit(data,opts);
    
    %save results
    %save([resdir '/' sessions{sess} '_LCM_results_MAP'],'results') %for quick MAP approximation
    save([resdir '/' sessions{sess} '_LCM_results_PF'],'results')   %for higher fidelity particle filtering
      
end
 
%-visualise model fit and output
%==========================================================================
% results = LCM_infer(X,opts)
% OUTPUTS:
%   results - structure wih the following fields:
%       .opts - options (missing fields set to defaults)
%       .V - [T x 1] US prediction
%       .post - [T x K] latent cause posterior, where post(t,k) is the
%            probability of latent cause k being active on trial t,
%            after observing the all the features. K (the number of
%            latent causes) is determined adaptively by the model.
opts.K=10; %maximum number of latent causes (default)
nAlpha=50; %number of alpha values evaluated (default)
res_long=[];

for sess=1:numel(sessions)    
    
    load([datadir '/' sessions{sess} '_raw_dataZ.mat'])
    load([resdir '/' sessions{sess} '_LCM_results_PF.mat'])
    
    %set up figures
    p1=figure; p2=figure; p3=figure;
    
    for s=1:numel(data)       
        %look at distribution of raw ratings data across trials
        figure(p1);
        subplot(ceil(sqrt(numel(data))),ceil(sqrt(numel(data))),s)
        histogram(data(s).CR,10)

        %extract some stuff from inside 
        LCs=zeros(numel(results(s).latents), opts.K);    %preallocate size for LCs as not all matrices same size) (K=10 is LCM default)
        LCs_t=zeros(length(data(s).CR), opts.K);         %preallocate size
        for p=1:numel(results(s).latents)
            Vs(p,:)=results(s).latents(p).V;             %extract series of stimulus (US prediction) values for each alpha value
            pred(p,:)=results(s).latents(p).CR;          %extract series of predicted CRs for each alpha value
            betas(p,:)=results(s).latents(p).b;          %extract beta parameter value for each alpha value
            %for overall latent cause estimates (average across trials, per alpha value
            LCs(p,1:length(mean(results(s).latents(p).post)))=mean(results(s).latents(p).post);    %extract posterior over latent causes for each alpha value
            %for probability-weighted latent cause estimates by trial
            LCs_t_p=zeros(length(data(s).CR), opts.K);
            LCs_t_p(:,1:length(mean(results(s).latents(p).post)))=results(s).latents(p).post.*results(s).P(p); %extract posterior over latent causes *per trial* multiplied by posterior prob of that alpha value
            LCs_t=(LCs_t+LCs_t_p);  %then just sum over all probability-weighted state estimates
        end
        
        %look at actual vs predicted CR (rating) on each trial:
        predCR=mean(pred.*results(s).P)*nAlpha;   %multiply by posterior P over alpha values and sum
        predV=mean(Vs.*results(s).P)*nAlpha; 
        figure(p2);
        subplot(ceil(sqrt(numel(data))),ceil(sqrt(numel(data))),s)
        scatter(data(s).CR,predCR)
        r=round(corrcoef(data(s).CR,predCR),3);
        title(['{\it r}= ' num2str(r(2))]) 
        Rs(s,1)=r(2);
        
        %extract posterior distribution over latent causes for this sub:
        LCpost_all(s,:)=mean(LCs.*results(s).P)*nAlpha;  %multiply by posterior P over alpha values and sum
        
        %look at posterior distributions over alpha values:
        figure(p3);
        subplot(ceil(sqrt(numel(data))),ceil(sqrt(numel(data))),s)
        plot(linspace(0,10,nAlpha), results(s).P)
        title(['a=' num2str(round(results(s).alpha,1)) '; lBF=' num2str(round(results(s).logBF,1))])
        ylim([0 0.15])
        
        %put results in long format for further analysis
        res_long_s(1:length(data(s).CR),1)=results(s).logBF;
        res_long_s(1:length(data(s).CR),2)=results(s).alpha;
        res_long_s(1:length(data(s).CR),3)=sum(betas.*results(s).P);
        res_long_s(:,4:opts.K+3)=LCs_t;
        res_long_s(:,opts.K+4)=predV';
        res_long_s(:,opts.K+5)=predCR';
        res_long = vertcat(res_long, res_long_s);
    end
    
    %save long format data+results
    data_res_long=horzcat(data_long, res_long);
    headers = {'subID','trial','CS','CR','CRz','US','RT',...
                'logBF','alpha','beta','LC1','LC2','LC3','LC4','LC5','LC6','LC7','LC8','LC9','LC10','predV','predCR'}; 
    headers = strjoin(headers, ',');
    fid = fopen([resdir '/' sessions{sess} '_data_results_long.csv'], 'w'); 
    fprintf(fid,'%s\n', headers);
    fclose(fid);
    dlmwrite([resdir '/' sessions{sess} '_data_results_long.csv'], data_res_long, '-append');

end

%plot mean distribution of latent causes (across all trials) across ppts
%1. averaged across all subjects/alpha values, weighted by posterior P:
p4=figure;
bar(1:1:opts.K,mean(LCpost_all))

%-random permutation of ratings data / model fit:
%=========================================================================
nperm=1000;     %number of permuated datasets to generate

for sess=1:numel(sessions)
    
    %load session raw real data
    load([resdir '/' sessions{sess} '_raw_data'])
    
    %extract actual ratings data across all subjects 
    for s=1:numel(data) 
        ratings_raw(s,:)=data(s).CR;
    end
    %reshape to single vector
    ratings_raw=reshape(ratings_raw,1,numel(data)*(n_t_exclRecall-n_t_block));
    %randomly permute ratings and store against real trial structure as randomly-resampled 'subjects'
    for p=1:nperm
        i=randperm(length(ratings_raw),50);
        data_perm(p).CR=zscore(ratings_raw(:,i))';
        data_perm(p).CS=data(1).CS;
        data_perm(p).US=data(1).CS;    
    end
    save([resdir '/' sessions{sess} '_data_permZ'],'data_perm')
    
    %fit LCM model to randomly permuted data:
    opts.a=1; opts.b=1; 
    opts.M=1;       %for MAP
    %opts.M = 500;  %for higher fidelity particle filtering
    [results] = LCM_fit(data_perm,opts);
    
    %save results
    save([resdir '/' sessions{sess} '_LCM_results_perm'],'results')      
end

%look at output
opts.K=10; %maximum number of latent causes (default)
nAlpha=50; %number of alpha values evaluated
for sess=1:numel(sessions)    
    
    load([resdir '/' sessions{sess} '_LCM_results_perm']) 
    
    for s=1:nperm      
        %extract some stuff from inside 
        LCs=zeros(numel(results(s).latents),opts.K);    %preallocate size for LCs as not all matrices same size) (K=10 is LCM default)
        for p=1:numel(results(s).latents)
            pred(p,:)=results(s).latents(p).CR;         %extract series of predicted CRs for each alpha value
            LCs(p,1:length(mean(results(s).latents(p).post)))=mean(results(s).latents(p).post);    %extract posterior over latent causes for each alpha value
        end
        
        %look at actual vs predicted CR (rating) on each trial:
        %1. averaged across all alpha values, weighted by posterior P:
        predCR=pred.*results(s).P;  %multiply by posterior P over alpha values
        predCR=mean(predCR)*nAlpha; %rescale by number of alpha values for visualisation
        r=round(corrcoef(data_perm(s).CR, predCR),3);
        Rperms(s,1)=r(2);
        
        %look at distribution of latent cause posteriors over subjects
        LCpost_all(s,:)=mean(LCs.*results(s).P)*nAlpha;
        LCpost_alphaMax(s,:)=LCs(results(s).P==max(results(s).P),:);  %for maximum a posterior alpha value, average latent cause estimates over all trials
        
    end

end

%plot random vs observed distribution of model fit (R values for actual vs model-predicted) ratings data)
%cdf
figure;
cdfplot(Rperms)
hold on
cdfplot(Rs)
%kernel density (smoothed)
figure;
ksdensity(Rperms)
hold on
ksdensity(Rs)

%permutation test for difference in distributions of R values
%here, we are using the permutation difference test by lkrol
%see: https://github.com/lrkrol/permutationTest
[p, observeddifference, effectsize] = permutationTest(Rs, Rperms, 10000, 'plotresult', 1);

