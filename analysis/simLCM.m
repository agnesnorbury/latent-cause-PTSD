% script to run latent cause modeling (LCM) described in Gershman, Blei and
% Niv (2010), https://github.com/sjgershm/LCM. See Gershman &
% Hartley, 2015, for example implementation in human Pavlovian fear
% conditioning (GSR) data.

clear; close all;

%-initialise file paths
%==========================================================================
homedir='/Users/agnesnorbury/Desktop/causal_gen/iCBT_gen';
%where to find LCM toolbox (for latent cause model fitting)
LCM_path=[homedir '/analysis/LCM']; addpath(LCM_path)
%where to find LCM results
resdir=[homedir '/results'];

%-specify sessions and relevant task info
%==========================================================================
sessions={'pilot'};
sess=1;

%-set simulation / recovery model options
%==========================================================================
opts.K=10;       %maximum number of latent causes (default)
opts.a=1;        %uniform prior across possible stimulus configurations
opts.b=1;        %uniform prior across possible stimulus configurations
opts.M = 500;    %particle filtering with M=500 particles

%-first, let's look at the distribution of posterior parameter estimates in
%-our sample
%==========================================================================
%load particle filter fit results
load([resdir '/' sessions{sess} '_LCM_results_PF.mat'])

%get posterior parameter estimates
alphas=[results.alpha];
betas=[results.beta];

%visualise distribution of alphas and fit simple pdf
figure;
hist(alphas); hold on;
alpha_pdf=fitdist(alphas', 'ExtremeValue');
xvals=0:0.1:10;
y=pdf(alpha_pdf,xvals);
plot(xvals, y*20);  %rescale y for easy vis
xlabel('posterior alpha'); ylabel('count');

%visualise distribution of betas anbd fit simple pdf
figure;
hist(betas); hold on;
beta_pdf=fitdist(betas', 'ExtremeValue');
xvals=-1:0.01:1.5;
y=pdf(beta_pdf,xvals);
plot(xvals, y*10);  %rescale y for easy vis
xlabel('posterior beta'); ylabel('count');

%-now let's try simulation and recovery analysis based on the empirical
%-distribution of these parameters in our sample
%==========================================================================
n_val=250;                                      %number of simulations to run

% Step 1: sample from observed pdfs
alpha_sim=random(alpha_pdf,n_val,1); 
alpha_sim(alpha_sim<0)=0;                       %set limit of alpha>=0
beta_sim=random(beta_pdf,length(alpha_sim),1);

% Step 1: simulate CRs based on real sequence of stimuli experienced by
% study ppts. Since all participants saw the same sequence of CSs and USs,
% we can just take stim config data from the first ppt:
load([resdir '/' sessions{sess} '_raw_dataZ_all'])  
stim.CS=data(1).CS;
stim.US=data(1).US;
data=stim;

for a=1:length(alpha_sim)
    %define alpha and beta
    alpha=alpha_sim(a);
    beta=beta_sim(a);
    %run simulation
    results_sim(a)=LCM_sim(alpha, beta, data, opts);
end

% Step 3: fit model to simulated data
[data_sim(1:length(alpha_sim)).CS]=deal(stim.CS);
[data_sim(1:length(alpha_sim)).US]=deal(stim.US);
for s=1:numel(alpha_sim)
    tmp3(s).CRz=zscore(results_sim(s).CR);
end
[data_sim(1:length(alpha_sim)).CR]=tmp3.CRz;
%fit
[sim_recov_results] = LCM_fit(data_sim, opts);
save([resdir '/' sessions{sess} '_LCM_sim_recover_results_PF'],'sim_recov_results', 'alpha_sim') 

%plot simulated alpha against recovered alpha
figure;
scatter(alpha_sim, [sim_recov_results.alpha]); hold on;
p=polyfit(alpha_sim', [sim_recov_results.alpha], 1); %linear trendline
f=polyval(p,[0:0.1:6]);
plot([0:0.1:6],f,'-');
xlabel('simulated alpha'); ylabel('recovered posterior alpha');
[r,p]=corrcoef(alpha_sim, [sim_recov_results.alpha]); %simple corr coef
aa_sr_r=diag(r,1)
aa_sr_p=diag(p,1)

%plot simulated alpha against recovered beta
figure;
scatter(alpha_sim, [sim_recov_results.beta]); hold on;
p=polyfit(alpha_sim', [sim_recov_results.beta], 1); %linear trendline
f=polyval(p,[0:0.1:6]);
plot([0:0.1:6],f,'-');
xlabel('simulated alpha'); ylabel('recovered posterior beta');
[r,p]=corrcoef(alpha_sim, [sim_recov_results.beta]); %simple corr coef
ab_sr_r=diag(r,1)
ab_sr_p=diag(p,1) 

%also regress simulated alpha values against all recovered (estimated)
%param values, to check for identifiability:
m1=fitlm([[sim_recov_results.alpha]',[sim_recov_results.beta]'], alpha_sim);
%"Of particular interest here is the relative amount of variance in each
%estimated parameter that can be explained by variations in simulated
%parameters"
an1=anova(m1)
rel_var1=an1.SumSq(1) / an1.SumSq(2)
%"any non-diagonal element in the matrix of regression coefficients signals
%a potential non-identifiability issue between the corresponding
%parameters"
m1.CoefficientCovariance

%-since we are interested in detecting logBF, lets enrich for alpha=0 and
%-look at logBF alpha=0 vs alpha>0
%==========================================================================
% Step 1: sample from observed pdfs
alpha_sim=random(alpha_pdf,n_val,1); 
alpha_sim(alpha_sim<0)=0;                          %set limit of alpha>=0
alpha_sim=vertcat(alpha_sim,ones(25,1)*0);         %enrich for alpha=0
beta_sim=random(beta_pdf,length(alpha_sim),1);

% Step 1: simulate CRs based on real data stim
data=stim;
for a=1:length(alpha_sim)
    %define alpha and beta
    alpha=alpha_sim(a);
    beta=beta_sim(a);
    %run simulation
    results_sim(a)=LCM_sim(alpha, beta, data, opts);
end

% Step 3: fit model to simulated data
[data_sim(1:length(alpha_sim)).CS]=deal(stim.CS);
[data_sim(1:length(alpha_sim)).US]=deal(stim.US);
for s=1:numel(alpha_sim)
    tmp3(s).CRz=zscore(results_sim(s).CR);
end
[data_sim(1:length(alpha_sim)).CR]=tmp3.CRz;
%fit
[sim_recov_results_a0] = LCM_fit(data_sim, opts);
save([resdir '/' sessions{sess} '_LCM_sim_recover_results_PF_a0enriched'],'sim_recov_results_a0', 'alpha_sim') 

%plot logBF (recovered) by simulated alpha group (> or = 0)
clear tmp;
tmp.alpha_Group=(alpha_sim>0);
tmp.logBF=[sim_recov_results_a0.logBF]';
y=[mean(tmp.logBF(tmp.alpha_Group==0)), mean(tmp.logBF(tmp.alpha_Group==1))];
sd=[std(tmp.logBF(tmp.alpha_Group==0)), std(tmp.logBF(tmp.alpha_Group==1))];
x=categorical({'alpha=0','alpha>0'});
figure;
bar(x,y); hold on;              %bar plot of recovered logBF values
er=errorbar(x,y,sd);            %add error bars   
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
ylim([-1.5,3.5]);
xlabel('simulated alpha'); ylabel('recovered model logBF[alpha>0]');
plot(x(1),tmp.logBF(tmp.alpha_Group==0),'ko');   %add raw data
plot(x(2),tmp.logBF(tmp.alpha_Group==1),'ko');   %add raw data

[h,p,ci,stats]=ttest2(tmp.logBF(tmp.alpha_Group==0),tmp.logBF(tmp.alpha_Group==1));   %Welch's 2-sample t-test (unequal variances)

%-finally, let's just try and simuluate some data based on different alphas
%==========================================================================
alpha_test=[0;0.001;0.002;0.003;0.20;0.30;1;3];   %alpha values to simulate
n_sims=25;                                        %number of repetitions per value
alpha_sim=repmat(alpha_test,n_sims,1);

for a=1:length(alpha_sim)
    %define alpha and beta 
    alpha=alpha_sim(a);
    beta=1;
    %define data
    data=stim;
    %run simulation
    results_sim(a)=LCM_sim_noNoise(alpha, beta, data, opts);
end
%save results
save([resdir '/' sessions{sess} '_LCM_ap_simulations'],'results_sim') 

%extract some behavioural summary data for visualisation:
n_t_block=10;               %number of trials per block
for s=1:length(alpha_sim)
    tmp(:,1:2)=stim.CS;
    tmp(:,3)=results_sim(s).CR;
    b_i=1;
    for b=1:length(tmp)/n_t_block
        tmp2=tmp(b_i:b_i+n_t_block-1,:);
        CSplus_rating_byBlock(s,b)=mean(tmp2(tmp2(:,1)==1,3));
        CSminus_rating_byBlock(s,b)=mean(tmp2(tmp2(:,2)==1,3)); 
        b_i=b_i+n_t_block;
    end       
end
