clear
load Data_XT_rollingmean.mat %data 

%%
% The model sum of squares given in the model structure.
model.ssfun = @f2;

%%
% All parameters are constrained to be positive. The initial
% concentrations are also unknown and are treated as extra parameters.
params1 = {

    
      %parameter
      %name,mean,min.max,mean,std
      {'sigma_a', 5.17,3,14,4.9,1}  %incubation period
%     {'report_w', 0.05,0,1} %report rate
%     {'gamma',10, 2,15}%infection period of asymptomatic infectious population
      {'beta', 1,0.01,5}  %transmission rate 
       {'Theta',0.4,0,1} % the proportion of contact tracing
       %{'Tau1',0.01,0,1} % the testing rate
       %{'Tau2',0.01,0,1} % the testing rate
       {'HLa0',1,0,Inf}
       {'HIa0',1,0,Inf}
       {'HLs0',1,0,Inf}
       {'HLp0',1,0,Inf}
       %{'P',0.25,0.17,0.82}
       %{'Recovery',14,10,21} % the days to recovery
%     {'Iw0', 5,1,11} %the estimated reported cases of first 5 days in wuhan
%     {'c1high',0.9,0,1} %control level1-1
%     {'c1mid',0.5,0,1} %control level1-2
%     {'c1low',0.1,0,1} %control level1-3 
%     {'c2nd',0.15, 0,1}%2nd control
%     {'c_gamma',0.50, 0,1}%gamma*c_gamma:infection period of symptomatic infectious population 

    };

%%
% We assume having at least some prior information on the
% repeatability of the observation and assign rather non informational
% prior for the residual variances of the observed states. The default
% prior distribution is sigma2 ~ invchisq(S20,N0), the inverse chi
% squared distribution (see for example Gelman et al.). The 3
% components (_A_, _Z_, _P_) all have separate variances.
model.S20 = [4];
model.N0  = [1];

City = 'Xingtai'

%%
% First generate an initial chain.
options.nsimu = 100000;
options.stats = 1;
[results, chain, s2chain,sschain]= mcmcrun(model,data,params1,options);



%regenerate chain to convergence
options.nsimu = 100000;
options.stats = 1;
[results2, chain2, s2chain2] = mcmcrun(model,data,params1,options,results);

save('results2.mat','results2')
save('chain2.mat','chain2')
save('s2chain2.mat','s2chain2')


%%
% Chain plots should reveal that the chain has converged and we can
% % use the results for estimation and predictive inference.
  
figure
mcmcplot(chain2,[],results2,'denspanel',2);
saveas(gcf,strcat(City,'_Chain2_posterior_probability'),'epsc')

figure
mcmcplot(chain2,[],results2); %,'pairs'
saveas(gcf,strcat(City,'_Chain2_MCMC'),'epsc')

%%
% Function |chainstats| calculates mean ans std from the chain and
% estimates the Monte Carlfigure
% the integrated autocorrelation time and |geweke| is a simple test
% for a null hypothesis that the chain has converged.

results2.sstype = 1; % needed for mcmcpred and sqrt transformation
tmp = chainstats(chain2,results2) %statistic results of parameter estimation


T = table(Names,tmp(:,1),tmp(:,2),tmp(:,3),tmp(:,4),tmp(:,5), 'VariableNames', { 'Parameter','Mean', 'std','MC_err','tau','geweke'} )
writetable(T, strcat(City,'_chain2_parameter_stats.csv'))



%%
% In order to use the |mcmcpred| function we need
% function |modelfun| with input arguments given as
% |modelfun(xdata,theta)|. We construct this as an anonymous function.
modelfun = @(d,th) f3(d(:,1),th,d);


%%o
% We sample 1000 parameter realizations from |chain| and |s2chain|
% and calculate the predictive plots.
nsample = 1000;
results2.sstype = 1;
out = mcmcpred(results2,chain2,s2chain2,data.xdata,modelfun,nsample);%data.ydata-->data

time=36;
% DHT
figure
subplot(2,1,1)
fillyy(1:time,out.obslims{1,1}{1,1}(3,:),out.obslims{1,1}{1,1}(1,:),[0.8 0.8 0.8]);
hold on 
plot(data.ydata(:,2),'.r');
plot(out.obslims{1,1}{1,1}(2,:));
hold off
ylabel('DHT')
%DHC
subplot(2,1,2)
fillyy(1:time,out.obslims{1,1}{1,2}(3,:),out.obslims{1,1}{1,2}(1,:),[0.8 0.8 0.8]);
hold on 
plot(data.ydata(:,3),'.r');
plot(out.obslims{1,1}{1,2}(2,:));
hold off
ylabel('DHC')

saveas(gcf,strcat(City,'_Fitting_wide'),'epsc')


DownCI = transpose(out.obslims{1,1}{1,1}(1,:))
Mean = transpose(out.obslims{1,1}{1,1}(2,:))
UpCI = transpose(out.obslims{1,1}{1,1}(3,:))

T = table(DownCI,Mean,UpCI, 'VariableNames', { 'DownCI', 'Mean','UpCI'} )
writetable(T, strcat(City,'_fitting_PopulationTesting_wide.csv'))

DownCI = transpose(out.obslims{1,1}{1,2}(1,:))
Mean = transpose(out.obslims{1,1}{1,2}(2,:))
UpCI = transpose(out.obslims{1,1}{1,2}(3,:))

T = table(DownCI,Mean,UpCI, 'VariableNames', { 'DownCI', 'Mean','UpCI'} )
writetable(T, strcat(City,'_fitting_ContactTracing_wide.csv'))




CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));

Mean = transpose(mean(chain2,1))
tmp=transpose(CIFcn(chain2,95))
UpCI = tmp(:,2)
DownCI = tmp(:,1)
T = table(Names,DownCI,Mean,UpCI, 'VariableNames', {'Parameter' ,'DownCI', 'Mean','UpCI'} )

writetable(T, strcat(City,'_fitting_parameter.csv'))




