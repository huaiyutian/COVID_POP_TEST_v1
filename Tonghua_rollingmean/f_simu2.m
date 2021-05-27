%theta 

load chain2.mat
load results2.mat
xdata=data.xdata;



theta=mean(chain2,1);

City = 'Tonghua'
PopSize = 2147400
theta
%single parameter
%Theta:Contact tracing rate


Theta=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,theta(3)];
res_Theta=zeros(size(data.ydata,1),length(Theta));
for i = 1:length(Theta)
    theta(3)=Theta(i);
    res_Theta(:,i)=fsimu_deterministic(size(data.ydata,1),theta,xdata,PopSize);%Ã¿ÈÕÐÂÔö²¡ÀýÇúÏß
end

T = table(res_Theta(:,1),res_Theta(:,2),res_Theta(:,3),res_Theta(:,4),res_Theta(:,5),res_Theta(:,6),res_Theta(:,7),res_Theta(:,8),res_Theta(:,9),res_Theta(:,10),res_Theta(:,11), 'VariableNames', { 'Theta=0','Theta=0.1', 'Theta=0.2','Theta=0.3','Theta=0.4','Theta=0.5','Theta=0.6','Theta=0.7','Theta=0.8','Theta=0.9',strcat('Theta=',sprintf('%.6f',theta(3)))} )
writetable(T, strcat(City,'_ContactTracingEffect_currentPopulationTesting.csv'))






%-------------------------------------------------------------------------------------

%-----Using the current settings from Tonghua, How many round of population testing

theta = mean(chain2,1);
 
ContactTracing = [0.2:0.02:0.5,theta(3),0]
Break = [2,1,0]
for i = 1:length(Break)
    for j = 1:length(ContactTracing)
        tmp = xlsread(strcat("population_testing_strategy_break",int2str(Break(i)),".xlsx"));
        theta(3) = ContactTracing(j);
        Data_pop_test = tmp(:,2:size(tmp,2));
        N = size(tmp,1);
        res_test=zeros(N,size(Data_pop_test,2));
        for k = 1:size(Data_pop_test,2) 
            tau=Data_pop_test(:,k);
            res_test(:,k)=fsimu_deterministic(N,theta,xdata,PopSize,tau);
        end
        writematrix(res_test,strcat(City,"_DailyCaseNum_UnderDifferentTestingStrategy_break",int2str(Break(i)),"_ContactTracing",sprintf('%.6f',ContactTracing(j)),".csv"))
    end
end


theta = mean(chain2,1);
 
ContactTracing = [0.52:0.02:0.65];
Break = [2];
for i = 1:length(Break)
    for j = 1:length(ContactTracing)
        tmp = xlsread(strcat("population_testing_strategy_break",int2str(Break(i)),".xlsx"));
        theta(3) = ContactTracing(j);
        Data_pop_test = tmp(:,2:size(tmp,2));
        N = size(tmp,1);
        res_test=zeros(N,size(Data_pop_test,2));
        for k = 1:size(Data_pop_test,2) 
            tau=Data_pop_test(:,k);
            res_test(:,k)=fsimu_deterministic(N,theta,xdata,PopSize,tau);
        end
        writematrix(res_test,strcat(City,"_DailyCaseNum_UnderDifferentTestingStrategy_break",int2str(Break(i)),"_ContactTracing",sprintf('%.6f',ContactTracing(j)),".csv"))
    end
end



theta = mean(chain2,1);
 
ContactTracing = [0.1:0.1:0.5];
Break = [0];
for i = 1:length(Break)
    for j = 1:length(ContactTracing)
        tmp = xlsread(strcat("population_testing_strategy_break",int2str(Break(i)),".xlsx"));
        theta(3) = ContactTracing(j);
        Data_pop_test = tmp(:,2:size(tmp,2));
        N = size(tmp,1);
        res_test=zeros(N,size(Data_pop_test,2));
        for k = 1:size(Data_pop_test,2) 
            tau=Data_pop_test(:,k);
            res_test(:,k)=fsimu_deterministic(N,theta,xdata,PopSize,tau);
        end
        writematrix(res_test,strcat(City,"_DailyCaseNum_UnderDifferentTestingStrategy_break",int2str(Break(i)),"_ContactTracing",sprintf('%.6f',ContactTracing(j)),".csv"))
    end
end

%-------------------------------------------------------------------------------------------------------------



%----relative transimmision effect----------------------------------------

tmp = mean(chain2,1);
Ratio = [0.1:0.1:2];
Baseline = tmp(2);

ContactTracing = [0.4];
Break = [2];
theta = mean(chain2,1);


for r = 1:length(Ratio)
    for i = 1:length(Break)
        for j = 1:length(ContactTracing)
            tmp = xlsread(strcat("population_testing_strategy_break",int2str(Break(i)),".xlsx"));
            theta(3) = ContactTracing(j);
            theta(2) = Baseline*Ratio(r);
            Data_pop_test = tmp(:,2:size(tmp,2));
            N = size(tmp,1);
            res_test=zeros(N,size(Data_pop_test,2));
            for k = 1:size(Data_pop_test,2) 
                tau=Data_pop_test(:,k);
                res_test(:,k)=fsimu_deterministic(N,theta,xdata,PopSize,tau);
            end
            writematrix(res_test,strcat(City,"_DailyCaseNum_UnderDifferentTestingStrategy_break",int2str(Break(i)),"_ContactTracing",sprintf('%.6f',ContactTracing(j)),"_Ratio",sprintf('%.6f',Ratio(r)),".csv"))
        end
    end
end

%------------------------------------------------------------------------------



%-----asympotomaic proportion effect--------------------------
theta = mean(chain2,1);
 
%ContactTracing = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,theta(3)];
%Break = [0,1,2,3,4,5];

ContactTracing = [theta(3)];
Break = [2];
Prob = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9];
for i = 1:length(Break)
    for j = 1:length(ContactTracing)
        tmp = xlsread(strcat("population_testing_strategy_break",int2str(Break(i)),".xlsx"));
        theta(3) = ContactTracing(j);
        Data_pop_test = tmp(:,2:size(tmp,2));
        N = size(tmp,1);
        for m =1:length(Prob)
            res_test=zeros(N,size(Data_pop_test,2));
            for k = 1:size(Data_pop_test,2) 
                tau=Data_pop_test(:,k);
                res_test(:,k)=fsimu_P_deterministic(N,theta,xdata,PopSize,tau,Prob(m));
            end
            writematrix(res_test,strcat(City,"_DailyCaseNum_UnderDifferentTestingStrategy_break",int2str(Break(i)),"_ContactTracing",sprintf('%.6f',ContactTracing(j)),"_AymProb",sprintf('%.6f',Prob(m)),".csv"))
        end
        
    end
end


%----------------------------------------------------------------------------------------





%----P.1------------------------

theta = mean(chain2,1);
 
ContactTracing = [0.2:0.02:0.5,theta(3),0];
Break = [2];

theta = mean(chain2,1);
tmp = mean(chain2,1);
theta(2) = tmp(2)*2;


for i = 1:length(Break)
    for j = 1:length(ContactTracing)
        tmp = xlsread(strcat("population_testing_strategy_break",int2str(Break(i)),".xlsx"));
        theta(3) = ContactTracing(j);
        Data_pop_test = tmp(:,2:size(tmp,2));
        N = size(tmp,1);
        res_test=zeros(N,size(Data_pop_test,2));
        for k = 1:size(Data_pop_test,2) 
            tau=Data_pop_test(:,k);
            res_test(:,k)=fsimu_deterministic(N,theta,xdata,PopSize,tau);
        end
        writematrix(res_test,strcat(City,"_DailyCaseNum_UnderDifferentTestingStrategy_break",int2str(Break(i)),"_ContactTracing",sprintf('%.6f',ContactTracing(j)),"_P.1.csv"))
    end
end

%--------------------------------------



%----B.1.1.7------------------------

theta = mean(chain2,1);
 
ContactTracing = [0.2:0.02:0.5,theta(3),0];
Break = [2];

theta = mean(chain2,1);
tmp = mean(chain2,1);
theta(2) = tmp(2)*1.5;


for i = 1:length(Break)
    for j = 1:length(ContactTracing)
        tmp = xlsread(strcat("population_testing_strategy_break",int2str(Break(i)),".xlsx"));
        theta(3) = ContactTracing(j);
        Data_pop_test = tmp(:,2:size(tmp,2));
        N = size(tmp,1);
        res_test=zeros(N,size(Data_pop_test,2));
        for k = 1:size(Data_pop_test,2) 
            tau=Data_pop_test(:,k);
            res_test(:,k)=fsimu_deterministic(N,theta,xdata,PopSize,tau);
        end
        writematrix(res_test,strcat(City,"_DailyCaseNum_UnderDifferentTestingStrategy_break",int2str(Break(i)),"_ContactTracing",sprintf('%.6f',ContactTracing(j)),"_B.1.1.7.csv"))
    end
end

%--------------------------------------


%---middle size outbreak---------------


load Data_TH.mat %data 
load chain2.mat;
load results2.mat
xdata=data.xdata;

theta=mean(chain2,1);

City = 'Tonghua';
PopSize = 2147400;
theta

theta(3) = 0;
N=356;
tau = zeros(N,1);
MiddleSize=5000;
%%%% HS(i),HLa(i),HIa(i),HLs(i),HLp(i),HI(i),HR(i),i
Middle_value = fsimu_MiddleSize_deterministic(N,theta,xdata,PopSize,tau,MiddleSize);


Inital_Middle = Middle_value(size(Middle_value,1),1:7);


theta = mean(chain2,1);
 
ContactTracing = [0.2:0.02:0.5,theta(3),0]
Break = [2]
for i = 1:length(Break)
    for j = 1:length(ContactTracing)
        tmp = xlsread(strcat("population_testing_strategy_break",int2str(Break(i)),".xlsx"));
        theta(3) = ContactTracing(j);
        Data_pop_test = tmp(:,2:size(tmp,2));
        N = size(tmp,1);
        res_test=zeros(N,size(Data_pop_test,2));
        for k = 1:size(Data_pop_test,2) 
            tau=Data_pop_test(:,k);
            res_test(:,k)=fsimu_runFromMiddle_deterministic(N,theta,xdata,PopSize,tau,Inital_Middle);
        end
        writematrix(res_test,strcat(City,"_DailyCaseNum_UnderDifferentTestingStrategy_break",int2str(Break(i)),"_ContactTracing",sprintf('%.6f',ContactTracing(j)),"_MiddleOutBreak.csv"))
    end
end


%----vaccine-----



%--------------Vaccine--------------------



%----P.1-----------------------------


theta=mean(chain2,1);
theta(2) = theta(2)*2;

ContactTracing = [0.35,0.40,0.45];
Break = [2];

Cross_immunity = 0.68;
Vac_vector = [0:0.05:1];

for l = 1:length(Vac_vector)
    Vac = Vac_vector(l);
    Vac_eff = Vac*Cross_immunity;
    HS0 = PopSize*(1-Vac_eff);
    HR0 = PopSize*Vac_eff;
    for i = 1:length(Break)
        for j = 1:length(ContactTracing)
            tmp = xlsread(strcat("population_testing_strategy_break",int2str(Break(i)),".xlsx"));
            theta(3) = ContactTracing(j);
            Data_pop_test = tmp(:,2:size(tmp,2));
            N = size(tmp,1);
            res_test=zeros(N,size(Data_pop_test,2));
            for k = 1:size(Data_pop_test,2) 
                tau=Data_pop_test(:,k);
                res_test(:,k)=fsimu_Vaccine_deterministic(N,theta,xdata,PopSize,tau,HS0,HR0);
            end
            writematrix(res_test,strcat(City,"_DailyCaseNum_UnderDifferentTestingStrategy_break",int2str(Break(i)),"_ContactTracing",sprintf('%.6f',ContactTracing(j)),"_Vaccine",sprintf('%.6f',Vac),"_P.1.csv"))
        end
    end
end


theta=mean(chain2,1);
theta(2) = theta(2)*2;

p=0.5;

ContactTracing = [0.35,0.40,0.45];
Break = [2];

Cross_immunity = 0.68;
Vac_vector = [0:0.05:1];

for l = 1:length(Vac_vector)
    Vac = Vac_vector(l);
    Vac_eff = Vac*Cross_immunity;
    HS0 = PopSize*(1-Vac_eff);
    HR0 = PopSize*Vac_eff;
    for i = 1:length(Break)
        for j = 1:length(ContactTracing)
            tmp = xlsread(strcat("population_testing_strategy_break",int2str(Break(i)),".xlsx"));
            theta(3) = ContactTracing(j);
            Data_pop_test = tmp(:,2:size(tmp,2));
            N = size(tmp,1);
            res_test=zeros(N,size(Data_pop_test,2));
            for k = 1:size(Data_pop_test,2) 
                tau=Data_pop_test(:,k);
                res_test(:,k)=fsimu_Vaccine_deterministic_asym(N,theta,xdata,PopSize,tau,HS0,HR0,p);
            end
            writematrix(res_test,strcat(City,"_DailyCaseNum_UnderDifferentTestingStrategy_break",int2str(Break(i)),"_ContactTracing",sprintf('%.6f',ContactTracing(j)),"_Vaccine",sprintf('%.6f',Vac),"_Asym",sprintf('%.6f',p),"_P.1.csv"))
        end
    end
end




%--------------------------------------------------------------------------------

%----------------reinfection with the same virus with primary infection------------


theta=mean(chain2,1);

ContactTracing = [0.35,0.40,0.45];
Break = [2];

Cross_immunity = 0.80;
Vac_vector = [0:0.05:1];

for l = 1:length(Vac_vector)
    Vac = Vac_vector(l);
    Vac_eff = Vac*Cross_immunity;
    HS0 = PopSize*(1-Vac_eff);
    HR0 = PopSize*Vac_eff;
    for i = 1:length(Break)
        for j = 1:length(ContactTracing)
            tmp = xlsread(strcat("population_testing_strategy_break",int2str(Break(i)),".xlsx"));
            theta(3) = ContactTracing(j);
            Data_pop_test = tmp(:,2:size(tmp,2));
            N = size(tmp,1);
            res_test=zeros(N,size(Data_pop_test,2));
            for k = 1:size(Data_pop_test,2) 
                tau=Data_pop_test(:,k);
                res_test(:,k)=fsimu_Vaccine_deterministic(N,theta,xdata,PopSize,tau,HS0,HR0);
            end
            writematrix(res_test,strcat(City,"_DailyCaseNum_UnderDifferentTestingStrategy_break",int2str(Break(i)),"_ContactTracing",sprintf('%.6f',ContactTracing(j)),"_Vaccine",sprintf('%.6f',Vac),"_PrimaryVirus.csv"))
        end
    end
end





theta=mean(chain2,1);
p=0.5;

ContactTracing = [0.35,0.40,0.45];
Break = [2];

Cross_immunity = 0.80;
Vac_vector = [0:0.05:1];

for l = 1:length(Vac_vector)
    Vac = Vac_vector(l);
    Vac_eff = Vac*Cross_immunity;
    HS0 = PopSize*(1-Vac_eff);
    HR0 = PopSize*Vac_eff;
    for i = 1:length(Break)
        for j = 1:length(ContactTracing)
            tmp = xlsread(strcat("population_testing_strategy_break",int2str(Break(i)),".xlsx"));
            theta(3) = ContactTracing(j);
            Data_pop_test = tmp(:,2:size(tmp,2));
            N = size(tmp,1);
            res_test=zeros(N,size(Data_pop_test,2));
            for k = 1:size(Data_pop_test,2) 
                tau=Data_pop_test(:,k);
                res_test(:,k)=fsimu_Vaccine_deterministic_asym(N,theta,xdata,PopSize,tau,HS0,HR0,p);
            end
            writematrix(res_test,strcat(City,"_DailyCaseNum_UnderDifferentTestingStrategy_break",int2str(Break(i)),"_ContactTracing",sprintf('%.6f',ContactTracing(j)),"_Vaccine",sprintf('%.6f',Vac),"_Asym",sprintf('%.6f',p),"_PrimaryVirus.csv"))
        end
    end
end




%--------------------------------------------------------------------

%----reinfection with B.1.1.7-----------------------------


theta=mean(chain2,1);
theta(2) = theta(2)*1.5;

ContactTracing = [0.35,0.40,0.45];
Break = [2];

Cross_immunity = 0.70;
Vac_vector = [0:0.05:1];

for l = 1:length(Vac_vector)
    Vac = Vac_vector(l);
    Vac_eff = Vac*Cross_immunity;
    HS0 = PopSize*(1-Vac_eff);
    HR0 = PopSize*Vac_eff;
    for i = 1:length(Break)
        for j = 1:length(ContactTracing)
            tmp = xlsread(strcat("population_testing_strategy_break",int2str(Break(i)),".xlsx"));
            theta(3) = ContactTracing(j);
            Data_pop_test = tmp(:,2:size(tmp,2));
            N = size(tmp,1);
            res_test=zeros(N,size(Data_pop_test,2));
            for k = 1:size(Data_pop_test,2) 
                tau=Data_pop_test(:,k);
                res_test(:,k)=fsimu_Vaccine_deterministic(N,theta,xdata,PopSize,tau,HS0,HR0);
            end
            writematrix(res_test,strcat(City,"_DailyCaseNum_UnderDifferentTestingStrategy_break",int2str(Break(i)),"_ContactTracing",sprintf('%.6f',ContactTracing(j)),"_Vaccine",sprintf('%.6f',Vac),"_B.1.1.7.csv"))
        end
    end
end



theta=mean(chain2,1);
theta(2) = theta(2)*1.5;
p=0.5;

ContactTracing = [0.35,0.40,0.45];
Break = [2];

Cross_immunity = 0.70;
Vac_vector = [0:0.05:1];

for l = 1:length(Vac_vector)
    Vac = Vac_vector(l);
    Vac_eff = Vac*Cross_immunity;
    HS0 = PopSize*(1-Vac_eff);
    HR0 = PopSize*Vac_eff;
    for i = 1:length(Break)
        for j = 1:length(ContactTracing)
            tmp = xlsread(strcat("population_testing_strategy_break",int2str(Break(i)),".xlsx"));
            theta(3) = ContactTracing(j);
            Data_pop_test = tmp(:,2:size(tmp,2));
            N = size(tmp,1);
            res_test=zeros(N,size(Data_pop_test,2));
            for k = 1:size(Data_pop_test,2) 
                tau=Data_pop_test(:,k);
                res_test(:,k)=fsimu_Vaccine_deterministic_asym(N,theta,xdata,PopSize,tau,HS0,HR0,p);
            end
            writematrix(res_test,strcat(City,"_DailyCaseNum_UnderDifferentTestingStrategy_break",int2str(Break(i)),"_ContactTracing",sprintf('%.6f',ContactTracing(j)),"_Vaccine",sprintf('%.6f',Vac),"_Asym",sprintf('%.6f',p),"_B.1.1.7.csv"))
        end
    end
end

%--------------------------------------------------------------------------------
