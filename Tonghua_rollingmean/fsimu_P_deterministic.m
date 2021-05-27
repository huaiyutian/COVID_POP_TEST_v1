
function ydot = fsimu(T,theta,xdata,parpop,par_tau,p)

t=1:T; %time period 2021/1.21-2021/2.7

Tau(1:3) = 0;
Tau(4:7) = 1.0/4; %the first round. Jan 15, 2021---Jan 18, 2021.
Tau(8) = 0;
Tau(9:11)=1.0/3; %the second round. the start date is Jan 20, 2021.
Tau(12:13)=0;
Tau(14:16)=1.0/3; % the third round. the start date is Jan 25, 2021.
Tau(17:18)=0;

Tau(19:27)=0; % 




beta=theta(2);%basic transimission rate

Epsilon_prime=1/theta(1);%Mean latent period (days)

Theta = theta(3);% the proportion of contact tracing 

Gamma=1/2.3;%latent period BUT presymptomatic 
Epsilon=1/2.9;%latent period of symptomatic individuals
Mu = 1/5; % recovery rate

%p=0.0346;%proportion of asymptomatic  
afa=1/3.85;%relative infectiousness of asymptomatic individuals 
beta_s=0.15;%proportion of presymptomatic transmission equals percent of infectioness [ref.Alberto,2020,nature human behavior]

% data

%population of Tonghua
 Ptonghua=parpop;

%Population testing strategy: 
if (nargin > 4)
        Tau(1:T) = par_tau;
end

HLa(1) = theta(4); %latent asymptomatic 
HIa(1) = theta(5); % infectious asymptomatic
HLs(1) = theta(6); %latent symptomatic 
HLp(1) = theta(7); %presymptomatic



HI(1) = 0; %infectious symptomatic
HS(1) = Ptonghua; %susceptible
HR(1) = 0; %removed infectious indivial
HN(1)=HS(1);%all population
HT(1) = 0; %the identified cases through the population-level-testing
HC(1) = 5; %the identified cases through the contact tracing.
DHC(1)=5;
DHT(1)=0;


Beta1=afa*beta;
Beta2=beta_s*beta;
Beta3=beta;


DHS(1)= HLa(1)+HIa(1)+HLs(1)+HLp(1)+HI(1);

Cutoff = 1;
for i = 2:T 
if HIa(i-1)<Cutoff
 HIa(i-1)=0;
end
if HI(i-1)<Cutoff
 HI(i-1)=0;
end
if HLp(i-1)<Cutoff
 HLp(i-1)=0;
end

Prob_infected = (Beta1.*HIa(i-1)+Beta3.*HI(i-1)+Beta2.*HLp(i-1))/HN(i-1);
Prob_infected = HS(i-1)*Prob_infected;
tmp_Exporue = Prob_infected;

%suspecitble population

HS(i)=HS(i-1)-tmp_Exporue;

%asymptomatic population

HLa(i)=HLa(i-1)+p*tmp_Exporue-(Tau(i-1)+Theta-Tau(i-1)*Theta)*HLa(i-1)-(1-Tau(i-1))*(1-Theta)*Epsilon_prime*HLa(i-1); %latent
HIa(i)=HIa(i-1)+(1-Tau(i-1))*(1-Theta)*Epsilon_prime*HLa(i-1)-(Tau(i-1)+Theta-Tau(i-1)*Theta).*HIa(i-1)-(1-Tau(i-1))*(1-Theta)*Mu*HIa(i-1);%infectious

%symptomatic population

HLs(i)=HLs(i-1)+(1-p)*tmp_Exporue-(Tau(i-1)+Theta-Tau(i-1)*Theta)*HLs(i-1)-(1-Tau(i-1))*(1-Theta)*Epsilon*HLs(i-1);%latent
HLp(i)=HLp(i-1)+(1-Tau(i-1))*(1-Theta)*Epsilon*HLs(i-1)-(Tau(i-1)+Theta-Tau(i-1)*Theta)*HLp(i-1)-(1-Tau(i-1))*(1-Theta)*Gamma*HLp(i-1);%presymtomatic
HI(i)=HI(i-1)+(1-Tau(i-1))*(1-Theta)*Gamma*HLp(i-1)-(Tau(i-1)+Theta-Tau(i-1)*Theta).*HI(i-1)-(1-Tau(i-1))*(1-Theta)*Mu*HI(i-1); %infectious


% removed population
HR(i)=HR(i-1)+(1-Tau(i-1))*(1-Theta)*Mu.*HI(i-1)+(1-Tau(i-1))*(1-Theta)*Mu.*HIa(i-1);
HT(i)=HT(i-1)+(1-Theta)*Tau(i-1)*(HLa(i-1)+HIa(i-1)+HLs(i-1)+HLp(i-1)+HI(i-1));
HC(i)=HC(i-1)+Theta*(HLa(i-1)+HIa(i-1)+HLs(i-1)+HLp(i-1)+HI(i-1));
%all population
HN(i)=HS(i)+HLa(i)+HIa(i)+HLs(i)+HLp(i)+HI(i)+HR(i)+HT(i)+HC(i); 
% the daily cases

DHS(i)=tmp_Exporue;


end

%ydot=[HC(:),HT(:)]; % the accumulated obervations
ydot=[DHS(:)]; % the daily obervations


