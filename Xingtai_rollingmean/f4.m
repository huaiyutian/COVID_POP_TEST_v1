function ydot = f4(t,theta,xdata)


T = 36;% the days of outbreak

t=1:T; %time period 2021/1.3-2021/2.27



Tau(1:3) = 0; % 2021.01.03-2021.01.05
Tau(4:5) = 1.0/2; % the first round. 2021.01.06-2020.01.07
Tau(6) = 0; % 2021.01.08
Tau(7:9) = 1.0/3; % the second round 2021.01.9-2021.01.11
Tau(10:14) = 0; % the second round 2021.01.12-2021.01.16

Tau(15:16) = 1.0/2; % the third round % 2021.01.17-2021.01.18

Tau(17) = 0; % the fourth round % 2021.01.19

Tau(18:19) = 1.0/2; % the fifth round % 2021.01.20-2021.01.21
Tau(20) = 0; % Jan 22
Tau(21:22) = 1.0/2; % Jan 23- 24
Tau(23:24) = 1.0/2; % Jan 25-26
Tau(25) = 0; % 2021.01.27

Tau(26:27) = 1.0/2; % Jan 28-29
Tau(28) = 0; % Jan 30
Tau(29:31) = 1.0/3; % Jan 31-Feb 2
Tau(32:33) = 1.0/2; % Feb 3-4
Tau(34) = 0; % Feb 05
Tau(35:36) = 1.0/2; % Feb 06-07





Epsilon_prime=1/theta(1);%Mean latent period (days)

beta=theta(2);%basic transimission rate


Theta = theta(3);% the proportion of contact tracing 

Gamma=1/2.3;%latent period BUT presymptomatic 
Epsilon=1/2.9;%latent period of symptomatic individuals
Mu = 1/5; % recovery rate % ref [Chadi M. Saad-Roy, Science,2021]

p=0.1125;%proportion of asymptomatic  
afa=1/3.85;%relative infectiousness of asymptomatic individuals 
beta_s=0.15;%proportion of presymptomatic transmission equals percent of infectioness [ref.Alberto,2020,nature human behavior]

% data
PXingtai=8013700; %population of Xingtai


HLa(1) = theta(4); %latent asymptomatic 
HIa(1) = theta(5); % infectious asymptomatic
HLs(1) = theta(6); %latent symptomatic 
HLp(1) = theta(7); %presymptomatic



HI(1) = 0; %infectious symptomatic
HS(1) = PXingtai; %susceptible
HR(1) = 0; %removed infectious indivial
HN(1)=HS(1);%all population
HT(1) = 1; %the identified cases through the population-level-testing
HC(1) = 2; %the identified cases through the contact tracing.
DHC(1)=2;
DHT(1)=1;


Beta1=afa*beta;
Beta2=beta_s*beta;
Beta3=beta;



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
tmp_Exporue = poissrnd(Prob_infected,1); % Poisson distribution

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
DHC(i) = Theta*(HLa(i-1)+HIa(i-1)+HLs(i-1)+HLp(i-1)+HI(i-1));
DHT(i)=(1-Theta)*Tau(i-1)*(HLa(i-1)+HIa(i-1)+HLs(i-1)+HLp(i-1)+HI(i-1));


end

%ydot=[HC(:),HT(:)]; % the accumulated obervations
ydot=[DHT(:),DHC(:)]; % the daily obervations


