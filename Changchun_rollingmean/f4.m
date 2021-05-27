function ydot = f4(t,theta,xdata)


T = 28;% the days of outbreak

t=1:T; %time period 2021/1.11-2021/2.7



Tau(1:4) = 0;% 20210111-20210114
Tau(5:6) = 1.0/2; % 20210115-20210116
Tau(7) = 0; % Jan 17
Tau(8:11)=1.0/4; %20210118-20210121.
Tau(12) = 0; % Jan 22
Tau(13:14) = 1.0/2; % 20210123-20210124
Tau(15:16) = 0; % Jan 25-26
Tau(17:18) = 1.0/2; % Jan 27-28
Tau(19) = 0; % Jan 29
Tau(20:21) = 1.0/4; % 20210130-20210131
Tau(22:28) = 0; % 





Epsilon_prime=1/theta(1);%Mean latent period (days)

beta=theta(2);%basic transimission rate


Theta = theta(3);% the proportion of contact tracing 

Gamma=1/2.3;%latent period BUT presymptomatic 
Epsilon=1/2.9;%latent period of symptomatic individuals
Mu = 1/5; % recovery rate % ref [Chadi M. Saad-Roy, Science,2021]

p=0.0187;%proportion of asymptomatic  
afa=1/3.85;%relative infectiousness of asymptomatic individuals
beta_s=0.15;%proportion of presymptomatic transmission equals percent of infectioness [ref.Alberto,2020,nature human behavior]

% data
PChangchun=8569200; %population of Changchun


HLa(1) = theta(4); %latent asymptomatic 
HIa(1) = theta(5); % infectious asymptomatic
HLs(1) = theta(6); %latent symptomatic 
HLp(1) = theta(7); %presymptomatic



HI(1) = 0; %infectious symptomatic
HS(1) = PChangchun; %susceptible
HR(1) = 0; %removed infectious indivial
HN(1)=HS(1);%all population
HT(1) = 0; %the identified cases through the population-level-testing
HC(1) = 3; %the identified cases through the contact tracing.
DHC(1)=3;
DHT(1)=0;


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


