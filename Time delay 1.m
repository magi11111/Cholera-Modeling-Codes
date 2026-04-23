%
% Cholera Model - Time Delay Analysis Part 1
% Author: Magdoleen Saad
% Purpose: Solve cholera model with time delay using dde23
%
% The time delay represents the incubation period
% (time between exposure and becoming infectious)
%
% This script solves the delayed differential equations (DDEs)
% and plots the results for Susceptible, Infected, Recovered, and Bacteria
%

lags=2.5;%1.261143044;


x0=[12346 1 0 0 ];
sol = dde23(@codelay,lags,@mydelayhist,[0,10000]);
t=linspace(0,0.8,10000);
y= deval(sol,t);
%plot(t,y);

subplot(2,2,1),plot(y(1,:),'-b');
grid on ;
xlabel('Time/days');
ylabel('Susceptible population ');

subplot(2,2,2),plot(y(2,:),'g');
grid on ; xlabel('Time/days');
ylabel('Infected  population');

subplot(2,2,3),plot(y(3,:),'r');
grid on ; xlabel('Time/days');
ylabel('Rrecovered population ');

subplot(2,2,4),plot(y(4,:),'-');
grid on ; xlabel('Time/days');
ylabel('V.colera concentration');


function dy = codelay(~,y,z) 

mu=1/(43.5*365);
k=2*10^6;   
gamma_1=5;
zeta_1=10;
zeta=1/30;
alpha=1.57*10^(-5);%3*10^(-5);%5.8991*10^(-5);
eta=0.011;%0.02;%0.2668;
N=12347;


ylag1= z(:,1);


dy = zeros(4,1);
dy(1) = mu*N-eta*(y(1)*(y(4))/(k+y(4)))-alpha*y(1)*y(2)-mu*y(1);
dy(2) = eta*(ylag1(1))*(ylag1(4))/(k+ylag1(4))+alpha*y(1)*y(2)-(mu+gamma_1)*y(2);
dy(3) = gamma_1*y(2)-mu*y(4);
dy(4) = zeta_1 *y(2)-zeta*y(4);
end

function y = mydelayhist(~)

y=[12346 1 0 0];

end
