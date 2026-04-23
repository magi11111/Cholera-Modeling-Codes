

lags =3;
   %0.001; %0.145624740;

x0=[12346 1 0 0 ];
t=linspace(0,0.8,10000);
sol = dde23(@codelay,lags,@mydelayhist,t,x0);
%t= linspace(0,260);
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
alpha=1.57*10^(-5);%5.8991*10^(-5);
eta=0.011;%0.2668;
N=12347;


ylag1= z(:,1);


dy = zeros(4,1);
dy(1) = mu*N-eta*(y(1)*(y(4))/(k+y(4)))-alpha*y(1)*y(2)-mu*y(1);
dy(2) = eta*(y(1))*(y(4))/(k+y(4))+alpha*ylag1(1)*ylag1(2)-(mu+gamma_1)*y(2);
dy(3) = gamma_1*y(2)-mu*y(4);
dy(4) = zeta_1 *y(2)-zeta*y(4);
end

function S = mydelayhist(~)

S=[12346 1 0 0];

end
