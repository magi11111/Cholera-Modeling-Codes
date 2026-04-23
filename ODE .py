
Cholera Model - ODE Solver (SEIRB)
Author: Magdoleen Saad
Purpose: Solve the differential equations of the cholera model

Model equations (SEIRB):
dS/dt = Lambda - eta*(S*B^2/(k+B^2)) - alpha*S*I + rho*R - mu*S
dE/dt = eta*(S*B^2/(k+B^2)) + alpha*S*I - (gamma+mu)*E
dI/dt = (gamma+mu)*E - (gamma_1+mu+d)*I
dR/dt = gamma_1*I - (mu+rho)*R
dB/dt = zeta_1*I - zeta*B

Compartments:
S = Susceptible  |  E = Exposed  |  I = Infected  |  R = Recovered  |  B = Bacteria

Parameters:
alpha = direct transmission  |  eta = environmental transmission
gamma = progression rate     |  gamma_1 = recovery rate
zeta_1 = bacterial shedding  |  zeta = bacterial decay
k = half-saturation constant |  mu = mortality rate


from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

Lambda =3.1726e-04
mu=2.9694e-05

rho=1/10**3
k=10**6

d= 1.0794
gamma=1/(1.4)

gamma_1=1/5
zeta_1=10

zeta=0.33

alpha=9.904985905854343e-07
eta=0.00030816050042329604

N=1464906

N0= 1339461
I0= 1
E0= 0
R0= 0
S0 =1339460
B0=1000
z0=[S0,E0,I0,R0,B0]

t=np.linspace(0,62,1000)

def col(z,t):
      
    S=z[0]
    E=z[1]
    I=z[2]
    R=z[3]
    B=z[4]
    
    dSdt=Lambda-eta*((S*(B**2))/(k+B**2))-alpha*S*I+rho*R-(mu)*S
    dEdt=eta*((S*(B**2))/(k+B**2))+alpha*S*I-(gamma+mu)*E
    dIdt=(gamma+mu)*E-(gamma_1+mu+d)*I
    dRdt=gamma_1*I-(mu+rho)*R
    dBdt=zeta_1 *I-zeta*B
    
    return [dSdt,dEdt,dIdt,dRdt,dBdt]

#print(col([S0,E0,I0,R0,B0], t)) 

z=odeint(col,z0,t)

S=z[:,0]
E=z[:,1]
I=z[:,2]
R=z[:,3]
B=z[:,4]

#plt.semilogy(t,S,'b',label='Suseptiaple')
plt.semilogy(t,E,'r--',label='Exposed')
plt.semilogy(t,I,'orange',label='Infected')
plt.semilogy(t,R,'purple',label='Recoverd')
#plt.semilogy(t,B,'yellow',label='Vibrio')
plt.ylabel('human')
plt.xlabel('Time/days')
plt.legend(['S(t)','E(t)','I(t)','R(t)','B(t)'])
plt.grid(True)
plt.show()
