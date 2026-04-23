
Cholera Model - Parameter Estimation (Alpha & Beta)
Author: Magdoleen Saad
Purpose: Estimate transmission parameters from real outbreak data

Model: SEIRB (Susceptible-Exposed-Infected-Recovered-Bacteria)
Estimated parameters:
- alpha: person-to-person transmission rate
- eta: waterborne transmission rate

Fixed parameters from literature:
- gamma = 1/1.4 (progression rate)
- gamma_1 = 1/5 (recovery rate)
- zeta_1 = 10 (bacterial shedding)
- zeta = 0.33 (bacterial decay)
- k = 1,000,000 (half-saturation constant)

Data: 10 time points (days 25, 43-51) with infected counts from 1 to 3826

Method: Non-linear least squares optimization (scipy.optimize.fmin)

Output: Estimated alpha, eta + plot comparing model vs real data


import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
import scipy.optimize

tf=[25,43,44,45,46,47,48,49,50,51]
Inf_pop=[1,2685,2904,3036,3145,3307 ,3359,3525,3708,3826]

plt.scatter(tf,Inf_pop,c='b')
plt.xlabel('Time')
plt.ylabel('Infected population')

Lambda = 4.351584656666667                     #43.498918764
mu=2.9705555555555556e-06    #2.9694e-05

rho=1/10**3
k=10**6

d=1.0793
gamma=1/1.4

gamma_1=1/5
zeta_1=10
zeta=0.33


N=1464906

N0= 1339461
I0= 1
E0= 0
R0= 0
S0 =1339460
B0=1000

def col(variables,t,params):
    
    S=variables[0] 
    E=variables[1] 
    I=variables[2]  
    R=variables[3]
    B=variables[4]
    
    
    alpha=params[0]
    eta=params[1]
    

    dSdt=Lambda-eta*((S*(B**2))/(k+B**2))-alpha*S*I+rho*R-mu*S
    dEdt=eta*((S*(B**2))/(k+B**2))+alpha*S*I-(gamma+mu)*E
    dIdt=(gamma+mu)*E-(gamma_1+mu+d)*I
    dRdt=gamma_1*I-(mu+rho)*R
    dBdt=zeta_1 *I-zeta*B
    
    return ([dSdt,dEdt,dIdt,dRdt,dBdt])


def  loss_function(params,tf,Inf_pop)  :
    
    
    y0=[S0,E0,I0,R0,B0] 
    
    output=odeint(col, y0, tf,args=(params,))
    
    loss=0

    for i in range(len(tf)):
        data_Inf=Inf_pop[i]
        model_Inf=output[i,2]
        
        res=(data_Inf-model_Inf)**2
        
        loss+=res
     
    return(loss)
    
    
params0=np.array([3.3508 /10**9,1.6699/10**6])


minimum=scipy.optimize.fmin(loss_function,params0,args=(tf,Inf_pop))

print(minimum)

alpha_fit=minimum[0]
eta_fit=minimum[1]

params=[alpha_fit,eta_fit]
print(params)

y0=[S0,E0,I0,R0,B0] 

output=odeint(col, y0, tf,args=(params,))
U=output[:,2]

plt.plot(tf[1:],U[1:],c='purple')
plt.xlim([40 ,55])
plt.ylim([2500 , 4000])

plt.xlabel('Time')
plt.ylabel('Infected population')
plt.grid(True)
plt.show()



       
