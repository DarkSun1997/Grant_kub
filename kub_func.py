import math
from scipy import integrate


def fi1(a,
        m,
        T,
        tau):
    return math.exp((-a) / m * (T - tau)) / m

def fi2(a,
        m,
        T,
        tau):
    return (1 - math.exp((-a) / m * (T - tau))) / a

def g1(dparametr_T,
       dparametr_0,
       a,
       m,
       T):
    return dparametr_T - dparametr_0 * math.exp((-a) / m * T)

def g2(parametr_T,
       parametr_0,
       dparametr_0,
       a,
       m,
       T):
    return parametr_T - parametr_0 - dparametr_0 * (1 - math.exp((-a) / m * T)) * m / a

def normal_function_L2(func,
                       a,
                       m,
                       time_start,
                       time_finish,
                       tau):
    func_integrate = lambda tau_integrate: func(a=a, m=m, T=time_finish, tau=tau_integrate) ** 2
    return math.sqrt(integrate.quad(func_integrate, time_start, time_finish)[0])

def psi1(a,
         m,
         time_start,
         time_finish,
         tau):
    return fi1(a=a, m=m, T=time_finish, tau=tau) / \
           normal_function_L2(func=fi1, a=a, m=m, time_start=time_start, time_finish=time_finish, tau=tau)




def a21(a,
        m,
        time_start,
        time_finish,
        tau):
    func_integrate = lambda tau_integrate: - fi2(a=a, m=m, T=time_finish, tau=tau_integrate) * \
                                           psi1(a=a, m=m, time_start=time_start, time_finish=time_finish, tau=tau_integrate)
    return integrate.quad(func_integrate, time_start, time_finish)[0]


def psi2_p(a, m, time_start, time_finish, tau):
    return a21(a=a, m=m, time_start=time_start, time_finish=time_finish, tau=tau)*\
           psi1(a=a, m=m, time_start=time_start, time_finish=time_finish, tau=tau) + \
           fi2(a=a, m=m, T=time_finish, tau=tau)



def z(a,m,T,tau):
    return tau * 3