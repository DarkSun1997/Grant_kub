import math
from scipy import integrate


def fi1(a, m, T, tau):
    return math.exp((-a) / m * (T - tau)) / m

def fi2(a, m, T, tau):
    return (1 - math.exp((-a) / m * (T - tau))) / a

def g1(dparametr_T, dparametr_0, a, m, T):
    return dparametr_T - dparametr_0 * math.exp((-a) / m * T)

def g2(parametr_T, parametr_0, dparametr_0, a, m, T):
    return parametr_T - parametr_0 - dparametr_0 * (1 - math.exp((-a) / m * T)) * m / a

def z(a,m,t,tau):
    return tau * 3

def normal_function_L2(func, a, m, time_start, time_finish, tau):
    func_integrate = lambda tau_integrate: func(a, m, time_finish, tau_integrate) ** 2
    return math.sqrt(integrate.quad(func_integrate, time_start, time_finish)[0])

def phi1(a, m, time_start, time_finish, tau):
    return fi1(a, m, time_finish, tau) / normal_function_L2(fi1, a, m, time_start, time_finish, tau)

def ksi1(a, m, time_start, time_finish, tau):
    return fi1(a, m, time_finish, tau) / normal_function_L2(fi1, a, m, time_start, time_finish, tau)


def a21( a, m, time_start, time_finish, tau):
    func_integrate = lambda tau_integrate: - fi2(a, m, time_finish, tau_integrate) * \
                                           ksi1(a, m, time_start, time_finish, tau_integrate)
    return integrate.quad(func_integrate, time_start, time_finish)[0]


def ksi2_p(a, m, time_start, time_finish, tau):
    return a21()*ksi1(a, m, time_start, time_finish, tau) + fi2(a, m, T, tau)