import json
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

def z(a, m, T, tau):
    print(tau)
    return tau

def norm(func, a, m, T, tau):
    func_integrate = lambda tau_integrate: func(a, m, T, tau_integrate)
    return integrate.quad(func_integrate**2, 0, T)[0]

def psi1(a, m, T, tau):
    fi1_value=fi1(a, m, T, tau)
    fi1_norm = norm(fi1, a, m, T, tau)
    return fi1_value / fi1_norm

def psi2(a, m, T, tau):
    psi2_hut_value = psi2_hut(a, m, T, tau)
    psi2_hut_norm = norm(psi2_hut,a , m, T, tau)
    return  psi2_hut_value/psi2_hut_norm

def psi2_hut(a, m, T, tau):
    fi2_integrate = lambda tau_integrate: fi2(a, m, T, tau_integrate)
    psi1_integrate = lambda tau_integrate: psi1(a, m, T, tau_integrate)
    a_21=- integrate.quad(fi2_integrate*psi1_integrate, 0, T)[0]
    return a_21*psi1(a, m, T, tau) + fi2(a, m, T, tau)


def normal_func_L2(func, a, m, T, time_start, time_finish):
    sum = 0
    step = (time_finish - time_start) / 100
    time = time_start
    while time_finish - time > 0.00001:
        sum = sum + (func(a, m, T, tau = time + step / 2) )  * step
        print(sum)
        time = time + step
    return sum


file_information_BLA = open('file_information_BLA.json')
Data = json.load(file_information_BLA)
print(Data['info_BLA'][0])

time_start = Data["time_start"]
time_finish = Data["time_finish"]
T = time_finish
using_condition = Data["using_condition"]
data_BLA = Data["info_BLA"]
print(data_BLA)
print(time_start)
print(time_finish)
print(normal_func_L2(z, 10, 10, T, time_start, time_finish))