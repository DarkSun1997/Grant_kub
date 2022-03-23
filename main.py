import json
import math
import kub_func
from scipy import integrate
import matplotlib.pyplot as plt
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


def norm(func, a, m, T):
    func_integrate = lambda tau_integrate: func(a, m, T, tau_integrate)
    return integrate.quad(func_integrate ** 2, 0, T)[0]


def psi1(a, m, T, tau):
    fi1_value = fi1(a, m, T, tau)
    fi1_norm = norm(fi1, a, m, T)
    return fi1_value / fi1_norm


def psi2(a, m, T, tau):
    psi2_hut_value = psi2_hut(a, m, T, tau)
    psi2_hut_norm = norm(psi2_hut, a, m, T)
    return psi2_hut_value / psi2_hut_norm


def psi2_hut(a, m, T, tau):
    fi2_integrate = lambda tau_integrate: fi2(a, m, T, tau_integrate)
    psi1_integrate = lambda tau_integrate: psi1(a, m, T, tau_integrate)
    a_21 = - integrate.quad(fi2_integrate * psi1_integrate, 0, T)[0]
    return a_21 * psi1(a, m, T, tau) + fi2(a, m, T, tau)


def B1(dparametr_T, dparametr_0, a, m, T):
    b_11 = 1 / norm(fi1, a, m, T)
    return b_11 * g1(dparametr_T, dparametr_0, a, m, T)


def B2(parametr_T, parametr_0, dparametr_T, dparametr_0, a, m, T):
    fi2_integrate = lambda tau_integrate: fi2(a, m, T, tau_integrate)
    psi1_integrate = lambda tau_integrate: psi1(a, m, T, tau_integrate)
    a_21 = - integrate.quad(fi2_integrate * psi1_integrate, 0, T)[0]
    b_21 = a_21 / (norm(psi2_hut, a, m, T) * norm(fi1, a, m, T))
    b_22 = 1 / norm(psi2_hut, a, m, T)
    return b_21 * g1(dparametr_T, dparametr_0, a, m, T) + b_22 * g2(parametr_T, parametr_0, dparametr_0, a, m, T)


def U(parametr_T, parametr_0, dparametr_T, dparametr_0, a, m, T, tau):
    return B1(dparametr_T, dparametr_0, a, m, T) * psi1(a, m, T, tau) + B2(parametr_T, parametr_0, dparametr_T,
                                                                           dparametr_0, a, m, T) * psi2(a, m, T, tau)


def normal_func_L2(func, a, m, T, time_start, time_finish):
    sum = 0
    step = (time_finish - time_start) / 100
    time = time_start
    while time_finish - time > 0.00001:
        sum = sum + (func(a, m, T, tau=time + step / 2)) * step
        print(sum)
        time = time + step
    return sum


def eler(q,
         w,
        time,
         time_step,
        parametr,
        a,
        m,
        time_start,
        time_finish,
        tau):
    k1 = w
    q1 = (kub_func.U(parametr, a, m, time_start, time_finish, time) - a * w) / m

    k2 = w + time_step/2
    q2 = (kub_func.U(parametr, a, m, time_start, time_finish, time + time_step/2) - a * (w + time_step/2 * q1)) / m

    k3 = w + time_step / 2
    q3 = (kub_func.U(parametr, a, m, time_start, time_finish, time+ time_step/2) - a * (w + time_step / 2 * q2)) / m

    k4 = w + time_step
    q4 = (kub_func.U(parametr, a, m, time_start, time_finish, time+ time_step) - a * (w + time_step * q3)) / m

    dq = q + time_step/6 * (k1 + 2*k2 + 2*k3 + k4)
    dw = w + time_step/6 * (q1 + 2*q2 + 2*q3 + q4)
    #dq = q + time_step * w
    #dw = w + time_step * ((kub_func.U(parametr, a, m, time_start, time_finish, time) - a * w) / m)
    return dq, dw


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


a = Data["info_BLA"][0]["a"]
m = Data["info_BLA"][0]["m"]
result_x = []
e = [Data["info_BLA"][0]["x"][0], Data["info_BLA"][0]["x"][1]]
result_x.append(e)
result_y = []
e_y = [Data["info_BLA"][0]["y"][0], Data["info_BLA"][0]["y"][1]]
result_y.append(e_y)
time_step = 0.1
time = time_start
gg = 0
while abs(time_finish - time) > 0.00000001:
    e = eler(result_x[len(result_x)-1][0], result_x[len(result_x)-1][1], time, time_step, Data["info_BLA"][0]["x"],a,m,time_start,time_finish,gg)
    e_y = eler(result_y[len(result_y)-1][0], result_y[len(result_y)-1][1], time, time_step, Data["info_BLA"][0]["y"],a,m,time_start,time_finish,gg)
    result_x.append(e)
    result_y.append(e_y)
    time = time + time_step

print(result_x)
print("------------")
print(result_y)

