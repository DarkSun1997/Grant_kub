import json
import math
import info_cache
import kub_func
import Runge_Kutta
from scipy import integrate
import matplotlib.pyplot as plt
import numpy as np


#функции нормализации
def normal_function_L2_2func(func1,
                             func2,
                       id,
                       axis,
                       tau):
    func_integrate = lambda tau_integrate: func1(id, axis, tau_integrate) * func2(id, axis, tau_integrate)
    return -(integrate.quad(func_integrate, info_cache.time_start, info_cache.time_finish)[0])




def step_func_g(q,
         w,
         time,
         time_step,
         id,
         axis,
         tau,
         ksi):
    dq = w
    dw = (ksi(id, axis, time) - info_cache.info_BLA[id]["a"] * w) / info_cache.info_BLA[id]["m"]
    return dq, dw


def Runge_kutta(q,
         w,
         time,
         time_step,
         id,
         axis,
         tau,
         ksi):
    k1 = step_func_g(q, w, time, time_step, id, axis, tau, ksi)
    k2 = step_func_g(q + time_step * k1[0] / 2, w + time_step * k1[1] / 2, time, time_step, id, axis, tau, ksi)
    k3 = step_func_g(q + time_step * k2[0] / 2, w + time_step * k2[1] / 2, time, time_step, id, axis, tau, ksi)
    k4 = step_func_g(q + time_step * k3[0], w + time_step * k3[1], time, time_step, id, axis, tau, ksi)

    dq = q + time_step / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
    dw = w + time_step / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])

    return dq, dw