import math
import info_cache
from scipy import integrate
# функции которые дал Евгений Павлович

#fi1
def fi1(id,
        axis,
        tau):
    return math.exp((-info_cache.info_BLA[id]["a"]) / info_cache.info_BLA[id]["m"] *
                    (info_cache.time_finish - tau)) / info_cache.info_BLA[id]["m"]


#fi2
def fi2(id,
        axis,
        tau):
    return (1 - math.exp((-info_cache.info_BLA[id]["a"]) / info_cache.info_BLA[id]["m"] *
                         (info_cache.time_finish - tau))) / info_cache.info_BLA[id]["a"]

#g1
def g1(id,
        axis,
        tau):
    if info_cache.info_BLA[id][axis]["g1"] == None :
        info_cache.info_BLA[id][axis]["g1"] = info_cache.info_BLA[id][axis]["dT"] - info_cache.info_BLA[id][axis]["d0"] * \
           math.exp((-info_cache.info_BLA[id]["a"]) / info_cache.info_BLA[id]["m"] * info_cache.time_finish)
    return info_cache.info_BLA[id][axis]["g1"]

#g2
def g2(id,
        axis,
        tau):
    if info_cache.info_BLA[id][axis]["g2"] == None :
        info_cache.info_BLA[id][axis]["g2"] = info_cache.info_BLA[id][axis]["T"] - info_cache.info_BLA[id][axis]["0"] \
            - info_cache.info_BLA[id][axis]["d0"] * (1 - math.exp((-info_cache.info_BLA[id]["a"]) /
            info_cache.info_BLA[id]["m"] * info_cache.time_finish)) * info_cache.info_BLA[id]["m"] / info_cache.info_BLA[id]["a"]
    return info_cache.info_BLA[id][axis]["g2"]

#функции нормализации
def normal_function_L2(func,
                       id,
                       axis,
                       tau):
    func_integrate = lambda tau_integrate: func(id, axis, tau_integrate) ** 2
    return math.sqrt(integrate.quad(func_integrate, info_cache.time_start, info_cache.time_finish)[0])


def psi1(id,
         axis,
         tau):
    if info_cache.info_BLA[id][axis]["fi1_L2_const"] == None:
        info_cache.info_BLA[id][axis]["fi1_L2_const"] = normal_function_L2(fi1, id, axis, tau)
    return fi1(id, axis, tau) / info_cache.info_BLA[id][axis]["fi1_L2_const"]

#a21
def a21(id,
         axis,
         tau):
    if info_cache.info_BLA[id][axis]["a21"] == None:
        func_integrate = lambda tau_integrate: - fi2(id, axis, tau=tau_integrate) * \
                                           psi1(id, axis, tau=tau_integrate)
        info_cache.info_BLA[id][axis]["a21"] = integrate.quad(func_integrate, info_cache.time_start, info_cache.time_finish)[0]
    return info_cache.info_BLA[id][axis]["a21"]

#psi2_p
def psi2_p(id,
            axis,
            tau):
    return a21(id, axis, tau) * psi1(id, axis, tau) + fi2(id, axis, tau)

#psi2
def psi2(id,
            axis,
            tau):
    if info_cache.info_BLA[id][axis]["psi2_p_L2_const"] == None:
        info_cache.info_BLA[id][axis]["psi2_p_L2_const"] = normal_function_L2(psi2_p, id, axis, tau)
    return psi2_p(id, axis, tau) / info_cache.info_BLA[id][axis]["psi2_p_L2_const"]

#b11
def b11(id,
            axis,
            tau):
    if info_cache.info_BLA[id][axis]["b11"] == None:
        info_cache.info_BLA[id][axis]["b11"] = 1 / normal_function_L2(fi1, id, axis, tau)
    return info_cache.info_BLA[id][axis]["b11"]

#b22
def b22(id,
            axis,
            tau):

    if info_cache.info_BLA[id][axis]["b22"] == None:
        info_cache.info_BLA[id][axis]["b22"] = 1 / normal_function_L2(psi2_p, id, axis, tau)
    return info_cache.info_BLA[id][axis]["b22"]

#b21
def b21(id,
            axis,
            tau):
    if info_cache.info_BLA[id][axis]["b21"] == None:
        info_cache.info_BLA[id][axis]["b21"] = a21(id, axis, tau) / \
                (normal_function_L2(psi2_p, id, axis, tau) * normal_function_L2(fi1, id, axis, tau))
    return info_cache.info_BLA[id][axis]["b21"]

#B1
def B1(id,
       axis,
       tau):
    return b11(id, axis, tau) * g1(id, axis, tau)

#B2
def B2(id,
       axis,
       tau):
    return b21(id, axis, tau) * g1(id, axis, tau) + b22(id, axis, tau) * g2(id, axis, tau)

#U
def U(id,
       axis,
       tau):
    return B1(id, axis, tau) * psi1(id, axis, tau) + B2(id, axis, tau) * psi2(id, axis, tau)


def U1(id,
       axis,
       tau):
    return B1(id, axis, tau) * psi1(id, axis, tau) + B2(id, axis, tau) * psi2(id, axis, tau) + -30 * psi3(id, axis, tau)


def print_info():
    print(info_cache.info_BLA)













#говно код чтобы разобраться что делать

def fi3(id,
        axis,
        tau):
    return tau ** 2

def psi3_p(id,
            axis,
            tau):
    return a31(id, axis, tau) * psi1(id, axis, tau) + a32(id, axis, tau) * psi2(id, axis, tau) + fi3(id, axis, tau)

def psi3(id,
            axis,
            tau):
    return psi3_p(id, axis, tau) / normal_function_L2(psi3_p, id, axis, tau)

def a31(id,
         axis,
         tau):
    func_integrate = lambda tau_integrate: - fi3(id, axis, tau=tau_integrate) * \
                                           psi1(id, axis, tau=tau_integrate)
    z = integrate.quad(func_integrate, info_cache.time_start, info_cache.time_finish)[0]
    return z


def a32(id,
         axis,
         tau):
    func_integrate = lambda tau_integrate: - fi3(id, axis, tau=tau_integrate) * \
                                           psi2(id, axis, tau=tau_integrate)
    z = integrate.quad(func_integrate, info_cache.time_start, info_cache.time_finish)[0]
    return z


def F(s,id,
         axis,
         tau, t):

    dg = s[1]
    dr = (psi3_p(id, axis, tau)-info_cache.info_BLA[id]["a"]*s[1])/info_cache.info_BLA[id]["m"]
    return [dg, dr]


def Det_c3():
    return