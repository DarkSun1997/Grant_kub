import math
import info_cache
from scipy import integrate



def fi1(id,
        axis,
        tau):
    return math.exp((-info_cache.info_BLA[id]["a"]) / info_cache.info_BLA[id]["m"] *
                    (info_cache.time_finish - tau)) / info_cache.info_BLA[id]["m"]


def fi2(id,
        axis,
        tau):
    return (1 - math.exp((-info_cache.info_BLA[id]["a"]) / info_cache.info_BLA[id]["m"] *
                         (info_cache.time_finish - tau))) / info_cache.info_BLA[id]["a"]


def g1(id,
        axis,
        tau):
    if info_cache.info_BLA[id][axis]["g1"] == None :
        info_cache.info_BLA[id][axis]["g1"] = info_cache.info_BLA[id][axis]["dT"] - info_cache.info_BLA[id][axis]["d0"] * \
           math.exp((-info_cache.info_BLA[id]["a"]) / info_cache.info_BLA[id]["m"] * info_cache.time_finish)
    return info_cache.info_BLA[id][axis]["g1"]


def g2(id,
        axis,
        tau):
    if info_cache.info_BLA[id][axis]["g2"] == None :
        info_cache.info_BLA[id][axis]["g2"] = info_cache.info_BLA[id][axis]["T"] - info_cache.info_BLA[id][axis]["0"] \
            - info_cache.info_BLA[id][axis]["d0"] * (1 - math.exp((-info_cache.info_BLA[id]["a"]) /
            info_cache.info_BLA[id]["m"] * info_cache.time_finish)) * info_cache.info_BLA[id]["m"] / info_cache.info_BLA[id]["a"]
    return info_cache.info_BLA[id][axis]["g2"]


def normal_function_L2(func,
                       id,
                       axis,
                       tau):
    func_integrate = lambda tau_integrate: func(id, axis, tau_integrate) ** 2
    return math.sqrt(integrate.quad(func_integrate, info_cache.time_start, info_cache.time_finish)[0])


def psi1(parametr,
            a,
            m,
            time_start,
            time_finish,
            tau):
    global fi1_L2_const
    if fi1_L2_const == None:
        fi1_L2_const = normal_function_L2(fi1, parametr, a, m, time_start, time_finish, tau)
    return fi1(parametr, a, m, time_start, time_finish, tau) / fi1_L2_const



def a21(parametr,
            a,
            m,
            time_start,
            time_finish,
            tau):

    func_integrate = lambda tau_integrate: - fi2(parametr, a, m, time_start, time_finish, tau=tau_integrate) * \
                                           psi1(parametr, a, m, time_start, time_finish, tau=tau_integrate)
    return integrate.quad(func_integrate, time_start, time_finish)[0]


def psi2_p(parametr,
            a,
            m,
            time_start,
            time_finish,
            tau):
    return a21(parametr, a, m, time_start, time_finish, tau) * \
           psi1(parametr, a, m, time_start, time_finish, tau) + \
           fi2(parametr, a, m, time_start, time_finish, tau)


def psi2(parametr,
            a,
            m,
            time_start,
            time_finish,
            tau):
    return psi2_p(parametr, a, m, time_start, time_finish, tau) / \
           normal_function_L2(psi2_p, parametr, a, m, time_start, time_finish, tau)


def b11(parametr,
            a,
            m,
            time_start,
            time_finish,
            tau):
    return 1 / normal_function_L2(fi1, parametr, a, m, time_start, time_finish, tau)


def b22(parametr,
            a,
            m,
            time_start,
            time_finish,
            tau):
    return 1 / normal_function_L2(psi2_p, parametr, a, m, time_start, time_finish, tau)


def b21(parametr,
            a,
            m,
            time_start,
            time_finish,
            tau):
    return a21(parametr, a, m, time_start, time_finish, tau) / \
           (normal_function_L2(psi2_p, parametr, a, m, time_start, time_finish, tau) * \
            normal_function_L2(fi1, parametr, a, m, time_start, time_finish, tau))


def B1(parametr,
        a,
        m,
        time_start,
        time_finish,
        tau):
    return b11(parametr, a, m, time_start, time_finish, tau) * \
           g1(parametr, a, m, time_start, time_finish, tau)


def B2(parametr,
       a,
       m,
       time_start,
       time_finish,
       tau):
    return b21(parametr, a, m, time_start, time_finish, tau) * \
           g1(parametr, a, m, time_start, time_finish, tau) + \
           b22(parametr, a, m, time_start, time_finish, tau) * \
           g2(parametr, a, m, time_start, time_finish, tau)


def U(parametr,
       a,
       m,
       time_start,
       time_finish,
       tau):
    return B1(parametr, a, m, time_start, time_finish, tau) * \
           psi1(parametr, a, m, time_start, time_finish, tau) + \
           B2(parametr, a, m, time_start, time_finish, tau) * \
           psi2(parametr, a, m, time_start, time_finish, tau)


def print_info():
    print(info_cache.info_BLA)