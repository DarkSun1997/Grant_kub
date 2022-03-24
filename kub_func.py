import math
import info_cache
from scipy import integrate



def fi1(parametr,
        a,
        m,
        time_start,
        time_finish,
        tau):
    return math.exp((-a) / m * (time_finish - tau)) / m


def fi2(parametr,
        a,
        m,
        time_start,
        time_finish,
        tau):
    return (1 - math.exp((-a) / m * (time_finish - tau))) / a


def g1(parametr,
        a,
        m,
        time_start,
        time_finish,
        tau):
    return parametr[3] - parametr[1] * math.exp((-a) / m * time_finish)


def g2(parametr,
        a,
        m,
        time_start,
        time_finish,
        tau):
    return parametr[2] - parametr[0] - parametr[1] * (1 - math.exp((-a) / m * time_finish)) * m / a


def normal_function_L2(func,
                        parametr,
                        a,
                        m,
                        time_start,
                        time_finish,
                        tau):
    func_integrate = lambda tau_integrate: func(parametr, a, m, time_start, time_finish, tau_integrate) ** 2
    return math.sqrt(integrate.quad(func_integrate, time_start, time_finish)[0])


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