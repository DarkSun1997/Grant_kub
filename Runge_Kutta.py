import info_cache
import kub_func


def eler(q,
         w,
         time,
         time_step,
         id,
         axis,
         tau):
    dq = q + time_step * w
    dw = w + time_step * ((kub_func.U(id, axis, time) - info_cache.info_BLA[id]["a"] * w) / info_cache.info_BLA[id]["m"])
    return dq, dw


def func(q,
         w,
         time,
         time_step,
         id,
         axis,
         tau):
    dq = w
    dw = (kub_func.U(id, axis, time) - info_cache.info_BLA[id]["a"] * w) / info_cache.info_BLA[id]["m"]
    return dq, dw


def Runge_kutta(q,
         w,
         time,
         time_step,
         id,
         axis,
         tau):
    k1 = func(q, w, time, time_step, id, axis, tau)
    k2 = func(q + time_step * k1[0] / 2, w + time_step * k1[1] / 2, time, time_step, id, axis, tau)
    k3 = func(q + time_step * k2[0] / 2, w + time_step * k2[1] / 2, time, time_step, id, axis, tau)
    k4 = func(q + time_step * k3[0], w + time_step * k3[1], time, time_step, id, axis, tau)

    new_q = q + time_step / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
    new_w = w + time_step / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])

    return new_q, new_w


def func1(q,
         w,
         time,
         time_step,
         id,
         axis,
         tau):
    dq = w
    dw = (kub_func.U1(id, axis, time) - info_cache.info_BLA[id]["a"] * w) / info_cache.info_BLA[id]["m"]
    return dq, dw


def Runge_kutta1(q,
         w,
         time,
         time_step,
         id,
         axis,
         tau):
    k1 = func1(q, w, time, time_step, id, axis, tau)
    k2 = func1(q + time_step * k1[0] / 2, w + time_step * k1[1] / 2, time, time_step, id, axis, tau)
    k3 = func1(q + time_step * k2[0] / 2, w + time_step * k2[1] / 2, time, time_step, id, axis, tau)
    k4 = func1(q + time_step * k3[0], w + time_step * k3[1], time, time_step, id, axis, tau)

    new_q = q + time_step / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
    new_w = w + time_step / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])

    return new_q, new_w