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


def func(pred_val,
         pred_dval,
         time,
         time_step,
         id,
         axis,
         tau,
         number):
    
    val = pred_val + time_step * pred_dval
    dval = pred_dval + time_step * (
                (kub_func.U(id, axis, time) - info_cache.info_BLA[id]["a"] * pred_dval) / info_cache.info_BLA[id]["m"])
    return val, dval