import json
import math
import info_cache
import kub_func
from scipy import integrate
import matplotlib.pyplot as plt

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
    dq = q + time_step * w
    dw = w + time_step * ((kub_func.U(parametr, a, m, time_start, time_finish, time) - a * w) / m)
    return dq, dw


def transmission_of_information():
    file_information_BLA = open('file_information_BLA.json')
    Data = json.load(file_information_BLA)
    info_cache.info_BLA = Data["info_BLA"]
    info_cache.time_start = Data["time_start"]
    info_cache.time_finish = Data["time_finish"]
    info_cache.using_condition = Data["using_condition"]
    info_cache.d = Data["d"]
    generate_information_calc_BLA()


def generate_information_calc_BLA():
    for i in range(len(info_cache.info_BLA)):
        #print(info_cache.info_BLA[i])
        info_cache.info_BLA[i]["x"]["a21"] = None
        info_cache.info_BLA[i]["y"]["a21"] = None
        info_cache.info_BLA[i]["z"]["a21"] = None
        #print(info_cache.info_BLA[i])



transmission_of_information()
print(info_cache.info_BLA)
"""
generate_information_calc_BLA(Data)


print(kub_func.information_calc_BLA)
print(len(Data['info_BLA']))

print(Data['info_BLA'][0])

time_start = Data["time_start"]
time_finish = Data["time_finish"]
T = time_finish
using_condition = Data["using_condition"]
data_BLA = Data["info_BLA"]
print(data_BLA)


a = Data["info_BLA"][0]["a"]
m = Data["info_BLA"][0]["m"]
result_x = []
e = [Data["info_BLA"][0]["x"][0], Data["info_BLA"][0]["x"][1]]
result_x.append(e)
print(e)
time_step = 0.1
time = time_start
gg = 0
while abs(time_finish - time) > 0.00000001:
    e = eler(result_x[len(result_x)-1][0], result_x[len(result_x)-1][1], time, time_step, Data["info_BLA"][0]["x"],a,m,time_start,time_finish,gg)
    result_x.append(e)
    time = time + time_step

print(result_x)
print(kub_func.fi1_L2_const)
"""