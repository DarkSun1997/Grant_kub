import json
import math
import info_cache
import kub_func
from scipy import integrate
import matplotlib.pyplot as plt
import numpy as np

def eler(q,
         w,
         time,
         time_step,
         id,
         axis,
         tau):
    dq = q + time_step * w
    dw = w + time_step * ((kub_func.U(id, axis, time) - info_cache.info_BLA[id]["a"] * w) / info_cache.info_BLA[id]["m"] )
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
        info_cache.info_BLA[i]["x"]["a21"] = None
        info_cache.info_BLA[i]["y"]["a21"] = None
        info_cache.info_BLA[i]["z"]["a21"] = None
        info_cache.info_BLA[i]["x"]["g1"] = None
        info_cache.info_BLA[i]["y"]["g1"] = None
        info_cache.info_BLA[i]["z"]["g1"] = None
        info_cache.info_BLA[i]["x"]["g2"] = None
        info_cache.info_BLA[i]["y"]["g2"] = None
        info_cache.info_BLA[i]["z"]["g2"] = None
        info_cache.info_BLA[i]["x"]["fi1_L2_const"] = None
        info_cache.info_BLA[i]["y"]["fi1_L2_const"] = None
        info_cache.info_BLA[i]["z"]["fi1_L2_const"] = None
        info_cache.info_BLA[i]["x"]["a21"] = None
        info_cache.info_BLA[i]["y"]["a21"] = None
        info_cache.info_BLA[i]["z"]["a21"] = None
        info_cache.info_BLA[i]["x"]["psi2_p_L2_const"] = None
        info_cache.info_BLA[i]["y"]["psi2_p_L2_const"] = None
        info_cache.info_BLA[i]["z"]["psi2_p_L2_const"] = None
        info_cache.info_BLA[i]["x"]["b11"] = None
        info_cache.info_BLA[i]["y"]["b11"] = None
        info_cache.info_BLA[i]["z"]["b11"] = None
        info_cache.info_BLA[i]["x"]["b22"] = None
        info_cache.info_BLA[i]["y"]["b22"] = None
        info_cache.info_BLA[i]["z"]["b22"] = None
        info_cache.info_BLA[i]["x"]["b21"] = None
        info_cache.info_BLA[i]["y"]["b21"] = None
        info_cache.info_BLA[i]["z"]["b21"] = None



transmission_of_information()
print(info_cache.info_BLA)


result_x = []
time_step = 0.001
time = info_cache.time_start
gg = 0
e = [info_cache.info_BLA[0]["x"]["0"],info_cache.info_BLA[0]["x"]["d0"]]
result_x.append(e)
print(e)
while abs(info_cache.time_finish - time) > 0.00000001:
    e = eler(result_x[len(result_x)-1][0], result_x[len(result_x)-1][1], time, time_step, 0, "x", gg)
    result_x.append(e)
    time = time + time_step

print(result_x)


result_y = []
time_step = 0.001
time = info_cache.time_start
gg = 0
e = [info_cache.info_BLA[0]["y"]["0"],info_cache.info_BLA[0]["y"]["d0"]]
result_y.append(e)
print(e)
while abs(info_cache.time_finish - time) > 0.00000001:
    e = eler(result_y[len(result_y)-1][0], result_y[len(result_y)-1][1], time, time_step, 0, "y", gg)
    result_y.append(e)
    time = time + time_step

print(result_y)

full_res = []
for i in range(len(result_x)):
    full_res.append(result_x[i][0])

full_res1_x = full_res
full_res = []
for i in range(len(result_y)):
    full_res.append(result_y[i][0])

full_res1_y = full_res
plt.plot(full_res1_x,full_res1_y)




result_x = []
time_step = 0.001
time = info_cache.time_start
gg = 0
e = [info_cache.info_BLA[1]["x"]["0"],info_cache.info_BLA[0]["x"]["d0"]]
result_x.append(e)
print(e)
while abs(info_cache.time_finish - time) > 0.00000001:
    e = eler(result_x[len(result_x)-1][0], result_x[len(result_x)-1][1], time, time_step, 1, "x", gg)
    result_x.append(e)
    time = time + time_step

print(result_x)


result_y = []
time_step = 0.001
time = info_cache.time_start
gg = 0
e = [info_cache.info_BLA[1]["y"]["0"],info_cache.info_BLA[0]["y"]["d0"]]
result_y.append(e)
print(e)
while abs(info_cache.time_finish - time) > 0.00000001:
    e = eler(result_y[len(result_y)-1][0], result_y[len(result_y)-1][1], time, time_step, 1, "y", gg)
    result_y.append(e)
    time = time + time_step

print(result_y)

full_res = []
for i in range(len(result_x)):
    full_res.append(result_x[i][0])

full_res1_x = full_res
full_res = []
for i in range(len(result_y)):
    full_res.append(result_y[i][0])

full_res1_y = full_res
plt.plot(full_res1_x,full_res1_y)


plt.show()
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