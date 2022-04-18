import json
import math
import info_cache
import kub_func
import Runge_Kutta
from scipy import integrate
import matplotlib.pyplot as plt
import numpy as np

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

result = []
for id in range(len(info_cache.info_BLA)):
    result.append({"x": [], "y": [], "z": []})
    e = [info_cache.info_BLA[id]["x"]["0"], info_cache.info_BLA[id]["x"]["d0"]]
    result[id]["x"].append(e)
    e = [info_cache.info_BLA[id]["y"]["0"], info_cache.info_BLA[id]["y"]["d0"]]
    result[id]["y"].append(e)
    e = [info_cache.info_BLA[id]["z"]["0"], info_cache.info_BLA[id]["z"]["d0"]]
    result[id]["z"].append(e)


time_step = 0.01
time = info_cache.time_start
gg = 0
name_axis = ["x", "y", "z"]
while abs(info_cache.time_finish - time) > 0.00000001:
    for id in range(len(info_cache.info_BLA)):
        if (id==1):

            for axis in name_axis:
                if(axis=="x"):
                    e = Runge_Kutta.Runge_kutta1(result[id][axis][len(result[id][axis]) - 1][0],
                                    result[id][axis][len(result[id][axis]) - 1][1], time, time_step, id, axis, gg)
                    result[id][axis].append(e)
                else:
                    e = Runge_Kutta.Runge_kutta(result[id][axis][len(result[id][axis]) - 1][0],
                                                result[id][axis][len(result[id][axis]) - 1][1], time, time_step, id,
                                                axis, gg)
                    result[id][axis].append(e)
        else:

            for axis in name_axis:
                e = Runge_Kutta.Runge_kutta(result[id][axis][len(result[id][axis]) - 1][0],
                                     result[id][axis][len(result[id][axis]) - 1][1], time, time_step, id, axis, gg)
                result[id][axis].append(e)
    time = time + time_step



for id in range(len(info_cache.info_BLA)):

    full_res = []
    for i in range(len(result[id]["x"])):
        full_res.append(result[id]["x"][i][0])

    full_res1_x = full_res
    full_res = []
    for i in range(len(result[id]["y"])):
        full_res.append(result[id]["y"][i][0])

    full_res1_y = full_res
    plt.plot(full_res1_x, full_res1_y)




plt.show()
