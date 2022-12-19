import json
import math
import info_cache
import kub_func
import Runge_Kutta
import Calc_function_g
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
    info_cache.axes_correction = Data["axes_correction"]
    generate_information_calc_BLA()


def generate_information_calc_BLA():
    for id in range(len(info_cache.info_BLA)):
        name_axis = ["x", "y", "z"]
        for axis in name_axis:
            info_cache.info_BLA[id][axis]["g1"] = None
            info_cache.info_BLA[id][axis]["g2"] = None
            info_cache.info_BLA[id][axis]["fi1_L2_const"] = None
            info_cache.info_BLA[id][axis]["psi2_p_L2_const"] = None


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
epsilon = 0.00000001
name_axis = ["x", "y", "z"]
for id in range(len(info_cache.info_BLA)):
    time = info_cache.time_start
    while abs(info_cache.time_finish - time) > epsilon:
        if info_cache.using_condition == False:
            for axis in name_axis:
                e = Runge_Kutta.Runge_kutta(result[id][axis][len(result[id][axis]) - 1][0],
                                            result[id][axis][len(result[id][axis]) - 1][1], time, time_step, id, axis,
                                            gg)
                result[id][axis].append(e)
        else:
            if (id > 0):
                for axis in name_axis:
                    if axis == info_cache.axes_correction:
                        e = Runge_Kutta.Runge_kutta1(result[id][axis][len(result[id][axis]) - 1][0],
                                                    result[id][axis][len(result[id][axis]) - 1][1], time, time_step, id,
                                                     axis, gg)
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


#Отрисовка результата
for id in range(len(info_cache.info_BLA)):
    full_res = []
    for i in range(len(result[id]["x"])):
        full_res.append(result[id]["x"][i][0])
    full_result_x = full_res

    full_res = []
    for i in range(len(result[id]["y"])):
        full_res.append(result[id]["y"][i][0])
    full_result_y = full_res

    plt.plot(full_result_x, full_result_y)

#plt.show()


#Некоторый мусор для просмотра промежуточных результатов
result_g = []
result_g.append([0, 0])
time = info_cache.time_start
print(time_step)
while abs(info_cache.time_finish - time) > epsilon:
    e = Calc_function_g.Runge_kutta(result_g[len(result_g) - 1][0],
                                                    result_g[len(result_g) - 1][1], time, time_step, 0,
                                                     "x", 0, kub_func.psi3)
    result_g.append(e)
    time = time + time_step
print(result_g)
