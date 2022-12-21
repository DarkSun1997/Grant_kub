import json
import math
import info_cache
import kub_func
import Runge_Kutta
import drawing_generation


def transmission_of_information():
    file_information_BLA = open('file_information_BLA.json')
    Data = json.load(file_information_BLA)
    info_cache.info_BLA = Data["info_BLA"]
    info_cache.time_start = Data["time_start"]
    info_cache.time_finish = Data["time_finish"]
    info_cache.using_condition = Data["using_condition"]
    info_cache.d = Data["d"]
    info_cache.axes_correction = Data["axes_correction"]
    info_cache.info_objects = Data["info_objects"]

transmission_of_information()

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
    info_cache.func_g.append([])
    time = info_cache.time_start
    while abs(info_cache.time_finish - time) > epsilon:
        for axis in name_axis:
            e = Runge_Kutta.Runge_kutta(result[id][axis][len(result[id][axis]) - 1][0],
                                        result[id][axis][len(result[id][axis]) - 1][1], time, time_step, id, axis,
                                        gg)
            result[id][axis].append(e)
        time = time + time_step
    if info_cache.using_condition == True:
        flag = True
        for id_elem in range(id):
            for h in range(len(result[id]['x'])):
                dist = math.sqrt((result[id]['x'][h][0] - result[id_elem]['x'][h][0]) ** 2 +
                                 (result[id]['y'][h][0] - result[id_elem]['y'][h][0]) ** 2 +
                                 (result[id]['z'][h][0] - result[id_elem]['z'][h][0]) ** 2)
                if(dist < info_cache.d):
                    flag = False

        for number in range(len(info_cache.info_objects)):
            for h in range(len(result[id]['x'])):
                for id_ob in range(len(info_cache.info_objects[number]['cylinder'])):
                    dist = math.sqrt((result[id]['x'][h][0] - info_cache.info_objects[number]['cylinder'][id_ob]["x"]) ** 2 +
                                     (result[id]['y'][h][0] - info_cache.info_objects[number]['cylinder'][id_ob]["y"]) ** 2)
                    if(dist - info_cache.info_objects[number]['cylinder'][id_ob]["radius"] < info_cache.d):
                        flag = False

        if flag == False:
            info_cache.func_g[id] = []
            time = info_cache.time_start
            result_g = []
            result_g.append((0, 0))
            for index in range(3):
                while abs(info_cache.time_finish - time) > epsilon:
                    e = Runge_Kutta.Runge_kutta_func_g(result_g[len(result_g) - 1][0],
                                                       result_g[len(result_g) - 1][1], time, time_step, 0,
                                                       "x", 0, kub_func.psi3)
                    result_g.append(e)
                    time = time + time_step
            info_cache.func_g[id].append(result_g)

            time = info_cache.time_start
            result_g = []
            result_g.append((0, 0))
            for index in range(3):
                while abs(info_cache.time_finish - time) > epsilon:
                    e = Runge_Kutta.Runge_kutta_func_g(result_g[len(result_g) - 1][0],
                                                       result_g[len(result_g) - 1][1], time, time_step, 0,
                                                       "x", 0, kub_func.psi4)
                    result_g.append(e)
                    time = time + time_step
            info_cache.func_g[id].append(result_g)

            result_g = []
            time = info_cache.time_start
            result_g.append((0, 0))
            for index in range(3):
                while abs(info_cache.time_finish - time) > epsilon:
                    e = Runge_Kutta.Runge_kutta_func_g(result_g[len(result_g) - 1][0],
                                                       result_g[len(result_g) - 1][1], time, time_step, 0,
                                                       "x", 0, kub_func.psi5)
                    result_g.append(e)
                    time = time + time_step
            info_cache.func_g[id].append(result_g)

            #решение с коррекцией
            time = info_cache.time_start
            result[id][info_cache.axes_correction] = []
            e = [info_cache.info_BLA[id][info_cache.axes_correction]["0"], info_cache.info_BLA[id][info_cache.axes_correction]["d0"]]
            result[id][info_cache.axes_correction].append(e)
            while abs(info_cache.time_finish - time) > epsilon:

                e = Runge_Kutta.Runge_kutta1(result[id][info_cache.axes_correction][len(result[id][info_cache.axes_correction]) - 1][0],
                                                result[id][info_cache.axes_correction][len(result[id][info_cache.axes_correction]) - 1][1],
                                                time, time_step, id,
                                                info_cache.axes_correction,
                                                gg)
                result[id][info_cache.axes_correction].append(e)
                time = time + time_step

#Отрисовка результата
drawing_generation.generate_ris(result)


