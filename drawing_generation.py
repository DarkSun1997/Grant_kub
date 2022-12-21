import math
import info_cache
import matplotlib.pyplot as plt


def generate_ris(result):
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

    if info_cache.using_condition == True:
        for number in range(len(info_cache.info_objects)):
            if info_cache.info_objects[number]['cylinder']:
                for id in range(len(info_cache.info_objects[number]['cylinder'])):
                    mas_point_x = []
                    mas_point_y = []
                    for ange in range(370):
                        mas_point_x.append(info_cache.info_objects[number]['cylinder'][id]["x"] +
                                           info_cache.info_objects[number]['cylinder'][id]["radius"] * math.cos(ange/180*3.1415))
                        mas_point_y.append(info_cache.info_objects[number]['cylinder'][id]["y"] +
                                           info_cache.info_objects[number]['cylinder'][id]["radius"] * math.sin(ange/180*3.1415))
                    plt.plot(mas_point_x, mas_point_y)

    plt.show()