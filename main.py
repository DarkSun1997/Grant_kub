import json
import math
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
    q = q + time_step * w
    w = w + time_step * ((kub_func.U(parametr, a, m, time_start, time_finish, time) - a * w) / m)
    return q, w



file_information_BLA = open('file_information_BLA.json')
Data = json.load(file_information_BLA)
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