import json
import math
import kub_func
from scipy import integrate





file_information_BLA = open('file_information_BLA.json')
Data = json.load(file_information_BLA)
print(Data['info_BLA'][0])

time_start = Data["time_start"]
time_finish = Data["time_finish"]
T = time_finish
using_condition = Data["using_condition"]
data_BLA = Data["info_BLA"]
ggg = kub_func.normal_function_L2(kub_func.z, 5.0, 10.0, -2.0, 2.0,1.0)
print(ggg)
print(data_BLA)
