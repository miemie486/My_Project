import numpy as np
a = np.array([1,2,3,4,5,6,7,8,9])
indices = np.array([0, 8, 4]) # 这是一个坐标列表，表示(0,0)，(1,1)，(2,2)三个位置的元素
values = np.array([10, 20, 30]) # 这是一个值列表，表示要替换的值
np.put(a, indices, values) # 根据坐标列表和值列表替换元素
print(a) # 输出结果
# [[10  2  3]
#  [ 4 20  6]
#  [ 7  8 30]]
