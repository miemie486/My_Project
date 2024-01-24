# 导入随机数模块
import random
import time
import numpy as np


# 定义参数
POP_SIZE = 10000 # 种群大小
GENES = [-1, 1] # 基因范围
CROSS_RATE = 0.5 # 交叉概率
MUTATE_RATE = 0.1 # 变异概率
MAX_GEN = 100 # 最大代数

# 定义目标函数
def f(x, y):
    return x**2 + y**2

# 定义个体类
class Individual:
    def __init__(self, x, y):
        self.x = x # x坐标
        self.y = y # y坐标
        self.fitness = f(x, y) # 适应度

    def __str__(self):
        return f"({self.x}, {self.y})"

# 初始化种群
def init_pop():
    pop = []
    for i in range(POP_SIZE):
        random.seed(time.perf_counter_ns())
        x = random.uniform(GENES[0], GENES[1]) # 随机生成x
        random.seed(time.perf_counter_ns())
        y = random.uniform(GENES[0], GENES[1]) # 随机生成y
        ind = Individual(x, y) # 创建个体
        pop.append(ind) # 加入种群
    return pop

# 选择操作
def selection(pop):
    # 按适应度从小到大排序
    pop.sort(key=lambda ind: ind.fitness)
    # 选择前一半的个体
    return pop[:POP_SIZE//2]

# 交叉操作
def crossover(pop):
    new_pop = []
    for i in range(POP_SIZE//2):
        # 随机选择两个父代
        random.seed(time.perf_counter_ns())
        p1 = random.choice(pop)
        random.seed(time.perf_counter_ns())
        p2 = random.choice(pop)
        # 以一定概率进行交叉
        random.seed(time.perf_counter_ns())
        if random.random() < CROSS_RATE:
            # 随机选择一个交叉点
            random.seed(time.perf_counter_ns())
            cross_point = random.randint(0, 1)
            # 交换基因
            if cross_point == 0:
                c1 = Individual(p1.x, p2.y)
                c2 = Individual(p2.x, p1.y)
            else:
                c1 = Individual(p2.x, p1.y)
                c2 = Individual(p1.x, p2.y)
        else:
            # 不交叉则保持不变
            c1 = Individual(p1.x, p1.y)
            c2 = Individual(p2.x, p2.y)
        # 加入新种群
        new_pop.append(c1)
        new_pop.append(c2)
    return new_pop

# 变异操作
def mutate(pop):
    for ind in pop:
        # 以一定概率进行变异
        random.seed(time.perf_counter_ns())
        if random.random() < MUTATE_RATE:
            # 随机选择一个变异点
            random.seed(time.perf_counter_ns())
            mutate_point = random.randint(0, 1)
            # 随机生成一个新基因
            random.seed(time.perf_counter_ns())
            new_gene = random.uniform(GENES[0], GENES[1])
            # 替换基因
            if mutate_point == 0:
                ind.x = new_gene
            else:
                ind.y = new_gene
            # 更新适应度
            ind.fitness = f(ind.x, ind.y)

# 主函数
def main():
    # 初始化种群
    pop = init_pop()
    # 进化指定代数
    for gen in range(MAX_GEN):
        # 选择
        pop = selection(pop)
        # 交叉
        pop = crossover(pop)
        # 变异
        mutate(pop)
        # 打印当前最优解
        print(f"Generation {gen+1}: {pop[0]}")
    # 打印最终结果
    print(f"Final solution: {pop[0]}")

# 调用主函数
if __name__ == "__main__":
    s=time.perf_counter()
    main()
    f=time.perf_counter()
    print(f-s)
