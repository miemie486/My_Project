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

# 初始化种群
def init_pop():
    pop = []
    for i in range(POP_SIZE):
        random.seed(time.perf_counter_ns())
        x = random.uniform(GENES[0], GENES[1]) # 随机生成x
        random.seed(time.perf_counter_ns())
        y = random.uniform(GENES[0], GENES[1]) # 随机生成y
        fitness = f(x, y) # 计算适应度
        pop.append((x, y, fitness)) # 加入种群，用元组表示个体
    return pop

# 选择操作
def selection(pop):
    # 按适应度从小到大排序
    pop.sort(key=lambda ind: ind[2])
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
                c1 = (p1[0], p2[1], f(p1[0], p2[1]))
                c2 = (p2[0], p1[1], f(p2[0], p1[1]))
            else:
                c1 = (p2[0], p1[1], f(p2[0], p1[1]))
                c2 = (p1[0], p2[1], f(p1[0], p2[1]))
        else:
            # 不交叉则保持不变
            c1 = (p1[0], p1[1], p1[2])
            c2 = (p2[0], p2[1], p2[2])
        # 加入新种群
        new_pop.append(c1)
        new_pop.append(c2)
    return new_pop

# 变异操作
def mutate(pop):
    for i in range(len(pop)):
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
                pop[i] = (new_gene, pop[i][1], f(new_gene, pop[i][1]))
            else:
                pop[i] = (pop[i][0], new_gene, f(pop[i][0], new_gene))
    return pop

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
        pop = mutate(pop)
        # 打印当前最优解
        print(f"Generation {gen+1}: ({pop[0][0]}, {pop[0][1]})")
    # 打印最终结果
    print(f"Final solution: ({pop[0][0]}, {pop[0][1]})")

# 调用主函数
if __name__ == "__main__":
    s=time.perf_counter()
    main()
    f=time.perf_counter()
    print(f-s)
