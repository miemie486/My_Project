# 导入NumPy模块
import numpy as np
import time
import numba as nb
# 定义参数

POP_SIZE = 10000 # 种群大小
GENES = [-1, 1] # 基因范围
CROSS_RATE = 0.5 # 交叉概率
MUTATE_RATE = 0.1 # 变异概率
MAX_GEN = 100 # 最大代数
DEMENSION = 10

# 定义目标函数
# Genes_arrary
@nb.jit
def f(Genes_arrary):
    fitness=0
    for i in Genes_arrary:
        fitness=fitness+i**2
    return fitness

# 初始化种群
def init_pop():
    Genes_arrary=np.zeros(DEMENSION)
    # 使用numpy.random.uniform来生成随机数
    for i in range(DEMENSION):
        Genes_arrary[i]=np.random.uniform(GENES[0], GENES[1], POP_SIZE)
    fitness = f(Genes_arrary) # 计算适应度
    # 使用numpy.array来表示个体和种群
    pop = np.array([x, y, fitness]).T # 转置为种群矩阵
    return pop

# 选择操作
def selection(pop):
    # 按适应度从小到大排序
    pop = pop[pop[:, 2].argsort()]
    # 选择前一半的个体
    return pop[:POP_SIZE//2]

# 交叉操作
def crossover(pop):
    new_pop = []
    for i in range(POP_SIZE//2):
        # 随机选择两个父代
        p1 = pop[np.random.randint(POP_SIZE//2)]
        p2 = pop[np.random.randint(POP_SIZE//2)]
        # 以一定概率进行交叉
        if np.random.random() < CROSS_RATE:
            # 随机选择一个交叉点
            cross_point = np.random.randint(0, 1)
            # 交换基因
            if cross_point == 0:
                c1 = np.array([p1[0], p2[1], f(p1[0], p2[1])])
                c2 = np.array([p2[0], p1[1], f(p2[0], p1[1])])
            else:
                c1 = np.array([p2[0], p1[1], f(p2[0], p1[1])])
                c2 = np.array([p1[0], p2[1], f(p1[0], p2[1])])
        else:
            # 不交叉则保持不变
            c1 = np.array([p1[0], p1[1], p1[2]])
            c2 = np.array([p2[0], p2[1], p2[2]])
        # 加入新种群
        new_pop.append(c1)
        new_pop.append(c2)
    # 转换为numpy.array
    new_pop = np.array(new_pop)
    return new_pop

# 变异操作
def mutate(pop):
    for i in range(len(pop)):
        # 以一定概率进行变异
        if np.random.random() < MUTATE_RATE:
            # 随机选择一个变异点
            mutate_point = np.random.randint(0, 1)
            # 随机生成一个新基因
            new_gene = np.random.uniform(GENES[0], GENES[1])
            # 替换基因
            if mutate_point == 0:
                pop[i] = np.array([new_gene, pop[i][1], f(new_gene, pop[i][1])])
            else:
                pop[i] = np.array([pop[i][0], new_gene, f(pop[i][0], new_gene)])
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
    s = time.perf_counter()
    main()
    f = time.perf_counter()
    print(f-s)
