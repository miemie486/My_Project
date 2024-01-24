# 导入NumPy模块
import numpy as np
import time
import numba as nb
# 定义参数

POP_SIZE = 4 # 种群大小
GENES = [-1, 1] # 基因范围
CROSS_RATE = 1 # 交叉概率
MUTATE_RATE = 1 # 变异概率
MAX_GEN = 100 # 最大代数
DEMENSION = 10

# 定义目标函数
# Genes_arrary
#@nb.jit(nopython=True)
def f(Genes_arrary):
    fitness=0
    for i in Genes_arrary:
        fitness=fitness+i**2
    return fitness

# 初始化种群
def init_pop():
    Genes_arrary=np.array([np.random.uniform(GENES[0], GENES[1], POP_SIZE)])
    # 使用numpy.random.uniform来生成随机数
    for i in range(DEMENSION-1):
        Genes_arrary=np.r_[Genes_arrary,[np.random.uniform(GENES[0], GENES[1], POP_SIZE)]]
    #Genes_arrary=np.array(Genes_arrary)
    fitness = f(Genes_arrary) # 计算适应度
    Genes_arrary=np.r_[Genes_arrary,[fitness]]
    print(fitness)
    # 使用numpy.array来表示个体和种群
    pop = np.array(Genes_arrary).T # 转置为种群矩阵
    return pop

# 选择操作
def selection(pop):
    # 按适应度从小到大排序
    pop = pop[pop[:, DEMENSION].argsort()]
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
            # 随机选择n个交叉点
            # cp = cross_point
            cp=np.random.randint(0,DEMENSION,size=np.random.randint(1, DEMENSION+1))
            c1=p1[:DEMENSION]
            c2=p2[:DEMENSION]
            for i in cp:
                c1[i]=p2[i]
                c2[i]=p1[i]
            c1=np.r_[c1,[f(c1)]]
            c2=np.r_[c2,[f(c2)]]
        else:
            # 不交叉则保持不变
            c1 = p1
            c2 = p2
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
            print(i)
            # 随机选择n个变异点
            mutate_point = np.random.randint(0,DEMENSION,size=np.random.randint(1, DEMENSION+1))
            print(mutate_point)
            # 随机生成n个新基因
            new_gene = np.random.uniform(GENES[0], GENES[1],len(mutate_point))
            print(new_gene)
            # 替换基因
            print(pop[i])
            print("--")
            np.put(pop[i],mutate_point,new_gene)
            print(pop[i])
    return pop

def main():
    # 初始化种群
    pop = init_pop()
    # 进化指定代数
    print(pop)
    print("---------------------------初始化")
    pop = selection(pop)
    print(pop)
    print("---------------------------选出前一半")
    pop = crossover(pop)
    print(pop)
    print("---------------------------杂交")
    pop = mutate(pop)
    print(pop)
    print("---------------------------变异")

if __name__ == "__main__":
    s = time.perf_counter()
    main()
    f = time.perf_counter()
    print(f-s)
