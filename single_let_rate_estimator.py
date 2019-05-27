from sys import argv
from math import pow
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

def calculate_single_let_rate(cell_num, drop_num):
    cell_num = float(cell_num) * 1000
    drop_num = float(drop_num) * 1000

    drop_prob = 1/drop_num

    no_drop_prob = 1 - drop_prob

    single_let_rate = cell_num * pow(no_drop_prob, cell_num-1) / (drop_num * (1 - pow(no_drop_prob, cell_num)) )

    return single_let_rate


def plot_exponential_dist(cell_num, drop_num):
    lambd = cell_num / drop_num / 1000
    scale = 1000 * drop_num / cell_num 
    x = np.arange(0,50000, 1)
    f = np.vectorize(lambda x : stats.expon.pdf(x, scale=scale), otypes=[np.float])
    print(stats.expon.cdf(1, scale=scale))
    y = f(x)
    #y = lambd * np.exp(-lambd * x)
    plt.plot(x,y)
    plt.title("Plot")
    plt.xlabel('x')
    plt.ylabel('P')
    plt.show()


if __name__ == "__main__":
    single_let_rate = calculate_single_let_rate(float(argv[1]), float(argv[2]))
    print(single_let_rate)
    #plot_exponential_dist(float(argv[1]), float(argv[2]) )

