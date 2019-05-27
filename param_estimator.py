from sys import argv
from math import pow
from math import log 

def cell_num_estimator(a_num, captured_drop_num, capture_rate):
    estimated_drop_num = captured_drop_num / capture_rate
    base = 1 - 1 / estimated_drop_num
    power = 1 - a_num / captured_drop_num
    print("power:", power)
    print("base:", base)
    cell_num = log(power, base)
    return cell_num


def drop_num_estimator(a_num, b_num, shared_num):
    #A_rate = a_num / A_num
    #no_drop_rate = pow(A_rate, 1 / A_num)
    #drop_num = 1 / (1 - no_drop_rate)
    drop_num = a_num * b_num / shared_num 
    return drop_num 


def compute_shared_num(drop_num, A_num, B_num):
#    no_drop_rate = (1 - 1 / drop_num)
#    A_rate = pow(no_drop_rate, B_num)
#    print("A_rate: ", A_rate)
#    shared_num = (1 - A_rate) *  A_num
    A_rate = compute_mix_rate(drop_num, B_num) 
    print("A_rate: ", A_rate)
    shared_num = A_rate *  A_num
    return shared_num


# Computes the rate of drops that have certain cells in them.
def compute_mix_rate(drop_num, cell_num):
    no_drop_rate = (1 - 1 / drop_num)
    cell_in_rate = 1 - pow(no_drop_rate, cell_num)
    return cell_in_rate


#def compute_SSM_rate(a_num, drop_num):
#    cell_num = cell_num_estimator(a_num, drop_num) 
#    SSM_rate = compute_SSM_rate_with_cell_num(cell_num, drop_num)
#    return SSM_rate
   

def compute_SSM_rate_with_cell_num(cell_num, drop_num):
    no_drop_rate = (1 - 1 / drop_num)
    singlet_drops = cell_num * pow(no_drop_rate, cell_num - 1)
    drop_with_cells = (1 - pow(no_drop_rate, cell_num)) * drop_num
    SSM_rate = 1 - singlet_drops / drop_with_cells
    return SSM_rate


if __name__ == "__main__":
    print(compute_shared_num(float(argv[1]), float(argv[2]), float(argv[3])))

