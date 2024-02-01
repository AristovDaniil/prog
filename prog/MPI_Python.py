import math
import time
from random import randint
from mpi4py import MPI
import numpy as np

comm = MPI.COMM_WORLD
rank_mpi = comm.Get_rank()
size_mpi = comm.Get_size()

SIZE = 256
STEP = 10

lattice = np.zeros((SIZE, SIZE))
segment = np.zeros((int(SIZE / size_mpi), SIZE))

w = [0.0] * 5
T_MAX = 4
T = T_MAX
M = 0
E = 0
ratio = 0
nmcs = 0
ecum = 0.0
e2cum = 0.0
mcum = 0.0
m2cum = 0


def calc_w():
    global T
    global w
    e4 = math.exp(-1)
    e8 = e4 * e4
    w[0] = w[4] = e8
    w[1] = w[3] = e4
    w[2] = 0


def load_data():
    global lattice
    #print("Start load data\n")
    fp = open("data.txt", "r")
    for i in range(SIZE):
        for j in range(SIZE):
            random_number = int(fp.readline())
            lattice[i][j] = ((random_number / float(0x7fff) >= 0.5) - 1) * 2 + 1
    fp.close()
    #print("End load data\n")


def init():
    #print("Start initialization")
    global M, E
    M = 0
    E = 0

    for i in range(SIZE):
        for j in range(SIZE):
            M += lattice[i][j]
            if i + 1 != SIZE:
                E += lattice[i][j] * lattice[i + 1][j]
            if j + 1 != SIZE:
                E += lattice[i][j] * lattice[i][j + 1]


def metropolis():
    global ratio, M, E, w, segment

    for count in range(int(SIZE / size_mpi) * SIZE):
        x = randint(0, int(SIZE / size_mpi) - 1)
        y = randint(0, SIZE - 1)

        sum_element = \
            segment[(x - 1 + int(SIZE / size_mpi)) % int(SIZE / size_mpi)][y] +\
            segment[(x + 1 + int(SIZE / size_mpi)) % int(SIZE / size_mpi)][y] +\
            segment[x][(y - 1 + SIZE) % SIZE] +\
            segment[x][(y + 1 + SIZE) % SIZE]
        if sum_element * segment[x][y] <= 0 or float(randint(0, 0x7fff)) / float(0x7fff) < w[int(sum_element / 2 + 2)]:
            segment[x][y] = -segment[x][y]
            ratio += 1
            M += 2 * segment[x][y]
            E -= 2 * segment[x][y] * sum_element


def step(step_monte_carlo=0):
    global nmcs, ecum, e2cum, mcum, m2cum, M, E, ratio

    for step_count in range(step_monte_carlo):
        metropolis()
        if rank_mpi == 0:
            #print("step ", step_count + 1)
            for num_range in range(1, size_mpi):
                buf = comm.recv(source=num_range, tag=3)
                E += buf[0]
                M += buf[1]
                ratio += buf[2]
            nmcs += 1
            ecum += E
            e2cum += E * E
            mcum += M
            m2cum += M * M
        else:
            comm.send([E, M, ratio], dest=0, tag=3)
            M = E = ratio = 0


def output_data():
    global E, M, ratio, ecum, e2cum, mcum, m2cum
    norm = 1 / float(nmcs * SIZE * SIZE)

    print("Cf =", ratio * norm)
    print("Avg energy =", ecum * norm)
    print("Avg energy2 =", e2cum * norm)


def test():
    norm = 1 / float(nmcs * SIZE * SIZE)
    test_q_prin = 0.624280
    test_energi_spin = -0.537854
    test_avg_energi_spin = 19061.834424
    q_prin = ratio * norm
    energi_spin = ecum * norm
    avg_energi_spin = e2cum * norm

    if 0.9 * math.fabs(test_q_prin) < math.fabs(q_prin) < 1.1 * math.fabs(test_q_prin) and \
            0.9 * math.fabs(test_energi_spin) < math.fabs(energi_spin) < 1.1 * math.fabs(test_energi_spin) and \
            0.9 * math.fabs(test_avg_energi_spin) < math.fabs(avg_energi_spin) < 1.1 * math.fabs(test_avg_energi_spin):
        print("Test passed\n")
    else:
        print("Test failed")
        print("q_prin: ", 0.9 * math.fabs(test_q_prin), math.fabs(q_prin), 1.1 * math.fabs(test_q_prin))
        print("energi_spin: ", 0.9 * math.fabs(test_energi_spin), math.fabs(energi_spin),
              1.1 * math.fabs(test_energi_spin))
        print("avg_energi_spin: ", 0.9 * math.fabs(test_avg_energi_spin), math.fabs(avg_energi_spin),
              1.1 * math.fabs(test_avg_energi_spin))


# разборка массивов в потоки
def segmentation():
    global segment
    if rank_mpi == 0:
        for num_range in range(1, size_mpi):  # кол-во сообщений
            for num_row in range(int(SIZE / size_mpi)):
                segment[num_row] = lattice[int(SIZE / size_mpi) * num_range + num_row]
            comm.send(segment, dest=num_range, tag=1)
        for num_row in range(int(SIZE / size_mpi)):
            segment[num_row] = lattice[num_row]
    else:
        segment = comm.recv(source=0, tag=1)


def make_segment():
    global lattice, segment
    if rank_mpi == 0:
        for i in range(1, size_mpi):  # кол-во сообщений
            segment = comm.recv(source=i, tag=2)
            for j in range(int(SIZE / size_mpi)):
                lattice[int(SIZE / size_mpi) * i + j] = segment[j]
    else:
        comm.send(segment, dest=0, tag=2)


start_time = 0.0
calc_w()
if rank_mpi == 0:
    print("Start program metropol MPI Python")
    start_time = time.time()
    load_data()
    init()
segmentation()
step(STEP)

if rank_mpi == 0:
    output_data()
    print("Program metropol MPI Python (thread ", size_mpi, ",size ", SIZE, ",step ", STEP, ") work ", (time.time() - start_time), " sec")
    test()
