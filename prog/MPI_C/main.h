#ifndef MAIN_H_
#define MAIN_H_

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <mpi.h>

#define SIZE 256
#define STEP 10

int lattice[SIZE][SIZE];
int **segment;
double w[5];
double M = 0.0, E = 0.0;

int ratio = 0;
size_t nmcs = 0;
double ecum = 0.0, e2cum = 0.0, mcum = 0.0, m2cum = 0.0;

int size_mpi;
int rank_mpi;

MPI_Status status;

struct data_bufer{
	double E;
	double M;
	int ratio;
};

void calcW();//вычисление вероятностей перехода
void init();//создание начальной конфигурации
void metropolis();//реализация алгоритма Метрополиса
void step(int step);//число шагов Монте-Карло
void output_data();//вычисление среднего на спин, вывод необходимых данных
void create_data();//создание файла с данными
void load_data();//загрузка файла с данными
void test();//проверка на корректность результата
int **alloc_2d_int(int rows, int cols);

#endif /* MAIN_H_ */
