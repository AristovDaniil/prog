#ifndef MAIN_H_
#define MAIN_H_

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <pthread.h>

#define SIZE 256
#define THREADS 4
#define STEP 10

short lattice[SIZE][SIZE];
double w[5];
double M, E;

int ratio = 0;
size_t nmcs = 0;
double ecum = 0.0, e2cum = 0.0, mcum = 0.0, m2cum = 0.0;

//структура для данных потока
typedef struct{
	int sigment;
	int ratio;
	double M;
	double E;
} pthrData;

void calcW();//вычисление вероятностей перехода
void init();//создание начальной конфигурации
void* metropolis(void* thread_data);//реализация алгоритма Метрополиса
void step(int step);//число шагов Монте-Карло
void outputData();//вычисление среднего на спин, вывод необходимых данных
void load_data();//загрузка файла с данными
void test();//проверка на корректность результата

#endif /* MAIN_H_ */
