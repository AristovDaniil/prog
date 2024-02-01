#ifndef MAIN_H_
#define MAIN_H_

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define SIZE 256
#define STEP 10

void calcW();//вычисление вероятностей перехода
void init();//создание начальной конфигурации
void metropolis();//реализация алгоритма Метрополиса
void step(int step);//число шагов Монте-Карло
void output_data();//вычисление среднего на спин, вывод необходимых данных
void create_data();//создание файла с данными
void load_data();//загрузка файла с данными
void test();//проверка на корректность результата

#endif /* MAIN_H_ */
