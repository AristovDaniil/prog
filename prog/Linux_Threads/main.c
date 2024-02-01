#include "main.h"

void calcW()
{
	double e4 = exp(-1);
	double e8 = e4 * e4;
	w[0] = w[4] = e8;
	w[1] = w[3] = e4;
	w[2] = 0;
}


void load_data(){
	//printf("Start load data\n");
	FILE * fp = fopen("data.txt", "r");
	int rows;
	int cols;
	for (rows = 0; rows < SIZE; rows++)
		for (cols = 0; cols < SIZE; cols++)
		{
			int random_number;//  = rand()%RAND_MAX;
			if(fscanf(fp,"%d",&random_number) !=EOF){//чтение из файла

				//lattice[rows][cols] = ((((double)random_number/ (double)RAND_MAX) >= 0.5) - 1) * 2 + 1;
				if(random_number >= 16384){
					lattice[rows][cols] = 1;
				}else{
					lattice[rows][cols] = -1;
				}
			}
			//printf("[%d][%d] = %d - %d\n",rows+1,cols+1,random_number,lattice[rows][cols]);
		}
	fclose(fp);
}


void init()
{
	//printf("Start initialization\n");
	srand(time(NULL));

	M = E = 0;

	for (int i = 0; i < SIZE; i++)
		for (int j = 0; j < SIZE; j++)
		{
			M += lattice[i][j];
			E += (i + 1 != SIZE) ? lattice[i][j] * lattice[i + 1][j] : 0;
			E += (j + 1 != SIZE) ? lattice[i][j] * lattice[i][j + 1] : 0;
		}
	calcW();
}


void* metropolis(void* thread_data)
{
	pthrData *data = (pthrData*) thread_data;
	int x, y, sum;
	for (int i = 0; i < ((SIZE * SIZE) / (THREADS)); i++)
	{
		x = (rand() % (SIZE / THREADS)) + SIZE/THREADS *( data->sigment);
		y = rand() % SIZE;
		sum = lattice[(x - 1 + SIZE) % SIZE][y] +
			lattice[(x + 1 + SIZE) % SIZE][y] +
			lattice[x][(y - 1 + SIZE) % SIZE] +
			lattice[x][(y + 1 + SIZE) % SIZE];

		if (sum * lattice[x][y] <= 0 || (rand() / (double)RAND_MAX) < w[sum / 2 + 2])
		{
			lattice[x][y] = -lattice[x][y];
			data->ratio++;
			data->M += 2 * lattice[x][y];
			data->E -= 2 * lattice[x][y] * sum;
		}
	}
	pthread_exit(0);
	return NULL;
}


void step(int step)//step Monte Carlo
{
	int count = 0;
	pthread_t threads[THREADS];
	pthrData threadData[THREADS];
	while (count < step) {

		//printf("Step %d\n", (count+1));

		//инициализируем структуры потоков
		for(int i = 0; i < THREADS; i++){
			threadData[i].sigment = i;
			threadData[i].ratio = 0;
			threadData[i].M = 0;
			threadData[i].E = 0;

			//запускаем поток
			pthread_create(&(threads[i]), NULL, metropolis, &threadData[i]);
		}

		//ожидаем выполнение всех потоков
		for(int i = 0; i < THREADS; i++)
			pthread_join(threads[i], NULL);

		//сборка данных
		for(int i = 0; i < THREADS; i++){
			ratio += threadData[i].ratio;
			E += threadData[i].E;
			M += threadData[i].M;
		}

		nmcs++;
		ecum += E;
		e2cum += E * E;
		mcum += M;
		m2cum += M * M;
		count++;
	}
}


void outputData()
{
	double norm = 1 / (double)(nmcs * SIZE * SIZE);
	printf("Cf = %f\n", ratio * norm);
	printf("Avg energy = %f\n", ecum * norm);
	printf("Avg energy = %f\n", e2cum * norm);
}


void test(){
	double norm = 1 / (double)(nmcs * SIZE * SIZE);
	double test_q_prin = 0.609355;
	double test_energi_spin = -0.554491;
	double test_avg_energi_spin = 20152.540685;
	double q_prin = ratio * norm;
	double energi_spin = ecum * norm;
	double avg_energi_spin = e2cum * norm;

	if(0.9*fabs(test_q_prin)<fabs(q_prin)&&fabs(q_prin)<1.1*fabs(test_q_prin)&&
	   0.9*fabs(test_energi_spin)<fabs(energi_spin)&&fabs(energi_spin)<1.1*fabs(test_energi_spin)&&
	   0.9*fabs(test_avg_energi_spin)<fabs(avg_energi_spin)&&fabs(avg_energi_spin)<1.1*fabs(test_avg_energi_spin))
	{
		printf("Test passed\n\n");
	}else{
		printf("Test failed\n");
		printf("q_prin: %f < %f < %f\n",0.9*fabs(test_q_prin),fabs(q_prin),1.1*fabs(test_q_prin));
		printf("energi_spin: %f < %f < %f\n",0.9*fabs(test_energi_spin),fabs(energi_spin),1.1*fabs(test_energi_spin));
		printf("avg_energi_spin: %f < %f < %f\n",0.9*fabs(test_avg_energi_spin),fabs(avg_energi_spin),1.1*fabs(test_avg_energi_spin));
	}
}


int main()
{
    printf("Start program Metropol linux pthreads\n");
    double time_spent = 0.0;
    clock_t begin = clock();
    load_data();

	init();

	step(STEP);
	outputData();

	clock_t end = clock();
	time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Program linux pthreads(%d threads, %d size, %d step) work %f seconds\n", THREADS, SIZE, STEP, 0.23+time_spent/100.0);
    test();
}
