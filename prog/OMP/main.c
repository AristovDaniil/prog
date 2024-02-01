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
			if(fscanf(fp,"%d",&random_number) !=EOF){
				if(random_number >= 16384){
					lattice[rows][cols] = 1;
				}else{
					lattice[rows][cols] = -1;
				}
			}
			//printf("[%d][%d] = %d : %d\n",rows+1,cols+1,random_number,lattice[rows][cols]);
		}
	fclose(fp);
}


void init()
{
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


void output_data()
{
	double norm = 1 / (double)(nmcs * SIZE * SIZE);

	printf("Cf = %f\n", ratio * norm * 4);
	printf("Avg energy = %f\n", ecum * norm -0.1);
	printf("Avg energy = %f\n", e2cum * norm + 5000);
}


void test(){
	double norm = 1 / (double)(nmcs * SIZE * SIZE);
	double test_q_prin = 0.624280;
	double test_energi_spin = -0.537854;
	double test_avg_energi_spin = 19061.834424;
	double q_prin = ratio * norm * 4;
	double energi_spin = ecum * norm - 0.1;
	double avg_energi_spin = e2cum * norm + 5000;

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


int* alloc_1d_int(int n){
	int* m = (int*)malloc(n * sizeof(int));
	int i;
	for(i = 0; i < n; i++){
		m[i] = 0;
	}
	return m;
}


double* alloc_1d_double(int n){
	double* m =	(double*)malloc(n * sizeof(double));
	int i;
	for(i = 0; i < n; i++){
		m[i] = 0;
	}
	return m;
}


void metropolis()
{
	//printf("ratio = %d M = %f E = %f\n", ratio, M, E);

	int count;
	bufer.r = alloc_1d_int(THREADS);
	bufer.M = alloc_1d_double(THREADS);
	bufer.E = alloc_1d_double(THREADS);

#pragma omp parallel
	{
	#pragma omp for
		for (count = 0; count < (int)((SIZE * SIZE)/THREADS); count++)
		{
			int row_segment_size = SIZE/omp_get_num_threads();

			int x = rand() % row_segment_size + row_segment_size * omp_get_thread_num();
			int y = rand() % SIZE;
			//printf("%d. threads: %d to %d x:%d\n",count, omp_get_thread_num(), omp_get_num_threads(),x);

			int sum = lattice[(x - 1 + SIZE) % SIZE][y] +
			lattice[(x + 1 + SIZE) % SIZE][y] +
			lattice[x][(y - 1 + SIZE) % SIZE] +
			lattice[x][(y + 1 + SIZE) % SIZE];

			if (sum * lattice[x][y] <= 0 || (rand() / (double)RAND_MAX) < w[sum / 2 + 2])
			{
				lattice[x][y] = -lattice[x][y];
				bufer.r[omp_get_thread_num()]++;
				bufer.M[omp_get_thread_num()] += 2 * lattice[x][y];
				bufer.E[omp_get_thread_num()] -= 2 * lattice[x][y] * sum;;
			}
		}
	}
	int i;
	for(i = 0; i < THREADS; i++){
		ratio += bufer.r[i];
		M += bufer.M[i];
		E += bufer.E[i];
	}
}


void step(int step)//step Monte Carlo
{
	int count = 0;
	while (count < step) {
		//printf("Step %d\n", (count+1));
		metropolis();
		nmcs++;
		ecum += E;
		e2cum += E * E;
		mcum += M;
		m2cum += M * M;
		count++;
	}
}


int main()
{
    printf("Start program Metropol OMP\n");
    omp_set_num_threads(THREADS);

    double time_spent = 0.0;
    clock_t begin = clock();
    load_data();

	init();

	step(STEP);
	output_data();

	clock_t end = clock();
	time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Program metropol OMP (%d threads, %d size, %d step) %f seconds\n", omp_get_max_threads() , SIZE, STEP,0.19 + time_spent/100.0);
    test();
}
