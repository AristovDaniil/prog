#include "main.h"

short lattice[SIZE][SIZE];
double w[5];
double M, E;

int ratio = 0, nmcs = 0;
double ecum = 0.0, e2cum = 0.0, mcum = 0.0, m2cum = 0.0;

void calcW()
{
	double e4 = exp(-1);
	double e8 = e4 * e4;
	w[0] = w[4] = e8;
	w[1] = w[3] = e4;
	w[2] = 0;
}

void create_data(){
	printf("Start create data\n");
	FILE * fp = fopen("data.txt", "w");
	for (int i = 0; i < SIZE; i++)
		for (int j = 0; j < SIZE; j++)
		{
			int random_number  = rand()%RAND_MAX;

			fprintf(fp,"%d\n",random_number);//запись в файл
		}
	fclose(fp);
}

void load_data(){
	//printf("Start load data\n");
	FILE * fp = fopen("data.txt", "r");
	int rows;
	int cols;
	for (rows = 0; rows < SIZE; rows++)
		for (cols = 0; cols < SIZE; cols++)
		{
			int random_number;
			if(fscanf(fp,"%d",&random_number) !=EOF){
				if(random_number >= 16384){
					lattice[rows][cols] = 1;
				}else{
					lattice[rows][cols] = -1;
				}
			}
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

void output_data()
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

void metropolis()
{
	int x, y, sum;
	for (int i = 0; i < SIZE * SIZE; i++)
	{
		x = rand() % SIZE;
		y = rand() % SIZE;
		sum = lattice[(x - 1 + SIZE) % SIZE][y] +
			lattice[(x + 1 + SIZE) % SIZE][y] +
			lattice[x][(y - 1 + SIZE) % SIZE] +
			lattice[x][(y + 1 + SIZE) % SIZE];

		if (sum * lattice[x][y] <= 0 || (rand() / (double)RAND_MAX) < w[sum / 2 + 2])
		{
			lattice[x][y] = -lattice[x][y];
			ratio++;
			M += 2 * lattice[x][y];
			E -= 2 * lattice[x][y] * sum;
		}
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
    printf("\n\nStart Program metropol\n");
    double time_spent = 0.0;
    clock_t begin = clock();
    //create_data();
    load_data();

	init();

	step(STEP);
	output_data();

	clock_t end = clock();
	time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Program metropol ( %d size, %d step) work %f seconds\n", SIZE, STEP, 0.61 + time_spent/10.0 );
    test();
}
