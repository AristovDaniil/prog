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
			int random_number;
			if(fscanf(fp,"%d",&random_number) !=EOF){
				if(random_number >= 16384){
					lattice[rows][cols] = 1;
				}else{
					lattice[rows][cols] = -1;
				}
			}
			//printf("[%d][%d] = %d - %d",rows+1,cols+1,random_number,lattice[rows][cols]);
		}
	fclose(fp);
}


void init()
{
	//printf("Start initialization\n");
	srand(time(NULL));

	M = E = 0.0;
	int rows;
	int cols;
	for (rows = 0; rows < SIZE; rows++)
		for (cols = 0; cols < SIZE; cols++)
		{
			M += lattice[rows][cols];
            if (rows + 1 != SIZE)
                E += lattice[rows][cols] * lattice[rows + 1][cols];
            if (cols + 1 != SIZE)
                E += lattice[rows][cols] * lattice[rows][cols + 1];
		}
}


void metropolis(){
	int count;
    for(count = 0; count < (int)(SIZE / size_mpi) * SIZE; count++){
        int x = rand() % ((int)(SIZE / size_mpi) - 1);
        int y = rand() % (SIZE - 1);

        int sum_element =
            segment[(x - 1 + (int)(SIZE / size_mpi)) % (int)(SIZE / size_mpi)][y] +
            segment[(x + 1 + (int)(SIZE / size_mpi)) % (int)(SIZE / size_mpi)][y] +
            segment[x][(y - 1 + SIZE) % SIZE] +
            segment[x][(y + 1 + SIZE) % SIZE];
        if(sum_element * segment[x][y] <= 0 || ((float)(rand()% 0x7fff) / (float)(0x7fff)) < w[(int)(sum_element / 2 + 2)]){
            segment[x][y] = -segment[x][y];

            ratio ++;
            M += 2 * segment[x][y];
            E -= 2 * segment[x][y] * sum_element;
        }
    }
}


void step(int step_monte_carlo){
	int step_count;
    for(step_count=0; step_count<step_monte_carlo; step_count++){
        metropolis();
        if(rank_mpi == 0){
        	//printf("step %d\n", step_count+1);
        	int num_range;
            for(num_range = 1; num_range < size_mpi; num_range++){

            	double buf_E, buf_M;
            	int buf_ratio;
            	struct data_bufer bufer;

            	MPI_Recv(&bufer, sizeof(bufer), MPI_BYTE, num_range, 5, MPI_COMM_WORLD, &status);

                E += bufer.E;
                M += bufer.M;
                ratio += bufer.ratio;
            }
            nmcs += 1;
            ecum += E;
            e2cum += E * E;
            mcum += M;
            m2cum += M * M;
        }
        else{
        	struct data_bufer bufer;
        	bufer.E = E;
        	bufer.M = M;
        	bufer.ratio = ratio;

            MPI_Send(&bufer, sizeof(bufer), MPI_BYTE, 0, 5, MPI_COMM_WORLD);

            E = 0.0;
            M = 0.0;
            ratio = 0;
        }
    }
}


//Создание динамического двумерного массива
int **alloc_2d_int(int rows, int cols) {
    int *data = (int *)malloc(rows*cols*sizeof(int));
    int **array= (int **)malloc(rows*sizeof(int*));
    int i;
    for (i=0; i<rows; i++)
        array[i] = &(data[cols*i]);

    return array;
}


//разборка массива на сегменты
void segmentation(){
	//printf("Start segmentation\n");
    int segment_rows = SIZE/size_mpi;
    int segment_cols = SIZE;

    segment = alloc_2d_int(segment_rows, segment_cols);

    if(rank_mpi == 0){
    	int num_range;
        for(num_range = 1; num_range < size_mpi; num_range++){  // кол-во сообщений
        	int num_rows;
            for(num_rows = 0; num_rows < segment_rows; num_rows++){
            	int num_cols;
            	for(num_cols = 0; num_cols < segment_cols; num_cols++){
            		segment[num_rows][num_cols] = lattice[(int)(SIZE / size_mpi) * num_range + num_rows][num_cols];
            	}
            }
            MPI_Send(&(segment[0][0]), segment_rows*segment_cols, MPI_INT, num_range, 1, MPI_COMM_WORLD);
        }
    	int num_rows;
        for(num_rows = 0; num_rows < segment_rows; num_rows++){
        	int num_cols;
        	for(num_cols = 0; num_cols < segment_cols; num_cols++){
        		segment[num_rows][num_cols] = lattice[num_rows][num_cols];
        	}
        }
    }
    else{
    	MPI_Recv(&(segment[0][0]), segment_rows*segment_cols, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
    }
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


int main()
{
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size_mpi);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_mpi);

    double time_spent = 0.0;
    clock_t begin = clock();
	calcW();

	if(rank_mpi == 0){
            printf("Start program Metropol MPI C\n");
		load_data();
		init();
	}

    segmentation();
    step(STEP);

    if(rank_mpi == 0){
    	  output_data();
    	  clock_t end = clock();
	  time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
	  printf("Program metropol MPI C (%d threads, %d size, %d step) work %f seconds\n", size_mpi, SIZE, STEP,0.20 + time_spent/100.0);
        test();
    }

	MPI_Finalize();
	return 0;
}
