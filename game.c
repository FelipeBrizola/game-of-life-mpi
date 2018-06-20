#include <stdio.h>
#include <omp.h> // compilado em macos. Se utilizar linux, troca aspas por <>

#define CELL(I, J) (field[size * (I) + (J)])
#define ALIVE(I, J) t[size * (I) + (J)] = 1
#define DEAD(I, J) t[size * (I) + (J)] = 0

int FIELD_SIZE = 0;
int FIELD_GEN = 0;

int count_alive(const char *field, int i, int j, int size) {
    int x, y, a = 0;

    for (x = i - 1; x <= (i + 1); x++) {
        for (y = j - 1; y <= (j + 1); y++) {
            if ((x == i) && (y == j))
                continue;
            if ((y < size) && (x < size) &&
                (x >= 0) && (y >= 0)) {
                a += CELL(x, y);
            }
        }
    }
    return a;
}

void evolve(const char *field, char *t, int size) {
    int i, j, alive, cs;

    #pragma omp parallel private (i, j, cs, alive)
    for (i = 0; i < size; i++) {
        
      #pragma omp for
        for (j = 0; j < size; j++) {
            alive = count_alive(field, i, j, size);
            cs = CELL(i, j);
            
            if (cs) {
                if ((alive > 3) || (alive < 2))
                    DEAD(i, j);
                else
                    ALIVE(i, j);
            }
            else {
                if (alive == 3)
                    ALIVE(i, j);
                else
                    DEAD(i, j);
            }
        }
    }
}


/* set das celulas i, j como ativas */
#define SCELL(I, J) field[FIELD_SIZE * (I) + (J)] = 1

void dump_field(const char *f, int size) {
    int i;
    
    for (i = 0; i < (size * size); i++) {
        if ((i % size) == 0)
            printf("\n");
        printf("%c", f[i] ? 'X' : '.');
    }
    printf("\n");
}

int main(int argc, char **argv) {
    int i;
    char *fa, *fb, *tt;
    double starttime, stoptime;
    int NUM_THREADS = 0;

    if (argc != 4) {
        fprintf(stderr, "error - the parameters are: gen, gen_size, num_threads\n");
        return 1;
    }
    
    if (sscanf (argv[1], "%i", &FIELD_GEN) != 1) {
        fprintf(stderr, "error - not an integer");
    }
    if (sscanf (argv[2], "%i", &FIELD_SIZE) != 1) {
        fprintf(stderr, "error - not an integer");
    }
    if (sscanf (argv[3], "%i", &NUM_THREADS) != 1) {
        fprintf(stderr, "error - not an integer");
    }
    
    // quantidade de threads
    omp_set_num_threads(NUM_THREADS);
    
    char field[FIELD_SIZE * FIELD_SIZE];
    char temp_field[FIELD_SIZE * FIELD_SIZE];
    int size = FIELD_SIZE * FIELD_SIZE;

    for (i = 0; i < size; i++)
        field[i] = 0;

    // Set de valores da matriz
    SCELL(10, 1);
    SCELL(10, 2);
    SCELL(11, 2);
    SCELL(2, 0);
    SCELL(2, 1);
    SCELL(2, 2);
    SCELL(2, 10);
    SCELL(3, 11);
    SCELL(3, 9);
    SCELL(1, 12);
    SCELL(2, 15);
    SCELL(3, 14);
    SCELL(3, 13);
    SCELL(1, 16);
    SCELL(20, 21);
    SCELL(19, 20);
    SCELL(18, 22);

    fa = field;
    fb = temp_field;
    
    fprintf(stdout, "Running with %d  threads... print first gen\n", omp_get_num_threads());
    dump_field(fa, FIELD_SIZE); // print da primeira execucao
    
    starttime = omp_get_wtime();

    for (i = 0; i < FIELD_GEN; i++) {
        evolve(fa, fb, FIELD_SIZE);
        tt = fb;
        fb = fa;
        fa = tt;
    }

    stoptime = omp_get_wtime();
    fprintf(stdout, "Done... print last gen\n");
    dump_field(fa, FIELD_SIZE); // print da ultima execucao
    printf("Tempo de execucao %3.2f segundos \n", stoptime - starttime);

    return 0;

}






