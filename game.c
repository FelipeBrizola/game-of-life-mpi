#include <stdio.h>
#include <mpi.h>
#define MASTER 0

#define CELL(I, J) (field[size * (I) + (J)])
// #define CELL(I, J) (field[size * (I) + (J)])
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

    for (i = 0; i < size; i++) {
        
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
    char *fa, *fb, *fc, *tt;

    int size, rank, tag, rc, N, generations, outPoints, s, z, k, dest;
    MPI_Status Stat;

    if (argc != 3) {
        fprintf(stderr, "error - the parameters are: gen, gen_size\n");
        return 1;
    }
    
    if (sscanf (argv[1], "%i", &FIELD_GEN) != 1) {
        fprintf(stderr, "error - not an integer");
    }
    if (sscanf (argv[2], "%i", &FIELD_SIZE) != 1) {
        fprintf(stderr, "error - not an integer");
    }
    N = FIELD_SIZE;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    char field[FIELD_SIZE * FIELD_SIZE];
    char temp_field[FIELD_SIZE * FIELD_SIZE];

    if (rank == MASTER) {


        s = N/size; // quantos pedacos

        for (i = 0; i < FIELD_SIZE * FIELD_SIZE; i++)
            field[i] = 0;

        // Set de valores da matriz
        SCELL(10, 1); SCELL(10, 2); SCELL(11, 2); SCELL(2, 0); SCELL(2, 1);
        SCELL(2, 2); SCELL(2, 10); SCELL(3, 11); SCELL(3, 9); SCELL(1, 12);
        SCELL(2, 15); SCELL(3, 14); SCELL(3, 13); SCELL(1, 16); SCELL(20, 21);
        SCELL(19, 20); SCELL(18, 22);

        fa = field;
        fb = temp_field;
        fc = field;

        dump_field(fa, FIELD_SIZE); // print da primeira execucao

        for (i = 0; i < FIELD_GEN; i++) {
            evolve(fa, fb, FIELD_SIZE);
            tt = fb;
            fb = fa;
            fa = tt;
        }
        


		//SENDING INITIAL INFORMATION (N, k, #generations, output points) TO EVERYONE
        int info[4];
		info[0] = N;
        info[1] = s;
        info[2] = generations;
        info[3] = outPoints;

        for (dest = 0; dest < size; dest++)
            MPI_Send(&info, 4, MPI_INT, dest, 1, MPI_COMM_WORLD); //send info

        int slice[N/size][N];
        for (z = 0; z < size; z++) {
			for (k = 0; k < s; k++) 
				for (int l=0; l<N; l++) 
					slice[k][l]=CELL(k+(z*s), l);	//cut a slice from the the board

			MPI_Send(&slice, N*s, MPI_INT, z, 1, MPI_COMM_WORLD);	//and send it
        }
        
    }
    
    //RECEIVED INITIAL INFORMATION
	int localinfo[4];		// local info for initial information
	MPI_Recv(&localinfo, 4, MPI_INT, 0, 1, MPI_COMM_WORLD, &Stat);	//receive info
	int myslice[localinfo[1]][localinfo[0]]; //my own slice of the board
	MPI_Recv(&myslice, localinfo[0]*localinfo[1], MPI_INT, 0, 1, MPI_COMM_WORLD, &Stat);	//receive slice
	N = localinfo[0];			
	s = localinfo[1];			
	generations=localinfo[2];	
	outPoints=localinfo[3];		
	
	int todown[N], toup[N], fromdown[N], fromup[N]; //arrays to send and to receive

    for(i = 0; i < FIELD_GEN; i++) {
        // all except for last send down
        if (rank!=size-1)  {
			for (int j = 0; j < N; j++)
                todown[j]=myslice[s-1][j];
			MPI_Send(&todown, N, MPI_INT, rank+1, 1, MPI_COMM_WORLD);

		} else {
            // last one generates empty stripe "from down"
			for (int k=0; k<N; k++)
                fromdown[k]=0;
        } 

        // all except for first receive from up
		if (rank!=0)
			MPI_Recv(&fromup, N, MPI_INT, rank-1, 1, MPI_COMM_WORLD, &Stat);	

		else {
            // first one generats empty line "from up"	
            for (int k=0; k<N; k++)
                fromup[k]=0;
        } 
	
        // all except for first send up
		if (rank != MASTER)  {
			for (int j=0; j<N; j++)
                toup[j]=myslice[0][j];
			MPI_Send(&toup, N, MPI_INT, rank-1, 1, MPI_COMM_WORLD);
		}

        // all except for last receive from down
		if (rank != size-1)
			MPI_Recv(&fromdown, N, MPI_INT, rank+1, 1, MPI_COMM_WORLD, &Stat);


        int sum = 0;

        //for each row
        //for each column
        for (int x=0; x<s; x++)  {	
			for (int y=0; y<N; y++)  {

                // count_alive(fc, x, y, 3)

                if (x==0 && y==0) //upper-left cell
					sum = count_alive(fc, x+1, y+1, 3) + count_alive(fc, x+1, y+1, 3) + count_alive(fc, 0, y+1, 3) +fromup[0] + fromup[1];

				else if (x==0 && y==N-1) //upper-right cell
                    sum = count_alive(fc, x, y-1, 3) + count_alive(fc, x+1, y-1, 3) + count_alive(fc, x+1, y+1, 3) + fromup[N-1]+fromup[N-2];

				else if (x==s-1 && y==0) //lower-left cell
                    sum = count_alive(fc, x, y+1, 3) + count_alive(fc, x-1, y+1, 3) + count_alive(fc, x-1, y, 3) + fromdown[0]+fromdown[1];

				else if (x==s-1 && y==N-1) //lower-right cell
                    sum = count_alive(fc, x-1, y, 3) + count_alive(fc, x-1, y-1, 3) + count_alive(fc, x, y-1, 3) + fromdown[N-1]+fromdown[N-2];

				else // not corner cells    
				{
					if (y==0) // leftmost line, not corner
                        sum  = count_alive(fc, x-1, y, 3) + count_alive(fc, x-1, y+1, 3) + count_alive(fc, x, y+1, 3) + count_alive(fc, x+1, y+1, 3) +count_alive(fc, x+1, y, 3);

					else if (y==N-1) //rightmost line, not corner
                        sum  = count_alive(fc, x-1, y, 3) + count_alive(fc, x-1, y-1, 3) + count_alive(fc, x, y-1, 3) + count_alive(fc, x+1, y-1, 3) +count_alive(fc, x+1, y, 3);

					else if (x==0) //uppermost line, not corner
                        sum  = count_alive(fc, x, y-1, 3) + count_alive(fc, x+1, y-1, 3) + count_alive(fc, x+1, y, 3) + count_alive(fc, x+1, y+1, 3) + fromup[y-1] + fromup[y]+fromup[y+1];

						
					else if (x==s-1) //lowermost line, not corner
                        sum = count_alive(fc, x-1, y-1, 3) + count_alive(fc, x-1, y, 3) + count_alive(fc, x-1, y+1, 3) + count_alive(fc, x, y+1, 3) +  count_alive(fc, x, y-1, 3) +fromdown[y-1]+fromdown[y]+fromdown[y+1];
                    
					else //general case, any cell within
                        sum = count_alive(fc, x-1, y-1, 3) + count_alive(fc, x-1, y, 3) + count_alive(fc, x-1, y+1, 3) + count_alive(fc, x, y+1, 3) +  count_alive(fc, x+1, y+1, 3) + count_alive(fc, x+1, y, 3) + count_alive(fc, x+1, y-1, 3) + count_alive(fc, x, y-1, 3);

				}
				
                evolve(fc, fb, 3);

				
                
            }
        }


        // evolve(fc, fb, size);

        // copy new slice onto myslice
		for (int x=0; x<s; x++)
			for (int y=0; y<N; y++)
                fa = &CELL(x,y);


        //s-th generation, send everything to node 0
        if (FIELD_GEN % outPoints==0) {
			if (rank==0)  {
                dump_field(fc, FIELD_SIZE); // print da ultima execucao
            }
        }

    }

    
    MPI_Finalize();    
    
    
    return 0;

}

// ./mpirun -np 1 hello 25 25
