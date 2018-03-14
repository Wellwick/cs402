#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include "datadef.h"

/* Given the boundary conditions defined by the flag matrix, update
 * the u and v velocities. Also enforce the boundary conditions at the
 * edges of the matrix.
 */
void apply_boundary_conditions(float **u, float **v, float **p, char **flag,
    int imax, int jmax, float ui, float vi, int rank, int size)
{
    int i, j;

    for (j=0; j<=jmax+1; j++) {
        if (rank == 0) {
			/* Fluid freely flows in from the west */
			u[0][j] = u[1][j];
			v[0][j] = v[1][j];
		}
		if (rank == size-1) {
			/* Fluid freely flows out to the east */
			u[imax][j] = u[imax-1][j];
			v[imax+1][j] = v[imax][j];
		}
    }

    for (i=0; i<=imax+1; i++) {
        /* The vertical velocity approaches 0 at the north and south
         * boundaries, but fluid flows freely in the horizontal direction */
        v[i][jmax] = 0.0;
        u[i][jmax+1] = u[i][jmax];

        v[i][0] = 0.0;
        u[i][0] = u[i][1];
    }

    /* Apply no-slip boundary conditions to cells that are adjacent to
     * internal obstacle cells. This forces the u and v velocity to
     * tend towards zero in these cells.
     */
	 
	// Unfortunately there is a dependency here, since each of the previous
	// nodes must have completed before this node starts
	int canGo = 1;
	if (rank != 0) {
		MPI_Status stat;
		MPI_Recv(&canGo, 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD, &stat);
	}
	 
    for (i=1; i<=imax; i++) {
        for (j=1; j<=jmax; j++) {
            if (flag[i][j] & B_NSEW) {
                switch (flag[i][j]) {
                    case B_N: 
                        v[i][j]   = 0.0;
                        u[i][j]   = -u[i][j+1];
                        u[i-1][j] = -u[i-1][j+1];
                        break;
                    case B_E: 
                        u[i][j]   = 0.0;
                        v[i][j]   = -v[i+1][j];
                        v[i][j-1] = -v[i+1][j-1];
                        break;
                    case B_S:
                        v[i][j-1] = 0.0;
                        u[i][j]   = -u[i][j-1];
                        u[i-1][j] = -u[i-1][j-1];
                        break;
                    case B_W: 
                        u[i-1][j] = 0.0;
                        v[i][j]   = -v[i-1][j];
                        v[i][j-1] = -v[i-1][j-1];
                        break;
                    case B_NE:
                        v[i][j]   = 0.0;
                        u[i][j]   = 0.0;
                        v[i][j-1] = -v[i+1][j-1];
                        u[i-1][j] = -u[i-1][j+1];
                        break;
                    case B_SE:
                        v[i][j-1] = 0.0;
                        u[i][j]   = 0.0;
                        v[i][j]   = -v[i+1][j];
                        u[i-1][j] = -u[i-1][j-1];
                        break;
                    case B_SW:
                        v[i][j-1] = 0.0;
                        u[i-1][j] = 0.0;
                        v[i][j]   = -v[i-1][j];
                        u[i][j]   = -u[i][j-1];
                        break;
                    case B_NW:
                        v[i][j]   = 0.0;
                        u[i-1][j] = 0.0;
                        v[i][j-1] = -v[i-1][j-1];
                        u[i][j]   = -u[i][j+1];
                        break;
                }
            }
        }
    }
	
	// Have to let the next node go!
	if (rank != size-1)
		MPI_Send(&canGo, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
	
	MPI_Status stat;
	
	MPI_Request sendLeft;
	MPI_Request recieveLeft;
	MPI_Request sendRight;
	MPI_Request recieveRight;
	
	float recieveData[(jmax+2)*3];
	
	if (rank != 0) {
		float sendData[(jmax+2)*3];
		for (i = 0; i < jmax+2; i++) {
			sendData[i] 			 = u[1][i];
			sendData[i+jmax+2] 		 = v[1][i];
			sendData[i+((jmax+2)*2)] = p[1][i];
		}
		// Need to send out u, v and p values to the left
		MPI_Isend(&sendData, (jmax+2)*3, MPI_FLOAT, rank-1, 3, MPI_COMM_WORLD, &sendLeft);
		
		// Need to receive from the left as well
		MPI_Irecv(&recieveData, (jmax+2)*3, MPI_FLOAT, rank-1, 4, MPI_COMM_WORLD, &recieveRight);
	} 
	if (rank != size-1) {
		float sendData[(jmax+2)*3];
		for (i = 0; i < jmax+2; i++) {
			sendData[i] 			 = u[imax][i];
			sendData[i+jmax+2] 		 = v[imax][i];
			sendData[i+((jmax+2)*2)] = p[imax][i];
		}
		// Need to send to the right
		MPI_Isend(&sendData, (jmax+2)*3, MPI_FLOAT, rank+1, 4, MPI_COMM_WORLD, &sendRight);
		
		// Need to receive from the right
		MPI_Irecv(&recieveData, (jmax+2)*3, MPI_FLOAT, rank+1, 3, MPI_COMM_WORLD, &recieveLeft);
	}

    /* Finally, fix the horizontal velocity at the  western edge to have
     * a continual flow of fluid into the simulation.
     */
	if (rank == 0) {
		v[0][0] = 2*vi-v[1][0];
		for (j=1;j<=jmax;j++) {
			u[0][j] = ui;
			v[0][j] = 2*vi-v[1][j];
		}
	}
	
	if (rank != 0) {
		MPI_Wait(&sendLeft, &stat);
		MPI_Wait(&recieveRight, &stat);
		for (i = 0; i < jmax+2; i++) {
			u[0][i] = recieveData[i];
			v[0][i] = recieveData[i+(jmax+2)];
			p[0][i] = recieveData[i+((jmax+2)*2)];
		}
	}
	if (rank != size-1) {
		MPI_Wait(&sendRight, &stat);
		MPI_Wait(&recieveLeft, &stat);
		for (i = 0; i < jmax+2; i++) {
			u[imax+1][i] = recieveData[i];
			v[imax+1][i] = recieveData[i+(jmax+2)];
			p[imax+1][i] = recieveData[i+((jmax+2)*2)];
		}
		
	}
}
