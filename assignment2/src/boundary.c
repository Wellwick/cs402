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
	
	MPI_Status stat;
	
	MPI_Request *sendULeft;
	MPI_Request *sendVLeft;
	MPI_Request *sendPLeft;
	MPI_Request *recieveULeft;
	MPI_Request *recieveVLeft;
	MPI_Request *recievePLeft;
	MPI_Request *sendURight;
	MPI_Request *sendVRight;
	MPI_Request *sendPRight;
	MPI_Request *recieveURight;
	MPI_Request *recieveVRight;
	MPI_Request *recievePRight;
	
	if (rank != 0) {
		// Need to send out u, v and p values to the left
		MPI_Isend(u[1], jmax+2, MPI_FLOAT, rank-1, 1, MPI_COMM_WORLD, sendULeft);
		MPI_Isend(v[1], jmax+2, MPI_FLOAT, rank-1, 2, MPI_COMM_WORLD, sendVLeft);
		MPI_Isend(p[1], jmax+2, MPI_FLOAT, rank-1, 3, MPI_COMM_WORLD, sendPLeft);
		
		// Need to receive from the left as well
		MPI_Irecv(u[0], jmax+2, MPI_FLOAT, rank-1, 4, MPI_COMM_WORLD, recieveURight);
		MPI_Irecv(v[0], jmax+2, MPI_FLOAT, rank-1, 5, MPI_COMM_WORLD, recieveVRight);
		MPI_Irecv(p[0], jmax+2, MPI_FLOAT, rank-1, 6, MPI_COMM_WORLD, recievePRight);
	} 
	if (rank != size-1) {
		// Need to send to the right
		MPI_Isend(u[imax], jmax+2, MPI_FLOAT, rank+1, 4, MPI_COMM_WORLD, sendURight);
		MPI_Isend(v[imax], jmax+2, MPI_FLOAT, rank+1, 5, MPI_COMM_WORLD, sendVRight);
		MPI_Isend(p[imax], jmax+2, MPI_FLOAT, rank+1, 6, MPI_COMM_WORLD, sendPRight);
		
		// Need to receive from the right
		MPI_Irecv(u[imax+1], jmax+2, MPI_FLOAT, rank+1, 1, MPI_COMM_WORLD, recieveULeft);
		MPI_Irecv(v[imax+1], jmax+2, MPI_FLOAT, rank+1, 2, MPI_COMM_WORLD, recieveVLeft);
		MPI_Irecv(p[imax+1], jmax+2, MPI_FLOAT, rank+1, 3, MPI_COMM_WORLD, recievePLeft);
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
		MPI_Wait(sendULeft, &stat);
		MPI_Wait(sendVLeft, &stat);
		MPI_Wait(sendPLeft, &stat);
		MPI_Wait(recieveURight, &stat);
		MPI_Wait(recieveVRight, &stat);
		MPI_Wait(recievePRight, &stat);
	}
	if (rank != size-1) {
		MPI_Wait(sendURight, &stat);
		MPI_Wait(sendVRight, &stat);
		MPI_Wait(sendPRight, &stat);
		MPI_Wait(recieveULeft, &stat);
		MPI_Wait(recieveVLeft, &stat);
		MPI_Wait(recievePLeft, &stat);
	}
}
