#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "datadef.h"
#include "init.h"

#define max(x,y) ((x)>(y)?(x):(y))
#define min(x,y) ((x)<(y)?(x):(y))

extern int *ileft, *iright;
extern int nprocs, proc;

/*
	Due to the complexity of MPI, there needs to be a deterministic selection
	of which process handles which section of the fluid simulation
*/

/* Computation of tentative velocity field (f, g) */
void compute_tentative_velocity(float **u, float **v, float **f, float **g,
    char **flag, int imax, int jmax, float del_t, float delx, float dely,
    float gamma, float Re)
{
    int  i, j;
    float du2dx, duvdy, duvdx, dv2dy, laplu, laplv;

    for (i=1; i<=imax-1; i++) {
        for (j=1; j<=jmax; j++) {
            /* only if both adjacent cells are fluid cells */
            if ((flag[i][j] & C_F) && (flag[i+1][j] & C_F)) {
                du2dx = ((u[i][j]+u[i+1][j])*(u[i][j]+u[i+1][j])+
                    gamma*fabs(u[i][j]+u[i+1][j])*(u[i][j]-u[i+1][j])-
                    (u[i-1][j]+u[i][j])*(u[i-1][j]+u[i][j])-
                    gamma*fabs(u[i-1][j]+u[i][j])*(u[i-1][j]-u[i][j]))
                    /(4.0*delx);
                duvdy = ((v[i][j]+v[i+1][j])*(u[i][j]+u[i][j+1])+
                    gamma*fabs(v[i][j]+v[i+1][j])*(u[i][j]-u[i][j+1])-
                    (v[i][j-1]+v[i+1][j-1])*(u[i][j-1]+u[i][j])-
                    gamma*fabs(v[i][j-1]+v[i+1][j-1])*(u[i][j-1]-u[i][j]))
                    /(4.0*dely);
                laplu = (u[i+1][j]-2.0*u[i][j]+u[i-1][j])/delx/delx+
                    (u[i][j+1]-2.0*u[i][j]+u[i][j-1])/dely/dely;
   
                f[i][j] = u[i][j]+del_t*(laplu/Re-du2dx-duvdy);
            } else {
                f[i][j] = u[i][j];
            }
        }
    }

    for (i=1; i<=imax; i++) {
        for (j=1; j<=jmax-1; j++) {
            /* only if both adjacent cells are fluid cells */
            if ((flag[i][j] & C_F) && (flag[i][j+1] & C_F)) {
                duvdx = ((u[i][j]+u[i][j+1])*(v[i][j]+v[i+1][j])+
                    gamma*fabs(u[i][j]+u[i][j+1])*(v[i][j]-v[i+1][j])-
                    (u[i-1][j]+u[i-1][j+1])*(v[i-1][j]+v[i][j])-
                    gamma*fabs(u[i-1][j]+u[i-1][j+1])*(v[i-1][j]-v[i][j]))
                    /(4.0*delx);
                dv2dy = ((v[i][j]+v[i][j+1])*(v[i][j]+v[i][j+1])+
                    gamma*fabs(v[i][j]+v[i][j+1])*(v[i][j]-v[i][j+1])-
                    (v[i][j-1]+v[i][j])*(v[i][j-1]+v[i][j])-
                    gamma*fabs(v[i][j-1]+v[i][j])*(v[i][j-1]-v[i][j]))
                    /(4.0*dely);

                laplv = (v[i+1][j]-2.0*v[i][j]+v[i-1][j])/delx/delx+
                    (v[i][j+1]-2.0*v[i][j]+v[i][j-1])/dely/dely;

                g[i][j] = v[i][j]+del_t*(laplv/Re-duvdx-dv2dy);
            } else {
                g[i][j] = v[i][j];
            }
        }
    }

    /* f & g at external boundaries */
    for (j=1; j<=jmax; j++) {
        f[0][j]    = u[0][j];
        f[imax][j] = u[imax][j];
    }
    for (i=1; i<=imax; i++) {
        g[i][0]    = v[i][0];
        g[i][jmax] = v[i][jmax];
    }
}


/* Calculate the right hand side of the pressure equation */
void compute_rhs(float **f, float **g, float **rhs, char **flag, int imax,
    int jmax, float del_t, float delx, float dely)
{
    int i, j;

    for (i=1;i<=imax;i++) {
        for (j=1;j<=jmax;j++) {
            if (flag[i][j] & C_F) {
                /* only for fluid and non-surface cells */
                rhs[i][j] = (
                             (f[i][j]-f[i-1][j])/delx +
                             (g[i][j]-g[i][j-1])/dely
                            ) / del_t;
            }
        }
    }
}


/* Red/Black SOR to solve the poisson equation */
int poisson(float **p, float **rhs, char **flag, int imax, int jmax,
    float delx, float dely, float eps, int itermax, float omega,
    float *res, int ifull, int rank, int size, int iStartPos)
{
    int i, j, iter;
    float add, beta_2, beta_mod;
    float p0 = 0.0;
    
    int rb; /* Red-black value. */
	
	/* Bear in mind there are three new variable
	 * rank - Which MPI process we are in
	 * size - How many MPI processes there are
	 * iStartPos - necessary for making sure the rb value is correct
	 */

    float rdx2 = 1.0/(delx*delx);
    float rdy2 = 1.0/(dely*dely);
    beta_2 = -omega/(2.0*(rdx2+rdy2));
	
    /* Calculate sum of squares */
    for (i = 1; i <= imax; i++) {
        for (j=1; j<=jmax; j++) {
            if (flag[i][j] & C_F) { p0 += p[i][j]*p[i][j]; }
        }
    }
	
	// MPI Reduce the calculated p0 value to the root node
	float pTot = 0.0;
	MPI_Reduce(&p0, &pTot, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
	
	MPI_Status stat;
	int tag = 0;
	
	if (rank == 0) {
		// Now redistribute the p0 value
		p0 = sqrt(pTot/ifull);
		if (p0 < 0.0001) { p0 = 1.0; }
		
		printf("Node 0 has reduced the value of p0 as %f\n", p0);
		
		// Need to send to the rest of the nodes
		int node = 0;
		for (node = 1; node < size; node++) {
			MPI_Send(&p0, 1, MPI_FLOAT, node, tag, MPI_COMM_WORLD);
		}
	} else {
		MPI_Recv(&p0, 1, MPI_FLOAT, 0, tag, MPI_COMM_WORLD, &stat);
		printf("Node %d has received the value of p0\n", rank);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
			
	MPI_Request requestLeftSend;
	MPI_Request requestLeftRecv;
	MPI_Request requestRightSend;
	MPI_Request requestRightRecv;
	

    /* Red/Black SOR-iteration */
    for (iter = 0; iter < itermax; iter++) {
        for (rb = 0; rb <= 1; rb++) {
            for (i = 1; i <= imax; i++) {
                for (j = 1; j <= jmax; j++) {
                    if ((i+j+iStartPos) % 2 != rb) { continue; }
                    if (flag[i][j] == (C_F | B_NSEW)) {
                        /* five point star for interior fluid cells */
                        p[i][j] = (1.-omega)*p[i][j] - 
                              beta_2*(
                                    (p[i+1][j]+p[i-1][j])*rdx2
                                  + (p[i][j+1]+p[i][j-1])*rdy2
                                  -  rhs[i][j]
                              );
                    } else if (flag[i][j] & C_F) { 
                        /* modified star near boundary */
                        beta_mod = -omega/((eps_E+eps_W)*rdx2+(eps_N+eps_S)*rdy2);
                        p[i][j] = (1.-omega)*p[i][j] -
                            beta_mod*(
                                  (eps_E*p[i+1][j]+eps_W*p[i-1][j])*rdx2
                                + (eps_N*p[i][j+1]+eps_S*p[i][j-1])*rdy2
                                - rhs[i][j]
                            );
                    }
                } /* end of j */
            } /* end of i */

			
			// This gets performed for both steps of the Red/Black SOR-iteration
			// Need to transfer p values for the next phase of red/black SOR iteration
			if (rank != 0) { // Need to transfer to the left
				// Transfer with tag value of 1, with a non blocking send
				MPI_Isend(p[1], jmax+2, MPI_FLOAT, rank-1, 1, MPI_COMM_WORLD, &requestLeftSend);
				MPI_Irecv(p[0], jmax+2, MPI_FLOAT, rank-1, 2, MPI_COMM_WORLD, &requestRightRecv);
			}
			if (rank != size-1) { // Need to transfer to the right
				MPI_Irecv(p[imax+1], jmax+2, MPI_FLOAT, rank+1, 1, MPI_COMM_WORLD, &requestLeftRecv);
				// Transfer with tag value of 2, non blocking send
				MPI_Isend(p[imax],   jmax+2, MPI_FLOAT, rank+1, 2, MPI_COMM_WORLD, &requestRightSend);
			}
			
			// Now make sure these are received!
			if (rank != 0)  {
				MPI_Wait(&requestRightRecv, &stat);
				MPI_Wait(&requestLeftSend, &stat);
			}
			if (rank != size-1) {
				MPI_Wait(&requestLeftRecv, &stat);
				MPI_Wait(&requestRightSend, &stat);
			}
			
			
        } /* end of rb */
		
        /* Partial computation of residual */
        *res = 0.0;
        for (i = 1; i <= imax; i++) {
            for (j = 1; j <= jmax; j++) {
                if (flag[i][j] & C_F) {
                    /* only fluid cells */
                    add = (eps_E*(p[i+1][j]-p[i][j]) - 
                        eps_W*(p[i][j]-p[i-1][j])) * rdx2  +
                        (eps_N*(p[i][j+1]-p[i][j]) -
                        eps_S*(p[i][j]-p[i][j-1])) * rdy2  -  rhs[i][j];
                    *res += add*add;
                }
            }
        }
		
		float resTot = 0.0;
		MPI_Reduce(res, &resTot, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
		
		if (rank == 0) {
			*res = resTot;
			*res = sqrt((*res)/ifull)/p0;
			int node;
			// Send back out to the rest of the nodes
			for (node = 1; node < size; node++) {
				MPI_Send(res, 1, MPI_FLOAT, node, tag, MPI_COMM_WORLD);
			}
		} else {
			MPI_Recv(res, 1, MPI_FLOAT, 0, tag, MPI_COMM_WORLD, &stat);
		}
		
		
        /* convergence? */
        if (*res<eps) break;
		
    } /* end of iter */
	
	printf("Node %d has completed the Red/Black SOR iteration\n", rank);

    return iter;
}


/* Update the velocity values based on the tentative
 * velocity values and the new pressure matrix
 */
void update_velocity(float **u, float **v, float **f, float **g, float **p,
    char **flag, int imax, int jmax, float del_t, float delx, float dely)
{
    int i, j;

    for (i=1; i<=imax-1; i++) {
        for (j=1; j<=jmax; j++) {
            /* only if both adjacent cells are fluid cells */
            if ((flag[i][j] & C_F) && (flag[i+1][j] & C_F)) {
                u[i][j] = f[i][j]-(p[i+1][j]-p[i][j])*del_t/delx;
            }
        }
    }
    for (i=1; i<=imax; i++) {
        for (j=1; j<=jmax-1; j++) {
            /* only if both adjacent cells are fluid cells */
            if ((flag[i][j] & C_F) && (flag[i][j+1] & C_F)) {
                v[i][j] = g[i][j]-(p[i][j+1]-p[i][j])*del_t/dely;
            }
        }
    }
}


/* Set the timestep size so that we satisfy the Courant-Friedrichs-Lewy
 * conditions (ie no particle moves more than one cell width in one
 * timestep). Otherwise the simulation becomes unstable.
 */
void set_timestep_interval(float *del_t, int imax, int jmax, float delx,
    float dely, float **u, float **v, float Re, float tau, int rank, int size)
{
    int i, j;
    float umax, vmax, deltu, deltv, deltRe; 

    /* del_t satisfying CFL conditions */
    if (tau >= 1.0e-10) { /* else no time stepsize control */
        umax = 1.0e-10;
        vmax = 1.0e-10; 
		float temp;
		// Need to MPI Reduce the umax and vmax values
        for (i=0; i<=imax+1; i++) {
            for (j=1; j<=jmax+1; j++) {
                umax = max(fabs(u[i][j]), umax);
            }
        }
		
		MPI_Reduce(&umax, &temp, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
		if (rank == 0) umax = temp;
        
		for (i=1; i<=imax+1; i++) {
            for (j=0; j<=jmax+1; j++) {
                vmax = max(fabs(v[i][j]), vmax);
            }
        }
		
		MPI_Reduce(&vmax, &temp, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
		if (rank == 0) {
			vmax = temp;

			deltu = delx/umax;
			deltv = dely/vmax; 
			deltRe = 1/(1/(delx*delx)+1/(dely*dely))*Re/2.0;

			if (deltu<deltv) {
				*del_t = min(deltu, deltRe);
			} else {
				*del_t = min(deltv, deltRe);
			}
			*del_t = tau * (*del_t); /* multiply by safety factor */
			
			int node;
			for (node = 1; node < size; node++) {
				MPI_Send(del_t, 1, MPI_FLOAT, node, 0, MPI_COMM_WORLD);
			}
		} else {
			MPI_Status stat;
			MPI_Recv(del_t, 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &stat);
		}
    }
}
