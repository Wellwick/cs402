#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <errno.h>
#include "alloc.h"
#include "boundary.h"
#include "datadef.h"
#include "init.h"
#include "simulation.h"

void write_bin(float **u, float **v, float **p, char **flag,
     int imax, int jmax, float xlength, float ylength, char *file);

int read_bin(float **u, float **v, float **p, char **flag,
    int imax, int jmax, float xlength, float ylength, char *file);

static void print_usage(void);
static void print_version(void);
static void print_help(void);

static char *progname;

int proc = 0;                       /* Rank of the current process */
int nprocs = 0;                /* Number of processes in communicator */

int *ileft, *iright;           /* Array bounds for each processor */

#define PACKAGE "karman"
#define VERSION "1.0"

/* Command line options */
static struct option long_opts[] = {
    { "del-t",   1, NULL, 'd' },
    { "help",    0, NULL, 'h' },
    { "imax",    1, NULL, 'x' },
    { "infile",  1, NULL, 'i' },
    { "jmax",    1, NULL, 'y' },
    { "outfile", 1, NULL, 'o' },
    { "t-end",   1, NULL, 't' },
    { "verbose", 1, NULL, 'v' },
    { "version", 1, NULL, 'V' },
    { 0,         0, 0,    0   } 
};
#define GETOPTS "d:hi:o:t:v:Vx:y:"

int main(int argc, char *argv[])
{
    int verbose = 1;          /* Verbosity level */
    float xlength = 22.0;     /* Width of simulated domain */
    float ylength = 4.1;      /* Height of simulated domain */
    int imax = 660;           /* Number of cells horizontally */
    int jmax = 120;           /* Number of cells vertically */

    char *infile;             /* Input raw initial conditions */
    char *outfile;            /* Output raw simulation results */

    float t_end = 2.1;        /* Simulation runtime */
    float del_t = 0.003;      /* Duration of each timestep */
    float tau = 0.5;          /* Safety factor for timestep control */

    int itermax = 100;        /* Maximum number of iterations in SOR */
    float eps = 0.001;        /* Stopping error threshold for SOR */
    float omega = 1.7;        /* Relaxation parameter for SOR */
    float gamma = 0.9;        /* Upwind differencing factor in PDE
                                 discretisation */

    float Re = 150.0;         /* Reynolds number */
    float ui = 1.0;           /* Initial X velocity */
    float vi = 0.0;           /* Initial Y velocity */

    float t, delx, dely;
    int  i, j, itersor = 0, ifluid = 0, ibound = 0;
    float res;
    float **u, **v, **p, **rhs, **f, **g;
    char  **flag;
    int init_case, iters = 0;
    int show_help = 0, show_usage = 0, show_version = 0;

    progname = argv[0];
    infile = strdup("karman.bin");
    outfile = strdup("karman.bin");
	
	int size, rank, tag = 0;
	
	MPI_Status stat; 

    /* all MPI programs start with MPI_Init */
    MPI_Init(&argc,&argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	double runningTime;
	if (rank == 0) runningTime = MPI_Wtime();

    int optc;
    while ((optc = getopt_long(argc, argv, GETOPTS, long_opts, NULL)) != -1) {
        switch (optc) {
            case 'h':
                show_help = 1;
                break;
            case 'V':
                show_version = 1;
                break;
            case 'v':
                verbose = atoi(optarg);
                break;
            case 'x':
                imax = atoi(optarg);
                break;
            case 'y':
                jmax = atoi(optarg);
                break;
            case 'i':
                free(infile);
                infile = strdup(optarg);
                break;
            case 'o':
                free(outfile);
                outfile = strdup(optarg);
                break;
            case 'd':
                del_t = atof(optarg);
                break;
            case 't':
                t_end = atof(optarg);
                break;
            default:
                show_usage = 1;
        }
    }
    if (show_usage || optind < argc) {
        print_usage();
        return 1;
    }
    
    if (show_version) {
        print_version();
        if (!show_help) {
            return 0;
        }
    }
    
    if (show_help) {
        print_help();
        return 0;
    }

    delx = xlength/imax;
    dely = ylength/jmax;

	/*
		There are two approaches to the start of this program.
		1. The root process (0) can read in the file, and communicate all of the necessary data to each node.
		2. Each process reads the file, and only retains information for the area they need.
			a. But can multiple MPI processes read the same file simultaneously. 
	*/
	
	/* This method uses vertical division between the seperate nodes */
	int imaxNode = (imax/size);
	
	// Due to a rounding down error with the int values, the first node
	// may need to make up the extra data space
	
	// look at that, this stays at imax for 1 node
	int imaxPrimary = imax -(imaxNode*(size-1)); 
			
	// Needed for moving the actual values to each of the nodes
	// Keeping this available so they can be written into at the end of the process
	// Can then be used to write to the file
	float **uTemp, **vTemp, **pTemp;
	char  **flagTemp;
	
	if (rank==0) {
		/* Allocate arrays */
		u    = alloc_floatmatrix(imaxPrimary+2, jmax+2);
		v    = alloc_floatmatrix(imaxPrimary+2, jmax+2);
		f    = alloc_floatmatrix(imaxPrimary+2, jmax+2);
		g    = alloc_floatmatrix(imaxPrimary+2, jmax+2);
		p    = alloc_floatmatrix(imaxPrimary+2, jmax+2);
		rhs  = alloc_floatmatrix(imaxPrimary+2, jmax+2); 
		flag = alloc_charmatrix(imaxPrimary+2, jmax+2);
		
		uTemp    = alloc_floatmatrix(imax+2, jmax+2);
		vTemp    = alloc_floatmatrix(imax+2, jmax+2);
		pTemp    = alloc_floatmatrix(imax+2, jmax+2);
		flagTemp = alloc_charmatrix(imax+2, jmax+2);
		
		// Used for sending to each node only their values
		float **uNode, **vNode, **pNode;
		char  **flagNode;
		
		uNode    = alloc_floatmatrix(imaxNode+2, jmax+2);
		vNode    = alloc_floatmatrix(imaxNode+2, jmax+2);
		pNode    = alloc_floatmatrix(imaxNode+2, jmax+2);
		flagNode = alloc_charmatrix(imaxNode+2, jmax+2);
		
		if (!u || !v || !f || !g || !p || !rhs || !flag 
		|| !uNode || !vNode || !pNode || !flagNode
		|| !uTemp || !vTemp || !pTemp || !flagTemp) {
			// Let the other nodes know that there we were not successful at allocating memory
			fprintf(stderr, "Couldn't allocate memory for matrices.\n");
			int fail[1] = {1};
			for (i = 1; i < size; i++) {
				// Sends the value 1 to each 
				MPI_Send(fail, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
			}
			MPI_Finalize(); 
			return 1;
		}

		/* Read in initial values from a file if it exists */
		init_case = read_bin(uTemp, vTemp, pTemp, flagTemp, imax, jmax, xlength, ylength, infile);
			
		if (init_case > 0) {
			/* Error while reading file */
			// Let the other nodes know that there we were not successful at reading
			int fail[1] = {1};
			for (i = 1; i < size; i++) {
				// Sends the value 1 to each 
				MPI_Send(fail, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
			}
			MPI_Finalize(); 
			return 1;
		}

		if (init_case < 0) {
			/* Set initial values if file doesn't exist */
			for (i=0;i<=imax+1;i++) {
				for (j=0;j<=jmax+1;j++) {
					uTemp[i][j] = ui;
					vTemp[i][j] = vi;
					pTemp[i][j] = 0.0;
				}
			}
			init_flag(flagTemp, imax, jmax, delx, dely, &ibound);
			// If we set size to 0, no communication will be done
			apply_boundary_conditions(uTemp, vTemp, pTemp, flagTemp, imax, jmax, ui, vi, rank, 0);
		}
		
		//printf("Root node has completed read. Starting handshake to %d nodes\n", size-1);
		
		// Now that these values have been set, we need to split them up for the
		// different nodes
		int node;
		int success[1] = {0}; // 0 always means things are okay!
		// Need confirmation from each other node that they succesfully allocated memory
		for (node = 1; node < size; node++) {
			MPI_Send(success, 1, MPI_INT, node, tag, MPI_COMM_WORLD);
		}
		printf("Root node has sent confirmation to all nodes\n");
		// Now get receives
		int failure = 0;
		for (node = 1; node < size; node++) {
			MPI_Recv(success, 1, MPI_INT, node, tag, MPI_COMM_WORLD, &stat);
			if (success[0] == 1) failure++;
		}
		// If we have reached end, with no failures, complete the handshake and get on with things!
		if (failure > 0) {
			success[0] = failure;
			printf("%d nodes have been unable to allocate resources\n", failure);
		}
		
		for (node = 1; node < size; node++) {
			MPI_Send(success, 1, MPI_INT, node, tag, MPI_COMM_WORLD);
		}
		
		if (failure > 0) {
			// Have already notified everyone that this broke
			MPI_Finalize();
			return 1;
		}
		
		//printf("Root node has finished the handshake with all nodes\n");
		
		
		for (node = 1; node<size; node++) {
			// now need to construct the array specifically for the ith MPI node
			for (i = 0; i <= imaxNode+1; i++) {
				int pos = imaxPrimary + ((node-1)*imaxNode) + i;
				for (j = 0; j <= jmax+1; j++) {
					uNode[i][j] 	= uTemp[pos][j];
					vNode[i][j] 	= vTemp[pos][j];
					pNode[i][j] 	= pTemp[pos][j];
					flagNode[i][j] 	= flagTemp[pos][j];
				}
			}
			// Send the four necessary arrays
			//printf("Root is now sending arrays to node %d\n",node);
			for (i=0; i <= imaxNode+1; i++) {
				MPI_Send(uNode[i], jmax+2, MPI_FLOAT, node, tag, MPI_COMM_WORLD);
				MPI_Send(vNode[i], jmax+2, MPI_FLOAT, node, tag, MPI_COMM_WORLD);
				MPI_Send(pNode[i], jmax+2, MPI_FLOAT, node, tag, MPI_COMM_WORLD);
				MPI_Send(flagNode[i], jmax+2, MPI_CHAR, node, tag, MPI_COMM_WORLD);
				// printf("Root has sent round %d of %d\n", i+1, imaxNode+2); // Debug
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		
		// Finally, fill in our own array
		for (i = 0; i <= imaxPrimary+1; i++) {
			for (j = 0; j <= jmax+1; j++) {
				u[i][j] 	= uTemp[i][j];
				v[i][j] 	= vTemp[i][j];
				p[i][j] 	= pTemp[i][j];
				flag[i][j] 	= flagTemp[i][j];
			}
		}
		
		free_matrix(uNode);
		free_matrix(vNode);
		free_matrix(pNode);
		free_matrix(flagNode);
		
	} else { // What to do if you are not the root node
		/* Allocate arrays */
		u    = alloc_floatmatrix(imaxNode+2, jmax+2);
		v    = alloc_floatmatrix(imaxNode+2, jmax+2);
		f    = alloc_floatmatrix(imaxNode+2, jmax+2);
		g    = alloc_floatmatrix(imaxNode+2, jmax+2);
		p    = alloc_floatmatrix(imaxNode+2, jmax+2);
		rhs  = alloc_floatmatrix(imaxNode+2, jmax+2); 
		flag = alloc_charmatrix(imaxNode+2, jmax+2);
		
		int success[1] = {0}; // 0 Means things are okay!
		
		// Need to wait on success message from root
		MPI_Recv(success, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &stat);
		if (success[0] > 0) {
			// Problem on the root, ending
			MPI_Finalize();
			return 1;
		}
		
		//printf("Node %d has started the handshake with the root node\n", rank);
		
		if (!u || !v || !f || !g || !p || !rhs || !flag) {
			// Let the other nodes know that there we were not successful at allocating memory
			fprintf(stderr, "Couldn't allocate memory for matrices on node %d.\n", rank);
			success[0] = 1;
			MPI_Send(success, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
			MPI_Recv(success, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &stat);
			MPI_Finalize();
			return 1;
		} else {
			MPI_Send(success, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
			MPI_Recv(success, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &stat);
			if (success[0] > 0) {
				MPI_Finalize();
				return 1;
			}
		}
		
		//printf("Node %d is still active having finished the handshake\n", rank);
		
		// Reaching this point means the handshake has been completed
		//need to loop through and get the columns separately!
		for (i=0; i <= imaxNode+1; i++) {
			MPI_Recv(u[i], jmax+2, MPI_FLOAT, 0, tag, MPI_COMM_WORLD, &stat);
			MPI_Recv(v[i], jmax+2, MPI_FLOAT, 0, tag, MPI_COMM_WORLD, &stat);
			MPI_Recv(p[i], jmax+2, MPI_FLOAT, 0, tag, MPI_COMM_WORLD, &stat);
			MPI_Recv(flag[i], jmax+2, MPI_CHAR, 0, tag, MPI_COMM_WORLD, &stat);
			// Debug line
			// printf("Node %d successfully received round %d of %d array values\n",rank, i+1, imaxNode+2);
		}
		
		MPI_Barrier(MPI_COMM_WORLD);
	}
	
	//printf("Node %d has completed the handshake and is ready to start processing\n",rank);

	int imaxLocal;
	// Make use of a local variable for the width of the calculated area
	if (rank == 0)
		imaxLocal = imaxPrimary;
	else
		imaxLocal = imaxNode;
	
	// TODO need to calculate what the iStartPos is for each of the ranks
	int iStartPos = 0;
	for (i = 0; i < rank; i++) {
		// Need to add the imaxPrimary val
		if (i == 0) iStartPos += imaxPrimary;
		else iStartPos += imaxNode;
	}
	
	double intervalTimer = 0.0;
	double tentVelTimer = 0.0;
	double rhsTimer = 0.0;
	double poissonTimer = 0.0;
	double velTimer = 0.0;
	double boundaryTimer = 0.0;
	double tempTimer;
	int poissonCalcs = 0;
	
	// Most of the MPI handling for this gets dealt with in simulation
	/* Main loop */
	for (t = 0.0; t < t_end; t += del_t, iters++) {
		if (rank == 0) tempTimer = MPI_Wtime();
		set_timestep_interval(&del_t, imaxLocal, jmax, delx, dely, u, v, Re, tau, rank, size);
		if (rank == 0) {
			tempTimer = MPI_Wtime() - tempTimer;
			intervalTimer += tempTimer;
		}

		ifluid = (imax * jmax) - ibound;
		
		if (rank == 0) tempTimer = MPI_Wtime();
		compute_tentative_velocity(u, v, f, g, flag, imaxLocal, jmax,
			del_t, delx, dely, gamma, Re, rank, size);
		if (rank == 0) {
			tempTimer = MPI_Wtime() - tempTimer;
			tentVelTimer += tempTimer;
		}
		
		
		if (rank == 0) tempTimer = MPI_Wtime();
		compute_rhs(f, g, rhs, flag, imaxLocal, jmax, del_t, delx, dely);
		if (rank == 0) {
			tempTimer = MPI_Wtime() - tempTimer;
			rhsTimer += tempTimer;
		}
		
		if (ifluid > 0) {
			if (rank == 0) tempTimer = MPI_Wtime();
			//printf("Node %d is performing the poisson calculation\n", rank);
			itersor = poisson(p, rhs, flag, imaxLocal, jmax, delx, dely,
						eps, itermax, omega, &res, ifluid, rank, size, iStartPos);
			if (rank == 0) {
				tempTimer = MPI_Wtime() - tempTimer;
				poissonTimer += tempTimer;
				poissonCalcs++;
			}
		} else {
			itersor = 0;
			//printf("Node %d has an ifluid value of %d\n",rank,ifluid);
		}

		if (proc == 0 && verbose > 1) {
			printf("%d t:%g, del_t:%g, SOR iters:%3d, res:%e, bcells:%d\n",
				iters, t+del_t, del_t, itersor, res, ibound);
		}
		
		if (rank == 0) tempTimer = MPI_Wtime();
		update_velocity(u, v, f, g, p, flag, imaxLocal, jmax, del_t, delx, dely, rank, size);
		if (rank == 0) {
			tempTimer = MPI_Wtime() - tempTimer;
			velTimer += tempTimer;
		}
		
		if (rank == 0) tempTimer = MPI_Wtime();
		apply_boundary_conditions(u, v, p, flag, imaxLocal, jmax, ui, vi, rank, size);
		if (rank == 0) {
			tempTimer = MPI_Wtime() - tempTimer;
			boundaryTimer += tempTimer;
		}
	} /* End of main loop */
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	if (rank == 0) {
		printf("Program took %d iterations across %d MPI nodes\n", iters, size);
		// Let's display some running time information
		printf("Average time for set_timestep_interval      : %f\n", intervalTimer / (double)iters);
		printf("Average time for compute_tentative_velocity : %f\n", tentVelTimer / (double)iters);
		printf("Average time for compute_rhs                : %f\n", rhsTimer / (double)iters);
		printf("Average time for poisson                    : %f\n", poissonTimer / (double)poissonCalcs);
		printf("Average time for update_velocity            : %f\n", velTimer / (double)iters);
		printf("Average time for apply_boundary_conditions  : %f\n", boundaryTimer / (double)iters);
	}
	//printf("Node %d has completed the main loop\n", rank);
	
	// Do a (maybe unnecessary) check so every node is at the same point
	MPI_Barrier(MPI_COMM_WORLD);
	
	int node;
	
	// Need to collate all of the data
	if (rank == 0) {
		// Make sure to write the root information to the larger temp array as well
		for (i=0; i <= imaxPrimary; i++) {
			for (j=0; j < jmax+2; j++) {			
				uTemp[i][j] = u[i][j];
				vTemp[i][j] = v[i][j];
				pTemp[i][j] = p[i][j];
				flagTemp[i][j] = flag[i][j];	
			}
		}
		for (node = 1; node < size; node++) {
			// if this is the last node, we need to transfer the far edge again
			int imaxLocal = imaxNode;
			if (node == size-1) imaxLocal++;
			for (i=1; i <= imaxLocal; i++) {
				int pos = imaxPrimary + ((node-1)*imaxNode) + i;
				MPI_Recv(uTemp[pos], jmax+2, MPI_FLOAT, node, tag, MPI_COMM_WORLD, &stat);
				MPI_Recv(vTemp[pos], jmax+2, MPI_FLOAT, node, tag, MPI_COMM_WORLD, &stat);
				MPI_Recv(pTemp[pos], jmax+2, MPI_FLOAT, node, tag, MPI_COMM_WORLD, &stat);
				MPI_Recv(flagTemp[pos], jmax+2, MPI_CHAR, node, tag, MPI_COMM_WORLD, &stat);
				// printf("Root has received round %d of %d\n", i, imaxNode+1); // Debug
			}
			/*
			 * There are additional values in the far right boundary of the final node
			 * however these are never modified, so they do not need to be retransferred
			 */
		}
	} else {
		// Send the information that is within out vertical bounds
		int imaxLocal = imaxNode;
		if (rank == size-1) imaxLocal++;
		for (i=1; i <= imaxLocal; i++) {
			MPI_Send(u[i], jmax+2, MPI_FLOAT, 0, tag, MPI_COMM_WORLD);
			MPI_Send(v[i], jmax+2, MPI_FLOAT, 0, tag, MPI_COMM_WORLD);
			MPI_Send(p[i], jmax+2, MPI_FLOAT, 0, tag, MPI_COMM_WORLD);
			MPI_Send(flag[i], jmax+2, MPI_CHAR, 0, tag, MPI_COMM_WORLD);
		}
		//printf("Node %d has send their arrays back to root\n", rank);
	}
  
	MPI_Barrier(MPI_COMM_WORLD);
	
	// Only the first node can write the file
    if (outfile != NULL && strcmp(outfile, "") != 0 && proc == 0 && rank == 0) {
        write_bin(uTemp, vTemp, pTemp, flagTemp, imax, jmax, xlength, ylength, outfile);
		printf("Root node has written the file\n");
    }

    free_matrix(u);
    free_matrix(v);
    free_matrix(f);
    free_matrix(g);
    free_matrix(p);
    free_matrix(rhs);
    free_matrix(flag);
	
	//printf("Node %d has finished\n", rank);
	
	// If this is the root node, make sure to free the temp arrays
	if (rank == 0) {
		free_matrix(uTemp);
		free_matrix(vTemp);
		free_matrix(pTemp);
		free_matrix(flagTemp);
	}
	
	if (rank == 0) {
		runningTime = MPI_Wtime() - runningTime;
		printf("Running time                                : %f\n", runningTime);
	}
	
    /* MPI Programs end with MPI Finalize; this is a weak synchronization point */
    MPI_Finalize(); 

    return 0;
}

/* Save the simulation state to a file */
void write_bin(float **u, float **v, float **p, char **flag,
    int imax, int jmax, float xlength, float ylength, char* file)
{
    int i;
    FILE *fp;

    fp = fopen(file, "wb"); 

    if (fp == NULL) {
        fprintf(stderr, "Could not open file '%s': %s\n", file,
            strerror(errno));
        return;
    }

    fwrite(&imax, sizeof(int), 1, fp);
    fwrite(&jmax, sizeof(int), 1, fp);
    fwrite(&xlength, sizeof(float), 1, fp);
    fwrite(&ylength, sizeof(float), 1, fp);

    for (i=0;i<imax+2;i++) {
        fwrite(u[i], sizeof(float), jmax+2, fp);
        fwrite(v[i], sizeof(float), jmax+2, fp);
        fwrite(p[i], sizeof(float), jmax+2, fp);
        fwrite(flag[i], sizeof(char), jmax+2, fp);
    }
    fclose(fp);
}

/* Read the simulation state from a file */
int read_bin(float **u, float **v, float **p, char **flag,
    int imax, int jmax, float xlength, float ylength, char* file)
{
    int i,j;
    FILE *fp;

    if (file == NULL) return -1;

    if ((fp = fopen(file, "rb")) == NULL) {
        fprintf(stderr, "Could not open file '%s': %s\n", file,
            strerror(errno));
        fprintf(stderr, "Generating default state instead.\n");
        return -1;
    }

    fread(&i, sizeof(int), 1, fp);
    fread(&j, sizeof(int), 1, fp);
    float xl, yl;
    fread(&xl, sizeof(float), 1, fp);
    fread(&yl, sizeof(float), 1, fp);

    if (i!=imax || j!=jmax) {
        fprintf(stderr, "Warning: imax/jmax have wrong values in %s\n", file);
        fprintf(stderr, "%s's imax = %d, jmax = %d\n", file, i, j);
        fprintf(stderr, "Program's imax = %d, jmax = %d\n", imax, jmax);
        return 1;
    }
    if (xl!=xlength || yl!=ylength) {
        fprintf(stderr, "Warning: xlength/ylength have wrong values in %s\n", file);
        fprintf(stderr, "%s's xlength = %g,  ylength = %g\n", file, xl, yl);
        fprintf(stderr, "Program's xlength = %g, ylength = %g\n", xlength,
            ylength);
        return 1;
    }

    for (i=0; i<imax+2; i++) {
        fread(u[i], sizeof(float), jmax+2, fp);
        fread(v[i], sizeof(float), jmax+2, fp);
        fread(p[i], sizeof(float), jmax+2, fp);
        fread(flag[i], sizeof(char), jmax+2, fp);
    }
    fclose(fp);
    return 0;
}

static void print_usage(void)
{
    fprintf(stderr, "Try '%s --help' for more information.\n", progname);
}

static void print_version(void)
{
    fprintf(stderr, "%s %s\n", PACKAGE, VERSION);
}

static void print_help(void)
{
    fprintf(stderr, "%s. A simple computational fluid dynamics tutorial.\n\n",
        PACKAGE);
    fprintf(stderr, "Usage: %s [OPTIONS]...\n\n", progname);
    fprintf(stderr, "  -h, --help            Print a summary of the options\n");
    fprintf(stderr, "  -V, --version         Print the version number\n");
    fprintf(stderr, "  -v, --verbose=LEVEL   Set the verbosity level. 0 is silent\n");
    fprintf(stderr, "  -x, --imax=IMAX       Set the number of interior cells in the X direction\n");
    fprintf(stderr, "  -y, --jmax=JMAX       Set the number of interior cells in the Y direction\n");
    fprintf(stderr, "  -t, --t-end=TEND      Set the simulation end time\n");
    fprintf(stderr, "  -d, --del-t=DELT      Set the simulation timestep size\n");
    fprintf(stderr, "  -i, --infile=FILE     Read the initial simulation state from this file\n");
    fprintf(stderr, "                        (default is 'karman.bin')\n");
    fprintf(stderr, "  -o, --outfile=FILE    Write the final simulation state to this file\n");
    fprintf(stderr, "                        (default is 'karman.bin')\n");
}