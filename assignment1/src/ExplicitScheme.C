#include "ExplicitScheme.h"

#include <iostream>
#include <omp.h>

#define POLY2(i, j, imin, jmin, ni) (((i) - (imin)) + (((j)-(jmin)) * (ni)))

static double diffuseTimer = 0.0;
static double resetTimer = 0.0;
static double boundaryTimer = 0.0;
static double stepCounter = 0.0;

ExplicitScheme::ExplicitScheme(const InputFile* input, Mesh* m) :
    mesh(m)
{
    int nx = mesh->getNx()[0];
    int ny = mesh->getNx()[1];
}

void ExplicitScheme::doAdvance(const double dt)
{
    stepCounter++;
    diffuse(dt);

    reset();

    updateBoundaries();
}

void ExplicitScheme::updateBoundaries()
{
    double startTime = omp_get_wtime();
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < 4; i++) {
        reflectBoundaries(i);
    }
    std::cout << "+\tBoundary loops time taken: " << omp_get_wtime()-startTime << std::endl;
    boundaryTimer += omp_get_wtime()-startTime;
}

void ExplicitScheme::init()
{
    updateBoundaries();
}

void ExplicitScheme::reset()
{
    double* u0 = mesh->getU0();
    double* u1 = mesh->getU1();
    int x_min = mesh->getMin()[0];
    int x_max = mesh->getMax()[0];
    int y_min = mesh->getMin()[1]; 
    int y_max = mesh->getMax()[1]; 

    int nx = mesh->getNx()[0]+2;

    
    /* 
     * Tried both collapse and simd
     * These didn't work
    */ 
    
    int chunk = (x_max-x_min+1)/omp_get_num_procs();
    //despite this being what static schedule should use, using this produces speedup sometimes
    //int j;
    double startTime = omp_get_wtime();
    #pragma omp parallel for schedule(static) shared(u0, u1, x_min, y_min, x_max, y_max, nx)
    for(int k = y_min-1; k <= y_max+1; k++) {
	//#pragma omp for schedule(static, chunk)
        for(int j = x_min-1; j <=  x_max+1; j++) {
            int i = POLY2(j,k,x_min-1,y_min-1,nx);
            u0[i] = u1[i];
        }
    }
    #pragma omp barrier
    std::cout << "+\tReset loop time taken: " << omp_get_wtime()-startTime << std::endl;
    resetTimer += omp_get_wtime()-startTime;
}

void ExplicitScheme::diffuse(double dt)
{
    double* u0 = mesh->getU0();
    double* u1 = mesh->getU1();
    int x_min = mesh->getMin()[0];
    int x_max = mesh->getMax()[0];
    int y_min = mesh->getMin()[1]; 
    int y_max = mesh->getMax()[1]; 
    double dx = mesh->getDx()[0];
    double dy = mesh->getDx()[1];

    int nx = mesh->getNx()[0]+2;

    double rx = dt/(dx*dx);
    double ry = dt/(dy*dy);
    
    //put this outside so we only calculate it once
    double rXY = (1.0-2.0*rx-2.0*ry);
    
    /* Definitely no to SIMD, despite the fact this is working with multiple data, single instruction
     * there is no speedup!
     */
    
    //int k;
    int chunk = (x_max-x_min+1)/omp_get_num_procs();
    //despite this being what static schedule should use, using this produces speedup sometimes
    //int j;
    double startTime = omp_get_wtime();
    #pragma omp parallel shared(u0, u1, nx, rx, ry)
    for(int k=y_min; k <= y_max; k++) {
	//faster to divide threads through inner loop due to spatial locality
	//but only for sufficiently large values of nx
	#pragma omp for schedule(static)
	//#pragma omp simd parallel private(k) schedule(dynamic) shared(u0, u1, nx, rx, ry)
        for(int j=x_min; j <= x_max; j++) {
	    
	    // Can we calculate these once so we don't need to do them every time?
	    // We can! But it's slower >:(
            int n1 = POLY2(j,k,x_min-1,y_min-1,nx);
	    //int n2 = n1-1;
            int n2 = POLY2(j-1,k,x_min-1,y_min-1,nx);
            //int n3 = n1+1;
	    int n3 = POLY2(j+1,k,x_min-1,y_min-1,nx);
	    //int n4 = n1-nx;
            int n4 = POLY2(j,k-1,x_min-1,y_min-1,nx);
	    //int n5 = n1+nx;
            int n5 = POLY2(j,k+1,x_min-1,y_min-1,nx);
	    
	    double v1 = rXY*u0[n1];
	    double v2 = rx*u0[n2];
	    double v3 = rx*u0[n3];
	    double v4 = ry*u0[n4];
	    double v5 = ry*u0[n5];

            u1[n1] = v1 + v2 + v3 + v4 + v5;
        }
    }
    #pragma omp barrier
    std::cout << "+\tDiffuse loop time taken: " << omp_get_wtime()-startTime << std::endl;
    diffuseTimer += omp_get_wtime()-startTime;
}

void ExplicitScheme::reflectBoundaries(int boundary_id)
{
    double* u0 = mesh->getU0();
    int x_min = mesh->getMin()[0];
    int x_max = mesh->getMax()[0];
    int y_min = mesh->getMin()[1]; 
    int y_max = mesh->getMax()[1]; 

    int nx = mesh->getNx()[0]+2;

    switch(boundary_id) {
        case 0: 
            /* top */
            {
		
		int j;
		#pragma omp parallel for private(j) schedule(static)
		for(j = x_min; j <= x_max; j++) {
                    int n1 = POLY2(j, y_max, x_min-1, y_min-1, nx);
                    int n2 = POLY2(j, y_max+1, x_min-1, y_min-1, nx);

                    u0[n2] = u0[n1];
                }
            } break;
        case 1:
            /* right */
            {
		
		int k;
		#pragma omp parallel for private(k) schedule(static)
		for(k = y_min; k <= y_max; k++) {
                    int n1 = POLY2(x_max, k, x_min-1, y_min-1, nx);
                    int n2 = POLY2(x_max+1, k, x_min-1, y_min-1, nx);

                    u0[n2] = u0[n1];
                }
            } break;
        case 2: 
            /* bottom */
            {
		int j;
		#pragma omp parallel for private(j) schedule(static)
		for(j = x_min; j <= x_max; j++) {
                    int n1 = POLY2(j, y_min, x_min-1, y_min-1, nx);
                    int n2 = POLY2(j, y_min-1, x_min-1, y_min-1, nx);

                    u0[n2] = u0[n1];
                }
            } break;
        case 3: 
            /* left */
            {
		int k;
		#pragma omp parallel for private(k) schedule(static)
		for(k = y_min; k <= y_max; k++) {
                    int n1 = POLY2(x_min, k, x_min-1, y_min-1, nx);
                    int n2 = POLY2(x_min-1, k, x_min-1, y_min-1, nx);

                    u0[n2] = u0[n1];
                }
            } break;
        default: std::cerr << "Error in reflectBoundaries(): unknown boundary id (" << boundary_id << ")" << std::endl;
    }
}

void ExplicitScheme::getAverages()
{
    std::cout << "\t\tDiffuse loop average time taken: " << diffuseTimer / stepCounter << std::endl;
    std::cout << "\t\tReset loop average time taken: " << resetTimer / stepCounter << std::endl;
    std::cout << "\t\tBoundary loop average time taken: " << boundaryTimer / (stepCounter+1) << std::endl;
}