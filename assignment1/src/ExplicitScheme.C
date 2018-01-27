#include "ExplicitScheme.h"

#include <iostream>

#define POLY2(i, j, imin, jmin, ni) (((i) - (imin)) + (((j)-(jmin)) * (ni)))

ExplicitScheme::ExplicitScheme(const InputFile* input, Mesh* m) :
    mesh(m)
{
    int nx = mesh->getNx()[0];
    int ny = mesh->getNx()[1];
}

void ExplicitScheme::doAdvance(const double dt)
{
    diffuse(dt);

    reset();

    updateBoundaries();
}

void ExplicitScheme::updateBoundaries()
{
    for (int i = 0; i < 4; i++) {
        reflectBoundaries(i);
    }
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

    for(int k = y_min-1; k <= y_max+1; k++) {
        for(int j = x_min-1; j <=  x_max+1; j++) {
            int i = POLY2(j,k,x_min-1,y_min-1,nx);
            u0[i] = u1[i];
        }
    }
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
    
    #pragma omp parallel for
    for(int k=y_min; k <= y_max; k++) {
	#pragma omp parallel for
        for(int j=x_min; j <= x_max; j++) {

            int n1 = POLY2(j,k,x_min-1,y_min-1,nx);
            int n2 = POLY2(j-1,k,x_min-1,y_min-1,nx);
            int n3 = POLY2(j+1,k,x_min-1,y_min-1,nx);
            int n4 = POLY2(j,k-1,x_min-1,y_min-1,nx);
            int n5 = POLY2(j,k+1,x_min-1,y_min-1,nx);

            u1[n1] = (1.0-2.0*rx-2.0*ry)*u0[n1] + rx*u0[n2] + rx*u0[n3]
                + ry*u0[n4] + ry*u0[n5];
        }
    }
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
                for(int j = x_min; j <= x_max; j++) {
                    int n1 = POLY2(j, y_max, x_min-1, y_min-1, nx);
                    int n2 = POLY2(j, y_max+1, x_min-1, y_min-1, nx);

                    u0[n2] = u0[n1];
                }
            } break;
        case 1:
            /* right */
            {
                for(int k = y_min; k <= y_max; k++) {
                    int n1 = POLY2(x_max, k, x_min-1, y_min-1, nx);
                    int n2 = POLY2(x_max+1, k, x_min-1, y_min-1, nx);

                    u0[n2] = u0[n1];
                }
            } break;
        case 2: 
            /* bottom */
            {
                for(int j = x_min; j <= x_max; j++) {
                    int n1 = POLY2(j, y_min, x_min-1, y_min-1, nx);
                    int n2 = POLY2(j, y_min-1, x_min-1, y_min-1, nx);

                    u0[n2] = u0[n1];
                }
            } break;
        case 3: 
            /* left */
            {
                for(int k = y_min; k <= y_max; k++) {
                    int n1 = POLY2(x_min, k, x_min-1, y_min-1, nx);
                    int n2 = POLY2(x_min-1, k, x_min-1, y_min-1, nx);

                    u0[n2] = u0[n1];
                }
            } break;
        default: std::cerr << "Error in reflectBoundaries(): unknown boundary id (" << boundary_id << ")" << std::endl;
    }
}
