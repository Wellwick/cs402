#ifndef DIFFUSION_MESH_H_
#define DIFFUSION_MESH_H_

#include "InputFile.h"

class Mesh {
    private:
        const InputFile* input;

        double* u1;
        double* u0;
        double* cellx;
        double* celly;

        double* min_coords;
        double* max_coords;

        int NDIM;

        int* n; 
        int* min;
        int* max;

        double* dx;

        /*
         * A mesh has four neighbours, and they are 
         * accessed in the following order:
         * - top
         * - right
         * - bottom
         * - left
         */
        int* neighbours;

        void allocate();
        bool allocated;
    public:
        Mesh(const InputFile* input);
                
        double* getU0();
        double* getU1();

        double* getDx();
        int* getNx();
        int* getMin();
        int* getMax();
        int getDim();

        double* getCellX();
        double* getCellY();

        int* getNeighbours();

        double getTotalTemperature();
};
#endif
