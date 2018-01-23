#include "Diffusion.h"

#include "ExplicitScheme.h"

#include <iostream>
#include <cstdlib>

Diffusion::Diffusion(const InputFile* input, Mesh* m) :
    mesh(m) 
{

    std::string scheme_str = input->getString("scheme", "explicit");

    if(scheme_str.compare("explicit") == 0) {
        scheme = new ExplicitScheme(input, mesh);
    } else {
        std::cerr << "Error: unknown scheme \"" << scheme_str << "\"" << std::endl;
        exit(1);
    }

    subregion = input->getDoubleList("subregion", std::vector<double>());

    if (subregion.size() != 0 && subregion.size() != 4) {
        std::cerr << "Error:  region must have 4 entries (xmin, ymin, xmax, ymax)" << std::endl;
        exit(1);
    }

    init();
}

Diffusion::~Diffusion()
{
    delete scheme;
}

void Diffusion::init()
{
    double* u0 = mesh->getU0();

    int x_max = mesh->getNx()[0];
    int y_max = mesh->getNx()[1];

    double* cellx = mesh->getCellX();
    double* celly = mesh->getCellY();

    int nx = x_max+2;

    if(!subregion.empty()) {
        for (int j = 0; j < y_max+2; j++) {
            for (int i = 0; i < x_max+2; i++) {
                if (celly[j] > subregion[1] && celly[j] <= subregion[3] &&
                        cellx[i] > subregion[0] && cellx[i] <= subregion[2]) {
                    u0[i+j*nx] = 10.0;
                } else {
                    u0[i+j*nx] = 0.0;
                }

            }
        }
    } else {
        for (int j = 0; j < y_max+2; j++) {
            for (int i = 0; i < x_max+2; i++) {
                u0[i+j*nx] = 0.0;
            }
        }
    }

    scheme->init();
}

void Diffusion::doCycle(const double dt)
{
    scheme->doAdvance(dt);
}
