#include "Driver.h"

#include <iostream>
#include <omp.h>

    Driver::Driver(const InputFile* input, const std::string& pname)
: problem_name(pname)
{

    std::cout << "+++++++++++++++++++++" << std::endl;
    std::cout << "  Running deqn v0.1  " << std::endl;
#ifdef DEBUG
    std::cout << "- input file: " << problem_name << std::endl;
#endif

    dt_max = input->getDouble("dt_max",  0.2);
    dt = input->getDouble("initial_dt", 0.2);

    t_start = input->getDouble("start_time", 0.0);
    t_end = input->getDouble("end_time", 10.0);

    vis_frequency = input->getInt("vis_frequency",-1);
    summary_frequency = input->getInt("summary_frequency", 1);

#ifdef DEBUG
    std::cout << "- dt_max: " << dt_max << std::endl;
    std::cout << "- initial_dt: " << dt << std::endl;
    std::cout << "- start_time: " << t_start << std::endl;
    std::cout << "- end_time: " << t_end << std::endl;
    std::cout << "- vis_frequency: " << vis_frequency << std::endl;
    std::cout << "- summary_frequency: " << summary_frequency << std::endl;
#endif
    std::cout << "+++++++++++++++++++++" << std::endl;
    std::cout << std::endl;

    mesh = new Mesh(input);

    diffusion = new Diffusion(input, mesh);
    writer = new VtkWriter(pname, mesh);

    /* Initial mesh dump */
    if(vis_frequency != -1)
        writer->write(0, 0.0);
}

Driver::~Driver() {
    delete mesh;
    delete diffusion;
    delete writer;
}

void Driver::run() {

    int step = 0;
    double t_current;
    double startTime = omp_get_wtime();

    /*	Can't parallelise this loop since there is a depenency between each iteration, with a call to the diffusion
    **	Could seperate the loop out, making a seperate diffusion for each iteration where a single cycle is performed
    **	but this could ruin the iterative nature of the execution.
    */
    for(t_current = t_start; t_current < t_end; t_current += dt) {
        step = t_current/dt + 1;

        std::cout << "+ step: " << step << ", dt:   " << dt << std::endl;

        diffusion->doCycle(dt);

        if(step % vis_frequency == 0 && vis_frequency != -1)
            writer->write(step, t_current);
        if(step % summary_frequency == 0 && summary_frequency != -1) {
            double temperature = mesh->getTotalTemperature();
            std::cout << "+\tcurrent total temperature: " << temperature << std::endl;
        }

    }

    if(step % vis_frequency != 0 && vis_frequency != -1)
        writer->write(step, t_current);

    std::cout << std::endl;
    std::cout << "+++++++++++++++++++++" << std::endl;
    std::cout << "   Run complete.   " << std::endl;
    std::cout << "+++++++++++++++++++++" << std::endl;
    std::cout << "Total time taken: " << omp_get_wtime()-startTime << std::endl;
}
