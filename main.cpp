#include <iostream>
#include "variable.h"
#include <fstream>

using namespace std;

int main()
{
    int t;

    arq='result.dat';
    std::ofstream out;
    out.open("arq", ios::out | ios::app);
    //out << 4x,"Beta",4x,"Energy",4x,"Specific Heat",4x,"Uniform Susceptibility",  4x,"Magnetization^2",4x,"Number of H-operators", 4x,"(Number of H-operators)^2";

    for (Nn=0; Nn<100; Nn++){
        for(t=0; t<t_steps;t++){
            beta = 1.0 / (temp_ini + dt * t);
            NH = 0;
            NH_max = 0;
            L = std::max(4, N / 4);

            Allocate(opstring[0:L - 1])
            Allocate(first_vertex_visitted[N], last_vertex_visitted[N])
            Allocate(vertex_link[0:4 * L - 1])

            opstring = 0;
            lattice();
            init();

            add_prob = 0.5 * beta * Nb;
            rem_prob = 1.0 / (0.5 * beta * Nb);

            termalization();
            qmc_steps();
            results();
            make_dist();
            free_memory();
        }
    }
    out.close();

    return 0;
}

void lattice()
{

}
