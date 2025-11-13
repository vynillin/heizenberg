#ifndef VARIABLE_H_INCLUDED
#define VARIABLE_H_INCLUDED

    //Critical temperature for a 3D Heisenberg model is Tc = 0.945d0.

    //lx=10, ly=1, lz=1 means a 1D chain lattice with size 10.
    // lx=10, ly=10, lz=1 means a 2D square lattice.
    // lx=10, ly=10, lz=10 means a 3D cubic lattice.

    constexpr int lx = 16;                     //Number of spins in x.
    constexpr int ly = 16;                     //Number of spins in y.
    constexpr int lz = 16;                     //Number of spins in z.
    constexpr int Nbins = 20;                  //Averages written to file after mcsteps.
    constexpr int mcsteps = 1e4;               // Monte Carlo steps.
    constexpr int term_steps = 1e4;            // Termalization steps.
    constexpr int N = lx * ly * lz;            // Total number of spins.
    constexpr double temp_ini = 1.05;          // Initial temperature.
    constexpr int t_steps = 1;                 // Number of temperature steps.
    constexpr double dt = -0.025;              // Size of temperature steps.

    int NH;                                    // Number of H-operator.
    int L;                                     // Maximum string size. (Don't change it)
    int d;                                     // Dimension. (Selected by the program)
    int Nb;                                    // Number o bonds between spins. (N = d*lx*ly*lz)
    int Nn;

    int *spin;                                 // Spin value.
    int *opstring;                             // Operator string ID.
    int **bound_spin_ID;                       // List of the spin sites connect between the bound_spin_ID(1, bound_idx) and bound_spin_ID(2, bound_idx)
    int *vertex_link;                          // List of vertex links.
    int *first_vertex_visitted;                // First operator on each site in linked vertex list.
    int *last_vertex_visitted;                 // Last operator on each site in linked vertex list.

    int NH_max;
    char arq;

    double beta;
    double rem_prob;
    double add_prob;
    double n_opH = 0.0, n_opH2 = 0.0;
    double ususc = 0.0, staggered = 0.0;

#endif // VARIABLE_H_INCLUDED
