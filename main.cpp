#include <iostream>
#include "variable.h"
#include <fstream>

using namespace std;

void lattice()
{
    int i, j, k, bound_idx, spin_idx

    // mod(x, L), in this case, returns 0 if x == Lx or returns x if x != Lx.

    if ((lx > 1 && ly > 1 && lz == 1)){

        print*, '2D with Periodic boundary '

        d = 2
        Nb = d * lx*ly
        Allocate(bound_spin_ID(2, Nb))
        for (i=0; i<lx; i++){
            for (j=0; j<ly; j++){
                spin_idx = 1 + i + j * lx

                bound_idx = 1 + i + j * lx

                bound_spin_ID(1, bound_idx) = spin_idx
                bound_spin_ID(2, bound_idx) = 1 + mod(i + 1, lx) + j * lx

                bound_spin_ID(1, bound_idx + lx * ly) = spin_idx
                bound_spin_ID(2, bound_idx + lx * ly) = 1 + i + mod(j + 1, ly) * lx
            }
        }
    }elseif ( (lx > 1) .and. (ly == 1 .and. lz == 1) ){

    print*, '1D with Periodic boundary '

    d = 1
    Nb = d * lx
    Allocate(bound_spin_ID(2, Nb))

    Do i = 0, lx - 1
       spin_idx = i + 1

       bound_idx = i + 1

       bound_spin_ID(1, bound_idx) = spin_idx
       bound_spin_ID(2, bound_idx) = 1 + mod(i + 1, lx)
    end do

    }elseif (lx > 1 .and. ly > 1 .and. lz > 1) then

    print*, '3D with Periodic boundary '

    d = 3
    Nb = d * lx * ly * lz
    Allocate(bound_spin_ID(2, Nb))

    Do i = 0, lx - 1
      Do j = 0, ly - 1
        Do k = 0, lz - 1

          spin_idx = 1 + i + j * lx + k * lx * ly

          bound_idx = 1 + i + j * lx + k * lx * ly

          bound_spin_ID(1, bound_idx) = spin_idx
          bound_spin_ID(2, bound_idx) = 1 + mod(i + 1, lx) + j * lx + k * lx * ly

          bound_spin_ID(1, bound_idx + lz * lx * ly) = spin_idx
          bound_spin_ID(2, bound_idx + lz * lx * ly) = 1 + i + mod(j + 1, ly) * lx + k * lx * ly

          bound_spin_ID(1, bound_idx + 2 * lz * lx * ly) = spin_idx
          bound_spin_ID(2, bound_idx + 2 * lz * lx * ly) = 1 + i + j * lx + mod(k + 1, lz) * lx * ly

        end do
      end do
    end do

    else

    print*, 'For a 1D lattice: lx > 0, ly==1 and lz==1'
    print*, 'For a 2D lattice: lx > 0, ly > 0 and lz==1'
    print*, 'For a 3D lattice: lx > 0, ly > 0 and lz > 0'
    stop
    end if
}

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
