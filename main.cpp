#include <iostream>
#include "variable.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <random>
#include <vector>
#include <cmath>
#include <limits>

std::mt19937 rng(std::random_device{}());
std::uniform_real_distribution<double> uniform(0.0, 1.0);

void lattice()
{

    if (bound_spin_ID) {
        delete[] bound_spin_ID[0];
        delete[] bound_spin_ID;
        bound_spin_ID = nullptr;
    }

    std::cout << "\nheisenberg_sse.cpp\n";

    int bound_idx, spin_idx;

    // mod(x, L), in this case, returns 0 if x == Lx or returns x if x != Lx.

    if ((lx > 1 && ly > 1) && (lz == 1)){

        std::cout << "2D with Periodic boundary\n";

        d = 2;
        Nb = d * lx*ly;

        bound_spin_ID = new int*[2];
        bound_spin_ID[0] = new int[2 * Nb]();  // всё одним куском памяти
        bound_spin_ID[1] = bound_spin_ID[0] + Nb;


        for (int i=0; i<lx; i++){
            for (int j=0; j<ly; j++){
                spin_idx = 1 + i + j * lx;
                bound_idx = 1 + i + j * lx;

                // x-связь
                bound_spin_ID[0][bound_idx - 1] = spin_idx;
                bound_spin_ID[1][bound_idx - 1] = 1 + ( (i + 1) % lx ) + j * lx;
                // y-связь
                bound_spin_ID[0][bound_idx + lx * ly - 1] = spin_idx;
                bound_spin_ID[1][bound_idx + lx * ly - 1] = 1 + i + ( (j + 1) % ly ) * lx;
            }
        }
    }else if ((lx > 1) && (ly == 1 && lz == 1)){

    std::cout << "1D with Periodic boundary\n";

    d = 1;
    Nb = d * lx;

    bound_spin_ID = new int*[2];
    bound_spin_ID[0] = new int[2 * Nb]();
    bound_spin_ID[1] = bound_spin_ID[0] + Nb;

    for (int i = 0; i < lx; ++i) {
            spin_idx = i + 1;
            bound_idx = i + 1;

            bound_spin_ID[0][bound_idx - 1] = spin_idx;
            bound_spin_ID[1][bound_idx - 1] = 1 + ( (i + 1) % lx );
        }
    } else if (lx > 1 && ly > 1 && lz > 1) {

    std::cout << "3D with Periodic boundary\n";

    d = 3;
    Nb = d * lx * ly * lz;

    bound_spin_ID = new int*[2];
    bound_spin_ID[0] = new int[2 * Nb]();
    bound_spin_ID[1] = bound_spin_ID[0] + Nb;

    for (int i = 0; i < lx; ++i) {
            for (int j = 0; j < ly; ++j) {
                for (int k = 0; k < lz; ++k) {
                    spin_idx = 1 + i + j * lx + k * lx * ly;
                    bound_idx = 1 + i + j * lx + k * lx * ly;

                    // x-связь
                    bound_spin_ID[0][bound_idx - 1] = spin_idx;
                    bound_spin_ID[1][bound_idx - 1] =
                                    1 + ((i + 1) % lx) + j * lx + k * lx * ly;
                    // y-связь
                    bound_spin_ID[0][bound_idx + lz * lx * ly - 1] = spin_idx;
                    bound_spin_ID[1][bound_idx + lz * lx * ly - 1] =
                                    1 + i + ((j + 1) % ly) * lx + k * lx * ly;
                    // z-связь
                    bound_spin_ID[0][bound_idx + 2 * lz * lx * ly - 1] = spin_idx;
                    bound_spin_ID[1][bound_idx + 2 * lz * lx * ly - 1] =
                                    1 + i + j * lx + ((k + 1) % lz) * lx * ly;
                }
            }
        }
    } else {
        std::cerr << "For a 1D lattice: lx > 0, ly==1 and lz==1\n";
        std::cerr << "For a 2D lattice: lx > 0, ly > 0 and lz==1\n";
        std::cerr << "For a 3D lattice: lx > 0, ly > 0 and lz > 0\n";
        exit(1);
    }
}

void init() {
    spin = new int[N];

    for (int i = 0; i < N; ++i) {
        int x = i % lx;
        int y = (i / lx) % ly;
        int z = i / (lx * ly);
        int exponent = x + y + z;
        spin[i] = (exponent % 2 == 0) ? 1 : -1;
    }

    // Вывод информации о решетке и параметрах
    std::cout << "\nLattice:\tlx=" << lx << "\tly=" << ly << "\tlz=" << lz << "\n";
    std::cout << "Number of spins:\t" << N << "\n";
    std::cout << "Number of bins:\t" << Nbins << "\n";
    std::cout << "Monte Carlo steps:\t" << std::scientific << mcsteps << "\n";
    std::cout << "Termalization steps:\t" << std::scientific << term_steps << "\n";
    std::cout << std::fixed << std::setprecision(4)
              << "Temperatures:\t" << 1.0 / beta << "\t->\t" << temp_ini + (t_steps - 1) * dt << "\n";
    std::cout << "N:\t" << Nn << "\n";

    // Открытие bin-файла и заголовок
    {
        std::ostringstream fname;
        fname << "bin_" << lx << "x" << ly << "x" << lz << "_T="
                << std::setprecision(4) << 1.0 / beta << ".dat";
        std::ofstream out_bin(fname.str());
        out_bin << "Energy\tSpecific_Heat\tUniform_Susceptibility\tMagnetization^2\tNumber_of_H-operators\t(Number_of_H-operators)^2\n";
        out_bin.close();
    }
    // Открытие raw-файла и заголовок
    {
        std::ostringstream fname;
        fname << "raw_" << lx << "x" << ly << "x" << lz << "_T="
                << std::setprecision(4) << 1.0 / beta << ".dat";
        std::ofstream out_raw(fname.str());
        out_raw << "Number_of_H-operators\n";
        out_raw.close();
    }
    // Открытие cut-off-файла и заголовок
    {
        std::ostringstream fname;
        fname << "cut-off_" << lx << "x" << ly << "x" << lz << "_T="
                    << std::setprecision(4) << 1.0 / beta << ".dat";
        std::ofstream out_cut(fname.str());
        out_cut << "Step\tCut-off\tNumber_of_H-operators\tMaximum_Number_of_H-operators\n";
        out_cut.close();
    }
}

void diagonalupdate()
{
    for (int p = 0; p < L; ++p) {
        int op = opstring[p];

        if (op == 0) {
            double ran = uniform(rng);
            int bound_idx = static_cast<int>(ran * Nb) + 1;

            if (spin[bound_spin_ID[0][bound_idx - 1]-1] != spin[bound_spin_ID[1][bound_idx - 1]-1]) {
                ran = uniform(rng);
                if (ran * (NH + 1 - Nn) <= (add_prob / (L - NH)) * (NH + 1)) {
                    opstring[p] = 2 * bound_idx;
                    NH += 1;
                }
            }
        }
        else if ((op % 2) == 0) { //диагональный оператор
            double ran = uniform(rng);
            if (ran * NH <= (L - NH + 1) * rem_prob * (NH - Nn)) {
                opstring[p] = 0;
                NH -= 1;
            }
        }
        else {
            // Офф-диагональный оператор
            int bound_idx = op / 2;
            spin[bound_spin_ID[0][bound_idx - 1]-1] = -spin[bound_spin_ID[0][bound_idx - 1]-1];
            spin[bound_spin_ID[1][bound_idx - 1]-1] = -spin[bound_spin_ID[1][bound_idx - 1]-1];
        }
    }
}

void linkvertices()
{
    for (int i = 0; i < N; ++i) {
        first_vertex_visitted[i] = -1;
        last_vertex_visitted[i] = -1;
    }

    for (int v0 = 0; v0 < 4 * L; v0 += 4) {
        int p = v0 / 4;
        int op = opstring[p];

        if (op != 0) {
            int bound_idx = op / 2;
            int spin_idx1 = bound_spin_ID[0][bound_idx - 1];
            int spin_idx2 = bound_spin_ID[1][bound_idx - 1];

            int last_vertex1 = last_vertex_visitted[spin_idx1 - 1];
            int last_vertex2 = last_vertex_visitted[spin_idx2 - 1];

            if (last_vertex1 != -1) {
                vertex_link[last_vertex1] = v0;
                vertex_link[v0] = last_vertex1;
            } else {
                first_vertex_visitted[spin_idx1 - 1] = v0;
            }

            if (last_vertex2 != -1) {
                vertex_link[last_vertex2] = v0 + 1;
                vertex_link[v0 + 1] = last_vertex2;
            } else {
                first_vertex_visitted[spin_idx2 - 1] = v0 + 1;
            }

            last_vertex_visitted[spin_idx1 - 1] = v0 + 2;
            last_vertex_visitted[spin_idx2 - 1] = v0 + 3;
        } else {
            // vertex_link(v0:v0+3) = -1
            for (int i = 0; i < 4; ++i) {
                vertex_link[v0 + i] = -1;
            }
        }
    }
    for (int spin_idx = 1; spin_idx <= N; ++spin_idx) {
        int first_vertex = first_vertex_visitted[spin_idx - 1];
        if (first_vertex != -1) {
            int last_vertex = last_vertex_visitted[spin_idx - 1];
            vertex_link[last_vertex] = first_vertex;
            vertex_link[first_vertex] = last_vertex;
        }
    }
}

void loopupdate()
{
    for (int v0 = 0; v0 < 4 * L; v0 += 2) {
        if (vertex_link[v0] < 0)
            continue; // Fortran: cycle

        int vertex_in = v0;
        double ran = uniform(rng);

        if (ran < 0.5) {
            while (true) {
                opstring[vertex_in / 4] ^= 1;          // Change the operator type.
                vertex_link[vertex_in] = -2;           // Mark as visitted
                int vertex_out = vertex_in ^ 1;        // Next vertex
                vertex_in = vertex_link[vertex_out];
                vertex_link[vertex_out] = -2;
                if (vertex_in == v0)
                    break; // Выход — петля замкнулась
            }
        } else {
            while (true) {
                vertex_link[vertex_in] = -1;
                int vertex_out = vertex_in ^ 1;
                vertex_in = vertex_link[vertex_out];
                vertex_link[vertex_out] = -1;
                if (vertex_in == v0)
                    break;
            }
        }
    }
    for (int spin_idx = 1; spin_idx <= N; ++spin_idx) {
        if (first_vertex_visitted[spin_idx - 1] != -1) {
            if (vertex_link[first_vertex_visitted[spin_idx - 1]] == -2)
                spin[spin_idx - 1] = -spin[spin_idx - 1];
        } else {
            double ran = uniform(rng);
            if (ran < 0.5)
                spin[spin_idx - 1] = -spin[spin_idx - 1];
        }
    }
}

void adjustcutoff(int step)
{
    int L_new = NH + NH / 3;

    if (L_new <= L){
        return;
    }
    // --- Расширяем opstring ---
    // Копируем старое содержимое
    int* copy_opstring = new int[L];
    for (int i = 0; i < L; ++i)
        copy_opstring[i] = opstring[i];

    delete[] opstring;
    opstring = new int[L_new];

    for (int i = 0; i < L; ++i)
        opstring[i] = copy_opstring[i];
    for (int i = L; i < L_new; ++i)
        opstring[i] = 0;

    delete[] copy_opstring;

    L = L_new;

    delete[] vertex_link;
    vertex_link = new int[4 * L];
}

void measure() {
    int stag_mag = 0;

    for (int i = 0; i < N; ++i) {
        int x = i % lx;
        int y = (i / lx) % ly;
        int z = i / (lx * ly);
        int exponent = x + y + z;
        stag_mag += spin[i] * ((exponent % 2 == 0) ? 1 : -1);
    }
    stag_mag /= 2;
    double stag_mag2 = 0.0;

    for (int i = 0; i < L; ++i) {
        int op = opstring[i];
        if ((op % 2) == 1) {
            int bound_idx = op / 2;
            int spin_idx1 = bound_spin_ID[0][bound_idx - 1];
            int spin_idx2 = bound_spin_ID[1][bound_idx - 1];

            spin[spin_idx1 - 1] = -spin[spin_idx1 - 1];
            spin[spin_idx2 - 1] = -spin[spin_idx2 - 1];

            int x = i % lx;
            int y = (i / lx) % ly;
            int z = i / (lx * ly);
            int exponent = x + y + z;

            stag_mag += 2 * spin[spin_idx1 - 1] * ((exponent % 2 == 0) ? 1 : -1);
        }
        stag_mag2 += static_cast<double>(stag_mag * stag_mag);
    }

    stag_mag2 /= static_cast<double>(L);
    staggered += stag_mag2;

    long total_spin = 0;
    for (int i = 0; i < N; ++i) {
        total_spin += spin[i];
    }

    ususc += std::pow(static_cast<double>(total_spin) / 2.0, 2.0);
    n_opH += static_cast<double>(NH);
    n_opH2 += NH * NH;

    std::ostringstream fname;
    fname << "raw_" << lx << "x" << ly << "x" << lz << "_T="
                << std::setprecision(4) << 1.0 / beta << ".dat";
    std::ofstream out_raw(fname.str(), std::ios::app);
    if (out_raw.is_open()) {
        out_raw << NH << "\t" << -(static_cast<double>(NH) / (beta * N) - 0.25 *
                        static_cast<double>(Nb) / static_cast<double>(N)) << "\n";
        out_raw.close();
    }
}

void write_results() {
    n_opH   /= static_cast<double>(mcsteps);
    n_opH2  /= static_cast<double>(mcsteps);
    ususc   /= static_cast<double>(mcsteps);
    staggered /= static_cast<double>(mcsteps);

    double ener = - (n_opH / (beta * N) - 0.25 * static_cast<double>(Nb) / static_cast<double>(N));
    double Cv = (n_opH2 - n_opH * n_opH - n_opH) / static_cast<double>(N);
    // double Cv = (n_opH**2 - n_opH2 - n_opH) / dble(N)

    staggered = 3.0 * staggered / (static_cast<double>(N) * static_cast<double>(N));
    ususc = beta * ususc / static_cast<double>(N);

    n_opH   /= static_cast<double>(N);
    n_opH2  /= static_cast<double>(N);

    std::ostringstream fname;
    fname << "bin_" << lx << "x" << ly << "x" << lz << "_T="
                << std::setprecision(4) << 1.0 / beta << ".dat";
    std::ofstream out_bin(fname.str(), std::ios::app);
    if (out_bin.is_open()) {
        out_bin << std::setprecision(10) << ener << '\t' << Cv << '\t' << ususc << '\t'
                                << staggered << '\t' << n_opH << '\t' << n_opH2 << '\n';
        out_bin.close();
    }

    n_opH = 0.0;
    n_opH2 = 0.0;
    ususc = 0.0;
    staggered = 0.0;
}

void termalization()
{
    for (int i = 1; i <= term_steps; ++i) {
        diagonalupdate();
        linkvertices();
        loopupdate();
        adjustcutoff(i);
    }

    std::cout << "Series cut-off: " << L << "\n";
    std::cout << "----------------------------\n\n";
}

void qmc_steps()
{
    for (int j = 1; j <= Nbins; ++j) {
        for (int i = 1; i <= mcsteps; ++i) {
            diagonalupdate();
            linkvertices();
            loopupdate();
            measure();
        }
        write_results();
    }

   /* if (out_bin.is_open()) {
        out_bin.close();
    }
    if (out_raw.is_open()) {
        out_raw.close();
    }
    if (out_cut.is_open()) {
        out_cut.close();
    }*/
}

void results() {
    std::vector<double> e(Nbins, 0);
    std::vector<double> c(Nbins, 0);
    std::vector<double> x(Nbins, 0);
    std::vector<double> stag(Nbins, 0);
    std::vector<double> n_op(Nbins, 0);
    std::vector<double> n_op2(Nbins, 0);

    std::ostringstream fname;
    fname << "bin_" << lx << "x" << ly << "x" << lz << "_T="
                << std::setprecision(4) << 1.0 / beta << ".dat";
    std::ifstream in(fname.str());
    if (!in.is_open()) return;

    std::string header;
    std::getline(in, header);

    for (int i = 0; i < Nbins; ++i) {
        in >> e[i] >> c[i] >> x[i] >> stag[i] >> n_op[i] >> n_op2[i];
    }
    in.close();

    double inv = 1.0 / static_cast<double>(Nbins);
    for (int i = 0; i < Nbins; ++i) {
        e[i]    *= inv;
        c[i]    *= inv;
        x[i]    *= inv;
        stag[i] *= inv;
        n_op[i] *= inv;
        n_op2[i] *= inv;
    }

    std::ostringstream fname1;
    fname1 << "results_" << lx << "x" << ly << "x" << lz << ".dat";
    std::ofstream out(fname1.str(), std::ios::app); // дозапись!
    if (out.is_open()) {
        out << std::setprecision(10)
            << beta << "\t"
            << std::accumulate(e.begin(), e.end(), 0.0)    << "\t"
            << std::accumulate(c.begin(), c.end(), 0.0)    << "\t"
            << std::accumulate(x.begin(), x.end(), 0.0)    << "\t"
            << std::accumulate(stag.begin(), stag.end(), 0.0) << "\t"
            << std::accumulate(n_op.begin(), n_op.end(), 0.0) << "\t"
            << std::accumulate(n_op2.begin(), n_op2.end(), 0.0) << "\n";
        out.close();
    }
}

void make_dist() {
    int total = Nbins * mcsteps;
    int* opn = new int[total];
    double* energy = new double[total];

    std::ostringstream raw_name;
    raw_name << "raw_" << lx << "x" << ly << "x" << lz << "_T=" << std::setprecision(4) << 1.0 / beta << ".dat";
    std::ifstream in(raw_name.str());
    if (!in.is_open()) {
        delete[] opn; delete[] energy;
        return;
    }

    std::string header;
    std::getline(in, header);

    for (int i = 0; i < total; ++i) {
        in >> opn[i] >> energy[i];
    }
    in.close();

    int op_min = std::numeric_limits<int>::max();
    int op_max = std::numeric_limits<int>::min();
    for (int i = 0; i < total; ++i) {
        if (opn[i] < op_min) op_min = opn[i];
        if (opn[i] > op_max) op_max = opn[i];
    }
    int n_dist_sz = op_max - op_min + 1;
    int* n_dist = new int[n_dist_sz]();
    for (int i = 0; i < total; ++i) {
        n_dist[opn[i] - op_min] += 1;
    }

    std::ostringstream dist_n_name;
    dist_n_name << "dist_n_" << lx << "x" << ly << "x" << lz << "_T=" << std::setprecision(4) << 1.0 / beta << "Nn=" << Nn << ".dat";
    std::ofstream out_dist_n(dist_n_name.str());
    if (out_dist_n.is_open()) {
        for (int i = 0; i < n_dist_sz; ++i)
            out_dist_n << (op_min + i) << "\t" << n_dist[i] << "\n";
        out_dist_n.close();
    }
    delete[] n_dist;

    double e_min = energy[0], e_max = energy[0];
    for (int i = 1; i < total; ++i) {
        if (energy[i] < e_min) e_min = energy[i];
        if (energy[i] > e_max) e_max = energy[i];
    }
    for (int i = 0; i < total; ++i)
        energy[i] -= e_min;

    e_max = energy[0];
    for (int i = 1; i < total; ++i) {
        if (energy[i] > e_max) e_max = energy[i];
    }

    double energy_bin = 1.0 / (beta * N);
    int e_bins = static_cast<int>(e_max / energy_bin) + 1;
    int* e_dist = new int[e_bins]();
    for (int i = 0; i < total; ++i) {
        int idx = static_cast<int>(energy[i] / energy_bin);
        if (idx >= 0 && idx < e_bins)
            e_dist[idx] += 1;
    }

    std::ostringstream dist_e_name;
    dist_e_name << "dist_e_" << lx << "x" << ly << "x" << lz << "_T=" << std::setprecision(4) << 1.0 / beta << "Nn=" << Nn << ".dat";
    std::ofstream out_dist_e(dist_e_name.str());
    if (out_dist_e.is_open()) {
        for (int i = 0; i < e_bins; ++i)
            out_dist_e << (i * energy_bin + e_min) << "\t" << e_dist[i] << "\n";
        out_dist_e.close();
    }
    delete[] e_dist;
    delete[] opn;
    delete[] energy;
}

void free_memory()
{
    delete[] spin;
    delete[] opstring;
    delete[] vertex_link;
    delete[] first_vertex_visitted;
    delete[] last_vertex_visitted;
    if (bound_spin_ID) {
        delete[] bound_spin_ID[0];
        delete[] bound_spin_ID;
        bound_spin_ID = nullptr;
    }
}

int main()
{
    std::ostringstream filename;
    filename << "results_" << lx << "x" << ly << "x" << lz << ".dat";
    std::ofstream outfile(filename.str(), std::ios::app);

    outfile << "Beta\tEnergy\tSpecific Heat\tUniform Susceptibility\tMagnetization^2\tNumber of H-operators\t(Number of H-operators)^2" << std::endl;
    for (Nn=0; Nn<=100; Nn++){
        for(int t=0; t<t_steps;t++){
            beta = 1.0 / (temp_ini + dt * t);
            NH = 0;
            NH_max = 0;
            L = std::max(4, N / 4);

            opstring = new int[L](); // () — инициализация нулями
            first_vertex_visitted = new int[N]();
            last_vertex_visitted = new int[N]();
            vertex_link = new int[4 * L]();

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
    outfile.close();

    return 0;
}
