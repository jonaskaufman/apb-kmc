#ifndef PARAMETERS_H
#define PARAMETERS_H

/// Boltzmann constant, in eV/K
#define KB 8.617333262145e-5

/// Barrier for hops in zeta minus ordering, in eV
//  TODO this is average, add in actual barriers
#define ZETA_MINUS_BARRIER 0.03

#define ZETA_MINUS_KINK_FORM 0.0297
#define ZETA_MINUS_KINK_MOVE_I 0.0324
#define ZETA_MINUS_KINK_MOVE_II 0.0331
#define ZETA_MINUS_REPULSION 0.0543

// TODO does this need to be a function?
double zeta_minus_boundary_energy(const int& spacing)
{
    switch (spacing)
    {
    case 1:
        return 0.0543;
    default:
        return 0;
    }
};

#endif
