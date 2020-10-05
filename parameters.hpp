#ifndef PARAMETERS_H
#define PARAMETERS_H

/// Boltzmann constant, in eV/K
#define KB 8.617333262145e-5

/// Barrier for hops in zeta minus ordering, in eV
//  TODO this is average, add in actual barriers
#define ZETA_MINUS_BARRIER 0.03

double zeta_minus_boundary_energy(const int& spacing)
{
    switch (spacing)
    {
    case 1:
        return 0.0543;
//    case 2:
//        return 0.01011;
//    case 3:
//        return 0.0042;
//    case 4:
//        return 0.00003999999996;
    default:
        return 0;
    }
};

#endif
