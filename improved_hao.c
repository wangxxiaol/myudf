#include "udf.h"
#include "mem.h"

#define Me = 6.3;
#define Re = 1.77e+6;

DEFINE_TRANS_RETHETA_C(hpyersonic_rethetac, c, t)
{
    real Re_thetat = C_RETHETA(c, t); // from fluent re_theta transport equation
    real T = C_T(c, t);               // local temperature
    real S = C_STRAIN_RATE_MAG(c, t); // strain rate magnitude
    real rho = C_R(c, t);             // density
    real y = C_WALL_DIST(c, t);       // distance to nearest wall
    real miu = C_MU_L(c, t);          // local viscosity
    real Re_v, Re_thetac, Re_thetac_baseline, K, D, Fr;

    if (Re_thetat > 1870)
        Re_thetac_baseline = Re_thetat - (593.11 + (Re_thetat - 1870.0) * 0.482);
    else
        Re_thetac_baseline = Re_thetat - ((396.035e-2) + (-120.656e-4) * Re_thetat + (868.230e-6) * pow(Re_thetat, 2) + (-696.506e-9) * pow(Re_thetat, 3) + (174.105e-12) * pow(Re_thetat, 4));

    Re_v = rho * y * y * S / miu;

    K = 2.1 - 0.0895 * Me;

    D = 300 + 0.0000275 * Re;

    Fr = exp(-5 * Re_v / (2 * D * D));

    Re_thetac = K * Re_thetac_baseline * Re_v / 2.193 * (Re_v - D * (1 - Fr));

    return Re_thetac;
}

DEFINE_TRANS_FLENGTH(hyper_Flength, c, t)
{
    return 20;
}
