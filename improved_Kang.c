#include "udf.h"
#include "mem.h"

DEFINE_TRANS_RETHETA_C(hpyersonic_rethetac, c, t)
{
    real Re_thetat = C_RETHETA(c, t);
    real Re_thetac;
    real Re_thetac_baseline;
    real f_Me;
    real g_Me;
    real R_Me;
    real f_R_Me;
    real Me = 6.3;
    real Tw = 300;
    real T0 = 4868.427;

    if (Re_thetat > 1870)
        Re_thetac_baseline = Re_thetat - (593.11 + (Re_thetat - 1870.0) * 0.482);
    else
        Re_thetac_baseline = Re_thetat - ((396.035e-2) + (-120.656e-4) * Re_thetat + (868.230e-6) * pow(Re_thetat, 2) + (-696.506e-9) * pow(Re_thetat, 3) + (174.105e-12) * pow(Re_thetat, 4));

    f_Me = 0.8807 - 0.20323 * Me + (4.265e-4) * Tw + 0.3647 * (Tw / T0) + 0.092512 * Me * Me - (7.684e-5) * Me * Tw + 0.17956 * Me * (Tw / T0) - (6.05e-4) * Tw * (Tw / T0) + 0.9843 * pow(Tw / T0, 2) + 0.4192 * Me * Me * (Tw / T0);

    g_Me = 0.09490935 * pow(Me, 4) - 2.8455282 * pow(Me, 3) + 31.349949 * pow(Me, 2) - 149.21067 * Me + 260.69627;

    R_Me = (tanh(9.1667 * (Me - 2.2)) + 1) / 2;

    f_R_Me = f_Me * R_Me + (1 - R_Me) * 2.193;

    Re_thetac = f_R_Me * Re_thetac_baseline / (2.193 * g_Me);

    return Re_thetac;
}

DEFINE_TRANS_FLENGTH(hyper_Flength, c, t)
{
    real Flength;
    real Flength_Baeline;
    real R_Me;
    real h_Me;
    real h_R_Me;
    real Me = 6.3;
    real Re_thetat = C_RETHETA(c, t);

    if (Re_thetat < 400)
        Flength_Baeline = 398.189e-1 - (119.270e-4) * Re_thetat - (132.567e-6) * pow(Re_thetat, 2);
    else if (Re_thetat < 596)
        Flength_Baeline = 263.404 - (123.939e-2) * Re_thetat + (194.548e-5) * pow(Re_thetat, 2) - (101.695e-8) * pow(Re_thetat, 3);
    else if (Re_thetat < 1200)
        Flength_Baeline = 0.5 - (Re_thetat - 596.0) * (3e-4);
    else
        Flength_Baeline = 0.3188;
    h_Me = 0.5 / Me;
    R_Me = (tanh(9.1667 * (Me - 2.2)) + 1) / 2;
    h_R_Me = h_Me * R_Me + (1 - R_Me);
    Flength = Flength_Baeline / h_R_Me;
    return Flength;
}
