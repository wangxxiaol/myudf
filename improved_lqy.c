#include "udf.h"
#include "mem.h"
#include <math.h>

#define T_w 300

DEFINE_TRANS_RETHETA_C(user_Re_thetac, c, t)
{
    real Re_thetac;
    real f_Re = C_UDSI(c, t, 0);
    real Re_thetat_new = C_RETHETA(c, t) / f_Re;
    if (Re_thetat_new > 1870)
        Re_thetac = (Re_thetat_new - (593.11 + (Re_thetat_new - 1870.0) * 0.482));
    else
        Re_thetac = (Re_thetat_new - ((396.035e-2) - (120.656e-4) * Re_thetat_new + (868.230e-6) * pow(Re_thetat_new, 2) - (696.506e-9) * pow(Re_thetat_new, 3) + (174.105e-12) * pow(Re_thetat_new, 4)));
    return f_Re * Re_thetac;
}

DEFINE_TRANS_FLENGTH(user_Flength, c, t)
{
    real Flength;
    real f_Re = C_UDSI(c, t, 0);
    // Message("f_Re is: %f\n", f_Re);
    real Re_thetat_new = C_RETHETA(c, t) / f_Re;
    if (Re_thetat_new < 400)
        Flength = 398.189e-1 - (119.270e-4) * (Re_thetat_new) - (132.567e-6) * pow(Re_thetat_new, 2);
    else if (Re_thetat_new < 596)
        Flength = 263.404 - (123.939e-2) * (Re_thetat_new) + (194.548e-5) * pow(Re_thetat_new, 2) - (101.695e-8) * pow(Re_thetat_new, 3);
    else if (Re_thetat_new < 1200)
        Flength = 0.5 - (Re_thetat_new - 596.0) * (3e-4);
    else
        Flength = 0.3188;
    return Flength;
}

DEFINE_TRANS_RETHETA_T(user_Re_thetat, c, t)
{
    real Re_thetat;
    real F_Tu, Tu;
    real T = C_T(c, t);
    real T_R, T_aw;
    real gamma = C_GAMMA(c, t); // specific heat ration
    real rho = C_R(c, t);       // density
    real T_Ke = C_K(c, t);      // turb. kinetic energy
    real U = sqrt(C_VMAG2(c, t));
    real Me; // Mach number
    real F_lambda, lambda_theta = 0.0, temp = 0.0;
    real K, niu;
    real DUDx, DUDy, DUDz, DUDs;
    real dudx = C_DUDX(c, t), dudy = C_DUDY(c, t), dudz = C_DUDZ(c, t);
    real dvdx = C_DVDX(c, t), dvdy = C_DVDY(c, t), dvdz = C_DVDZ(c, t);
    real dwdx = C_DWDX(c, t), dwdy = C_DWDY(c, t), dwdz = C_DWDZ(c, t);
    real u = C_U(c, t), v = C_V(c, t), w = C_W(c, t);
    real c0, c1, c2, c_fc;

    Me = min(max(U / sqrt(gamma * 287 * T), 0.4), 8);

    DUDx = 0.5 * pow(U * U, -0.5) * (2.0 * u * dudx + 2.0 * v * dvdx + 2.0 * w * dwdx);
    DUDy = 0.5 * pow(U * U, -0.5) * (2.0 * u * dudy + 2.0 * v * dvdy + 2.0 * w * dwdy);
    DUDz = 0.5 * pow(U * U, -0.5) * (2.0 * u * dudz + 2.0 * v * dvdz + 2.0 * w * dwdz);

    DUDs = (u / U) * DUDx + (v / U) * DUDy / +(w / U) * DUDz;

    niu = pow(T / 288.15, 1.5) * (288.15 + 110.4) / (T + 110.4) * (1.7894e-5) / rho;

    K = niu / (U * U) * DUDs;

    Tu = max(100.0 * sqrt(2.0 * T_Ke / 3.0) / U, 0.027);

    for (int i = 0; i < 10; i++)
    {
        lambda_theta = temp;

        if (Tu > 1.3)
            F_Tu = 331.5 * pow(Tu - 0.5658, -0.671);
        else
            F_Tu = 1173.51 - 589.428 * Tu + 0.2196 / (Tu * Tu);

        if (lambda_theta > 0)
            F_lambda = 1.0 + 0.275 * (1.0 - exp(-35.0 * lambda_theta)) * exp(-Tu / 0.5);
        else
            F_lambda = 1.0 - (-12.986 * lambda_theta - 123.66 * pow(lambda_theta, 2) - 405.689 * pow(lambda_theta, 3)) * exp(-pow(Tu / 1.5, 1.5));

        Re_thetat = max(F_Tu * F_lambda, 20.0);

        temp = min(max(Re_thetat * Re_thetat * K, -0.1), 0.1);
    }

    c0 = -0.00141089 * pow(Me, 3) - 0.00467533 * pow(Me, 2) - 0.0270837 * Me + 0.00576259;

    c1 = 0.00298137 * pow(Me, 3) + 0.0103366 * pow(Me, 2) + 0.0453367 * Me + 1.02002;

    c2 = -0.00884236 * pow(Me, 3) + 0.0864964 * pow(Me, 2) - 0.323869 * Me - 0.404892;

    T_aw = T * (1 + 0.85 * ((gamma - 1) / 2) * Me * Me); // use local temperature and Mach number

    T_R = 0.5 * T + 0.22 * T_aw + 0.28 * T;

    c_fc = pow(10, (c2 * pow(log10(T_R / T), 2) + c1 * log10(T_R / T) + c0));

    return Re_thetat * c_fc;
}