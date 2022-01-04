#include "udf.h"
#include "mem.h"

DEFINE_TRANS_RETHETA_C(user_Re_thetac, c, t)
{
    real Re_thetac;
    real Re_thetat = C_RETHETA(c, t);
    if (Re_thetat > 1870)
        Re_thetac = Re_thetat - (593.11 + (Re_thetat - 1870.0) * 0.482);
    else
        Re_thetac = Re_thetat - ((396.035e-2) - (120.656e-4) * Re_thetat + (868.230e-6) * pow(Re_thetat, 2) - (696.506e-9) * pow(Re_thetat, 3) + (174.105e-12) * pow(Re_thetat, 4));
    return Re_thetac;
}

DEFINE_TRANS_FLENGTH(user_Flength, c, t)
{
    real Flength;
    real Re_thetat = C_RETHETA(c, t);
    if (Re_thetat < 400)
        Flength = 398.189e-1 - (119.270e-4) * Re_thetat - (132.567e-6) * pow(Re_thetat, 2);
    else if (Re_thetat < 596)
        Flength = 263.404 - (123.939e-2) * Re_thetat + (194.548e-5) * pow(Re_thetat, 2) - (101.695e-8) * pow(Re_thetat, 3);
    else if (Re_thetat < 1200)
        Flength = 0.5 - (Re_thetat - 596.0) * (3e-4);
    else
        Flength = 0.3188;
    return Flength;
}

DEFINE_TRANS_RETHETA_T(user_Re_thetat, c, t)
{
    real Re_thetat;
    real F_Tu, Tu;
    real T = C_T(c, t);
    real rou = C_R(c, t);
    real T_Ke = C_K(c, t);
    real U = sqrt(C_VMAG2(c, t));
    real F_lambda, lambda_theta, temp;
    real K, miu;
    real DUDx, DUDy, DUDz, DUDs;
    real dudx = C_DUDX(c, t);
    real dudy = C_DUDY(c, t);
    //real dudz = C_DUDZ(c, t);
    real dvdx = C_DVDX(c, t);
    real dvdy = C_DVDY(c, t);
    //real dvdz = C_DVDZ(c, t);
    //real dwdx = C_DWDX(c, t);
    //real dwdy = C_DWDY(c, t);
    //real dwdz = C_DWDZ(c, t);
    real u = C_U(c, t);
    real v = C_V(c, t);
    //real w = C_W(c, t);
    real G_Ma;
    real Ma = 6.3;

    DUDx = 0.5 * pow(U * U, -0.5) * (2.0 * u * dudx + 2.0 * v * dvdx /*+ 2.0 * w * dwdx*/);
    DUDy = 0.5 * pow(U * U, -0.5) * (2.0 * u * dudy + 2.0 * v * dvdy /*+ 2.0 * w * dwdy*/);
    //DUDz = 0.5 * pow(U * U, -0.5) * (2.0 * u * dudz + 2.0 * v * dvdz + 2.0*w * dwdz);
    DUDs = (u / U) * DUDx + (v / U) * DUDy /*+ (w / U) * DUDz*/;

    miu = pow(T / 288.15, 1.5) * (288.15 + 110.4) / (T + 110.4) * (1.7894e-5);
    K = miu / (U * U) * DUDs;

    Tu = max(100.0 * sqrt(2.0 * T_Ke / 3.0) / U, 0.027);

    lambda_theta = 0.0;
    temp = 0.0;

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
    return Re_thetat;
}