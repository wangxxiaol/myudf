#include "udf.h"

#define sigma_Re 2.0 // constant number
#define c_Re 0.03    // constant number
#define c_e2 50      // constant number

#define T_R 270.1044   // reference temperature
#define miu_R 1.701e-5 // reference viscosity
#define gamma 1.4      // specific heat ration

DEFINE_UDS_FLUX(lqy_uds_flux, f, t, i)
{
    cell_t c0, c1 = -1;
    Thread *t0, *t1 = NULL;
    real NV_VEC(psi), NV_VEC(A), flux = 0.0;
    c0 = F_C0(f, t);
    t0 = F_C0_THREAD(f, t);
    /* defining psi in terms of velocity field */
    // NV_D(psi, =, F_U(f, t), F_V(f, t), F_W(f, t));
    // NV_S(psi, *=, F_R(f, t)); /* multiplying density to get psi vector */
    // F_AREA(A, f, t); /* face normal vector returned from F_AREA */
    // return NV_DOT(psi, A); /* dot product of the two returned */
    // if (BOUNDARY_FACE_THREAD_P(t))
    //{
    // real rho;
    // if (NULLP(THREAD_STORAGE(t, SV_DENSITY)))
    //     rho = F_R(f, t);
    // else
    // rho = C_R(c0, t0);
    // NV_DS(psi, =, C_U(c0, t0), C_V(c0, t0), C_W(c0, t0), *, rho);
    // flux = NV_DOT(psi, A);
    // }
    // else
    //{
    // c1 = F_C1(f, t);
    // t1 = F_C1_THREAD(f, t);
    NV_DS(psi, =, C_U(c0, t0), C_V(c0, t0), C_W(c0, t0), *, C_R(c0, t0));
    // NV_DS(psi, =, C_U(c1, t1), C_V(c1, t1), C_W(c1, t1), *, C_R(c1, t1));
    flux = NV_DOT(psi, A);
    //}
    return flux;
}

DEFINE_UDS_UNSTEADY(lqy_uds_unsteady, c, t, i, apu, su)
{
    real physical_dt, vol, rho, fre_old;
    physical_dt = RP_Get_Real("physical-time-step");
    vol = C_VOLUME(c, t);
    rho = C_R_M1(c, t);
    *apu = -rho * vol / physical_dt; /*implicit part */
    fre_old = C_STORAGE_R(c, t, SV_UDSI_M1(i));
    *su = rho * vol * fre_old / physical_dt; /*explicit part*/
}

DEFINE_DIFFUSIVITY(lqy_diffusivity, c, t, i)
{
    return sigma_Re * (C_MU_L(c, t) + C_MU_T(c, t));
}

DEFINE_SOURCE(lqy_source, c, t, dS, equ)
{
    real source;
    real y = C_WALL_DIST(c, t);            // distance to nearest wall
    real rho = C_R(c, t);                  // local density
    real T = C_T(c, t);                    // local temperature
    real U = sqrt(C_VMAG2(c, t));          // velocity magnitude
    real miu = C_MU_L(c, t);               // local viscosity
    real f_Re = C_UDSI(c, t, 0);           // from scalar equation
    real Re_thetat = C_RETHETA(c, t);      // from fluent re_theta transport equation
    real intermittency = C_INTERMIT(c, t); // from fluent intermittency transport equation
    real omega = C_O(c, t);                // specific dissipation rate
    // real gamma = C_GAMMA(c, t);       // specific heat ration
    real time_scale;
    real f_ReL;
    real F_thetat, F_wake;
    real Re_omega;
    real delta_comp; // boundary-layer thickness
    real vorticity;  // vorticity magnitude
    real delta_BL, theta_BL;
    real a, b;

    /*****calculate time_scale---t*****/
    time_scale = 500 * miu / (rho * U * U);

    /*****calculate f_ReL*****/
    f_ReL = T_R * miu_R / (T * miu); // use local temperature and viscosity

    /*****calculate Re_omega*****/
    Re_omega = rho * omega * y * y / miu;

    /*****calculate F_wake*****/
    F_wake = exp(-pow(Re_omega / (1e+5), 2));

    /*****calculate vorticity magnitude*****/
#if RP_3D /*3D*/
    real wx = C_DWDY(c, t) - C_DVDZ(c, t), wy = C_DUDZ(c, t) - C_DWDX(c, t), wz = C_DVDX(c, t) - C_DUDY(c, t);
    vorticity = sqrt(wx * wx + wy * wy + wz * wz);
#else /*2D*/
    vorticity = C_DVDX(c, t) - C_DUDY(c, t);
#endif
    /*****calculate boundary-layer thickness*****/
    theta_BL = Re_thetat * miu / rho / U;

    delta_BL = f_Re * 7.5 * theta_BL; // improved by f_Re

    delta_comp = 50.0 * vorticity * y * delta_BL / U;

    /*****calculate F_thetat*****/
    a = F_wake * exp(-pow(y / delta_comp, 4));

    b = 1.0 - pow((intermittency - 1 / c_e2) / (1 - 1 / c_e2), 2);

    F_thetat = min(max(a, b), 1.0);

    /*****calculate source*****/
    source = c_Re * rho / time_scale * (f_ReL - f_Re) * (1.0 - F_thetat);

    if (a > b)
        if (a > 1)
            dS[equ] = -c_Re * rho / time_scale * (1.0 - F_thetat);
        else
            dS[equ] = -c_Re * rho / time_scale * ((1.0 - F_thetat) + (f_ReL - f_Re) * a * 4 * pow(y / (50.0 * vorticity * y * 7.5 * theta_BL / U), 4) * pow(f_Re, -5));
    else
        dS[equ] = -c_Re * rho / time_scale * (1.0 - F_thetat);

    return source;
}

DEFINE_TRANS_RETHETA_C(lqy_Re_thetac, c, t)
{
    real Re_thetac;
    real f_Re = C_UDSI(c, t, 0);
    real Re_thetat_new = C_RETHETA(c, t) / f_Re;

    if (Re_thetat_new > 1870)
        Re_thetac = Re_thetat_new - (593.11 + (Re_thetat_new - 1870.0) * 0.482);
    else
        Re_thetac = Re_thetat_new - ((396.035e-2) - (120.656e-4) * Re_thetat_new + (868.230e-6) * pow(Re_thetat_new, 2) - (696.506e-9) * pow(Re_thetat_new, 3) + (174.105e-12) * pow(Re_thetat_new, 4));
    return f_Re * Re_thetac;
}

DEFINE_TRANS_FLENGTH(lqy_Flength, c, t)
{
    real Flength;
    real f_Re = C_UDSI(c, t, 0);
    real Re_thetat_new = C_RETHETA(c, t) / f_Re;

    if (Re_thetat_new < 400)
        Flength = (398.189e-1) - (119.270e-4) * (Re_thetat_new) - (132.567e-6) * pow(Re_thetat_new, 2);
    else if (Re_thetat_new < 596)
        Flength = 263.404 - (123.939e-2) * (Re_thetat_new) + (194.548e-5) * pow(Re_thetat_new, 2) - (101.695e-8) * pow(Re_thetat_new, 3);
    else if (Re_thetat_new < 1200)
        Flength = 0.5 - (Re_thetat_new - 596.0) * (3e-4);
    else
        Flength = 0.3188;
    return Flength;
}

DEFINE_TRANS_RETHETA_T(lqy_Re_thetat, c, t)
{
    real Re_thetat;
    real F_Tu, Tu;
    real T = C_T(c, t); // local temperature
    // real gamma = C_GAMMA(c, t);   // specific heat ration
    real rho = C_R(c, t);         // density
    real T_Ke = C_K(c, t);        // turbulence kinetic energy
    real U = sqrt(C_VMAG2(c, t)); // velocity magnitude
    real Me;                      // local Mach number
    real F_lambda, lambda_theta = 0.0, temp = 0.0;
    real K, niu;
    real DUDx, DUDy, DUDz, DUDs;
    real c0, c1, c2, c_fc;

#if RP_3D /*3D*/
    DUDx = 0.5 * pow(U * U, -0.5) * (2.0 * C_U(c, t) * C_DUDX(c, t) + 2.0 * C_V(c, t) * C_DVDX(c, t) + 2.0 * C_W(c, t) * C_DWDX(c, t));

    DUDy = 0.5 * pow(U * U, -0.5) * (2.0 * C_U(c, t) * C_DUDY(c, t) + 2.0 * C_V(c, t) * C_DVDY(c, t) + 2.0 * C_W(c, t) * C_DWDY(c, t));

    DUDz = 0.5 * pow(U * U, -0.5) * (2.0 * C_U(c, t) * C_DUDZ(c, t) + 2.0 * C_V(c, t) * C_DVDZ(c, t) + 2.0 * C_W(c, t) * C_DWDZ(c, t));

    DUDs = (C_U(c, t) / U) * DUDx + (C_V(c, t) / U) * DUDy + (C_W(c, t) / U) * DUDz;
#else /*2D*/
    DUDx = 0.5 * pow(U * U, -0.5) * (2.0 * C_U(c, t) * C_DUDX(c, t) + 2.0 * C_V(c, t) * C_DVDX(c, t));

    DUDy = 0.5 * pow(U * U, -0.5) * (2.0 * C_U(c, t) * C_DUDY(c, t) + 2.0 * C_V(c, t) * C_DVDY(c, t));

    DUDs = (C_U(c, t) / U) * DUDx + (C_V(c, t) / U) * DUDy;
#endif

    niu = C_MU_L(c, t) / rho;

    K = niu / (U * U) * DUDs;

    Tu = max(100.0 * sqrt(2.0 * T_Ke / 3.0) / U, 0.027);

    /*****calculate Re_thetat*****/
    for (int i = 0; i < 10; i++)
    {
        lambda_theta = temp;
        /*****calculate F_Tu*****/
        if (Tu > 1.3)
            F_Tu = 331.5 * pow(Tu - 0.5658, -0.671);
        else
            F_Tu = 1173.51 - 589.428 * Tu + 0.2196 / (Tu * Tu);
        /*****calculate lambda_theta*****/
        if (lambda_theta > 0)
            F_lambda = 1.0 + 0.275 * (1.0 - exp(-35.0 * lambda_theta)) * exp(-Tu / 0.5);
        else
            F_lambda = 1.0 - (-12.986 * lambda_theta - 123.66 * pow(lambda_theta, 2) - 405.689 * pow(lambda_theta, 3)) * exp(-pow(Tu / 1.5, 1.5));

        Re_thetat = max(F_Tu * F_lambda, 20.0);

        temp = min(max(-0.1, Re_thetat * Re_thetat * K), 0.1);
    }
    /*****calculate local Mach number*****/
    Me = max(0.4, U / sqrt(gamma * 287 * T)); // use local temperature

    c0 = -0.00141089 * pow(Me, 3) - 0.00467533 * pow(Me, 2) - 0.0270837 * Me + 0.00576259; // use local temperature

    c1 = 0.00298137 * pow(Me, 3) + 0.0103366 * pow(Me, 2) + 0.0453367 * Me + 1.02002; // use local temperature

    c2 = -0.00884236 * pow(Me, 3) + 0.0864964 * pow(Me, 2) - 0.323869 * Me - 0.404892; // use local temperature

    /*****correction function*****/
    c_fc = pow(10, (c2 * pow(log10(T_R / T), 2) + c1 * log10(T_R / T) + c0)); // use local temperature

    return Re_thetat * c_fc;
}