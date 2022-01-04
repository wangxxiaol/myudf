#include "udf.h"

#define sigma_Re 2.0 // constant number
#define c_Re 0.03    // constant number
#define c_e2 50      // constant number
#define Me 6.3       // freesream Mach number
#define T_w 300      // iosthermal wall temperature
#define T_e 570      // freesream temperature

DEFINE_UDS_FLUX(my_uds_flux, f, t, i)
{
    cell_t c0, c1 = -1;
    Thread *t0, *t1 = NULL;
    real NV_VEC(psi), NV_VEC(A), flux = 0.0;
    c0 = F_C0(f, t);
    t0 = F_C0_THREAD(f, t);
    /* defining psi in terms of velocity field */
    // NV_D(psi, =, F_U(f, t), F_V(f, t), F_W(f, t));
    // NV_S(psi, *=, F_R(f, t)); /* multiplying density to get psi vector */
    F_AREA(A, f, t); /* face normal vector returned from F_AREA */
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

DEFINE_UDS_UNSTEADY(my_uds_unsteady, c, t, i, apu, su)
{
    real physical_dt, vol, rho, phi_old;
    physical_dt = RP_Get_Real("physical-time-step");
    vol = C_VOLUME(c, t);
    rho = C_R_M1(c, t);
    *apu = -rho * vol / physical_dt; /*implicit part */
    phi_old = C_STORAGE_R(c, t, SV_UDSI_M1(i));
    *su = rho * vol * phi_old / physical_dt; /*explicit part*/
}

DEFINE_DIFFUSIVITY(my_diffusivity, c, t, i)
{
    real T = C_T(c, t);
    real miu;
    real miu_t = C_MU_T(c, t);
    miu = pow(T / 288.15, 1.5) * (288.15 + 110.4) / (T + 110.4) * (1.7894e-5);
    return sigma_Re * (miu + miu_t);
}

DEFINE_SOURCE(my_source, c, t, dS, equ)
{
    real x[ND_ND]; /* this will hold the position vector */
    C_CENTROID(x, c, t);
    real source;
    real y = x[1];
    real rho = C_R(c, t);
    real T = C_T(c, t);
    real T_R, T_aw;               // reference temperature, Adiabatic wall temperature
    real U = sqrt(C_VMAG2(c, t)); // velocity
    real time_scale;
    real miu, miu_R; // viscosity, reference viscosity
    real f_ReL;
    real f_Re = C_UDSI(c, t, 0);      // from scalar equation
    real Re_thetat = C_RETHETA(c, t); // from fluent
    real gamma = C_GAMMA(c, t);       // specific heat ration
    real intermittency = C_INTERMIT(c, t);
    real F_thetat, F_wake;
    real omega = C_O(c, t); // specific dissipation rate
    real Re_omega;
    real delta_comp; // boundary-layer thickness
    real dudx = C_DUDX(c, t), dudy = C_DUDY(c, t), dudz = C_DUDZ(c, t);
    real dvdx = C_DVDX(c, t), dvdy = C_DVDY(c, t), dvdz = C_DVDZ(c, t);
    real dwdx = C_DWDX(c, t), dwdy = C_DWDY(c, t), dwdz = C_DWDZ(c, t);
    real vorticity[3][3] = {0, 0.5 * (dudy - dvdx), 0.5 * (dudz - dwdx), 0.5 * (dvdx - dudy), 0, 0.5 * (dvdz - dwdy), 0.5 * (dwdx - dudz), 0.5 * (dwdy - dvdz), 0};
    real temp[3][3];
    real sum1, sum2, sum;
    real abs_vorticity;
    real theta_BL;

    /*****calculate reference temperature*****/
    T_aw = T_e * (1 + 0.85 * ((gamma - 1) / 2) * Me * Me); // use freesream temperature
    T_R = 0.05 * T_w + 0.22 * T_aw + 0.28 * T_e;           // use freesream temperature

    /*****calculate local viscosity, reference viscosity*****/
    miu = pow(T / 288.15, 1.5) * (288.15 + 110.4) / (T + 110.4) * (1.7894e-5);       // use local temperature
    miu_R = pow(T_R / 288.15, 1.5) * (288.15 + 110.4) / (T_R + 110.4) * (1.7894e-5); // use reference temperature

    /*****calculate time_scale---t*****/
    time_scale = 500 * miu / (rho * U * U);

    /*****calculate f_ReL*****/
    f_ReL = T_R * miu_R / (T * miu); // use local temperature and viscosity

    /*****calculate Re_omega*****/
    Re_omega = rho * omega * y * y / miu;

    /*****calculate F_wake*****/
    F_wake = exp(-pow(Re_omega / (1e+5), 2));

    /*****calculate absolute value of vorticity*****/
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            temp[i][j] = 0;
            for (int k = 0; k < 3; k++)
                temp[i][j] += vorticity[i][k] * vorticity[k][j];
        }
    }
    sum1 = temp[0][2] * temp[1][0] * temp[2][1] + temp[0][1] * temp[1][2] * temp[2][0] + temp[0][0] * temp[1][1] * temp[2][2];

    sum2 = temp[2][2] * temp[0][1] * temp[1][0] + temp[2][1] * temp[1][2] * temp[0][0] + temp[0][2] * temp[1][1] * temp[2][0];

    sum = sum1 - sum2;

    abs_vorticity = sqrt(2 * sum);

    /*****calculate absolute value of vorticity*****/
    theta_BL = f_Re * 7.5 * Re_thetat * miu / rho / U;

    delta_comp = 50.0 * abs_vorticity * y * theta_BL / U;

    /*****calculate F_thetat*****/
    F_thetat = min(max(F_wake * exp(-pow(y / delta_comp, 4)), 1.0 - pow((intermittency - 1 / c_e2) / (1 - 1 / c_e2), 2)), 1.0);

    /*****calculate source*****/
    source = c_Re * rho / time_scale * (f_ReL - f_Re) * (1.0 - F_thetat);

    dS[equ] = 0.0;
}