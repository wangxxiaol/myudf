#include "udf.h"

DEFINE_TRANS_RETHETA_C(user_Re_thetac, c, t)
{
    real Re_thetac;
    real Re_thetat = C_RETHETA(c, t);
    if(Re_thetat > 1870)
        Re_thetac = Re_thetat-(593.11+(Re_thetat- 1870.0)*0.482);
    else
        Re_thetac = Re_thetat-((396.035e-2)- (120.656e-4)*Re_thetat+ (868.230e-6)*pow(Re_thetat, 2)- (696.506e-9)*pow(Re_thetat, 3)+ (174.105e-12)*pow(Re_thetat, 4));
    return Re_thetac;
}

DEFINE_TRANS_FLENGTH(user_Flength, c, t)
{
    real Flength;
    real Re_thetat = C_RETHETA(c, t);
    if(Re_thetat < 400)
        Flength = 398.189e-1-(119.270e-4)*Re_thetat-(132.567e-6)*pow(Re_thetat, 2);
    else if(Re_thetat < 596)
        Flength = 263.404-(123.939e-2)*Re_thetat+(194.548e-5)*pow(Re_thetat, 2)-(101.695e-8)*pow(Re_thetat, 3);
    else if(Re_thetat < 1200)
        Flength = 0.5-(Re_thetat-596.0)*(3e-4);
    else
        Flength = 0.3188;
    return Flength;
}