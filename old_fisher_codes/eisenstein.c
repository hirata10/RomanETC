#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define POW4(x) ((x)*(x)*(x)*(x))
#define SQR(x) ((x)*(x))

#define float double
float TFsound_horizon(float omega0hh, float f_baryon, float Tcmb)
/* Return the sound horizon in Mpc (not h^-1 Mpc) */
{
    float z_drag_b1, z_drag_b2, z_drag, sound_horizon;
    float omhh, obhh, theta_cmb, z_equality, k_equality, R_drag, R_equality;
    omhh = omega0hh;
    obhh = omhh*f_baryon;
    theta_cmb = Tcmb/2.7;
    z_equality = 2.50e4*omhh/POW4(theta_cmb);  /* Really 1+z */
    k_equality = 0.0746*omhh/SQR(theta_cmb);

    z_drag_b1 = 0.313*pow(omhh,-0.419)*(1+0.607*pow(omhh,0.674));
    z_drag_b2 = 0.238*pow(omhh,0.223);
    z_drag = 1291*pow(omhh,0.251)/(1+0.659*pow(omhh,0.828))*
                (1+z_drag_b1*pow(obhh,z_drag_b2));

    R_drag = 31.5*obhh/POW4(theta_cmb)*(1000/(1+z_drag));
    R_equality = 31.5*obhh/POW4(theta_cmb)*(1000/z_equality);

    sound_horizon = 2./3./k_equality*sqrt(6./R_equality)*
            log((sqrt(1+R_drag)+sqrt(R_drag+R_equality))/(1+sqrt(R_equality)));
    return sound_horizon;
}
#undef float
