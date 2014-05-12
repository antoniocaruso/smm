#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdio.h>
#include <math.h>
#include "smm.h"

/* Compile with gcc -std=c99 test.c smm.c -lgsl -lm */

#define M_PI  3.1415926535897932384626433832795029L  /* pi */
#define SIZE 1000   // max number of sensors.

double x[SIZE];
double y[SIZE];

int deployed = 0;

typedef struct { double min, avg, max; } Pair;

Pair stats(double v[], int n)
{
    Pair p;
    p.min = p.max = p.avg = v[0];
    for(int i=1; i<n; i++) {
        if (v[i] < p.min) p.min = v[i];
        if (v[i] > p.max) p.max = v[i];
        p.avg += v[i];
    }
    p.avg /= n;
    return p;
}

/* SmmTest:
   * Init the library for Meandering Jet
   * Deploy the sensors
   * Trace the fastest along the x axis...(the flow).
*/

int main(int argc, char *argv[])
{
    gsl_rng *r = gsl_rng_alloc(gsl_rng_ranlxs1);
    int seed = 1214314;
    gsl_rng_set(r, seed);

    const double c = 0.12;
    const double k = 2.0*M_PI/7.5;
    const double T = (2*M_PI)/(c*k);
    double randomshift = gsl_rng_uniform(r)*T;

    printf("# seed = %d, phase = %f\n", seed, randomshift);

    const int Streamfunction = SMM_MNDRJET;
    double params[SMM_MXPRMS]  = { 1.2, 0.12, 2.0*M_PI/7.5, 0.4, 0.3, randomshift };
    
    // Passo di Integrazione
    const double SmmTimeStep = 0.01; 
    // 1 t = 0.03 Days.
    const double DimensionalTime = 0.03;    
    // Durata in secondi del singolo passo = 25.92
    const double MicroStepTime = SmmTimeStep*24*60*60 * DimensionalTime;    

    SMM_init(SIZE, Streamfunction, params, SMM_RK2, SmmTimeStep);
    
    // initial deployment 200 sensors in the domain [0,4] x [-2,2] km

    int n = 200;
    for (int i=0; i<n; i++) {
        x[i] = gsl_rng_uniform(r) * 4.0;
        y[i] = gsl_rng_uniform(r) * 4.0 - 2.0;
    }
    deployed += n;
    SMM_deploy_nodes(n);

    // run for one day -> 24 h x 60x60 s.
    // trace the fastest particle;

    printf("# Time - min(X) - max(X)\n");
    Pair xp = stats(x,n);
    printf("%7.2f\t %.2f\t %.2f  %.2f\n", 0.0, 
            xp.min*1000, xp.avg*1000, xp.max*1000); 

    for (int i=0; i<100; i++) {
        SMM_Move_Sensors(x,y);
        Pair xp = stats(x,n);
        // time in seconds..
        double time = SMM_get_current_time() * (MicroStepTime/SmmTimeStep);
        printf("%7.2f\t %.2f\t %.2f  %.2f\n", 
                time, xp.min*1000, xp.avg*1000, xp.max*1000); 
    }
}
