
/********************************************************
Simple Mobility Model library

Note: GNU gsl <http://www.gnu.org/s/gsl/> is required.

   SmmLib code
   Copyright 2005-2011 Francesco Paparella, Antonio Caruso
   
   This file is part of SmmLib

   SmmLib is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   SmmLib is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with SmmLib.  If not, see <http://www.gnu.org/licenses/>.

*********************************************************/ 
/*
This library advects sensors according to a time-dependent, 
one-degree-of-freedom Hamiltonian (the streamfunction). 

NOTE: Error-checking is almost non-existent; 
you must know what you are doing!
*/


/* ChangeLog 24/07/2013: * Antonio *
    Added phaseshift parameter for Meandering Jet:
    The parameter is passed by the caller, in PsiPrms[5].
    The simulation time, is advanced before head, by phaseshift. Variable
    current_time is advanced in INIT, the function that reports the time
    is corrected accordingly. The caller must generate this shift in the
    intervarl [0,T] with T = 2PI/ck.
  -----------
*/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


#include "smm.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Private variables */

static double *PsiPrms;
static int Ns;
static double dt;
static double current_time;
static double *u;
static double *v;
static double *Xmed;
static double *Ymed;
static gsl_rng* r;
static double ur, vr;
static double U = 1.296;        // 0.5;   // 0.5 m/s
static double lambda = 0.015;   // (1/lambda) * 0.03g = 2g -> lambda = 0.015
static double expfact, sqrt_expfact2; // precalculated values for OU
static int Sf;
static int deployed = 0;

static double phaseshift = 0;

/* Function Prototipes */

void alternating_channels(double *posX, double *posY, double *u, double *v, double t);
void blinking_wavenumbers(double *posX, double *posY, double *u, double *v, double t);
void pulsating_vortex(double *posX, double *posY, double *u, double *v, double t);
void steady_vortex(double *posX, double *posY, double *u, double *v, double t);
void meandering_jet(double *posX, double *posY, double *u, double *v, double t);
void random_walk(double *posX, double *posY, double *u, double *v, double t);
void taylor_dispersion_jet(double *posX, double *posY, double *u, double *v, double t);
void rk2(double *posX, double *posY);


Streamfunction_Funptr SMM_Velocities;
Integrator_Funptr SMM_Move_Sensors;

double SMM_get_current_time(void) { return current_time - phaseshift; }

/* Deploy sensors: 
   Add sensors to the simulation, increase global variable deployed.
   Precondition: nodes > 0
   Post-condition: deployed += nodes if (nodes+deployed < Ns)
*/ 

void SMM_deploy_nodes(int nodes)
{
	if (deployed + nodes <= Ns)
        deployed += nodes;
	else {
		printf("Too many sensor deployed\n");
		exit(1);
	}
}


/* Initialization Function */

void SMM_init(int Nsensors, int Streamfunction, double *UserPrms, 
              int Integrator, double TimeStep)
{
    Ns = Nsensors;
    dt = TimeStep;
    current_time = 0;
    deployed = 0;
    Sf = Streamfunction;

    assert(UserPrms); // must be non NULL    
	if (PsiPrms == NULL)
	    PsiPrms = (double *)malloc(SMM_MXPRMS *sizeof(double));
    for (int i=0; i<SMM_MXPRMS; i++) 
        PsiPrms[i] = UserPrms[i];

    // new random number generator using gsl

	if (r == NULL)
	    r = gsl_rng_alloc(gsl_rng_ranlux);  
    const int seed = 1;
    gsl_rng_set(r, seed);

    /* We allocate space for the Ns sensors, users must call
       deploy_sensors to set the number of sensor to move 
    */
	
	if (u==NULL && v==NULL) {
	    u    = (double *)malloc((size_t)Ns*sizeof(double));
    	v    = (double *)malloc((size_t)Ns*sizeof(double));
    	Xmed = (double *)malloc((size_t)Ns*sizeof(double));
    	Ymed = (double *)malloc((size_t)Ns*sizeof(double));
	}
  
    switch(Streamfunction) {
    case SMM_ALTCHAN:   SMM_Velocities = &alternating_channels; break;
    case SMM_BLNKWAW:   SMM_Velocities = &blinking_wavenumbers; break;
    case SMM_PULSVOR:   SMM_Velocities = &pulsating_vortex; break;
    case SMM_STDYVOR:   SMM_Velocities = &steady_vortex; break;
	case SMM_MNDRJET:
		phaseshift = PsiPrms[5];	// Time shift. The caller provides it.
		current_time += phaseshift;
		SMM_Velocities = &meandering_jet;
	break;
	case SMM_RNDWALK:   SMM_Velocities = &random_walk; break; 	
	case SMM_RANDOM_MNDRJET: {
        int seed = (int)PsiPrms[5];
        gsl_rng_set (r, seed);
        U       = PsiPrms[6];
        lambda  = PsiPrms[7];
        ur      = U * gsl_ran_ugaussian(r);
        vr      = U * gsl_ran_ugaussian(r);
        expfact = exp(-lambda*dt);
        sqrt_expfact2 = U*sqrt(1.0-expfact*expfact);
    }
	SMM_Velocities = &meandering_jet; break;   
    case SMM_TAYLJET:   {
		/*
		The parameters are:
	    R    = PsiPrms[0];       RMS displacement of random walkers
    	Xmin = PsiPrms[1];       Reflecting barrier for random walk
	    Xmax = PsiPrms[2];       Reflecting barrier for random walk
	    Ymin = PsiPrms[3];       Reflecting barrier for random walk
	    Ymax = PsiPrms[4];       Reflecting barrier for random walk
	    seed = (int)PsiPrms[5];  Random number generator seed
	    U    = PsiPrms[6];       Max jet speed
		*/
        int seed = (int)PsiPrms[5];
        gsl_rng_set (r, seed);
        SMM_Velocities = &taylor_dispersion_jet; 
		break;
    }
    default:
        fprintf(stderr,"\nSMM ERROR: Unknown Streamfunction!\n\n");
        exit(1);
    }

    if (Integrator == SMM_RK2) {
        SMM_Move_Sensors = &rk2;
    } else {
        fprintf(stderr,"\nSMM ERROR: Unknown Integrator!\n\n");
        exit(2);
    }
}


/* Runge-Kutta II order integrator (midpoint method) */
/* The coordinate of the tracers in posX and posY given in input are
   overwritten with the new positions at current_time+dt. Also,
   current_time is incremented by dt.
*/
void rk2(double *posX, double *posY)
{
    /* Half step */
    SMM_Velocities(posX, posY, u, v, current_time);
    for (int i=0; i<deployed; i++) {
        Xmed[i] = posX[i] + 0.5*dt*u[i];
        Ymed[i] = posY[i] + 0.5*dt*v[i];
    }
    /* Full step */
    SMM_Velocities(Xmed, Ymed, u, v, current_time+0.5*dt);
    for (int i=0; i<deployed; i++) {
        posX[i] = posX[i] + dt*u[i];
        posY[i] = posY[i] + dt*v[i];
    }

    /* Zero-Mean Ornstein Uhlenbeck Process 
       i.e. let X(t) a standard brownian motion we have for a>0
       book def:  V(t) = e^(-at/2)X(e^(at))
       u_t+1 = u_t * e^(-at) + U * sqrt( 1 - e^(-2at))* X(t);
    */
    if (Sf == SMM_RANDOM_MNDRJET)
    {
        ur = ur * expfact + sqrt_expfact2*gsl_ran_ugaussian(r);
        vr = vr * expfact + sqrt_expfact2*gsl_ran_ugaussian(r);
        for (int i=0; i<deployed; i++) {
            posX[i] += dt*ur;
            posY[i] += dt*vr;
        }
    }   

    current_time += dt;
}

/* This implements the streamfunction (from A. Bower, J. Phys. Oceanography, 
21 (1991) 173; Cencini, et. al J. Phys. Ocean., 29 (1999) 2578)
psi = 1- tanh((y - B(t)*sin(k*(x-ct)))/sqrt(1 + B(t)^2*k^2*cos^2(k*(x-c*t))))
con B(t) = A + ep*cos(om*t)
A  is PsiPrms[0]  (suggested: 0.1)
c  is PsiPrms[1]
k  is PsiPrms[2]
om is PsiPrms[3]
ep is PsiPrms[4]
u = -dpsi/dy
v =  dpsi/dx
*/
void meandering_jet(double *posX, double *posY, double *u, double *v, double t)
{
    double A  = PsiPrms[0];
    double c  = PsiPrms[1];
    double k  = PsiPrms[2];
    double om = PsiPrms[3];
    double ep = PsiPrms[4];
    double sinxct, cosxct, sqrtcos, sech, sech2;
    double B = A + ep*cos(om*t);

    for (int i=0; i<deployed; i++) {
        sinxct  = sin(k*(posX[i] - c*t));
        cosxct  = cos(k*(posX[i] - c*t));
        sqrtcos = sqrt(B*B*k*k*cosxct*cosxct + 1);
        sech    = 1./cosh((posY[i] - B*sinxct)/sqrtcos);
        sech2   = sech*sech;
        u[i] = sech2/sqrtcos;
        v[i] = (B*k*cosxct/sqrtcos - B*B*k*k*k*cosxct*sinxct*(posY[i] - B*sinxct) /
            (sqrtcos*sqrtcos*sqrtcos) ) * sech2;
    }
}

/* This implements a brownian motion with gaussian jumps having zero mean and variance R^2.
   Boundary conditions are reflecting ones.
   R is PsiPrms[0] 
   The domain is a rectangle between Xmin, Xmax, Ymin, Ymax (respectively PsiPrms[1] to [4]) 
*/
void random_walk(double *posX, double *posY, double *u, double *v, double t)
{
    double R    = PsiPrms[0];
    double Xmin = PsiPrms[1];
    double Xmax = PsiPrms[2];
    double Ymin = PsiPrms[3];
    double Ymax = PsiPrms[4];
    double x, y;

    for (int i=0; i<deployed; i++) {
        x = posX[i] + R*gsl_ran_ugaussian(r);
        y = posY[i] + R*gsl_ran_ugaussian(r);
        if (x<Xmin) x = 2.*Xmin-x;
        if (x>Xmax) x = 2.*Xmax-x;
        if (y<Ymin) y = 2.*Ymin-y;
        if (y>Ymax) y = 2.*Ymax-y;
        posX[i] = x;
        posY[i] = y;
    }
}


/* This implements the streamfunction
psi = A*[sin(om*t)sin(pi*k*x) - cos(om*t)cos(pi*k*y)]sin(pi*x)sin(pi*y)/pi
A  is PsiPrms[0]  (suggested: 0.0225)
om is PsiPrms[1]  (suggested: 0.5)
k  is PsiPrms[2]  (suggested: 8; should be integer)
u = -dpsi/dy
v =  dpsi/dx
*/
void alternating_channels(double *posX, double *posY, 
			  double *u, double *v, double t){
  int i;
  double st, ct, SS, LS;
  double A  = PsiPrms[0];
  double om = PsiPrms[1];
  double k  = PsiPrms[2];

  st = sin(om*t);
  ct = cos(om*t);
  for(i=0; i<deployed; i++){
    SS = st*sin(M_PI*k*posX[i]) - ct*cos(M_PI*k*posY[i]);
    LS = sin(M_PI*posX[i])*sin(M_PI*posY[i]);
    u[i] = -A*(k*ct*sin(M_PI*k*posY[i])*LS 
	       + SS*sin(M_PI*posX[i])*cos(M_PI*posY[i]));
    v[i] =  A*(k*st*cos(M_PI*k*posX[i])*LS 
	       + SS*cos(M_PI*posX[i])*sin(M_PI*posY[i]));
  }
}


/* This implements the streamfunction
psi = A*[sin(om*t)sin(pi*k1*x)sin(pi*k1*y) + cos(om*t)sin(pi*k2*x)sin(pi*k2*y)]/pi
A  is PsiPrms[0]  (suggested: 0.0123)
om is PsiPrms[1]  (suggested: 0.5 to 1.5)
k1 is PsiPrms[2]  (suggested: 7; should be integer)
k2 is PsiPrms[3]  (suggested: 8; should be integer)
u = -dpsi/dy
v =  dpsi/dx
*/
void blinking_wavenumbers(double *posX, double *posY, 
			  double *u, double *v, double t){
  int i;
  double st, ct;
  double A  = PsiPrms[0];
  double om = PsiPrms[1];
  double k1 = PsiPrms[2];
  double k2 = PsiPrms[3];
  
  st = sin(om*t);
  ct = cos(om*t);
  for(i=0; i<deployed; i++){
    u[i] = -A*(k1*st*sin(M_PI*k1*posX[i])*cos(M_PI*k1*posY[i]) 
	  + k2*ct*sin(M_PI*k2*posX[i])*cos(M_PI*k2*posY[i]));
    v[i] =  A*(k1*st*cos(M_PI*k1*posX[i])*sin(M_PI*k1*posY[i]) 
	  + k2*ct*cos(M_PI*k2*posX[i])*sin(M_PI*k2*posY[i]));
  }
}


/* This implements the streamfunction
psi = A*[1+ecc*(sin(om*t)sin(pi*x) + cos(om*t)sin(pi*y))]sin(pi*x)sin(pi*y)/pi
A   is PsiPrms[0]  (suggested: 0.1)
om  is PsiPrms[1]  (suggested: 0.2)
ecc is PsiPrms[1]  (suggested: 0.5)
u = -dpsi/dy
v =  dpsi/dx
*/
void pulsating_vortex(double *posX, double *posY, 
		      double *u, double *v, double t){
  int i;
  double st, ct, sx, sy, cx, cy;
  double A   = PsiPrms[0];
  double om  = PsiPrms[1];
  double ecc = PsiPrms[2];
  
  st = sin(om*t);
  ct = cos(om*t);
  for(i=0; i<deployed; i++) {
    sx = sin(M_PI*posX[i]);
    cx = cos(M_PI*posX[i]);
    sy = sin(M_PI*posY[i]);
    cy = cos(M_PI*posY[i]);
    u[i] = -A*(2.*ecc*ct*sx*cy*sy + (ecc*st*sx*sx+sx)*cy);
    v[i] =  A*(ecc*ct*cx*sy*sy + (2.*ecc*st*cx*sx+cx)*sy);
  }
}

/* This implements the streamfunction
psi = A*sin(pi*x)sin(pi*y)/pi
A is PsiPrms[0]  (suggested: 0.1)
u = -dpsi/dy
v =  dpsi/dx
*/
void steady_vortex(double *posX, double *posY, 
		     double *u, double *v, double t){
  int i;
  double A = PsiPrms[0];
  
  for(i=0; i<deployed; i++) {
    u[i] = -A*sin(M_PI*posX[i])*cos(M_PI*posY[i]);
    v[i] =  A*cos(M_PI*posX[i])*sin(M_PI*posY[i]);
  }
}


/* This implements the streamfunction of 'meandering_jet' without the meanders */
/* that is:                                                                    */
/*  psi = -U*tanh(y)                                                           */
void taylor_dispersion_jet(double *posX, double *posY, 
			   double *u, double *v, double t){
    double U = PsiPrms[6];
    double sech, sech2;

    /* There is a brownian motion superimposed onto the deterministic motion of the jet */
    random_walk(posX, posY, u, v, t);
    for (int i=0; i<deployed; i++) {
        sech  = 1./cosh(posY[i]);
        sech2 = sech*sech;
        u[i]  = U*sech2;
	v[i]  = 0.;
    }

}
