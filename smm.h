/* SmmLib header
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
 */


#ifndef SMM
#define SMM

/* Max number of parameters for the streamfunction */
#define SMM_MXPRMS 10

/* Available Streamfunctions */

#define SMM_ALTCHAN   1
#define SMM_BLNKWAW   2
#define SMM_PULSVOR   3
#define SMM_STDYVOR   4
#define SMM_MNDRJET   5
#define SMM_RANDOM_MNDRJET 6
#define SMM_TAYLJET   7
#define SMM_RNDWALK   100
/* Available Integrators*/

#define SMM_RK2       1


/* Function Prototipes */

typedef void (*Streamfunction_Funptr)(double *posX, double *posY, double *u, double *v, double t);
typedef void (*Integrator_Funptr)(double *posX, double *posY);
extern Integrator_Funptr SMM_Move_Sensors;

void SMM_init(int Nsensors, int Streamfunction, double *UserPrms, int Integrator, double TimeStep);
void SMM_deploy_nodes(int n);
double SMM_get_current_time(void);

#endif 
