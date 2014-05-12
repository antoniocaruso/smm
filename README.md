

Simple Jet Meandering Model (SMM)
===========================

SmmLib:

This library advects sensors according to a time-dependent, one-degree-of-freedom 
Hamiltonian (the streamfunction). In smm.c several different streamfunctions are 
implemented. The streamfunction used in [1] is SMM_MNDRJET. 

The first call to the library must be SMM_init, it specify the streamfunction used, the MAXIMUM number
of sensor that will be deployed (eventually in different times) with calls to SMM_deploy_nodes, 
the parameters of the streamfunction, the integrator (SMM_RK2), and the timestep size.

The physical positions of sensor nodes are stored in two double arrays managed by the client. These
arrays must be passed to SMM_Move_Sensors to move the nodes according to the velocity field defined
by the streamfunction, advancing the simulation time by the specified timestep.

The initial position of the nodes is defined by the client. Before any call to SMM_Move_Sensors, nodes
must be deployed in the appropriate domain (the domain is specified with the streamfunction) with a 
call to:

void SMM_deploy_nodes(int n);

The separation between the maximum number of sensors specified in SMM_init and the explicit call to
SMM_deploy_nodes can be used to implements multiple sensor deployments (if required).

The current time internal to the mobility model can be inspected calling:

double SMM_get_current_time(void);

The time scale depends on the streamfunction, the initial paramenters and the timestep.

For [1] with params  = { A = 1.2, c = 0.12, k = 2.0*M_PI/7.5, 0.4, epsilon = 0.3, phase }, 

phase is a random phaseshift in the interval [ 0, T] with T = 2*M_PI/ck, 
this is used to change the relative position of the first deployment
with respect to the streamfunction. 

With  timestep = 0.01 a single call of SMM_Move_Sensors is equivalent 
to 25.92s = timestep *24*60*60 * 0.03

see [1] for further infomation.

------------------

[1] A. Caruso, F. Paparella, Luiz Vieira, Melike Erol , Mario Gerla, 
    The Meandering Current Mobility Model and its impact on Underwater 
	Mobile Sensor Networks, INFOCOM 2008, Phoenix, AZ, USA.
	
