#!/bin/csh
#____________________________ ATOMS 
setenv NATOMS                           1000    # Number of atoms in the problem
setenv type1                            Si      # Atomic types - add a line for each type
# setenv type 2    C

#_____________________________ ART

# The event_type determines the nature of the event. Use either 'NEW' for generating new events
#                    'REFINE_SADDLE' when further converging a saddle point
#                    'REFINE_AND_RELAX', to refine at the saddle and check
#
# To refine a state near saddle:
#       1. assign the file containing the configuration to REFCONFIG,  
#       2. remove restart.dat
#       3. set the new converging parameters and set REFINE_SADDLE or REFINE_AND_RELAX and launch the code
#
# Note that if a file "restart.dat" exists, a "new" simulation will
# continue the event that was being generated when the code was
# stopped. This is particularly useful with ab initio calculations. 
setenv EVENT_TYPE                     REFINE_AND_RELAX  
                                           
# ENERGY_CALC defines the potential used, use either SWP (Stillinger-Weber) or DFT 
setenv ENERGY_CALC                    SWP  

# The temperature used for Metropolis accept/reject in eV. If the
# temperature is negative, we never move and sampled from the minimum.
setenv Temperature                   0.4  

# Set the maximum number of events in the simulation
setenv Max_Number_Events             1 
setenv Increment_Size                0.09   # Overall scale for the increment moves (in Ang)

#_____________________________ Initial deformation

# The initial deformation can be made around a single atom (local)  or by deforming randomly the whole box. 
# If a local event is chosen, the the radius of the deformation must be indicated. 

setenv Initial_Step_Size             0.05   # Size of initial displacement, (in Ang)
setenv Type_of_Events               local   # Activation: global, local
setenv Radius_Initial_Deformation     3.5   # Cutoff for local and list_local (in Ang)

# It is also possible to select the atom at the center of the event
# and also, often useful for ab-initio calculations, whether a general
# small random deformation is applied to the box in order to
# facilitate convergence.

# setenv Central_Atom                     1   # Number of the atom around which the initial move takes place
# setenv sym_break_dist               0.001   # Breaks the symmetry of the crystal by randomly displacing
                                            # all atoms by this distance

#_____________________________ HARMONIC WELL

# After the initial deformation, the system is pushed along this
# direction with a small perpendicular relaxation to prevent crashing
# against other atoms.

# This number, "Max_Pep_Moves_Basin" is typically between 1 and 3 and
# must adjusted as a function of the sysystem's properties.

setenv Basin_Factor                   2.1   # Factor multiplying Increment_Size for leaving the basin
setenv Max_Perp_Moves_Basin             3   # Maximum number of perpendicular steps leaving basin

# Because Lanczos call are expensive, we wait a few steps before starting them.
setenv Min_Number_KSteps                3   # Min. number of ksteps before calling lanczos 
setenv Max_Iter_Basin                  12   # Maximum number of iteraction for leaving the basin (kter)

# The threshold must be adapted to each system but one digit precision is generally sufficient. 
setenv Eigenvalue_Threshold         -0.5   # Eigenvalue threshold for leaving basin


#_____________________________ ACTIVATION 

# In this phase, the system is pushed up along the direction of
# negative curvature while the energy is minimized in the
# perpendicular direction. Iterations are stopped then the total force
# is below a given threshold.

# The important parameter, here, is the force threshold for deciding
# that a saddle point has been reached. For ab initio calculations,
# this parameter can be set to 0.20 or 0.25 eV/A. For empirical
# potential, 0.10, for example.

setenv Exit_Force_Threshold           0.10   # Threshold for convergence at saddle point (eV/A)

# These two parameters can be kept unchanged. 
setenv Activation_MaxIter             800   # Maximum number of iteraction for reaching the saddle point
setenv Max_Perp_Moves_Activ            20   # Maximum number of perpendicular steps during activation

# The push over the saddle point to minimize into a new
# state. Typically between 0.15 and 0.3 depending on the structure of
# the energy landscape.
setenv Prefactor_Push_Over_Saddle     0.15   # Fraction of displacement over the saddle

# These last two paramters are used to try to minimize the number of
# force calls in ab initio calculations. Because they presuppose
# specific properties for saddles, they can bias the results. 
#
# By default, they are set small enough to have not impact. 
setenv delta_threshold               0.01   # if delta_e < delta_thr .and. delr < delr_thr, in the 
setenv delr_threshold                0.01   # convergence stage, then we kill the event (eV and Ang)


#_____________________________ LANCZOS 

# Overall, the parameters for the Lanczos steps are rather sturdy and
# do not need to be adapted except for the small displacement used for the numerical derivation of the force.

# While delta_disp_lanczos can be set to 0.001 A or less for empirical
# potential, it is around 0.01 for ab initio, the exact value
# depending on the precision on the force obtained in the calculations
setenv delta_disp_Lanczos            0.01   # Step of the numerical derivative of forces in lanczos (Ang)

setenv Lanczos_of_minimum         .False.   # Calculation of the Hessian for each minimum
setenv Lanczos_SCLoop                   1   # of iterations in the lanczos Self consistent loop
setenv Lanczos_collinear             0.73   # Exit criteria of the lanczos Self consistent loop
setenv Number_Lanczos_Vectors_H        13   # Number of vectors included in lanczos procedure in the Harmonic well
setenv Number_Lanczos_Vectors_C        12   # Number of vectors included in lanczos procedure in convergence


#_____________________________ DIIS
# We do not use DIIS anymore. So leave this section unchanged.
setenv Use_DIIS                   .False.   # Use DIIS for the final convergence to saddle
setenv Iterative                  .False.   # Iterative use of Lanczos & DIIS
setenv Inflection                     100   # Number of Lanczos steps after an inflection in the eigenvalue
setenv DIIS_Force_Threshold           0.4   # Force threshold for call DIIS
setenv DIIS_Memory                      5   # Number of vectors kepts in memory for algorithm
setenv DIIS_Step_size                0.02   # prefactor multiplying forces
setenv FACTOR_DIIS                    4.0   # times Increment_Size, max allowed diis step size
setenv MAX_DIIS                       100   # max diis iterations per call
setenv DIIS_Check_Eigenvector      .True.   # Check that the final state is indeed a saddle

#______________________ INPUT/OUTPUT 
# Not much to change, here, example for two flags particularly useful with
# an ab initio simulation where each force calculation counts.

# The first flag, Write_restart_file, ensures that after each iteraction,
# all positions and details are stored to relaunch from it, whether during
# minimization or activation. When the restart file exists, the code
# automatically uses it (in mode NEW), to continue where the simulation was
# interrupted.   
setenv Write_restart_file          .True.   # It is useful only for ab-initio

# It is also possible to save the configuration as every iteraction step,
# for further convergence, for example. While this rapidely becomes
# out-of-hand with empirical potential, it is often useful with more heavy
# calculations. 
setenv Save_Conf_Int               .False.   # Save the configuration at every step?

# The file names themselves can be left unchanged.
setenv FILECOUNTER            filecounter   # File tracking  the file (event) number - facultative
setenv LOGFILE                   log.file   # General output for message
setenv EVENTSLIST             events.list   # list of events with success or failure
setenv RESTART                restart.dat  # current data for restarting event
setenv Write_xyz                   .True.   # Writes min. and sad. configurations in .xyz format.
setenv REFCONFIG             ExcitedState   # Reference configuration for refine saddle. Without ext.
setenv INITDIR                    initdir   # Configuration for GUESS_direction. Without ext.


###### RUN THE SIMULATION #######################################################################
../../source/artdft
