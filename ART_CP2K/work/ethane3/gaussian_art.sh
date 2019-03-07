#!/bin/csh

#Sets the environment variables for the molecules in question

setenv Temperature 0.4      # Temperature in kcal/mol, if negative always reject the event  0.4 >> 0.1
setenv ENERGY_CALC GAU      # Choose: DFT, SWP (Stillinger-Weber) or GAU
setenv EVENT_TYPE  NEW      # Either 'NEW', 'REFINE_SADDLE' when further converging a saddle point
                            # Or "REFINE_AND_RELAX", to refine at the saddle and check the final minimum

setenv NATOMS     8         # Number of atoms in the problem

setenv Prefactor_Push_Over_Saddle 0.3 # The prefactor for pushing over the saddle point, fraction of distance from initial minimum to saddle point 0.15 (default) 0.3 (siesta)

######  ART ##############################################################################################

setenv Max_Number_Events                      1   # Maximum number of events
setenv Type_of_Events                    global   # Initial move for events - global or local list_local
setenv Radius_Initial_Deformation           3.0   # Cutoff for local-move (in angstroems)  3> 2 NOT NEEDED FOR "global"
setenv Central_Atom                           1   # Number of the atom # around which the initial move takes place NOT NEEDED FOR "global"

setenv Eigenvalue_Threshold                -0.2   # Eigenvalue threshold for leaving basin (units eV/angstrom)
setenv Exit_Force_Threshold                 0.1   # Threshold for convergence at saddle point

setenv Increment_Size                       0.1   # Overall scale for the increment moves in activation 0.06 >> 0.01
setenv Initial_Step_Size                    0.1   # Size of initial displacement, in A 0.1 >>> 0.001(default) 0.8 from Siesta

setenv sym_break_dist                       0.2   # Breaks the symmetry of the  0.1 >>> 0.3 (Siesta) 0(default)
                                                  # crystal by randomly displacing all atoms by this distance

setenv Max_Perp_Moves_Basin                   3   # Maximum number of perpendicular steps leaving basin  2>>3
setenv Min_Number_KSteps                      3   # Min. number of ksteps before calling lanczos   2>>3
setenv Basin_Factor                         2.1   # Factor multiplying Increment_Size for leaving the basin  3.0 >>2.1

#setenv Lanczos_of_minimum               .true.   # Calculation of the Hessian for each minimum
setenv Lanczos_SCLoop                         5   # Number of iterations in the lanczos Self consistent loop default 20 >> 1 (default)
setenv Activation_MaxIter                   400   # 400
setenv delta_threshold                      4.0   # Energy threshold during Lanczos

setenv Max_Perp_Moves_Activ                   5   # Maximum number of perpendicular steps during activation 12 (default) 8 (siesta)
setenv Force_Threshold_Perp_Rel             0.5   # Threshold for perpendicular relaxation 0.5 (Siesta)
setenv Max_Iter_Basin                        20   # Maximum number of iteraction for leaving the basin (kter) 100 (default) 200(Siesta)

setenv Write_xyz                         .true.   # Generates .xyz files


###### Automation parameters #########

setenv Strategy_of_Search                     0   # '0' --> Proceed randomly (as it was doing already), '1' --> Follow the vector, '2' --> Avoid the vector
setenv Odds_follow_avoid		      0   # Odds that ART will either 'follow' or 'avoid'
setenv Shared_History_Filename              ../shared_log/ethane3.log # Name the shared history file
setenv Shared_History_To_Read               ../shared_log/xxx.log # Shared vector file to read from
setenv nmin_read			    0 # nmin read from vector file
setenv nsad_read			    0 # nsad read from vector file
setenv natoms_read			    0 # natoms read from vector file
setenv natoms_correspond		    0 # natoms that correspond

############### Input              #######################################################################

setenv FILECOUNTER      filecounter               # File tracking  the file (event) number - facultative
setenv REFCONFIG        refconfig.dat             # Reference configuration (actual local minimum)

############### Output             #######################################################################

setenv LOGFILE             log.file               # General output for message
setenv EVENTSLIST       events.list               # list of events with success or failure
setenv RESTART          restart.dat               # current data.sh for restarting event

############### Run the simulation #######################################################################

set art_location = "../../source/artgaussian"
$art_location


