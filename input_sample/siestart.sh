#!/bin/csh 

#bharat added STARTS
setenv ENERGY_CALC GAU     #choose: DFT or SWP (Stillinger-Weber) or GAU
#setenv ENERGY_CALC DFT     #choose: DFT or SWP (Stillinger-Weber) or GAU
# bharat added ENDS


setenv Temperature 0.4   # Temperature in kcal/mol, if negative always reject the event bharat 0.4 >> 0.1
setenv EVENT_TYPE  NEW     # Either 'NEW', 'REFINE_SADDLE' when further converging a saddle point
                           # Or "REFINE_AND_RELAX", to refine at the saddle
			   # and check the final minimum

#setenv LABEL_SIESTA SIGA_INT
setenv NATOMS     8    # Number of atoms in the problem bharat 217 > 2

setenv Prefactor_Push_Over_Saddle 0.3 # The prefactor for pushing over the saddle point, fraction of distance from initial minimum to saddle point 0.15 (default) 0.3 (siesta)

setenv type1      C       # Name of the atomic types - used for writing out xmol-type file
setenv type2      O    
setenv type3      N
setenv type4      H     
 

######  ART ##############################################################################################

setenv Max_Number_Events                  10   # Maximum number of events
#setenv Type_of_Events                    list_local   # Initial move for events - global or local list_local
setenv Type_of_Events                    global
#setenv Radius_Initial_Deformation          3.0   #Cutoff for local-move (in angstroems) bharat 3> 2 NOT NEEDED FOR "global"
#setenv Central_Atom                        30   # Number of the atom # around which the initial move takes place NOT NEEDED FOR "global"

setenv Eigenvalue_Threshold             -0.2   # Eigenvalue threshold for leaving basin
setenv Exit_Force_Threshold              0.1   # Threshold for convergence at saddle point

setenv Increment_Size                   0.1  # Overall scale for the increment moves in activation 0.06 >> 0.01
setenv Initial_Step_Size                0.1    # Size of initial displacement, in A 0.1 >>> 0.001(default) 0.8 from Siesta

setenv sym_break_dist                    0.2   # Breaks the symmetry of the bharat 0.1 >>> 0.3 (Siesta) 0(default)
                                               # crystal by randomly displacing all atoms by this distance

setenv Max_Perp_Moves_Basin                3   # Maximum number of perpendicular steps leaving basin bharat 2>>3
setenv Min_Number_KSteps                   3   # Min. number of ksteps before calling lanczos  bharat 2>>3
setenv Basin_Factor                      2.1   # Factor multiplying Increment_Size for leaving the basin bharat 3.0 >>2.1

# bharat starts
#setenv Lanczos_of_minimum             .true.   # Calculation of the Hessian for each minimum
setenv Lanczos_SCLoop                     5    # Number of iterations in the lanczos Self consistent loop default 20 >> 1 (default) 
setenv Activation_MaxIter                 400 # 400
setenv delta_threshold                    4.0  # Energy threshold during Lanczos 
# bharat ends

setenv Max_Perp_Moves_Activ                8   # Maximum number of perpendicular steps during activation 12 (default) 8 (siesta)
setenv Force_Threshold_Perp_Rel          0.5   # Threshold for perpendicular relaxation 0.5 (Siesta)
setenv Max_Iter_Basin                    20   # Maximum number of iteraction for leaving the basin (kter) 100 (default) 200(Siesta)

################## Direction inversion in iterative subspace
#setenv Use_DIIS                       .false.   # Use DIIS for the final convergence to saddle
#setenv DIIS_Force_Threshold            0.10    # Force threshold for convergence bharat 0.25 >>0.10
#setenv DIIS_Memory                       5     # Number of vectors kepts in memory for algorithm
#setenv DIIS_MaxIter                    50      # Maximum number of iteractions for the DIIS scheme bharat 200 >>50
#setenv DIIS_Check_Eigenvector         .false.   # Check that the final state is indeed a saddle
#setenv DIIS_Step_Size                 0.005     # Step size for the position

############### Input              #######################################################################
setenv FILECOUNTER      filecounter            # File tracking  the file (event) number - facultative
setenv REFCONFIG        refconfig.dat             # Reference configuration (actual local minimum)
############### Output             #######################################################################
setenv LOGFILE             log.file            # General output for message
setenv EVENTSLIST       events.list            # list of events with success or failure
setenv RESTART          restart.dat            # current data for restarting event

############### Run the simulation #######################################################################

# echo this is test1
#./siestart_old
#/global/scratch/sharmabh/art/ART/art_dft_v3.1-distr-2/art_dft_v3.1-distr/source/artdft
#/global/scratch/sharmabh/art/ART/art_dft_v3.1-distr-2/art_dft_v3.1-distr/source/gaussian_art/artdft

# #rhf/3-21g nosymm opt
#/global/scratch/sharmabh/art/ART/art_dft_v3.1-distr-2/art_dft_v3.1-distr/source/gaussian_art_v1/artdft

# add extra threshold energy threshold
#/global/scratch/sharmabh/art/ART/art_dft_v3.1-distr-2/art_dft_v3.1-distr/source/gaussian_art_v4_rand/artdft
#

# v6 opt
#/global/scratch/sharmabh/art/ART/art_dft_v3.1-distr-2/art_dft_v3.1-distr/source/gaussian_art_v6_opt/artdft

#v7 opt
#/global/scratch/sharmabh/art/ART/art_dft_v3.1-distr-2/art_dft_v3.1-distr/source/gaussian_art_v7/artdft

#PSI v7
/global/scratch/jwhatley/ART_project/gaussian_art_v8/artdft
#
# for negative forces opt=5steps
#/global/scratch/sharmabh/art/ART/art_dft_v3.1-distr-2/art_dft_v3.1-distr/source/gaussian_art_v3_Force_Negative/artdft

#using fire minimization on gaussian
#/global/scratch/sharmabh/art/ART/art_dft_v3.1-distr-2/art_dft_v3.1-distr/source/gaussian_art_v2_GAU_fire/artdft
# echo this is test2
