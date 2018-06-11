program align

        !! This program calculates a transformation matrix that minimizes the RMSD between two sets of coordinates.
        !The first step of this transformation involves a translation. This is simply done by re-centering the molecules so that their
        !centroids coincide with the origin of the coordinate system.
        !Next, a covariance matrix is calculate using the re-centered molecules.
        !Then, the singular value decomposition (SVD) of this covariance matrix is calculated, which yields the factors of this
        !covariance, viz U, S and V (Here U and V are the left- and right- singular vectors of the covariance matrix.)
        ! According to the Kabsch algorithm, the product of U and V is the optimal rotation matrix that minimizes the RMSD between
        ! the two sets of vectors.

        implicit none

        real(kind=8),dimension(:), allocatable, target :: vec_A
        real(kind=8),dimension(:), allocatable, target :: vec_B
        real(kind=8),dimension(:), allocatable, target :: vec_A_moved
        real(kind=8),dimension(:), allocatable, target :: vec_B_moved
        
        real(kind=8), dimension(:), pointer :: vec_Ax, vec_Ay, vec_Az
        real(kind=8), dimension(:), pointer :: vec_Bx, vec_By, vec_Bz
        real(kind=8), dimension(:), pointer :: vec_A_movedx, vec_A_movedy, vec_A_movedz
        real(kind=8), dimension(:), pointer :: vec_B_movedx, vec_B_movedy, vec_B_movedz

        character(len=120), parameter :: VECA = "ethane_A.xyz"
        integer, parameter :: A = 21
        character(len=120), parameter :: VECB = "ethane_B.xyz" 
        integer, parameter :: B = 22
        character(len=120), parameter :: VECT = "ethane_transformed.xyz" 
        integer, parameter :: T = 23
        character(len=120), parameter :: RLOG = "rmsd.log" 
        integer, parameter :: RMSD = 24
        
        real(kind=8) :: cov(3,3), U(3,3), S(3), VT(3,3), R(3,3), work(20)
       
        ! Here cov is the covariance matrix, U and VT are the left- and
        !right- singular vectors of the cov, (VT is actually the transpose of V),
        !S is one of the factors of the cov. R is the rotation matrix, U.VT.
        ! Work is one of the arguments required by the LAPACK library subroutine that calculates the SVD
        !Documentation at: 
        !http://www.netlib.org/lapack/explore-html/d1/d7e/group__double_g_esing_ga84fdf22a62b12ff364621e4713ce02f2.html#ga84fdf22a62b12ff364621e4713ce02f2

        real(kind=8),dimension(:,:), allocatable :: vec_transformed
        character(len=40),dimension(:), allocatable :: typat
        character(len=40) :: molecule_A, molecule_B
        

        integer :: i,natoms,info
        character(len=20) :: string_natoms
        real :: sum_Ax,sum_Ay,sum_Az,cx_A,cy_A,cz_A
        real :: sum_Bx,sum_By,sum_Bz,cx_B,cy_B,cz_B
        real :: sum_rmsd_before, rmsd_before
        real :: sum_rmsd_after, rmsd_after, det_VT, det_U

        sum_Ax = 0.0
        sum_Ay = 0.0
        sum_Az = 0.0
        
        sum_Bx = 0.0
        sum_By = 0.0
        sum_Bz = 0.0

        open(A, file=VECA, action = 'read', status = 'unknown') ! Reading the first set of coordinates

        read(A,*) natoms
        read(A,*) molecule_A
       
        allocate(typat(natoms))

        allocate(vec_A(3*natoms))
        vec_Ax => vec_A(1:natoms)
        vec_Ay => vec_A(natoms+1:2*natoms)
        vec_Az => vec_A(2*natoms+1:3*natoms)

        allocate(vec_A_moved(3*natoms))
        vec_A_movedx => vec_A_moved(1:natoms)
        vec_A_movedy => vec_A_moved(natoms+1:2*natoms)
        vec_A_movedz => vec_A_moved(2*natoms+1:3*natoms)
        
        allocate(vec_B(3*natoms))
        vec_Bx => vec_B(1:natoms)
        vec_By => vec_B(natoms+1:2*natoms)
        vec_Bz => vec_B(2*natoms+1:3*natoms)

        allocate(vec_B_moved(3*natoms))
        vec_B_movedx => vec_B_moved(1:natoms)
        vec_B_movedy => vec_B_moved(natoms+1:2*natoms)
        vec_B_movedz => vec_B_moved(2*natoms+1:3*natoms)

        allocate(vec_transformed(natoms,3))
        

        do i = 1,natoms

                read(A,*) typat(i), vec_Ax(i), vec_Ay(i), vec_Az(i)
        
        enddo

        close(A)
        
        open(B, file=VECB, action = 'read', status = 'unknown') !Reading the second set of coordinates
        read(B,*) 
        read(B,*) molecule_B


        do i = 1,natoms

                read(B,*) typat(i), vec_Bx(i), vec_By(i), vec_Bz(i)

        enddo

        close(B)
        
        do i=1,natoms

        sum_rmsd_before = sum_rmsd_before + ((vec_Ax(i)-vec_Bx(i))**2 + (vec_Ay(i)-vec_By(i))**2 + (vec_Az(i)-vec_Bz(i))**2)  

        enddo

        rmsd_before = sqrt(sum_rmsd_before/natoms) !Calculating RMSD before structural alignment

        do i=1,natoms

        sum_Ax = sum_Ax + vec_Ax(i)
        sum_Ay = sum_Ay + vec_Ay(i)
        sum_Az = sum_Az + vec_Az(i)


        sum_Bx = sum_Bx + vec_Bx(i)
        sum_By = sum_By + vec_By(i)
        sum_Bz = sum_Bz + vec_Bz(i)
        
        enddo

        !Calculating the centroids for all the points

        cx_A = sum_Ax/natoms
        cy_A = sum_Ay/natoms
        cz_A = sum_Az/natoms

        cx_B = sum_Bx/natoms
        cy_B = sum_By/natoms
        cz_B = sum_Bz/natoms

       
        !Recentering the structures so that their centroids coincide with the origin of the coordinate system.

        do i = 1,natoms

                vec_A_movedx(i) = vec_Ax(i) - cx_A
                vec_A_movedy(i) = vec_Ay(i) - cy_A
                vec_A_movedz(i) = vec_Az(i) - cz_A
       
                vec_B_movedx(i) = vec_Bx(i) - cx_B
                vec_B_movedy(i) = vec_By(i) - cy_B
                vec_B_movedz(i) = vec_Bz(i) - cz_B
        
        enddo

        !Calculating the covariance matrix         

        cov = matmul(transpose(reshape((vec_A_moved),(/natoms,3/))),reshape((vec_B_moved),(/natoms,3/)))

        !Invoking the LAPACK library subroutine "dgesvd" that calculates the SVD of the covariance matrix

        call dgesvd('S','S',3,3,cov,3,S,U,3,VT,3,work,20,info)

        det_U = U(1,1)*(U(2,2)*U(3,3) - U(2,3)*U(3,2)) - U(1,2)*(U(2,1)*U(3,3) - U(2,3)*U(3,1)) + U(3,1)*(U(2,1)*U(3,2) - U(2,2)*U(3,1))
        det_VT = VT(1,1)*(VT(2,2)*VT(3,3) - VT(2,3)*VT(3,2)) - VT(1,2)*(VT(2,1)*VT(3,3) - VT(2,3)*VT(3,1)) + VT(3,1)*(VT(2,1)*VT(3,2) - VT(2,2)*VT(3,1))
        
        if(det_U*det_VT .LT. 0) then       ! A small correction to ensure a right-handed coordinate system

                VT(:,3) = -VT(:,3)              

        endif

        R = (matmul(U,VT))  !This is the optimal rotation matrix
 

        vec_transformed = matmul(reshape((vec_A_moved),(/natoms,3/)),R)
        
        
        open(T, file=VECT, action = 'write', status = 'unknown')
        
        write(string_natoms,'(i20)') natoms
        write(T,*) adjustl(string_natoms)
        write(T,'(A20)') trim(molecule_A) // trim('_transformed')

        do i = 1,natoms

                write(T,'((A1,2X), 3(f10.6,2X))') typat(i), vec_transformed(i,1), vec_transformed(i,2), vec_transformed(i,3)
        
        enddo

        close(T)       

        do i=1,natoms

                sum_rmsd_after = sum_rmsd_after + ((vec_transformed(i,1) - vec_Bx(i))**2 + (vec_transformed(i,2) - vec_By(i))**2 + (vec_transformed(i,3) - vec_Bz(i))**2)

        enddo

        rmsd_after = sqrt(sum_rmsd_after/natoms) !Calculating RMSD after structural alignment
        
        open(RMSD, file=RLOG, action ='write', status = 'unknown')
        write(RMSD,*) "This is the RMSD before alignment: ", rmsd_before
        write(RMSD,*) new_line('a')
        write(RMSD,*) "This is the RMSD after alignment: ", rmsd_after
        close(RMSD)

        deallocate(typat)
        deallocate(vec_A)
        deallocate(vec_A_moved)
        deallocate(vec_B)
        deallocate(vec_B_moved)
        deallocate(vec_transformed)

end program align

