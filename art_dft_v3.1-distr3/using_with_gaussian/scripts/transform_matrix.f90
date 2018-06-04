program move

        implicit none 
        real(kind=8), dimension(4,1) :: vec
        real(kind=8), dimension(4,1) :: vec_moved
        real(kind=8), dimension(4,4) :: rx
        real(kind=8), dimension(4,4) :: ry
        real(kind=8), dimension(4,4) :: rz
        real(kind=8), dimension(4,4) :: txyz
        real(kind=8), dimension(4,4) :: T
        real(kind=8), parameter :: pi = 3.141592653
        real(kind=8), parameter :: deg_to_rad = pi/180.0
        real :: thetax, thetay, thetaz
        real :: radx, rady, radz
        real :: tx, ty, tz
        real :: cx, cy, cz 
        real :: sx, sy, sz
        real :: vec1, vec2, vec3
        integer :: i,j
        
        write(*,*)"Enter the components of the vector you want to move; Vx, Vy, Vz: "
        read(*,*) vec1, vec2, vec3

        write(*,*)"Enter the angles (in degrees) you want to rotate the vector by; Rx, Ry, Rz: "
        read(*,*) thetax, thetay, thetaz
        
        write(*,*)"Enter the distances you want to translate the vector by; tx, ty, tz: "
        read(*,*) tx, ty, tz

        radx = thetax * deg_to_rad

        rady = thetay * deg_to_rad

        radz = thetaz * deg_to_rad

        cx = cos(radx)
        cy = cos(rady)
        cz = cos(radz)

        sx = sin(radx)
        sy = sin(rady)
        sz = sin(radz) 

        rx = reshape((/1.0,0.0,0.0,0.0,0.0,cx,sx,0.0,0.0,-sx,cx,0.0,0.0,0.0,0.0,1.0/),(/4,4/))

        ry = reshape((/cy,0.0,-sy,0.0,0.0,1.0,0.0,0.0,sy,0.0,cy,0.0,0.0,0.0,0.0,1.0/),(/4,4/))

        rz = reshape((/cz,sz,0.0,0.0,-sz,cz,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0/),(/4,4/))

        txyz = reshape((/1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,tx,ty,tz,1.0/),(/4,4/))

        vec = reshape((/vec1,vec2,vec3,1.0/),(/4,1/))

        T= matmul(matmul(matmul(rx, ry), rz),txyz)

        vec_moved = matmul(T, vec)

        write(*,*)"This is the moved vector: "
        do i = 1,3
                do j = 1,1
                write(*,'(f10.6)') vec_moved(i,j)
                enddo
        enddo

end program move

