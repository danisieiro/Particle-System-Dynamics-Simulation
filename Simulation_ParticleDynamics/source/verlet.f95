

subroutine verlet (np,dt,rx,ry,rz,vx,vy,vz,ax,ay,az,epot,ecin,etot,dfiv,d2fiv)

    use def_precision

    implicit none

        
    integer (kind = entero), parameter :: npmax=500
    integer (kind=entero) :: np
    real (kind = doblep) :: dt,dt12,dt2,epot,ecin,dfiv,d2fiv,etot
    real (kind = doblep), dimension(npmax) :: rx,ry,rz,vx,vy,vz,ax,ay,az
    integer (kind = entero) :: i,j
    

    dt12 = dt/2.d00
    dt2 = dt * dt12


    do i=1,npmax
        rx(i) = rx(i) + vx(i) * dt + ax(i) * dt2
        ry(i) = ry(i) + vy(i) * dt + ay(i) * dt2
        rz(i) = rz(i) + vz(i) * dt + az(i) * dt2

        vx(i) = vx(i) + ax(i) * dt12
        vy(i) = vy(i) + ay(i) * dt12
        vz(i) = vz(i) + az(i) * dt12
        
    end do


    call potlj (np,rx,ry,rz,ax,ay,az,epot,dfiv,d2fiv)

    do j=1,npmax
        vx(j) = vx(j) + ax(j) * dt12
        vy(j) = vy(j) + ay(j) * dt12
        vz(j) = vz(j) + az(j) * dt12
        
    end do

    ecin = 1.d00/2.d00 * (sum(vx*vx + vy*vy + vz*vz))
    etot=ecin+epot
    
    return

    end subroutine verlet