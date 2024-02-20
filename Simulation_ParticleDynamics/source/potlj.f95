!*****************************************************************

!*****************************************************************

        subroutine potlj (np,rx,ry,rz,ax,ay,az,epot,dfiv,d2fiv)


        use def_precision

        implicit none


!   pi = numero pi
!   pl = lado de la caja
!   npmax = numero de particulas maxima (parametro), en este caso 500
!   rc = radio de corte, a partir del cual no se consideran interacciones entre pares (mitad del lado de la caja)
!   rc2 = radio de corte al cuadrado
!
!   rx,rx,rz = array de posiciones en las coordenadas x,y,z
!   vx,vy,vz = array de velocidades en las direcciones x,y,z
!   ax,ay,az = array de aceleraciones en las direcciones x,y,z
!
!   rijx,rijy,rijz = distancia entre dos pares en las coordenadas x,y,z
!
!   epot = energía potencial total
!
!   indices de bucles: i, j
!
!   dfiv = derivada del potencial con respecto al volumen
!   d2fiv = derivada segunda del potencial con respecto al volumen
!
!   corr_ener = correccion a la energia
!


!   DECLARAMOS LAS VARIABLES

        integer (kind = entero) :: np,pl
        integer (kind = entero), parameter :: npmax = 500
        real (kind = doblep), parameter :: pi=3.141592653589
        real (kind = doblep), dimension(npmax) :: rx,ry,rz,ax,ay,az
        real (kind = doblep), intent(out) :: epot,dfiv,d2fiv
        integer (kind = entero) :: i,j
        real (kind = doblep) :: factor,xnp,corr_sum_rvp,corr_sum_r2vpp,sum_r2vpp,sum_rvp
        real (kind = doblep) :: vol,rrx,rry,rrz,rijx,rijy,rijz,dis2,fmod,a2,a6,a12,aux,pli,rc,rc2,corr_ener
    
!   INTRODUCIR LOS VALORES NECESARIOS

        xnp = dble(npmax)
        pl = 10
        vol = pl * pl * pl
        pli=1.d00/dble(pl)
        rc = dble(pl)/2.d00
        rc2 = rc * rc
        epot = 0.d00
        sum_rvp = 0.d00
        sum_r2vpp = 0.d00
        ax = 0.d00
        ay = 0.d00
        az = 0.d00
        factor=pi*xnp*xnp/(vol*rc**3)
        corr_ener=8.d00*factor*(1.d00/(3.d00*rc**6)-1.d00)/3.d00
        corr_sum_rvp=16.d00*factor*(-2.d00/(3.d00*rc**6)+1.d00)
        corr_sum_r2vpp=16.d00*factor*(26.d00/(3.d00*rc**6)-7.d00)


!   EVALUAMOS LA DISTANCIA ENTRE PARES

        do i=1, npmax
            rrx = rx(i)
            rry = ry(i)
            rrz = rz(i)

            do j=i+1, npmax
                rijx = rrx - rx(j)
                rijy = rry - ry(j)
                rijz = rrz - rz(j)

                rijx = rijx - dble(pl)*dnint(rijx * pli)
                rijy = rijy - dble(pl)*dnint(rijy * pli)
                rijz = rijz - dble(pl)*dnint(rijz * pli)
                
                dis2 = rijx*rijx + rijy*rijy + rijz*rijz
                
        !   SOLO CALCULAMOS LA INTERACCION ENTRE PARES SEPARADOS POR MENOS DE rc
                if (dis2<=rc2) then
                    a2 = 1.d00/dis2
                    a6 = a2*a2*a2
                    a12 = a6*a6
                    epot = epot + a12 - a6
                    aux = -2.d00 * a12 + a6
                    sum_rvp = sum_rvp + aux
                    sum_r2vpp = sum_r2vpp + 26.d00 * a12 - 7.d00 * a6
                    fmod = -aux * a2
                    ax(i) = ax(i) + fmod * rijx
                    ay(i) = ay(i) + fmod * rijy
                    az(i) = az(i) + fmod * rijz
          
                    ax(j) = ax(j) - fmod * rijx
                    ay(j) = ay(j) - fmod * rijy
                    az(j) = az(j) - fmod * rijz
                endif
            end do
        end do
        

!   AÑADIMOS LOS FACTORES GLOBALES PARA OBTENER LOS RESULTADOS
      
    epot = 4.d00 * epot + corr_ener
    sum_rvp = 24.d00 * sum_rvp + corr_sum_rvp
    sum_r2vpp = 24.d00 * sum_r2vpp + corr_sum_r2vpp
    ax = 24.d00 * ax
    ay = 24.d00 * ay
    az = 24.d00 * az


    dfiv = sum_rvp/(3.d00*vol)
    d2fiv = (sum_r2vpp - 2.d00 * sum_rvp)/(9.d00*vol*vol)

    return

    end subroutine potlj