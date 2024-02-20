!
!   SIMULACION NVE DE UN SISTEMA DE 500 PARTICULAS QUE INTERACTUAN ENTRE PARES MEDIANTE POTENCIAL
!   L-J. EL PROGRAMA REQUIERE DE UN SISTEMA INICIAL DE 500 PARTICULAS EN UNA CAJA DE LADO L Y
!   CONSIDERANDO UN RADIO DE CORTE DE L/2, CALCULA ENERGÍAS POTENCIAL Y CINÉTICA FIJANDO UN VALOR
!   PARA LA ENERGÍA TOTAL. CALCULA ADEMAS LAS DERIVADAS DEL POTENCIAL
!   
!   EL PROGRAMA HACE EVOLUCIONAR EN EL TIEMPO EL SISTEMA A PASOS DE TIEMPO DELTAT=0.0001 A LO LARGO
!   DE 500.000 PASOS MEDIANTE EL ALGORITMO DE VERLET.
!   
!   FINALMENTE, OBTIENE VALORES MEDIOS DE LAS ENERGÍAS, DE LA INVERSA DE LA ENERGÍA CINÉTICA,
!   DE LAS DERIVADAS DEL POTENCIAL Y DEL PRODUCTO DE LA INVERSA DE LA ENERGÍA CINÉTICA CON LA
!   DERIVADA DEL POTENCIAL PARA PODER DESCRIBIR LA DINÁMICA DEL SISTEMA
!   
program dinamica
        


        use def_precision
        implicit none
        
        integer (kind = entero), parameter :: npmax = 500
        integer (kind = entero) :: np,kk,ktotal,kpaso,kcuenta,pl,i,simulacion
        real (kind=doblep), dimension(npmax) :: rx,ry,rz,vx,vy,vz,ax,ay,az
        real (kind=doblep) :: xnp,epot,ecin,etot,xpasos,pli,rc,rc2,vol,dt
        real (kind=doblep) :: grad_lib,temp,aux1,aux2,capcv,cespv
        real (kind=doblep) :: etot_media,epot_media,ecin_media,alfa_ef,alfa_e,compr_s,gama
        real (kind=doblep) :: presion,fiv2_ecin_med,fiv_ecin_med,ecin_inv,d2fiv_med,dfiv_med,ecin_inv_media
        real (kind=doblep) :: d2fiv,dfiv,dens
        real (kind=doblep), dimension(500000) :: lista_et,lista_ep,lista_ec,lista_ecinv,lista_dfiv,lista_d2fiv
        real (kind=doblep) :: error_et,error_ep,error_ec,error_ecinv,error_dfiv,error_d2fiv,error_fiv_ecin
        real (kind=doblep) :: error_fiv2_ecin
        real (kind=doblep) :: error_t,error_p,error_cv,error_ae,error_aef,error_g,error_k
        real (kind=doblep) :: aux_cv,aux_ae,aux_k,aux_aef
        real (kind=doblep), dimension(10) :: lista_t,lista_p,lista_cv,lista_k,lista_g,lista_ae,lista_aef
        real (kind=doblep) :: et,ep,ec,ek,eg,eae,eaef
        real (kind=doblep) :: mt,mp,mc,mk,mg,mae,maef
        character(len=50) :: fname,gname1,gname2,gname3,gname4,gname5


!LEER DATOS IGUAL QUE CON EQUILIBRACION

!ktotal ahora es 500000
!leer 2 veces gname1 porque en una de ellas sobreescribes, o igual es porque 

9000 format(a25)
    write(*,*) 'Archivo de datos simulacion (parametros.txt)??'
    read(*,9000) fname



    do simulacion = 1, 10
        
        open(10,file=fname)
            read(10,*) np
            read(10,*) pl
            read(10,*) pli
            read(10,*) rc
            read(10,*) rc2
            read(10,*) vol
            read(10,*) dens
            read(10,*) ktotal
            read(10,*) kpaso
            read(10,*) dt!(500000 pasos, grabamos cada 100, y deltat de 0.0001)
            read(10,*) etot
            read(10,*) ecin
            read(10,*) epot
            read(10,*) dfiv
            read(10,*) d2fiv
            read(10,9000) gname1
            read(10,9000) gname2
            read(10,9000) gname3
            read(10,9000) gname4
            read(10,9000) gname5
        close(10)

        xnp=dble(np)
        xpasos=dble(ktotal)
        grad_lib=3.d00*xnp-3.d00
        aux1=1.d00-2.d00/grad_lib
        aux2=grad_lib/2.d00-1.d00

        open(15,file=gname1,form='unformatted')
            read(15) rx,ry,rz,vx,vy,vz,ax,ay,az
        close(15)

        
        
        ecin_media=0.d00
        ecin_inv_media=0.d00
        epot_media=0.d00
        etot_media=0.d00
        dfiv_med=0.d00
        d2fiv_med=0.d00
        fiv_ecin_med = 0.d00
        fiv2_ecin_med = 0.d00

        kcuenta=0


        do kk=1,ktotal
            call verlet(np,dt,rx,ry,rz,vx,vy,vz,ax,ay,az,epot,ecin,etot,dfiv,d2fiv)
            ecin_inv=1.d00/ecin

            ecin_media=ecin_media+ecin
            ecin_inv_media=ecin_inv_media+ecin_inv
            epot_media=epot_media+epot
            etot_media=etot_media+etot
            dfiv_med=dfiv_med+dfiv
            d2fiv_med=d2fiv_med + d2fiv
            fiv_ecin_med=fiv_ecin_med+dfiv*ecin_inv
            fiv2_ecin_med=fiv2_ecin_med+dfiv*dfiv*ecin_inv

            lista_et(kk) = etot
            lista_ep(kk) = epot
            lista_ec(kk) = ecin
            lista_ecinv(kk) = ecin_inv
            lista_dfiv(kk) = dfiv
            lista_d2fiv(kk) = d2fiv


                if (mod(kk,kpaso)==0) then
                    kcuenta = kcuenta + 1
                    write(*,*) kcuenta*kpaso,'/5.000.000 realizado. Simulacion: ',simulacion,'. Et= ',etot
                end if          
        end do
            



        write(*,*) 'realizados ',ktotal,' pasos de la simulacion ',simulacion

        !Calculo los valores medios de esta simulacion
        ecin_media = ecin_media/xpasos
        ecin_inv_media=ecin_inv_media/xpasos
        epot_media=epot_media/xpasos
        etot_media=etot_media/xpasos
        dfiv_med=dfiv_med/xpasos
        d2fiv_med=d2fiv_med/xpasos
        fiv_ecin_med=fiv_ecin_med/xpasos
        fiv2_ecin_med=fiv2_ecin_med/xpasos
        
        error_ec = 0.d00
        error_ecinv = 0.d00
        error_ep = 0.d00
        error_et = 0.d00
        error_dfiv = 0.d00
        error_d2fiv = 0.d00
        error_fiv_ecin = 0.d00
        error_fiv2_ecin = 0.d00
        
        !Calculo los errores
        do i=1,ktotal
            error_ec = error_ec + (ecin_media - lista_ec(i))**2
            error_ecinv = error_ecinv + (ecin_inv_media - lista_ecinv(i))**2
            error_ep = error_ep + (epot_media - lista_ep(i))**2
            error_et = error_et + (etot_media - lista_et(i))**2
            error_dfiv = error_dfiv + (dfiv_med - lista_dfiv(i))**2
            error_d2fiv = error_d2fiv + (d2fiv_med - lista_d2fiv(i))**2
            error_fiv_ecin = error_fiv_ecin + (fiv_ecin_med - lista_dfiv(i)*lista_ecinv(i))**2
            error_fiv2_ecin = error_fiv2_ecin + (fiv2_ecin_med - lista_dfiv(i)*lista_dfiv(i)*lista_ecinv(i))**2
        end do

        error_ec = dsqrt(error_ec/dble(ktotal))
        error_ecinv = dsqrt(error_ecinv/dble(ktotal))
        error_ep = dsqrt(error_ep/dble(ktotal))
        error_et = dsqrt(error_et/dble(ktotal))
        error_dfiv = dsqrt(error_dfiv/dble(ktotal))
        error_d2fiv = dsqrt(error_d2fiv/dble(ktotal))
        error_fiv_ecin = dsqrt(error_fiv_ecin/dble(ktotal))
        error_fiv2_ecin = dsqrt(error_fiv2_ecin/dble(ktotal))


        !Escribimos sobre el inicial, sera nuestro punto de partida para la siguiente simulacion
        open(70,file=gname1,form='unformatted')
            write(70) rx,ry,rz,vx,vy,vz,ax,ay,az
        close(70)

        ! calculamos magnitudes con error mediante propagaciones
        temp = 2.d00*ecin_media/grad_lib
        error_t = 2.d00*error_ec/grad_lib

        presion = xnp*temp/vol-dfiv_med
        error_p = dsqrt((xnp*error_t/vol)**2 + error_dfiv*error_dfiv)

        aux_cv = (1.d00-aux1*ecin_media*ecin_inv_media)
        error_cv = dsqrt((error_ec*ecin_inv_media)**2 + (error_ecinv*ecin_media)**2)*abs(1.d00 - aux1)
        capcv = 1.d00/aux_cv
        error_cv = capcv**2 * error_cv
        cespv = capcv/xnp

        gama = 1.d00/cespv+vol*aux2*(dfiv_med*ecin_inv_media-fiv_ecin_med)
        error_g = dsqrt(error_cv*error_cv*(xnp/(capcv*capcv))**2 &
                  + vol*vol*aux2*aux2*(ecin_inv_media*ecin_inv_media*error_dfiv*error_dfiv &
                  + dfiv_med*dfiv_med*error_ecinv*error_ecinv + error_fiv_ecin*error_fiv_ecin))

        aux_k =xnp*temp*(1.d00 + 2.d00*gama-xnp/capcv)/vol + vol*d2fiv_med
        aux_k = aux_k - vol*aux2*(fiv2_ecin_med-2.d00*dfiv_med*fiv_ecin_med+dfiv_med*dfiv_med*ecin_inv_media)
        compr_s = 1.d00/aux_k
        error_k = dsqrt((xnp/vol)**2 * (error_t**2 * (1.d00+2.d00*gama-xnp/capcv)**2 &
                  + temp**2 * (4.d00*error_g**2 + (xnp/capcv/capcv *error_cv)**2)) &
                  + (vol * error_d2fiv)**2 + (vol*aux2)**2 * (error_fiv2_ecin**2 &
                  + 4.d00*(error_dfiv*fiv_ecin_med)**2+4.d00*(error_fiv_ecin*dfiv_med)**2 &
                  + (error_dfiv*2.d00*dfiv_med*ecin_inv_media)**2 + (error_ecinv*dfiv_med**2)**2))
        error_k = 1.d00/(aux_k**2) * error_k


        aux_ae = presion*vol/capcv-gama*temp
        alfa_e = 1.d00/aux_ae
        error_ae = dsqrt((vol*error_p/capcv)**2 + (presion*vol*error_cv/capcv/capcv)**2 + (temp*error_g)**2 + (error_t*gama)**2)
        error_ae = 1.d00/(aux_ae**2) * error_ae
        aux_aef= vol*(aux1*ecin_media*fiv_ecin_med-dfiv_med)
        alfa_ef = 1.d00/aux_aef
        error_aef = vol*dsqrt((aux1*ecin_media*error_fiv_ecin)**2 + (aux1*fiv_ecin_med*error_ec)**2 + error_dfiv**2)
        error_aef = 1.d00/(aux_aef**2) * error_aef
        !!Estas dos tienen que dar lo mismo. Si no, es que alguna de estas, o las anteriores, está mal

        !GRABAMOS LOS RESULTADOS

        open(50,file=gname5,access='append')
            write(50,*)'**** RESUMEN DE LA SIMULACION ',simulacion,' ****'
            write(50,*) 'etot_media=',etot_media,' error_etot_media=',error_et
            write(50,*) 'ecin_media=',ecin_media,' error_ecin_med =',error_ec
            WRITE(50,*) 'epot_media=',epot_media,' error_epot_media=',error_ep
            write(50,*) 'fiv_media=',dfiv_med,' error_dfiv_med=',error_dfiv
            write(50,*) 'd2fiv_media=',d2fiv_med,' error_d2fiv = ',error_d2fiv
            write(50,*) 'fiv_ecin_med=',fiv_ecin_med,' error_fiv_ecinv = ',error_fiv_ecin,
            write(50,*) 'fiv2_ecin_med=',fiv2_ecin_med,' error_fiv2_ecinv=',error_fiv2_ecin
            write(50,*) 'ecinv_med = ',ecin_inv_media,' error_ecinv_med =',error_ecinv
            write(50,*) 'T=',temp,'error_t = ',error_t
            write(50,*) 'P=',presion,' error_p = ',error_p
            write(50,*) 'cv =',capcv,' error_cv = ',error_cv
            write(50,*) 'compr_s=', compr_s,'error k =',error_k
            write(50,*) 'gama=',gama,' error_g=',error_g
            write(50,*) 'alfa_e=',alfa_e,' error_ae=',error_ae
            write(50,*) 'alfa_ef=',alfa_ef,'error_aef= ',error_aef
        close(50)
        
        !Añadimos los resultados a una lista
        lista_t(simulacion) = temp
        lista_p(simulacion) = presion
        lista_cv(simulacion) = capcv
        lista_k(simulacion) = compr_s
        lista_g(simulacion) = gama
        lista_ae(simulacion) = alfa_e
        lista_aef(simulacion) = alfa_ef
        
        !Grabamos la nueva configuración inicial
        open(60,file=fname)
            write(60,*) np
            write(60,*) pl
            write(60,*) pli
            write(60,*) rc
            write(60,*) rc2
            write(60,*) vol
            write(60,*) dens
            write(60,*) ktotal
            write(60,*) kpaso
            write(60,*) dt!
            write(60,*) etot
            write(60,*) ecin
            write(60,*) epot
            write(60,*) dfiv
            write(60,*) d2fiv
            write(60,9000) gname1
            write(60,9000) gname2
            write(60,9000) gname3
            write(60,9000) gname4
            write(60,9000) gname5
        close(60)

        !Siguiente simulacion        
    end do
    
    !Calculamos la media de los resultados
    mt = sum(lista_t)/10.d00
    mp = sum(lista_p)/10.d00
    mc = sum(lista_cv)/10.d00
    mk = sum(lista_k)/10.d00
    mg = sum(lista_g)/10.d00
    mae = sum(lista_ae)/10.d00
    maef = sum(lista_aef)/10.d00

    !Calculamos el error de la media
    et=0.d00
    ep=0.d00
    ec=0.d00
    ek=0.d00
    eg=0.d00
    eae=0.d00
    eaef=0.d00
    
    do i = 1,10
        et = et + (mt - lista_t(i))**2
        ep = ep + (mp - lista_p(i))**2
        ec = ec + (mc - lista_cv(i))**2
        ek = ek + (mk - lista_k(i))**2
        eg = eg + (mg - lista_g(i))**2
        eae = eae + (mae - lista_ae(i))**2
        eaef = eaef + (maef - lista_aef(i))**2
    end do

    et = dsqrt(et/10.d00)
    ep = dsqrt(ep/10.d00)
    ec = dsqrt(ec/10.d00)
    ek = dsqrt(ek/10.d00)
    eg = dsqrt(eg/10.d00)
    eae = dsqrt(eae/10.d00)
    eaef = dsqrt(eaef/10.d00)

    !Lo añadimos al archivo de resultados
    open(90,file=gname5,access='append')
        write(90,*) '**** MEDIA DE LAS 10 SIMULACIONES ****'
        write(90,*) 'Media T=',mt,' error T=',et
        write(90,*) 'Media P=',mp,' error P=',ep
        write(90,*) 'Media Cv=',mc,' error Cv=',ec
        write(90,*) 'Media k=',mk,' error k=',ek
        write(90,*) 'Media g=',mg,' error g=',eg
        write(90,*) 'Media ae=',mae,' error ae=',eae
        write(90,*) 'Media aef=',maef,' error aef=',eaef
    close(90)
        
end program