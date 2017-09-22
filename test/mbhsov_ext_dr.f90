
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! This program tests separation of variables as a means of 
! evaluating the far field due to a collection of modified
! biharmonic sources (quadrupole sources with random strengths
! and directions)
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


  program mbhrouts_dr

    implicit real *8 (a-h,o-z)
    real *8, allocatable, dimension(:) :: wsave 
    complex *16, allocatable, dimension(:) :: mbhmp,ymp
    complex *16, allocatable, dimension(:) :: mbhmp3,ymp3
    complex *16, allocatable, dimension(:) :: lmp2,ymp2
    complex *16, allocatable, dimension(:) :: lmp4,ymp4
    real *8, allocatable, dimension(:) :: utarg,utargex
    real *8, allocatable, dimension(:) :: utarg2, utarg3, utarg4
    real *8, allocatable, dimension(:) :: u, up
    real *8, allocatable, dimension(:,:) :: gradtarg,gradtargex
    real *8, allocatable, dimension(:,:) :: hesstarg,hesstargex
    real *8, allocatable, dimension(:,:) :: gradtargloc,hesstargloc
    real *8, allocatable, dimension(:,:) :: gradtarg2,hesstarg2
    real *8, allocatable, dimension(:,:) :: gradtarg3,hesstarg3
    real *8, allocatable, dimension(:,:) :: gradtarg4,hesstarg4
    real *8, allocatable, dimension(:,:) :: gradplot,hessplot
    real *8, allocatable, dimension(:) :: potplot, dmask
    complex *16, allocatable, dimension(:,:) :: gradtargqp,hesstargqp
    complex *16, allocatable, dimension(:,:) :: gradtargqplap, & 
         hesstargqplap
    real *8, allocatable, dimension(:,:) :: zptstarg, zptssrc, zcirc
    real *8, allocatable, dimension(:,:) :: znrmtarg, znrmsrc
    real *8, allocatable, dimension(:,:) :: quadvec, dir1, dir2, dipvec
    real *8, allocatable, dimension(:) :: quadstr, dipstr, charge
    complex *16, allocatable, dimension(:) :: zdipstr, zcharge
    real *8, allocatable, dimension(:,:) :: target, source, targplot
    real *8, allocatable, dimension(:,:,:) :: amatsmbh, amatsnaive
    real *8, allocatable, dimension(:) :: condsmbh, condsnaive
    real *8 :: amattemp(2,2), utemp(2,2), vtemp(2,2), wtemp(2)
    real *8 :: ztarg(2),center(2),pot,grad(2),hess(3),lambda,zdiff(2), zarb(2)

    complex *16, allocatable, dimension(:) :: zquadstr, &
         potqp, potqplap
    real *8, allocatable :: potqp1(:), gradqp1(:,:), hessqp1(:,:), potqp2(:)
    complex *16 :: eye, z, ztemp1, ztemp2

    complex *16 gradtemp, hesstemp

    real *8 :: work(10000)

    integer nterms, npts, ifpot, ifgrad, ifhess
    parameter (iseed = 281+3308004)
    data eye /(0.0d0,1.0d0)/

    logical ifrel, matu, matv

    ! set parameters

    npts = 100
    nterms = (npts+1)/2
    nwsave = 2*npts+15
    l = (npts+1)/2

    ! scale

    ! target domain = [-scale,scale] x [-scale,scale] box
    ! outside of circle of radius 0.5 scale

    ! sources located within circle of radius frac*0.5*scale

    ! nterms needed for machine precision depends on frac

    scale = 1.0d0
    frac = 0.5d0
    ns = 100
    nt = 100

    ngrid = 300
    ntplot = ngrid*ngrid

    ! open files for writing

    open(FILE='ext_targ.txt',UNIT=22,STATUS='REPLACE')
    open(FILE='ext_src.txt',UNIT=23,STATUS='REPLACE')
    open(FILE='ext_example.txt',UNIT=24,STATUS='REPLACE')
    open(FILE='ext_funs.txt',UNIT=30,STATUS='REPLACE')        
    open(FILE='ext_conds.txt',UNIT=28,STATUS='REPLACE')            
    open(FILE='ext_pot.txt',UNIT=25,STATUS='REPLACE')
    open(FILE='ext_grad.txt',UNIT=26,STATUS='REPLACE')
    open(FILE='ext_hess.txt',UNIT=27,STATUS='REPLACE')    
    
    ! allocate storage

    allocate(amatsmbh(2,2,0:nterms),amatsnaive(2,2,0:nterms), &
         condsmbh(0:nterms),condsnaive(0:nterms))
    allocate(u(npts),up(npts))
    allocate(wsave(nwsave))
    allocate(target(2,nt),source(2,ns),zptssrc(2,npts),zptstarg(2,npts))
    allocate(targplot(2,ntplot),potplot(ntplot),gradplot(2,ntplot), &
         hessplot(3,ntplot),dmask(ntplot))
    allocate(znrmtarg(2,npts),znrmsrc(2,npts))
    allocate(zquadstr(ns),quadstr(ns),quadvec(3,ns),dir1(2,ns),dir2(2,ns))
    allocate(dipstr(ns),dipvec(2,ns),charge(ns),zcharge(ns))
    allocate(zdipstr(ns))
    allocate(utargex(nt),gradtargex(2,nt),utarg(nt),gradtarg(2,nt))
    allocate(hesstarg(3,nt),hesstargex(3,nt))
    allocate(gradtargloc(2,nt),hesstargloc(3,nt))
    allocate(mbhmp(0:npts),ymp(0:npts))
    allocate(mbhmp3(0:npts),ymp3(0:npts))    
    allocate(ymp2(0:npts),lmp2(0:npts))
    allocate(ymp4(0:npts),lmp4(0:npts))
    allocate(gradtarg2(2,nt),hesstarg2(3,nt))
    allocate(gradtarg3(2,nt),hesstarg3(3,nt))
    allocate(gradtarg4(2,nt),hesstarg4(3,nt))
    allocate(gradtargqp(2,nt),hesstargqp(3,nt))
    allocate(utarg4(nt),utarg3(nt))
    allocate(gradtargqplap(2,nt),hesstargqplap(3,nt))
    allocate(potqp(nt),potqplap(nt),utarg2(nt))
    allocate(zcirc(2,npts))
    allocate(potqp1(nt),gradqp1(2,nt),hessqp1(3,nt),potqp2(nt))
    call prini(6,13)

    temp = hkrand(iseed)
    pi = 4.0d0*atan(1.0d0)
    call dffti(npts,wsave)

    ! type of sources 

    ifcharge = 0
    ifdipole = 0
    ifquad = 1
    ifoct = 0

    ! arbitrary offset, trying to avoid possible subnormal 
    ! shenanigans
    
    zarb(1) = hkrand(0)
    zarb(2) = hkrand(0)

    zarb = (/ 0.0d0, 0.0d0 /)

    ! set up some targets in the box [-scale,scale]^2 
    ! minus the circle of radius 0.5*scale

    ztarg(1) = 0.0d0
    ztarg(2) = 0.0d0

    do i = 1,nt
       do ii = 1,100
          target(1,i) = ztarg(1) - scale + 2*scale*hkrand(0)
          target(2,i) = ztarg(2) - scale + 2*scale*hkrand(0)
          if (target(1,i)**2 + target(2,i)**2 .ge. (0.5d0*scale)**2) then
             exit
          endif
       enddo
    enddo

    do i = 1,nt
       target(1:2,i) = target(1:2,i) + zarb
    enddo

    do i = 1,nt
       write(22,*) target(1,i),   target(2,i)
    enddo

    ! points for plotting function

    h = 2.0d0*scale/ngrid

    do i = 1,ngrid
       do j = 1,ngrid
          ii = (i-1)*ngrid + j
          xx = ztarg(1) -scale - 0.5d0*h + i*h
          yy = ztarg(2) -scale - 0.5d0*h + j*h
          targplot(1,ii) = xx
          targplot(2,ii) = yy

          if (xx**2 + yy**2 .le. (0.5d0*scale)**2) then
             dmask(ii) = 1.0d200
          else
             dmask(ii) = 1.0d0
          endif

       enddo
    enddo

    do i = 1,ntplot
       targplot(1:2,i) = targplot(1:2,i) + zarb
    enddo

    ztarg = ztarg+zarb

    ! set up some quadrupole sources in the disc
    ! of radius 0.5d0*frac*scale

    center(1) = 0.0d0
    center(2) = 0.0d0

    do i = 1,ns
       ! locations
       th = 2.0d0*pi*hkrand(0)
       rr = 0.5d0*frac*scale*hkrand(0)
       source(1,i) = center(1) + dcos(th)*rr
       source(2,i) = center(2) + dsin(th)*rr
       ! charge
       charge(i) = -1.0d0+2.0d0*hkrand(0)
       zcharge(i) = charge(i)
       
       ! dipole
       ! direction

       a = -0.5d0+hkrand(0)
       b = -0.5d0+hkrand(0)
       dab = dsqrt(a**2+b**2)
       a = a/dab
       b = b/dab
       dipvec(1,i) = a
       dipvec(2,i) = b

       dipstr(i) = hkrand(0)

       zdipstr(i) = dipstr(i)*(dipvec(1,i)+eye*dipvec(2,i))       

       ! directions
       a = -0.5d0+hkrand(0)
       b = -0.5d0+hkrand(0)
       c = -0.5d0+hkrand(0)
       d = -0.5d0+hkrand(0)
       dab = dsqrt(a**2+b**2)
       dcd = dsqrt(c**2+d**2)
       a = a/dab
       b = b/dab
       c = c/dcd
       d = d/dcd
       dir1(1,i) = a
       dir1(2,i) = b
       dir2(1,i) = c
       dir2(2,i) = d
       ! convert to xx, xy, yy
       quadvec(1,i) = a*c
       quadvec(2,i) = a*d + b*c
       quadvec(3,i) = b*d
       ! strengths
       quadstr(i) = hkrand(0)
       zquadstr(i) = quadstr(i)
    enddo

    center = center + zarb
    do i = 1,ns
       source(1:2,i) = source(1:2,i) + zarb
    enddo

    do i = 1,ns
       write(23,*) source(1,i),   source(2,i),   dipvec(1,i), &
              dipvec(2,i),   dir1(1,i), &
              dir1(2,i),   dir2(1,i),   dir2(2,i)
    enddo

    lambda = 2.0d0**(-24)

    ifpot = 1
    ifgrad = 1
    ifhess = 1
    
    do i = 1,ntplot
       potplot(i) = 0.0d0
       gradplot(1:2,i) = (/ 0.0d0, 0.0d0 /)
       hessplot(1:3,i) = (/ 0.0d0, 0.0d0, 0.0d0 /)
       do j = 1,ns
          if (ifcharge .eq. 1) then 
             call modbhgreen(lambda,targplot(1,i),source(1,j), &
                  ifpot,pot,ifgrad,grad,ifhess,hess)
             potplot(i) = potplot(i) + charge(j)*pot
             gradplot(1:2,i) = gradplot(1:2,i) + charge(j)*grad
             hessplot(1:3,i) = hessplot(1:3,i) + charge(j)*hess
          endif
          if (ifdipole .eq. 1) then 
             call modbhgreend1(lambda,targplot(1,i),source(1,j), &
                  ifpot,pot,ifgrad,grad,ifhess,hess,dipvec(1,j))
             potplot(i) = potplot(i) + dipstr(j)*pot
             gradplot(1:2,i) = gradplot(1:2,i) + dipstr(j)*grad
             hessplot(1:3,i) = hessplot(1:3,i) + dipstr(j)*hess
          endif
          if (ifquad .eq. 1) then
             call modbhgreend2(lambda,targplot(1,i),source(1,j), &
                  ifpot,pot,ifgrad,grad,ifhess,hess,dir1(1,j), &
                  dir2(1,j))
             potplot(i) = potplot(i) + quadstr(j)*pot
             gradplot(1:2,i) = gradplot(1:2,i) + quadstr(j)*grad
             hessplot(1:3,i) = hessplot(1:3,i) + quadstr(j)*hess
          endif
       enddo
    enddo


    do i = 1,ntplot
       write(24,*) potplot(i),   gradplot(1,i),   gradplot(2,i), &
              hessplot(1,i),   hessplot(2,i),   hessplot(3,i), &
              dmask(i),   lambda,  targplot(1,i)-zarb(1), &
             targplot(2,i)-zarb(2)
    enddo


    do iiii = -24,8
  
       do jjjj = 1,10

          lambda = 2.0d0**(iiii) + hkrand(0)*2.0d0**(iiii)

          rscale = min(scale*0.5d0*lambda,1.0d0)
          rscalelap = min(scale*0.5d0,1.0d0)


          rad = 0.5d0*scale

          ! getting exact values

          ifpot = 1
          ifgrad = 1
          ifhess = 1

          do i = 1,nt
             utargex(i) = 0.0d0
             gradtargex(1:2,i) = (/ 0.0d0, 0.0d0 /)
             hesstargex(1:3,i) = (/ 0.0d0, 0.0d0, 0.0d0 /)
             do j = 1,ns
                if (ifcharge .eq. 1) then 
                   call modbhgreen(lambda,target(1,i),source(1,j), &
                        ifpot,pot,ifgrad,grad,ifhess,hess)
                   utargex(i) = utargex(i) + charge(j)*pot
                   gradtargex(1:2,i) = gradtargex(1:2,i) + charge(j)*grad
                   hesstargex(1:3,i) = hesstargex(1:3,i) + charge(j)*hess
                endif
                if (ifdipole .eq. 1) then 
                   call modbhgreend1(lambda,target(1,i),source(1,j), &
                        ifpot,pot,ifgrad,grad,ifhess,hess,dipvec(1,j))
                   utargex(i) = utargex(i) + dipstr(j)*pot
                   gradtargex(1:2,i) = gradtargex(1:2,i) + dipstr(j)*grad
                   hesstargex(1:3,i) = hesstargex(1:3,i) + dipstr(j)*hess
                endif
                if (ifquad .eq. 1) then
                   call modbhgreend2(lambda,target(1,i),source(1,j), &
                        ifpot,pot,ifgrad,grad,ifhess,hess,dir1(1,j), &
                        dir2(1,j))
                   utargex(i) = utargex(i) + quadstr(j)*pot
                   gradtargex(1:2,i) = gradtargex(1:2,i) + quadstr(j)*grad
                   hesstargex(1:3,i) = hesstargex(1:3,i) + quadstr(j)*hess
                endif
             enddo
          enddo

          ! compute multipole coefficients using new and old functions

          call mbh2dsov_circpts(zcirc,npts)

          call mbh2dsov_circvals(u,up,source,ifcharge,charge, &
               ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec, &
               ifoct,octstr,octvec,ns,center,rad,lambda,npts,zcirc)

          do i = 1,npts
             write(30,*) u(i),   up(i)
          enddo
          
          call mbh2dsov_dtocmp(lambda,rad,rscale,u,up,npts,mbhmp, &
               ymp,wsave,work)

          call mbh2dformmp_all(ier,lambda,rscale,source,ifcharge,charge, &
               ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec,ifoct,octstr, &
               octvec,ns,center,nterms,mbhmp3,ymp3)

          call mbh2dsov_matsmp(lambda,rad,rscale,npts,amatsmbh,work)

          call mbh2dsov_dtocmp_naive(lambda,rad,rscale,rscalelap,u,up,npts,lmp2, &
               ymp2,wsave,work)

          call mbh2dsov_matsmp_naive(lambda,rad,rscale,rscalelap, &
               npts,amatsnaive,work)

          do i = 0,nterms
             ymp4(i) = 0.0d0
             lmp4(i) = 0.0d0
          enddo

          if (ifcharge .eq. 1) then
             call y2dformmp_add(ier,lambda,rscale,source,charge,ns, &
                  center,nterms,ymp4)
             call l2dformmp_add(ier,rscalelap,source,zcharge,ns, &
                  center,nterms,lmp4)
          endif
          if (ifdipole .eq. 1) then
             call y2dformmp_dp_add(ier,lambda,rscale,source,dipstr,dipvec, &
                  ns,center,nterms,ymp4)
             call l2dformmp_dp_add(ier,rscalelap,source,zdipstr,ns, &
                  center,nterms,lmp4)
          endif
          if (ifquad .eq. 1) then
             call y2dformmp_qp_add(ier,lambda,rscale,source,zquadstr,quadvec, &
                  ns,center,nterms,ymp4)
             call l2dformmp_qp_add(ier,rscalelap,source,zquadstr,quadvec,ns, &
                  center,nterms,lmp4)
          endif


!          do i = 0,nterms-1
!             write(*,*) i
!             write(*,*) ymp4(i)/(lambda**2)/ymp2(i), lmp4(i)/(lambda**2)/lmp2(i)
!          enddo
          
          ! obtain condition numbers of scaled matrices 

          do i = 0,l-1
             amattemp(1:2,1:2) = amatsmbh(1:2,1:2,i)
             na = 2
             call dnormcols(amattemp,na,na,na)

             matu= .false.
             matv= .false.
             call svd(na,na,na,amattemp,wtemp,matu,utemp, &
                  matv,vtemp,ierr,work)
             if (ierr .ne. 0) then
                write(*,*) 'ERROR IN SVD ROUTINE', ierr
                write(*,*) 'ABORT'
                stop
             endif


             
             condsmbh(i) = max(wtemp(1)/wtemp(2),wtemp(2)/wtemp(1))
             
             amattemp(1:2,1:2) = amatsnaive(1:2,1:2,i)

             na = 2
             call dnormcols(amattemp,na,na,na)
             matu= .false.
             matv= .false.
             call svd(na,na,na,amattemp,wtemp,matu,utemp, &
                  matv,vtemp,ierr,work)
             if (ierr .ne. 0) then
                write(*,*) 'ERROR IN SVD ROUTINE', ierr
                write(*,*) 'ABORT'
                stop
             endif
             
             condsnaive(i) = max(wtemp(1)/wtemp(2),wtemp(2)/wtemp(1))


             
          enddo

          do i = 0,l-1
             write(28,*) l,  lambda,  condsmbh(i),  condsnaive(i), scale
          enddo

          ! evaluate field at targets using new functions

          call mbh2dmpevalall(lambda,rscale,center,mbhmp,ymp, &
               nterms,target,nt,ifpot,utarg,ifgrad,gradtarg, &
               ifhess,hesstarg)

          call mbh2dmpevalall(lambda,rscale,center,mbhmp3,ymp3, &
               nterms,target,nt,ifpot,utarg3,ifgrad,gradtarg3, &
               ifhess,hesstarg3)

          ! evaluate field at targets using the difference of naive functions

          ifpot = 1
          ifgrad = 1
          ifhess = 1

          do i = 1,nt
             potqp1(i) = 0.0d0
             gradqp1(1:2,i) = (/ 0.0d0, 0.0d0 /)
             hessqp1(1:3,i) = (/ 0.0d0, 0.0d0, 0.0d0 /)
          enddo

          call y2dmpevalall(lambda,rscale,center,ymp2,nterms,target,nt, &
               ifpot,potqp1,ifgrad,gradqp1,ifhess,hessqp1)

          do i = 1,nt
             potqp(i) = potqp1(i)
             gradtargqp(1:2,i) = gradqp1(1:2,i)
             hesstargqp(1:3,i) = hessqp1(1:3,i)
          enddo

          do i = 1,nt
             potqplap(i) = 0.0d0
             gradtargqplap(1:2,i) = (/ 0.0d0, 0.0d0 /)
             hesstargqplap(1:3,i) = (/ 0.0d0, 0.0d0, 0.0d0 /)
          enddo

          call l2dmpevalall(rscalelap,center,lmp2,nterms,target,nt, &
               ifpot,potqplap,ifgrad,gradtargqplap,ifhess,hesstargqplap)

          do i = 1,nt
             utarg2(i) = dreal(potqp(i)-potqplap(i)/(2.0d0*pi))
             gradtarg2(1,i) = dreal(gradtargqp(1,i) &
                  -gradtargqplap(1,i)/(2.0d0*pi))
             gradtarg2(2,i) = dreal(gradtargqp(2,i) &
                  -gradtargqplap(2,i)/(2.0d0*pi))
             hesstarg2(1,i) = dreal(hesstargqp(1,i) &
                  -hesstargqplap(1,i)/(2.0d0*pi))
             hesstarg2(2,i) = dreal(hesstargqp(2,i) &
                  -hesstargqplap(2,i)/(2.0d0*pi))
             hesstarg2(3,i) = dreal(hesstargqp(3,i) &
                  -hesstargqplap(3,i)/(2.0d0*pi))
          enddo

          ! evaluate field at targets using the difference of naive functions
          ! "exact" multipole coefficient values

          ifpot = 1
          ifgrad = 1
          ifhess = 1

          do i = 1,nt
             potqp1(i) = 0.0d0
             gradqp1(1:2,i) = (/ 0.0d0, 0.0d0 /)
             hessqp1(1:3,i) = (/ 0.0d0, 0.0d0, 0.0d0 /)
          enddo

          call y2dmpevalall(lambda,rscale,center,ymp4,nterms,target,nt, &
               ifpot,potqp1,ifgrad,gradqp1,ifhess,hessqp1)

          do i = 1,nt
             potqp(i) = potqp1(i)
             gradtargqp(1:2,i) = gradqp1(1:2,i)
             hesstargqp(1:3,i) = hessqp1(1:3,i)
          enddo

          do i = 1,nt
             potqplap(i) = 0.0d0
             gradtargqplap(1:2,i) = (/ 0.0d0, 0.0d0 /)
             hesstargqplap(1:3,i) = (/ 0.0d0, 0.0d0, 0.0d0 /)
          enddo

          call l2dmpevalall(rscalelap,center,lmp4,nterms,target,nt, &
               ifpot,potqplap,ifgrad,gradtargqplap,ifhess,hesstargqplap)

          do i = 1,nt
             utarg4(i) = dreal(potqp(i)-potqplap(i)/(2.0d0*pi))/(lambda**2)
             gradtarg4(1,i) = dreal(gradtargqp(1,i) &
                  -gradtargqplap(1,i)/(2.0d0*pi))/(lambda**2)
             gradtarg4(2,i) = dreal(gradtargqp(2,i) &
                  -gradtargqplap(2,i)/(2.0d0*pi))/(lambda**2)
             hesstarg4(1,i) = dreal(hesstargqp(1,i) &
                  -hesstargqplap(1,i)/(2.0d0*pi))/(lambda**2)
             hesstarg4(2,i) = dreal(hesstargqp(2,i) &
                  -hesstargqplap(2,i)/(2.0d0*pi))/(lambda**2)
             hesstarg4(3,i) = dreal(hesstargqp(3,i) &
                  -hesstargqplap(3,i)/(2.0d0*pi))/(lambda**2)
          enddo

          ifprint = 1
          
          if (ifprint .eq. 1) then
             
             ifrel = .true.

             write(*,*) 'LAMBDA    '
             write(*,*) lambda
             write(*,*) 'SIZE OF SOURCE BOX '
             write(*,*) scale
             write(*,*) 'ERROR FOR POTENTIAL WITH:' 
             write(*,*) '(1) NEW OUTGOING BASIS FUNCTIONS '
             write(*,*) '(2) DIFFERENCE OF YUKAWA AND LAPLACE '
             write(*,*) '(3) EXACT COEFFS NEW '
             write(*,*) '(4) EXACT COEFFS DIFF '

             err_newmp = rmsfun1(utarg,utargex,nt,ifrel)
             err_oldmp = rmsfun1(utarg2,utargex,nt,ifrel)
             err3 = rmsfun1(utarg3,utargex,nt,ifrel)
             err4 = rmsfun1(utarg4,utargex,nt,ifrel)
             write(*,*) err_newmp, err_oldmp, err3, err4

             write(25,*) lambda,   err_newmp,   err_oldmp, &
                    err3,   err4, scale

             write(*,*) 'ERROR FOR GRAD WITH:' 
             write(*,*) '(1) NEW OUTGOING BASIS FUNCTIONS '
             write(*,*) '(2) DIFFERENCE OF YUKAWA AND LAPLACE '
             write(*,*) '(3) EXACT COEFFS NEW '
             write(*,*) '(4) EXACT COEFFS DIFF '

             nt2 = 2*nt
             err_newmp = rmsfun1(gradtarg,gradtargex,nt2,ifrel)
             err_oldmp = rmsfun1(gradtarg2,gradtargex,nt2,ifrel)
             err3 = rmsfun1(gradtarg3,gradtargex,nt2,ifrel)
             err4 = rmsfun1(gradtarg4,gradtargex,nt2,ifrel)
             write(*,*) err_newmp, err_oldmp, err3, err4

             write(26,*) lambda,   err_newmp,   err_oldmp, &
                    err3,   err4, scale

             write(*,*) 'ERROR FOR HESS WITH:' 
             write(*,*) '(1) NEW OUTGOING BASIS FUNCTIONS '
             write(*,*) '(2) DIFFERENCE OF YUKAWA AND LAPLACE '
             write(*,*) '(3) EXACT COEFFS NEW '
             write(*,*) '(4) EXACT COEFFS DIFF '

             nt3 = 3*nt
             err_newmp = rmsfun1(hesstarg,hesstargex,nt3,ifrel)
             err_oldmp = rmsfun1(hesstarg2,hesstargex,nt3,ifrel)
             err3 = rmsfun1(hesstarg3,hesstargex,nt3,ifrel)
             err4 = rmsfun1(hesstarg4,hesstargex,nt3,ifrel)
             write(*,*) err_newmp, err_oldmp, err3, err4

             write(27,*) lambda,   err_newmp,   err_oldmp, &
                    err3,   err4, scale

          endif

       enddo

    enddo


    stop
  end program mbhrouts_dr


  real *8 function rmsfun1(u,uexact,n,ifrel)

    implicit none
    real *8 u(*), uexact(*)
    integer n
    logical ifrel

    real *8 err, tot
    integer i

    err = 0.0d0
    tot = 0.0d0

    do i = 1,n
       err = err + (u(i)-uexact(i))**2
       tot = tot + (uexact(i))**2
    enddo

    if (ifrel) then
       rmsfun1 = dsqrt(err/tot)
    else
       rmsfun1 = dsqrt(err/n)
    endif
    
  end function rmsfun1

  subroutine dnormcols(amat,ma,m,n)

    implicit none
    real *8 amat(ma,n)
    integer ma, m, n
    
    integer i, j
    real *8 dtemp

    do j = 1,n
       dtemp = 0.0d0
       do i = 1,m
          dtemp = dtemp + amat(i,j)**2
       enddo
       dtemp = dsqrt(dtemp)
       if (abs(dtemp).gt. 0.0d0) then
          do i = 1,m
             amat(i,j) = amat(i,j)/dtemp
          enddo
       end if
    enddo

  end subroutine dnormcols
