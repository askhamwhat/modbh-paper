
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
    real *8, allocatable, dimension(:) :: utarg,utargex
    complex *16, allocatable, dimension(:) :: mbhloc,lloc
    complex *16, allocatable, dimension(:) :: mbhloc2,lloc2
    real *8, allocatable, dimension(:) :: utargmp, utargloc
    real *8, allocatable, dimension(:) :: u, up
    real *8, allocatable, dimension(:) :: diffs,ders,kvec,kder
    real *8, allocatable, dimension(:,:) :: gradtarg,gradtargex
    real *8, allocatable, dimension(:,:) :: hesstarg,hesstargex
    real *8, allocatable, dimension(:,:) :: gradtargloc,hesstargloc
    real *8, allocatable, dimension(:,:) :: gradtargmp,hesstargmp
    complex *16, allocatable, dimension(:,:) :: gradtargqp,hesstargqp
    complex *16, allocatable, dimension(:,:) :: gradtargqplap, & 
         hesstargqplap
    real *8, allocatable, dimension(:,:) :: zptstarg, zptssrc, zcirc
    real *8, allocatable, dimension(:,:) :: znrmtarg, znrmsrc
    real *8, allocatable, dimension(:,:) :: quadvec, dir1, dir2
    real *8, allocatable, dimension(:) :: quadstr
    real *8, allocatable, dimension(:,:) :: target, source
    real *8 :: ztarg(2),center(2),pot,grad(2),hess(3),lambda,zdiff(2), zarb(2)

    complex *16, allocatable, dimension(:) :: hvec,hder,zquadstr,mpoleqp, &
         mpoleqplap, potqp, potqplap, ympoleqp, ylocqp
    real *8, allocatable :: potqp1(:), gradqp1(:,:), hessqp1(:,:), potqp2(:)
    complex *16 :: eye, z, ztemp1, ztemp2, zk

    complex *16 gradtemp, hesstemp

    real *8 :: work(10000)

    integer nterms, npts, ifpot, ifgrad, ifhess
    parameter (iseed = 281+3308004)
    data eye /(0.0d0,1.0d0)/

    ! set parameters
  
    lambda = 1.0d-3
    npts = 80
    nterms = npts/2
    nwsave = 2*npts+15
    l = (npts+1)/2

    zk = eye*lambda

    ! scale (source points are within a [-.5*scale,.5*scale]^2 box)

    scale = 1.0d0
    rscale = min(scale*1.01d0*lambda,1.0d0)
    rscalelap = min(scale*1.01d0,1.0d0)
    ns = 20
    nt = 20


    ! allocate storage

    allocate(u(npts),up(npts))
    allocate(wsave(nwsave))
    allocate(target(2,nt),source(2,ns),zptssrc(2,npts),zptstarg(2,npts))
    allocate(znrmtarg(2,npts),znrmsrc(2,npts))
    allocate(zquadstr(ns),quadstr(ns),quadvec(3,ns),dir1(2,ns),dir2(2,ns))
    allocate(utargex(nt),gradtargex(2,nt),utarg(nt),gradtarg(2,nt))
    allocate(hesstarg(3,nt),hesstargex(3,nt))
    allocate(gradtargloc(2,nt),hesstargloc(3,nt))
    allocate(mbhmp(0:npts),ymp(0:npts),diffs(0:l),ders(0:l))
    allocate(mbhloc(0:npts),lloc(0:npts))
    allocate(mbhloc2(0:npts),lloc2(0:npts))    
    allocate(kvec(0:l+5),kder(0:l+5),hvec(0:l+5),hder(0:l+5))
    allocate(mpoleqp(-nterms:nterms),mpoleqplap(0:nterms))
    allocate(ympoleqp(0:nterms),ylocqp(0:nterms))
    allocate(gradtargmp(2,nt),hesstargmp(3,nt))
    allocate(gradtargqp(2,nt),hesstargqp(3,nt))
    allocate(gradtargqplap(2,nt),hesstargqplap(3,nt))
    allocate(potqp(nt),potqplap(nt),utargmp(nt),utargloc(nt))
    allocate(zcirc(2,npts))
    allocate(potqp1(nt),gradqp1(2,nt),hessqp1(3,nt),potqp2(nt))
    call prini(6,13)
    temp = hkrand(iseed)

    pi = 4.0d0*atan(1.0d0)

    call dffti(npts,wsave)

    zarb(1) = 4.0d0
    zarb(2) = -6.0d0

    ! set up some targets in the well separated box 
    ! [1.5*scale,2.5*scale] x [-.5*scale,.5*scale]

    ztarg(1) = 2.0d0
    ztarg(2) = 0.0d0
    
    do i = 1,nt
       target(1,i) = ztarg(1) -0.5d0 + hkrand(0)
       target(2,i) = ztarg(2) -0.5d0 + hkrand(0)
    enddo

    ztarg(1) = 0.0d0 + ztarg(1)*scale
    ztarg(2) = ztarg(2)*scale

    do i = 1,nt
       target(1,i) = 0.0d0+scale*target(1,i)
       target(2,i) = scale*target(2,i)
    enddo


    ztarg = ztarg+zarb

    do i = 1,nt
       target(1:2,i) = target(1:2,i) + zarb
    enddo

    do i = 1,nt
       write(22,*) target(1,i), target(2,i)
    enddo
    
    ! set up some quadrupole sources in the center box 
    ! [-.5*scale,.5*scale] x [-.5*scale,.5*scale]
    ! with random directions and strengths

    center(1) = 0.0d0
    center(2) = 0.0d0
    
    do i = 1,ns
       ! locations
       source(1,i) = center(1) + scale*(-0.5d0 + hkrand(0))
       source(2,i) = center(2) + scale*(-0.5d0 + hkrand(0))
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

    ! make all the directions the same

    do i = 1,ns
       write(23,*) source(1,i), source(2,i), dir1(1,i), &
            dir1(2,i), dir2(1,i), dir2(2,i)
    enddo

    rad = (1.5d0-dsqrt(2.0d0)/2.0d0+0.5d0-.001d0)*scale
    radtarg = dsqrt(2.0d0)/2.0d0*scale

    ! getting exact values
    
    ifpot = 1
    ifgrad = 1
    ifhess = 1

    do i = 1,nt
       utargex(i) = 0.0d0
       gradtargex(1:2,i) = (/ 0.0d0, 0.0d0 /)
       hesstargex(1:3,i) = (/ 0.0d0, 0.0d0, 0.0d0 /)
       do j = 1,ns
          call modbhgreend2(lambda,target(1,i),source(1,j), &
               ifpot,pot,ifgrad,grad,ifhess,hess,dir1(1,j), &
               dir2(1,j))
          utargex(i) = utargex(i) + quadstr(j)*pot
          gradtargex(1:2,i) = gradtargex(1:2,i) + quadstr(j)*grad
          hesstargex(1:3,i) = hesstargex(1:3,i) + quadstr(j)*hess
       enddo
    enddo

    ! compute expansion coefficients for outgoing and incoming

    call mbh2dsov_circpts(zcirc,npts)

    ncalls = 1
    time1 = second()

    ifcharge = 0
    ifdipole = 0
    ifquad = 1
    ifoct = 0
    
    call mbh2dsov_circvals(u,up,source,ifcharge,charge, &
         ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec, &
         ifoct,octstr,octvec,ns,center,rad,lambda,npts,zcirc)
    
    call mbh2dsov_dtocmp(lambda,rad,rscale,u,up,npts,mbhmp, &
         ymp,wsave,work)
    
    call mbh2dsov_convert_coeffs(mbhmp,ymp,npts)
    
    time2 = second()

    write(*,*) 'time for ', ns*ncalls, 'sources ', time2-time1
    write(*,*) 'time per source ', (time2-time1)/(ns*ncalls)

!    call prin2('mbhmp *',mbhmp,nterms)
!    call prin2('ymp *',ymp,nterms)


    call mbh2dsov_circvals(u,up,source,ifcharge,charge, &
         ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec, &
         ifoct,octstr,octvec,ns,ztarg,radtarg,lambda,npts,zcirc)

    call mbh2dsov_dtocta(lambda,radtarg,rscale,u,up,npts, &
         mbhloc,lloc,wsave,work)

!    call prin2('mbhloc *',mbhloc,nterms)
!    call prin2('lloc *',lloc,nterms)

    ifpot = 1
    ifgrad = 1
    ifhess = 1

    ! evaluate expansions

    call mbh2dmpevalall(lambda,rscale,center,mbhmp,ymp, &
         nterms,target,nt,ifpot,utarg,ifgrad,gradtarg, &
         ifhess,hesstarg)
    
    call mbh2dtaevalall(lambda,rscale,ztarg,mbhloc,lloc, &
         nterms,target,nt,ifpot,utargloc,ifgrad,gradtargloc, &
         ifhess,hesstargloc)
    
    ! form multipole expansion due to quadrupoles 
    ! (modified Helmholtz)

    call y2dformmp_qp(ier,lambda,rscale,source,zquadstr,quadvec, &
         ns,center,nterms,ympoleqp)

    call h2dformta_qp(ier,zk,rscale,source,zquadstr,quadvec, &
         ns,ztarg,nterms,mpoleqp)

    call y2dformta_qp(ier,lambda,rscale,source,zquadstr,quadvec, &
         ns,ztarg,nterms,ylocqp)


    write(*,*) mpoleqp(0)
    write(*,*) ylocqp(0)
    do i = 1,nterms
       write(*,*) i
       write(*,*) mpoleqp(-i)
       write(*,*) mpoleqp(i)
       write(*,*) ylocqp(i)
    enddo
    
    call l2dformmp_qp(ier,rscalelap,source,zquadstr,quadvec, &
         ns,center,nterms,mpoleqplap)

    ifpot = 1
    ifgrad = 1
    ifhess = 1

    !call h2dmpevalall(zk,rscale,center,mpoleqp,nterms,target,nt, &
    !     ifpot,potqp,ifgrad,gradtargqp,ifhess,hesstargqp)

    !call y2dmpevalall(lambda,rscale,center,ympoleqp,nterms, &
    !     target,nt,ifpot,potqp1,ifgrad,gradqp1,ifhess,hessqp1)

    call y2dtaevalall(lambda,rscale,ztarg,ylocqp,nterms, &
         target,nt,ifpot,potqp1,ifgrad,gradqp1,ifhess,hessqp1)

    do i = 1,nt
       potqp(i) = potqp1(i)
       gradtargqp(1:2,i) = gradqp1(1:2,i)
       hesstargqp(1:3,i) = hessqp1(1:3,i)
    enddo

    
    call l2dmpevalall(rscalelap,center,mpoleqplap,nterms,target,nt, &
         ifpot,potqplap,ifgrad,gradtargqplap,ifhess,hesstargqplap)

    do i = 1,nt
       utargmp(i) = dreal(potqp(i)-potqplap(i)/(2.0d0*pi))/(lambda**2)
       gradtargmp(1,i) = dreal(gradtargqp(1,i) &
            -gradtargqplap(1,i)/(2.0d0*pi))/lambda**2
       gradtargmp(2,i) = dreal(gradtargqp(2,i) &
            -gradtargqplap(2,i)/(2.0d0*pi))/lambda**2
       hesstargmp(1,i) = dreal(hesstargqp(1,i) &
            -hesstargqplap(1,i)/(2.0d0*pi))/lambda**2
       hesstargmp(2,i) = dreal(hesstargqp(2,i) &
            -hesstargqplap(2,i)/(2.0d0*pi))/lambda**2
       hesstargmp(3,i) = dreal(hesstargqp(3,i) &
            -hesstargqplap(3,i)/(2.0d0*pi))/lambda**2
    enddo

    write(*,*) 'LAMBDA    '
    write(*,*) lambda
    write(*,*) 'SIZE OF SOURCE BOX '
    write(*,*) scale
    write(*,*) 'ERROR FOR POTENTIAL WITH:' 
    write(*,*) '(1) NEW OUTGOING BASIS FUNCTIONS '
    write(*,*) '(2) NEW LOCAL BASIS FUNCTIONS '
    write(*,*) '(3) DIFFERENCE OF YUKAWA AND LAPLACE '

    do i = 1,nt
       err_newmp = dabs(utarg(i)-utargex(i))
       err_newta = dabs(utargloc(i)-utargex(i))
       err_oldmp = dabs(utargmp(i)-utargex(i))
       dnorm = dabs(utargex(i))
       write(*,*) err_newmp/dnorm, err_newta/dnorm, err_oldmp/dnorm
    enddo

    write(*,*) 'ERROR FOR GRAD WITH:' 
    write(*,*) '(1) NEW OUTGOING BASIS FUNCTIONS '
    write(*,*) '(2) NEW LOCAL BASIS FUNCTIONS '
    write(*,*) '(3) DIFFERENCE OF YUKAWA AND LAPLACE '

    do i = 1,nt
       err_newmp = dsqrt((gradtarg(1,i)-gradtargex(1,i))**2 + &
            (gradtarg(2,i)-gradtargex(2,i))**2)
       err_newta = dsqrt((gradtargloc(1,i)-gradtargex(1,i))**2 + &
            (gradtargloc(2,i)-gradtargex(2,i))**2)
       err_oldmp = dsqrt((gradtargmp(1,i)-gradtargex(1,i))**2 + &
            (gradtargmp(2,i)-gradtargex(2,i))**2)
       dnorm = dsqrt(gradtargex(1,i)**2+gradtargex(2,i)**2)
       write(*,*) err_newmp/dnorm, err_newta/dnorm, err_oldmp/dnorm            
    enddo
    

    write(*,*) 'ERROR FOR HESS WITH:' 
    write(*,*) '(1) NEW OUTGOING BASIS FUNCTIONS '
    write(*,*) '(2) NEW LOCAL BASIS FUNCTIONS '
    write(*,*) '(3) DIFFERENCE OF YUKAWA AND LAPLACE '

    do i = 1,nt
       err_newmp = dsqrt((hesstarg(1,i)-hesstargex(1,i))**2 + &
            (hesstarg(2,i)-hesstargex(2,i))**2 + &
            (hesstarg(3,i)-hesstargex(3,i))**2)
       err_newta = dsqrt((hesstargloc(1,i)-hesstargex(1,i))**2 + &
            (hesstargloc(2,i)-hesstargex(2,i))**2 + &
            (hesstargloc(3,i)-hesstargex(3,i))**2)
       err_oldmp = dsqrt((hesstargmp(1,i)-hesstargex(1,i))**2 + &
            (hesstargmp(2,i)-hesstargex(2,i))**2 + &
            (hesstargmp(3,i)-hesstargex(3,i))**2)
       dnorm = dsqrt(hesstargex(1,i)**2+hesstargex(2,i)**2 + &
            hesstargex(3,i)**2)
       write(*,*) err_newmp/dnorm, err_newta/dnorm, err_oldmp/dnorm, dnorm            
    enddo

    call prin2('mbhloc *',mbhloc,npts)
    call prin2('lloc *',lloc,npts)

    stop
  end program mbhrouts_dr
