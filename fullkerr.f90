program wrapper
! gfortran amodules.f90 super_kerr.f90
    implicit none
    integer ne,i, ifl,jmax, nr, np
    parameter (nr=300) 
    parameter (np=300)
    parameter (ne=500,jmax=8)
    real Emax,Emin,ear(0:ne),photar(ne),E,dE
    double precision D_kpc, norm_full_kerr, rgrid(10000), tgrid(10000)
    real fk_param(8)
    integer :: num_args, ix
    character(len=12), dimension(:), allocatable :: args

    
    num_args = command_argument_count()
    allocate(args(num_args))  ! I've omitted checking the return status of the allocation 
  
    do ix = 1, num_args
        call get_command_argument(ix,args(ix))
    end do
  
    read( args(1), '(f10.0)' )  fk_param(1)   !a     " "     spin
    read( args(2), '(f10.0)' )  fk_param(2)   !inc   deg     inclination
    read( args(3), '(f10.0)' )  fk_param(3)   !M     Msun    Black hole mass
    read( args(4), '(f10.0)' )  fk_param(4)   !mdot  Edd     Accrertion rate   
    read( args(5), '(f10.0)' )  fk_param(5)   !delta_J 
    ! read( args(6), '(f10.0)' )  fk_param(6)   !alp 
    read( args(6), '(f10.0)' )  fk_param(6)   !f1
    read( args(7), '(f10.0)' )  fk_param(7)   !f2 
    read( args(8), '(f10.0)' )  fk_param(8)   !chi 
    read( args(9), '(f10.0)' ) D_kpc         !source distance in kilo-parsec
    ! Important note: I am defining the Eddington mass accretion rate to be
    ! \dot M_edd = L_edd/c^2 = 1.26e31/c^2 * (M_bh/M_sun) [kg/s]. 
    ! or explicitly 
    ! \dot M_edd = 0.1402 (M_bh/M_sun) [10^15 kg/s]. 

  ! Set energy grid
    Emax  = 75.0
    Emin  = 0.05
    do i = 0,ne
      ear(i) = Emin * (Emax/Emin)**(real(i)/real(ne))
    end do

    ! call return_temp(fk_param, rgrid, tgrid)
    ! do i = 1, 10000
    !   write(103,*) rgrid(i), tgrid(i)
    ! end do 

    call fullkerr(ear,ne,fk_param,ifl,photar)
    norm_full_kerr = 1/D_kpc**2.0 


    
  ! Write out model output
    open(99, file = 'fk_output.dat')
    do i = 1,ne
       E  = 0.5 * ( ear(i) + ear(i-1) )
       dE =         ear(i) - ear(i-1)
       write(99,*)E, norm_full_kerr*E**2*photar(i)/dE
    end do
    close(99) 
  end program wrapper
      
! Can't use an include for compiling within XSPEC
! include 'amodules.f90'

  subroutine return_temp(param, rgrid, tgrid, tcgrid, sgrid)
        ! Calculates disk temperature profile and writes it to fort.99
    implicit none
    integer i
    real param(8)
    double precision a,inc,m,mdot,pi,rin,rout,mu0
    double precision Ca,Cb,Cg,xa,xb,xg,j0,xI,Tscale
    double precision delta_j ! ISCO boundary condition 
    double precision kT_fk, kT_I 
    double precision ris, rho,  dr, disco, csI, cs_I
    double precision e1, e2, ef, e_prime, e3, e4, e_dot
    double precision p1, pf
    double precision rh1, rh2, rhf, rh_prime, rh3, rh4, rh_dot
    double precision r_isco, s_isco, eps, the_func, rg, k_es, kff0 
    double precision dx, dt, f1, f2, chi, alp 
    double precision rgrid(10000), tgrid(10000), sgrid(10000)
    double precision tcgrid(10000), rhgrid(10000), kffgrid(10000)
    double precision tau_star_grid(10000)


    pi  = acos(-1.d0)
  
  ! Parameters
    a       = dble( param(1) )                !Spin parameter
    inc     = dble( param(2) ) * pi / 180.d0  !Inclination (degrees)
    m       = dble( param(3) )                !Black hole mass (solar)
    mdot    = dble( param(4) )                !Mass accretion rate (Eddington)
    delta_j = dble( param(5) )                !ISCO stress parameter
    alp     = 0.1!dble( param(6) )                !Stress parameter
    f1      = dble( param(6) ) 
    f2      = dble( param(7) ) 
    chi     = dble( param(8) ) 
    ! Important note: I am defining the Eddington mass accretion rate to be
    ! \dot M_edd = L_edd/c^2 = 1.26E31/c^2 * (M_bh/M_sun) [kg/s]. 
    ! or explicitly 
    ! \dot M_edd = 0.1402 (M_bh/M_sun) [10^15 kg/s]. 
    
  ! Derived and hardwired quantities
    rin     = 1 + (1 - a**2.0)!event horizon
    rho     = 1 + (1 - a**2.0)!event horizon
    ris     = disco(a)!isco radius
    rout    = 5e3
    mu0     = cos(inc)
    mdot    = mdot * m * 0.1402!Convert mdot from Eddington ratio to units of 10^{15} kg/s.
    ! See above for comment on Eddington mass accretion rate definition.  

    rg =  6.67e-11 * m * 2.0e30 / (3e8**2.0)
    k_es = 0.034!m^2/kg
    kff0 = 6.4e22 * 1e-4 * 1e3  * 1e3 * 1e-6!m^2/kg
 
  ! Initialise temperature calculation
    call Tconstants(a,Ca,Cb,Cg,xa,xb,xg,j0,xI,Tscale)
    cs_I = csI(param)
    kT_I = kT_fk(ris,a,Ca,Cb,Cg,xa,xb,xg,j0,xI,Tscale,mdot,m,delta_j,cs_I)
    eps = (alp**0.5 * cs_I) * (3.0/2.0)**0.5 * ris**0.5 
    r_isco = ris * rg
    s_isco = mdot * 1e15 / (2 * pi * r_isco * alp**0.5 * cs_I * 3e8)

    !Assumed Mdot is in units of 10^{15} kg/s
    dr = (ris - rho)/10000.
    rgrid(1) = ris * rg
    tgrid(1) = kT_I * (3.0 * k_es * s_isco / 8.0)**0.25
    sgrid(1) = s_isco

    ! write(89, *) kT_I, cs_I, eps, s_isco, r_isco, tgrid(1) 

    dx = 0.00001
    dt = 0.00001

    do i = 2, 10000
      rgrid(i) = rgrid(i-1) - dr * rg

      e1 = ef(rgrid(i-1), tgrid(i-1), m, r_isco, s_isco, eps)
      e2 = ef(rgrid(i-1)-dx, tgrid(i-1), m, r_isco, s_isco, eps)
      e_prime = (e1 - e2)/dx
      e3 = ef(rgrid(i-1), tgrid(i-1), m, r_isco, s_isco, eps)
      e4 = ef(rgrid(i-1), tgrid(i-1)-dt, m, r_isco, s_isco, eps)
      e_dot = (e3 - e4)/dt

      ! write(88, *) e1, e2, e_prime, e3, e4, e_dot
      
      p1 = pf(rgrid(i-1), tgrid(i-1), m, r_isco, s_isco, eps)

      ! write(86, *) e1, p1
      
      rh1 = rhf(rgrid(i-1), tgrid(i-1), m, r_isco, s_isco, eps)
      rh2 = rhf(rgrid(i-1)-dx, tgrid(i-1), m, r_isco, s_isco, eps)
      rh_prime = (rh1 - rh2)/dx
      rh3 = rhf(rgrid(i-1), tgrid(i-1), m, r_isco, s_isco, eps)
      rh4 = rhf(rgrid(i-1), tgrid(i-1)-dt, m, r_isco, s_isco, eps)
      rh_dot = (rh3 - rh4)/dt
      
      ! write(87, *) rh1, rh2, rh_prime, rh3, rh4, rh_dot

      the_func = (e_prime - (e1+p1)/rh1 * rh_prime)
      the_func = the_func/(e_dot - (e1+p1)/rh1 * rh_dot)

      tgrid(i) = tgrid(i-1) + the_func * dr * rg

      sgrid(i) = s_isco * (r_isco/rgrid(i)) * 1/(1.0+1/eps*(-1.0 + r_isco/rgrid(i))**1.5)

      ! write(90,*) rgrid(i)/rg, tgrid(i) / (3.0 * k_es * sgrid(i) / 8.0)**0.25
    end do 

    do i = 1, 10000
      tcgrid(i) = tgrid(i)
      rhgrid(i) = rhf(rgrid(i), tcgrid(i), m, r_isco, s_isco, eps)
      kffgrid(i) = kff0 * rhgrid(i) * (tcgrid(i)/8.6173E-8)**(-3.5)
      tgrid(i) = tgrid(i) / (3.0 * (k_es + kffgrid(i)) * sgrid(i) / 8.0)**0.25
      tau_star_grid(i) = sgrid(i) * (3 * (k_es + kffgrid(i)) * kffgrid(i)) ** 0.5
      rgrid(i) = rgrid(i)/rg 
      ! write (45, *) rgrid(i), kffgrid(i), tcgrid(i), rhgrid(i)
    end do
    

    return 
    end subroutine return_temp

    
function csI(param) 
    real param(8)
    double precision a,inc,m,mdot,pi,rin,rout,mu0
    double precision Ca,Cb,Cg,xa,xb,xg,j0,xI,Tscale
    double precision delta_j ! ISCO boundary condition 
    double precision kT_fk 
    double precision ris, rho, disco, T_isco
    double precision r_isco, G, M_bh, c, kb, O_sb, k_es, mu, mp, M_sun, alp   
    double precision const_A, const_B, M_dot, csI, f_x, f_xp

    pi  = acos(-1.d0)
  
  ! Parameters
    a       = dble( param(1) )                !Spin parameter
    inc     = dble( param(2) ) * pi / 180.d0  !Inclination (degrees)
    m       = dble( param(3) )                !Black hole mass (solar)
    mdot    = dble( param(4) )                !Mass accretion rate (Eddington)
    delta_j = dble( param(5) )                !ISCO stress parameter
    alp     = 0.1!dble( param(6) )                !ISCO stress parameter
    ! Important note: I am defining the Eddington mass accretion rate to be
    ! \dot M_edd = L_edd/c^2 = 1.26E31/c^2 * (M_bh/M_sun) [kg/s]. 
    ! or explicitly 
    ! \dot M_edd = 0.1402 (M_bh/M_sun) [10^15 kg/s]. 
    
  ! Derived and hardwired quantities
    rin     = 1 + (1 - a**2.0)!event horizon
    rho     = 1 + (1 - a**2.0)!event horizon
    ris     = disco(a)!isco radius
    rout    = 5e3
    mu0     = cos(inc)
    mdot    = mdot * m * 0.1402!Convert mdot from Eddington ratio to units of 10^{15} kg/s.
    ! See above for comment on Eddington mass accretion rate definition.  

    G = 6.6743e-11
    c = 2.9979e8 
    M_sun = 1.988e30 
    k_es = 0.034 
    kb = 1.38e-23 
    mu = 0.615 
    mp = 1.672e-27 
    O_sb = 5.67e-8 
    M_bh = m * M_sun
    M_dot = mdot * 1e15 
    

    ! Initialise temperature calculation
    call Tconstants(a,Ca,Cb,Cg,xa,xb,xg,j0,xI,Tscale)
    T_isco = kT_fk(ris+0.0001,a,Ca,Cb,Cg,xa,xb,xg,j0,xI,Tscale,mdot,m,delta_j,delta_j)/8.6173E-8!Kelvin
    r_isco = ris * G * M_bh / c**2.0 
    
    const_A = (3.0*k_es/(8.0))**0.25 * kb * T_isco / (mu * mp) * (M_dot/(2*pi*r_isco*alp**0.5))**0.25
    const_B = (rI**3.0/(2.0*G*M_bh))**0.5 * O_sb*T_isco**4.0/(2.0*c) * k_es

    csI = c
    do i = 1, 100
       f_x = ((csI)**(9.0/4.0) - const_A - const_B * (csI)**(5.0/4.0))
       f_xp = (9.0/4.0 * (csI)**(5.0/4.0) - 5.0/4.0 * const_B * (csI)**(1.0/4.0))
       csI = csI  - f_x/f_xp! Newton-Raphson
    end do 

    csI = csI / c

    return 

  end function csI

function ef(r, t, M, rI, sI, eps)
  double precision ef, r, t, M, rI, sI, eps
  ef = 10316877524046.0*t**4.0 
  ef =ef+9.6699684969198e+15*(M*rI**3.0*sI**2.0*t/(r**6.0*(1.0+1.0/eps*(-1.0+rI/r)**1.5)**2.0) + 1.26474776863348e-7*t**8.0)**0.5 
  return 
end function ef

function pf(r, t, M, rI, sI, eps)!! t in keV, r in meters, M in M_sun, sI in sI (ironically), eps usual. 
  double precision pf, r, t, M, rI, sI, eps
  pf = 2292639449788.01*t**4.0  
  pf =pf+6.4466456646132e+15*(M*rI**3.0*sI**2.0*t/(r**6.0*(1.0+1.0/eps*(-1.0 + rI/r)**1.5)**2.0) + 1.26474776863348e-7*t**8.0)**0.5 
  return 
end function pf

function rhf(r, t, M, rI, sI, eps)
  double precision rhf, r, t, M, rI, sI, eps
  rhf = 2292639449788.01*t**4.0 
  rhf=rhf+6.4466456646132e+15*(M*rI**3.0*sI**2.0*t/(r**6.0*(1.0+1.0/eps*(-1.0 + rI/r)**1.5)**2.0) + 1.26474776863348e-7*t**8.0)**0.5
  rhf=2.665332e+20*M*rI**3.0*sI**2.0/(r**6.0*(1.0+1.0/eps*(-1.0 + rI/r)**1.5)**2.0*rhf) 
return 
end function rhf


!=======================================================================
  subroutine fullkerr(ear,ne,param,ifl,photar)
    ! Calculates observed disk spectrum
      use internal_grids
      implicit none
      integer ne,ifl,i,j,k
      real ear(0:ne),param(8),photar(ne)
      double precision a,inc,m,mdot,pi,rin,rout,mu0
      double precision Ca,Cb,Cg,xa,xb,xg,j0,xI,Tscale
      double precision delta_j ! ISCO boundary condition 
      double precision rnmin,rnmax,d,rfunc,disco,mudisk,re,kT_fk, fc, T0
      double precision alpha(nro,nphi),beta(nro,nphi),dOmega(nro,nphi)
      double precision alphan(nro,nphi),betan(nro,nphi),dOmegan(nro,nphi)
      double precision g,dlgfac,dlgfac_isco,kT,phie, ris, rho, csI, cs_I
      double precision rgrid(10000), tgrid(10000), tcgrid(10000), sgrid(10000)
      double precision ri_min, ri_max, dri, rg, sI, kT_I
      double precision f1, f2, chi, alp
      real dNdE_! NAN killer. 
      real mybbody,E,dE,kTcol,dNbydE(nec),Eem,dEem
      integer t_r(nro, nphi), ihi, ilo
      logical needtrace
      pi  = acos(-1.d0)
      ifl = 1
    
    ! Parameters
      a       = dble( param(1) )                !Spin parameter
      inc     = dble( param(2) ) * pi / 180.d0  !Inclination (degrees)
      m       = dble( param(3) )                !Black hole mass (solar)
      mdot    = dble( param(4) )                !Mass accretion rate (Eddington)
      delta_j = dble( param(5) )                !ISCO stress parameter
      alp     = 0.1!dble( param(6) )                !Stress parameter
      f1      = dble( param(6) ) 
      f2      = dble( param(7) ) 
      chi     = dble( param(8) ) 
      ! Important note: I am defining the Eddington mass accretion rate to be
      ! \dot M_edd = L_edd/c^2 = 1.26E31/c^2 * (M_bh/M_sun) [kg/s]. 
      ! or explicitly 
      ! \dot M_edd = 0.1402 (M_bh/M_sun) [10^15 kg/s]. 
      
    ! Derived and hardwired quantities
      rin     = 1 + (1 - a**2.0)!event horizon
      rho     = 1 + (1 - a**2.0)!event horizon
      ris     = disco(a)!isco radius
      rout    = 1000
      mu0     = cos(inc)
      mdot    = mdot * m * 0.1402!Convert mdot from Eddington ratio to units of 10^{15} kg/s.
      ! See above for comment on Eddington mass accretion rate definition.  
      
    ! Initialize
      if( firstcall )then
         firstcall = .false.
         !Define coarse internal energy grid
         Emax  = 75.0
         Emin  = 0.05
         do i = 0,nec
           earc(i) = Emin * (Emax/Emin)**(real(i)/real(nec))
         end do
         !Assign impossible initial values to previous parameters
         ! Note that only impossible parameter is now cos(inc). 
         aprev   = 10.d0
         mu0prev = 10.d0
      end if
    
    ! Initialise temperature calculation
      call Tconstants(a,Ca,Cb,Cg,xa,xb,xg,j0,xI,Tscale) 
      cs_I = csI(param)
      !Assumed Mdot is in units of 10^{15} kg/s
      
    ! Set up full GR grid
      rnmax = 300d0                           !Sets outer boundary of full GR grid
      rnmin = rfunc(a,mu0)                    !Sets inner boundary of full GR grid
      call impactgrid(rnmin,rnmax,mu0,nro,nphi,alpha,beta,dOmega)
      d     = max( 1.0d4 , 2.0d2 * rnmax**2 ) !Sensible distance to BH  
    
    ! Set up `straight lines' grid
      rnmax = rout                            !Sets outer boundary of Newtonian grid
      rnmin = 300.d0                          !Sets inner boundary of Newtonian grid
      call impactgrid(rnmin,rnmax,mu0,nro,nphi,alphan,betan,dOmegan)
      
    ! Do the ray tracing in full GR
      needtrace = .false.
      if( abs( a - aprev ) .gt. tiny(a) ) needtrace = .true.
      if( abs(mu0 - mu0prev) .gt. tiny(mu0) ) needtrace = .true.
      mudisk = 0.d0       !razor thin disk
      if( needtrace )then
         call dGRtrace(nro,nphi,alpha,beta,mu0,a,rin,rout,mudisk,d,pem1,re1, t_r)
      end if
      aprev = a
      mu0prev = mu0

      call return_temp(param, rgrid, tgrid, tcgrid, sgrid)

      ! ! Set up radial grid
      ri_min = ris
      ri_max = rho
      dri    = ( ri_max - ri_min ) / real(10000-1)
      rg = 6.67e-11 * m * 2e30/(3e8**2.0)
      sI = sgrid(1)
      kT_I = kT_fk(ris+0.0001,a,Ca,Cb,Cg,xa,xb,xg,j0,xI,Tscale,mdot,m,delta_j,cs_I)

      
      ! Loop through inner relativistic grid
      dNbydE = 0.0
      do j = 1,nphi
        do i = 1,nro
          if( pem1(i,j) .gt. 0.d0 )then
            re = re1(i,j)
            if( re .gt. rin .and. re .le. rout )then
              !Calculate g-factor
              g = 0.0
              if( re .gt. ris .and. re .le. rout )then
                g = dlgfac( a,mu0,alpha(i,j),re )
                !Calclulate temperature
                kT = kT_fk(re,a,Ca,Cb,Cg,xa,xb,xg,j0,xI,Tscale,mdot,m,delta_j,cs_I)
                !Calculate colour-temperature
                T0 = kT / 8.6173E-8 
                fc = f1!fcol_DEA(T0, re, a, m, mdot)
                kTcol = kT * fc

              end if 
              
              if( re .gt. rho+0.000001 .and. re .le. ris) then 
                g = dlgfac_isco( a, mu0, alpha(i,j), beta(i, j), re, t_r(i, j))
                ihi = ceiling( (re-ri_min)/dri ) + 1
                ihi = max( ihi , 2 )
                ilo = ihi - 1
                kT  = tgrid(ilo) + ( tgrid(ihi) - tgrid(ilo) ) * ( re - rgrid(ilo) ) / dri

                fc = f2 * (ris/re)**chi
                kTcol = kT * fc

              end if 
                  do k = 1,nec
                     E          = 0.5 * ( earc(k) + earc(k-1) )
                     dE         =       ( earc(k) - earc(k-1) )
                     Eem        =  E / g
                     dEem       = dE / g
                     dNdE_      = g**3 * kT**4 * mybbody(kTcol,Eem,dEem) * dOmega(i,j) / dE
                     if (dNdE_ .eq. dNdE_) then 
                      dNbydE(k) = dNbydE(k) + dNdE_
                     else
                      dNbydE(k) = dNbydE(k)
                     end if 
                  end do
               end if
           end if
         end do
      end do
    
    ! Loop through outer Newtonian grid
      do j = 1,nphi
         do i = 1,nro
            call drandphithick(alphan(i,j),betan(i,j),mu0,mudisk,re,phie)
            if( re .gt. rin .and. re .le. rout )then
               !Calclulate temperature
               kT = kT_fk(re,a,Ca,Cb,Cg,xa,xb,xg,j0,xI,Tscale,mdot,m,delta_j,cs_I)
               !Calculate colour-temperature
               T0 = kT / 8.6173E-8 
               fc = f1
               kTcol = kT * fc
              !Calculate g-factor
               g = dlgfac( a,mu0,alpha(i,j),re )
             do k = 1,nec
                  E         = 0.5 * ( earc(k) + earc(k-1) )
                  dE        = earc(k) - earc(k-1)
                  Eem        =  E / g
                  dEem       = dE / g
                  dNdE_      = g**3 * kT**4 * mybbody(kTcol,Eem,dEem) * dOmega(i,j) / dE
                  if (dNdE_ .eq. dNdE_) then 
                   dNbydE(k) = dNbydE(k) + dNdE_
                  else
                   dNbydE(k) = dNbydE(k)
                  end if 
               end do
            end if
         end do
      end do
    
    ! Rebin onto input grid
      call myinterp(nec,earc,dNbydE,ne,ear,photar)
    
    ! Multiply by dE
      do i = 1,ne
         E  = 0.5 * ( ear(i) + ear(i-1) )
         dE = ear(i) - ear(i-1)
         photar(i) = photar(i) * dE
      end do
      
    ! Re-normalise so that norm = [ 1 / Dkpc ]^2
      photar = 0.467842078 * m**2 * photar
      
      return
    end subroutine fullkerr
    !=======================================================================
    


!-----------------------------------------------------------------------
function kT_fk(re,a,Ca,Cb,Cg,xa,xb,xg,j0,xI,Tscale,mdot,m,delta_j,cs_I)
  ! mdot in units of 10^{15} kg s^{-1}
  ! m in units of solar
    implicit none
    double precision kT_fk,re,a,Ca,Cb,Cg,xa,xb,xg,j0,xI,Tscale
    double precision mdot,m, eps_inv, rI
    double precision x,br1,br2
    double precision delta_j, kT_I, cs_I
  
    x    = sqrt(re)
    rI   = xI**2.0

    if (x .ge. xI ) then   
      br1 = dble((1.0 - (3.0*a)/(2.0*x)*log(x) + Ca/x*log(x-xa) + Cb/x*log(x-xb) + Cg/x*log(x-xg) - (1-delta_j)*j0*xI/x)**0.25)
      br2 = dble( (1.0/(1.0 - 3.0/(x**2) + 2.0*a/(x**3)))**0.25 )
      kT_fk = Tscale * (mdot**0.25) / (m**0.5) * re**(-0.75) * br1 * br2
    else
      br1 = dble((1.0 - (3.0*a)/(2.0*xI)*log(xI) + Ca/xI*log(xI-xa) + Cb/xI*log(xI-xb) + Cg/xI*log(xI-xg) - (1-delta_j)*j0)**0.25)
      br2 = dble( (1.0/(1.0 - 3.0/(xI**2) + 2.0*a/(xI**3)))**0.25 )
      kT_I = Tscale * (mdot**0.25) / (m**0.5) * xI**(-1.5) * br1 * br2
      eps_inv = 1/(0.1**0.5 * cs_I * (3.0/2.0)**0.5 * xI)
      kT_fk = kT_I * (rI/re)**(17.0/28.0) * (eps_inv * ((rI/re) - 1)**1.5 + 1.0)**(-1.0/28.0)
      ! kT_fk = kT_I * (rI/re)**(5.0/4.0) * (eps_inv * ((rI/re) - 1)**1.5 + 1.0)**(-1.0/4.0)
      ! write(98, *) cs_I
    end if 
  

    return
  end function kT_fk  
  !-----------------------------------------------------------------------
      
      
!-----------------------------------------------------------------------
subroutine Tconstants(a,Ca,Cb,Cg,xa,xb,xg,j0,xI,Tscale)
  ! Defines physical parameters and derived quantities for relativistic
  ! temperature calculation
  ! mdot in units of 10^{15} kg s^{-1}
  ! m in units of solar
  ! Input: a,rin
  ! Output: xa,xb,xg,j0,xI,Tscale,chi,psi,c_a,r_b,i_b
    implicit none
    double precision a
    double precision Ca,Cb,Cg,xa,xb,xg,j0,xI,Tscale
    double precision rI,disco,pi
    pi = acos(-1.d0)
    Tscale  = 8.31697777 !Tempaerature scale, in keV, for a one solar mass BH and an 10^15 kg/s accretion rate)
    if (abs(a) .lt. 1.0) then
      xa = dble( 2.0*cos(1.0/3.0 * acos(-a) ) )
      xb = dble( 2.0*cos(1.0/3.0 * acos(-a) - 2.0*pi/3.0 ) )
      xg = dble( 2.0*cos(1.0/3.0 * acos(-a) - 4.0*pi/3.0 ) )
      
      Ca = dble( (2.0*xa - a*(1.0 + xa*xa))/(2.0*(1.0 - xa*xa)) )
      Cb = dble( (2.0*xb - a*(1.0 + xb*xb))/(2.0*(1.0 - xb*xb)) )
      Cg = dble( (2.0*xg - a*(1.0 + xg*xg))/(2.0*(1.0 - xg*xg)) )
      
      rI = disco(a) ! Finds the ISCO in gravitational units
      xI = dble( sqrt(rI) )
      
      j0=dble( (1.0 - (3.0*a)/(2.0*xI) * log(xI) + Ca/xI * log(xI-xa) + Cb/xI * log(xI-xb) + Cg/xI * log(xI-xg)) )
      ! Vanishing stress angular momentum constant (in natural units).     
    end if
    if (a .eq. 1.0) then
      rI = disco(a)
      xI = dble( sqrt(rI))
      j0 = 1.0 + 3.0/(2.0*xI) * log(xI + 2.0) - 3.0/(2.0*xI)*log(xI)
      ! Vanishing stress angular momentum constant (in natural units). 
    end if 
    if (a .eq. -1.0) then
      rI = disco(a)
      xI = dble( sqrt(rI))
      j0 = 1.0 - 3.0/(2.0*xI) * log(xI - 2.0) + 3.0/(2.0*xI)*log(xI)
      ! Vanishing stress angular momentum constant (in natural units). 
    end if 
    return
  end subroutine Tconstants
  !-----------------------------------------------------------------------
  
  
    

!=======================================================================
subroutine raytrace_grid(alpha, beta, param, g_fac, rs, t_r)
  ! Calculates observed disk spectrum
    use internal_grids
    implicit none
    integer ifl, i, j
    real param(2)
    double precision a,inc,pi,rin,rout,mu0, rho, ris, rt_out
    double precision rnmin,rnmax,d,rfunc,disco,mudisk,re
    double precision alpha(nro,nphi),beta(nro,nphi),dOmega(nro,nphi), g_fac(nro, nphi)
    double precision alphan(nro,nphi),betan(nro,nphi),dOmegan(nro,nphi)
    double precision dlgfac, rs(nro, nphi), dlgfac_isco
    integer t_r(nro, nphi)
    logical needtrace
    
    pi  = acos(-1.d0)
    ifl = 1
  
  ! Parameters
    a       = dble( param(1) )                !Spin parameter
    inc     = dble( param(2) ) * pi / 180.d0  !Inclination (degrees)
  
    ! Derived and hardwired quantities
    rin     = 1 + (1 - a**2.0)!event horizon
    rho     = 1 + (1 - a**2.0)!event horizon
    ris     = disco(a)!isco radius
    rout    = 35
    rt_out  = 15!
    mu0     = cos(inc)
    
  
  ! Set up full GR grid
    rnmax = min( 300d0, rout )              !Sets outer boundary of full GR grid
    rnmin = rfunc(a,mu0)                    !Sets inner boundary of full GR grid
    call impactgrid(rnmin,rnmax,mu0,nro,nphi,alpha,beta,dOmega)
    d     = max( 1.0d4 , 2.0d2 * rnmax**2 ) !Sensible distance to BH  
  
  ! Set up `straight lines' grid
    rnmax = rout                            !Sets outer boundary of Newtonian grid
    rnmin = 300.d0                          !Sets inner boundary of Newtonian grid
    call impactgrid(rnmin,rnmax,mu0,nro,nphi,alphan,betan,dOmegan)
    
  ! Do the ray tracing in full GR
    needtrace = .false.
    if( abs( a - aprev ) .gt. tiny(a) ) needtrace = .true.
    if( abs(mu0 - mu0prev) .gt. tiny(mu0) ) needtrace = .true.
    mudisk = 0.0d0       !razor thin disk
    if( needtrace )then
       call dGRtrace(nro,nphi,alpha,beta,mu0,a,rin,rout,mudisk,d,pem1,re1,t_r)
    end if

    do j = 1,nphi
      do i = 1,nro
        if( pem1(i,j) .gt. 0.d0 )then
           re = re1(i,j)
            if( re .gt. ris .and. re .le. rt_out )then
               g_fac(i, j) = dlgfac( a,mu0,alpha(i,j),re )
               rs(i, j) = re
            end if 
            if( re .gt. rho+0.00001 .and. re .le. ris) then 
              g_fac(i, j) = dlgfac_isco( a, mu0, alpha(i,j), beta(i, j), re, t_r(i, j))
              rs(i, j) = re
            end if 
        end if
      end do
    end do 


 return 

end subroutine raytrace_grid



!-----------------------------------------------------------------------
function dlgfac(a,mu0,alpha,r)
!c Calculates g-factor for a disk in the BH equatorial plane
  implicit none
  double precision dlgfac,a,mu0,alpha,r
  double precision sin0,omega,Delta,Sigma2,gtt,gtp,gpp
  sin0   = sqrt( 1.0 - mu0**2 )
  omega  = 1. / (r**1.5+a)
  Delta  = r**2 - 2*r + a**2
  Sigma2 = (r**2+a**2)**2 - a**2 * Delta
  gtt    = 4*a**2/Sigma2 - r**2*Delta/Sigma2
  gtp    = -2*a/r
  gpp    = Sigma2/r**2
  dlgfac = sqrt( -gtt - 2*omega*gtp - omega**2.*gpp )
  dlgfac = dlgfac / ( 1.+omega*alpha*sin0 )
  return
end function dlgfac
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
function dlgfac_isco(a, mu0, alpha, beta, r, t_r)
!c Calculates g-factor for a disk in the BH equatorial plane
  implicit none
  double precision dlgfac_isco,a,mu0, alpha, r, disco, beta
  double precision sin0,Delta
  double precision kr,vt,vr, vp
  double precision ris, eis, jis
  double precision lam, q 
  integer t_r

  ris = disco(a)
  vr = -(2./(3.*ris))**0.5 * (ris/r - 1.0)**1.5
  eis = (1. - 2./(3.*ris))**0.5
  jis = 2. * 3.**0.5 * (1 - 2.*a/(3.*ris**0.5))

  Delta  = r**2. - 2.*r + a**2.0

  vp = (2*eis*a + jis*(r-2))/(r*(r**2.0-2.*r+a**2.0))
  vt = (eis*(r**3.0 + r * a**2.0 + 2.0 * a**2.0) - 2.0*jis*a)/(r*(r**2.0-2.*r+a**2.0))

  sin0   = sqrt( 1.0 - mu0**2 )
  lam = -alpha * sin0

  q = (beta**2.0 - a**2.0 * mu0**2.0 + alpha**2.0 * mu0**2.0)
  
  kr=(r**4.-(q+lam**2.-a**2.)*r**2.+2.*r*(q+(lam-a)**2.) - a**2.0*q)/Delta**2.0

  dlgfac_isco=1/(+vt - (-1.0)**t_r * sqrt(kr)*vr - vp*lam)
  
  return
end function dlgfac_isco
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine dGRtrace(nro,nphi,alpha,beta,mu0,spin,rmin,rout,mudisk,d,pem1,re1,t_r)
! Traces rays in the Kerr metric for a camera defined by the impact
! parameters at infinity: alpha(nro,nphi) and beta(nro,nphi).
! Traces back to a disk defined by mudisk = cos(theta_disk), where
! theta_disk is the angle between the vertical and the disk surface.
! i.e. tan( theta_disk ) = 1 / (h/r)
! OUTPUT:
! pem1(nro,nphi)
! pem > 1: there is a solution
! pem = -1 photon goes to infinity without hitting disk surface
! pem = -2 photon falls into horizon without hitting disk surface
! re1(nro,nphi)      radius that the geodesic hits the disc
  use blcoordinate     ! This is a YNOGK module
  implicit none
  integer nro,nphi,i,j
  double precision alpha(nro,nphi),beta(nro,nphi),mu0,spin,rmin,rout,mudisk,d
  double precision pem1(nro,nphi),re1(nro,nphi)
  integer t_r(nro, nphi)
  double precision cos0,sin0,scal,velocity(3),f1234(4),lambda,q
  double precision pem,re,mucros,phie,taudo,sigmacros
  integer tr1,tr2      
  cos0  = mu0
  sin0  = sqrt(1.0-cos0**2)
  scal     = 1.d0
  velocity = 0.d0
  re1      = 0.0
  do i = 1,nro
    do j = 1,NPHI
      call lambdaq(-alpha(i,j),-beta(i,j),d,sin0,cos0,spin,scal,velocity,f1234,lambda,q)
      pem = Pemdisk(f1234,lambda,q,sin0,cos0,spin,d,scal,mudisk,rout,rmin)  !Can try rin instead of rmin to save an if statement
      pem1(i,j) = pem      
      !pem > 1 means there is a solution
      !pem < 1 means there is no solution
      if( pem .gt. 0.0d0 )then
        call YNOGK(pem,f1234,lambda,q,sin0,cos0,spin,d,scal,re,mucros,phie,taudo,sigmacros, tr1, tr2)
        re1(i,j)    = re! Should also check mu_cross 
        t_r(i, j) = tr1
      end if
    end do
  end do
  return
end subroutine dGRtrace
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine impactgrid(rnmin,rnmax,mu0,nro,nphi,alpha,beta,dOmega)
! Calculates a grid of impact parameters
! INPUT:
! rnmin        Sets inner edge of impact parameter grid
! rnmax        Sets outer edge of impact parameter grid
! mu0          Sets `eccentricity' of the grid
! nro          Number of steps in radial impact parameter (b)
! nphi         Number of steps in azimuthal impact parameter (phi)
! OUTPUT:
! alpha(nro,nphi)   Horizontal impact parameter
! beta(nro,nphi)    Vertical impact parameter
! dOmega(nro,nphi)  dalpha*dbeta
  implicit none
  integer nro,nphi,i,j
  double precision rnmin,rnmax,mu0,alpha(nro,nphi),beta(nro,nphi)
  double precision dOmega(nro,nphi),mueff,pi,rar(0:nro),dlogr,rn(nro)
  double precision logr,phin
  pi     = acos(-1.d0)

  mueff = max( mu0 , 0.3d0 )
  
  rar(0) = rnmin
  dlogr  = log10( rnmax/rnmin ) / dble(nro)
  do i = 1,NRO
    logr = log10(rnmin) + dble(i) * dlogr
    rar(i)    = 10.d0**logr
    rn(i)     = 0.5 * ( rar(i) + rar(i-1) )
    do j = 1,nphi
       domega(i,j) = rn(i) * ( rar(i) - rar(i-1) ) * mueff * 2.d0 * pi / dble(nphi)
       phin       = (j-0.5) * 2.d0 * pi / dble(nphi) 
       alpha(i,j) = rn(i)  * sin(phin)
       beta(i,j)  = rn(i) * cos(phin) * mueff
    end do
  end do
  
  return
end subroutine impactgrid
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine myinterp(nfx,farx,Gfx,nf,far,Gf)
! Interpolates the function Gfx from the grid farx(0:nfx) to the
! function Gf on the grid far(0:nf)
  implicit none
  integer nfx,nf
  real farx(0:nfx),Gfx(nfx),far(0:nf),Gf(nf)
  integer ix,j
  real fx(nfx),f,fxhi,Gxhi,fxlo,Gxlo
! Define grid of central input frequencies
  do ix = 1,nfx
      fx(ix) = 0.5 * ( farx(ix) + farx(ix-1) )
  end do
! Run through grid of central output frequencies
  ix = 1
  do j = 1,nf
      !Find the input grid frequencies either side of the current
      !output grid frequency
      f = 0.5 * ( far(j) + far(j-1) )
      do while( fx(ix) .lt. f .and. ix .lt. nfx )
        ix = ix + 1
      end do
      ix = max( 2 , ix )
      fxhi = fx(ix)
      Gxhi = Gfx(ix)
      ix = ix - 1
      fxlo = fx(ix)
      Gxlo = Gfx(ix)
      !Interpolate
      Gf(j) = Gxlo + ( Gxhi - Gxlo ) * ( f - fxlo ) / ( fxhi - fxlo )
  end do
  return
end subroutine myinterp
!-----------------------------------------------------------------------
  

!-----------------------------------------------------------------------
function dISCO(a)
  !ISCO in Rg 
  implicit none
  double precision a,dISCO,z1,z2
  if(abs(a).ge.1.0)then
    z1 = -(a**2.0 - 1.0)**(1.0/3.0)
    z1 = z1 * ( (1.0+abs(a))**(1.0/3.0)-(abs(a)-1.0)**(1.0/3.0))+1.0
  else
    z1 = ( 1.0 - a**2.0 )**(1.0/3.0)
    z1 = z1 * ( (1.0+a)**(1.0/3.0)+(1.0-a)**(1.0/3.0))+1.0
  end if 
	  
  z2 = sqrt( 3.0 * a**2.0 + z1**2.0 )
  if(a.ge.0.0)then
    dISCO = 3.0 + z2 - sqrt( (3.0-z1) * (3.0 + z1 + 2.0*z2) )
  else
    dISCO = 3.0 + z2 + sqrt( (3.0-z1) * (3.0 + z1 + 2.0*z2) )
  end if
  return
end function dISCO
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
function mybbody(kT,E,dE)
! Blackbody function in terms of number of photons with energy
! between E-dE/2 and E+dE/2. i.e. This is photar!
! Function is normalized such that the integrated energy flux is 1
! i.e. sum E * photar(E) = 1
  implicit none
  real mybbody,E,dE,kT
  real pi,fac,f
  pi   = acos(-1.0)
  fac  = E/kT
  if(fac .lt. 1e-3)then
    f = E * kT   !Using a Taylor expansion
  else
    if (fac .lt. 70.) then
      f = E**2 / ( exp(fac) - 1.0 ) 
    else
      f = E**2 * exp(-fac)
    end  if
  end if
  mybbody = f * (15.0/pi**4) / kT**4 * dE
  return
end function mybbody
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
function fcol(T)
! fcol = colour-temperature correction
! T    = true temperature in keV
! Done et al. 2012 model
  implicit none
  double precision fcol,T
  if( T .lt. 2.585d-3 )then
      fcol = 1.d0
  else if( T .lt. 8.617e-3 )then
      fcol = ( T / 2.585d-3 )**0.833
  else
      fcol = ( 72.d0 / T )**(1.d0/9.d0)
  end if
  return
end function fcol
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
function fcol_DEA(T, r, a, m, mdot)
! fcol = colour-temperature correction
! T    = true temperature in keV
! Davis and El-Abd model
  implicit none
  double precision fcol_dea, T, r, a, AAAA
  double precision zeta, wrp, k_es, kb, alp, O_sb, mu, mp , m0
  double precision u0, omega_prime, G, rg, ur, m, sigma, pi, mdot, Q, Rz, up, u_0

  pi  = acos(-1.d0)

  G = 6.67e-11
  rg = G * m * 2e30 / 3e8**2.0
  O_sb = 5.67e-8 

  alp = 0.1 
  k_es = 0.034
  mu = 0.67 
  mp = 1.67e-27 
  kb = 1.38e-23 

  u0 = (1 + a* r**(-1.5))/(1 - 3./r + 2.*a*r**(-1.5))**0.5 
  omega_prime = 3./2. * (G*m*2e30/rg**5)**0.5 * r**(-2.5) / ((1 + a*r**(-1.5))**2.0)

  AAAA = 3 * k_es * kb**4.0 * alp ** 4.0 / (16.0 * O_sb * (mu * mp)**4.0 )
  
  zeta = 2 * O_sb * T**4.0 * r * rg / (u0**2.0 * omega_prime)
  Wrp = (A * r**2.0 * rg**2.0 * u0**3.0 * zeta ** 2.0 * omega_prime)**0.2

  ur =  (Wrp * r * rg) / (zeta * u0)

  sigma = mdot * 1e15 / (2 * pi * r * rg * ur)
  m0 = 0.5 * 0.1 * sigma

  up =  (r)**0.5 * (1 + a**2.0/r**2.0 - 2.0*a*r**(-1.5))/(1 - 3./r + 2.*a*r**(-1.5))**0.5
  u_0 =  - (1 - 2./r+ a * r**(-1.5))/(1 - 3./r + 2.*a*r**(-1.5))**0.5

  Rz = (up**2.0 + a**2.0 * (1 - u_0**2.0))/(r)

  Q = G * m * 2e30 / (r * rg)**3.0 * Rz

  fcol_DEA = 1.74 + 1.06 * (log10(T) - 7.0) - 0.14 * (log10(Q) - 7.0) - 0.07 * (log10(m0) - 5.0)
  
  return
end function fcol_DEA
!-----------------------------------------------------------------------
  

!-----------------------------------------------------------------------
function rfunc(a,mu0)
! Sets minimum rn to use for impact parameter grid depending on mu0
! This is just an analytic function based on empirical calculations:
! I simply set a=0.998, went through the full range of mu0, and then
! calculated the lowest rn value for which there was a disk crossing.
! The function used here makes sure the calculated rnmin is always
! slightly lower than the one required.
!
! Lack of event horizon makes r_min much smaller for a > 1.  Hard coded
! Simple lower bound for now. 
!
  implicit none
  double precision rfunc,mu0,a
  if (a .gt. 1.0) then 
    rfunc = 0.01   
  else if( a .gt. 0.8 )then
    rfunc = 1.5d0 + 0.5d0 * mu0**5.5d0
    rfunc = min( rfunc , -0.1d0 + 5.6d0*mu0 )
    rfunc = max( 0.1d0 , rfunc )
  else
    rfunc = 3.0d0 + 0.5d0 * mu0**5.5d0
    rfunc = min( rfunc , -0.2d0 + 10.0d0*mu0 )
    rfunc = max( 0.1d0 , rfunc )
  end if
  rfunc = rfunc-0.5
  end function rfunc
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine getrgrid(rnmin,rnmax,mueff,nro,nphi,rn,domega)
! Calculates an r-grid that will be used to define impact parameters
  implicit none
  integer nro,nphi,i
  double precision rnmin,rnmax,mueff,rn(nro),domega(nro)
  double precision rar(0:nro),dlogr,logr,pi
  pi     = acos(-1.d0)
  rar(0) = rnmin
  dlogr  = log10( rnmax/rnmin ) / dble(nro)
  do i = 1,NRO
    logr = log10(rnmin) + dble(i) * dlogr
    rar(i)    = 10.d0**logr
    rn(i)     = 0.5 * ( rar(i) + rar(i-1) )
    domega(i) = rn(i) * ( rar(i) - rar(i-1) ) * mueff * 2.d0 * pi / dble(nphi)
  end do
  return
end subroutine getrgrid
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
      subroutine drandphithick(alpha,beta,cosi,costheta,r,phi)
!
! A disk with an arbitrary thickness
! The angle between the normal to the midplane and the disk surface is theta
! The inclination angle is i
      implicit none
      double precision alpha,beta,cosi,sini,r,phi
      double precision pi,costheta,sintheta,x,a,b,c,det
      double precision mu,sinphi
!      double precision muplus,muminus,ra,rb,rab,xplus1,xminus1,xplus2,xminus2
      pi = acos(-1.d0)
      sintheta = sqrt( 1.d0 - costheta**2 )
      sini     = sqrt( 1.d0 - cosi**2 )
      x        = alpha / beta
      if( abs(alpha) .lt. abs(tiny(alpha)) .and. abs(beta) .lt. abs(tiny(beta))  )then
        mu = 0.d0
        r  = 0.d0
      else if( abs(beta) .lt. abs(tiny(beta)) )then
        mu     = sini*costheta/(cosi*sintheta)
        sinphi = sign( 1.d0 , alpha ) * sqrt( 1.d0 - mu**2 )
        r      = alpha / ( sintheta * sinphi )
      else if( abs(alpha) .lt. abs(tiny(alpha)) )then
        mu     = 1.d0
        sinphi = 0.d0
        r      = beta / ( sini*costheta - cosi*sintheta )
      else
        a      = sintheta**2 + x**2*cosi**2*sintheta**2
        b      = -2*x**2*sini*cosi*sintheta*costheta
        c      = x**2*sini**2*costheta**2-sintheta**2
        det    = b**2 - 4.d0 * a * c
        if( det .lt. 0.d0 ) write(*,*)"determinant <0!!!"
        if( beta .gt. 0.d0 )then
          mu     = ( -b + sqrt( det ) ) / ( 2.d0 * a )
        else
          mu     = ( -b - sqrt( det ) ) / ( 2.d0 * a )
        end if
        sinphi = sign( 1.d0 , alpha ) * sqrt( 1.d0 - mu**2 )
        r      = alpha / ( sintheta * sinphi )
      end if
      phi = atan2( sinphi , mu )
      return
      end subroutine drandphithick
!-----------------------------------------------------------------------

