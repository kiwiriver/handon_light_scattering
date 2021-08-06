!>\author Meng Gao
!>\date 01/11/2011
!>\details public functions and subroutines for polygons

MODULE m_Public
  INTEGER, PARAMETER :: DP = KIND(1.0D0)
  REAL(DP), PARAMETER :: PI=3.141592653589793238462643383279502884197_dp
  REAL(DP), PARAMETER :: tiny1=0.0000000001_dp, weight_cutoff=0.00000001_dp  
  !tiny1=10**(-10),weight_cutoff=10**(-8)
  Real(dp), Dimension(:,:), Allocatable::shape
  Integer :: Nvertex
  Real(dp),Dimension(4,2) :: TheseSurfaces
  Real(dp) :: mr0, mi0 !refractive indexxs
  Real(dp) :: mr, mi !adjusted refractive indexxs
  Logical :: l_inside
  Real(dp), Dimension(2) :: weighti, weights 
  Integer :: Nray, kray, krayhit, Nrayhit
  Integer :: ncol
CONTAINS

  !> \brief random number genetator
  !> \author Numerical Recipe, Chapter 7 p272
  FUNCTION RAN2(idum)
    INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
    REAL RAN2,AM,EPS,RNMX
    PARAMETER(IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
         IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,            &
         IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
    INTEGER idum2,j,k,iv(NTAB),iy
    SAVE iv,iy,idum2
    DATA idum2/123456789/,iv/NTAB*0/,iy/0/

    IF(idum .LE. 0) THEN
       idum = MAX(-idum,1)
       idum2 = idum
       DO j = NTAB+8,1,-1
          k = idum/IQ1
          idum = IA1*(idum-k*IQ1) - k*IR1
          IF(idum .LT. 0) idum = idum + IM1
          IF(j .LE. NTAB) iv(j) = idum
       END DO
       iy = iv(1)
    END IF

    k = idum/IQ1
    idum = IA1*(idum-k*IQ1) - k*IR1
    IF(idum .LT. 0) idum = idum + IM1
    k = idum2/IQ2
    idum2 = IA2*(idum2-k*IQ2) - k*IR2
    IF(idum2 .LT. 0) idum2 = idum2 + IM2
    j = 1 + iy/NDIV
    iy = iv(j) - idum2
    iv(j) = idum
    IF(iy .LT. 1) iy = iy + IMM1
    RAN2 = MIN(AM*iy,RNMX)
    RETURN
  END FUNCTION RAN2
  !-------------------------------------------

  Subroutine AdjustIndex(l_inside, ni, n)
    Real(dp), Dimension(2), Intent(IN) :: ni, n
    Logical, Intent(IN) :: l_inside
    Real(dp) :: cosi,cost ! save ni.n
    !--full formulae-----------
    !    mr=Sqrt(mr0**2-mi0**2+(1-nin**2)+Sqrt((mr0**2-mi0**2-1+nin**2)**2+4*mr0**2*mi0**2))/Sqrt(2.0_dp)
    !    mi=mr0*mi0/mr
    !--approximation when mi0 is small----
    mr=mr0
    mi=mi0
  End Subroutine AdjustIndex

  Subroutine reflection(ni,n,f,nr)
    Real(dp), Dimension(2), Intent(IN) :: ni, n, f 
    Real(dp), Dimension(2), Intent(Out) :: nr
    Real(dp) :: nin, norm ! nin:save ni.n
    nin=n(1)*ni(1)+n(2)*ni(2)
    nr(1)=ni(1)-2*nin*n(1)
    nr(2)=ni(2)-2*nin*n(2)
    norm=nr(1)**2+nr(2)**2
    nr(1)=nr(1)/norm
    nr(2)=nr(2)/norm
  End Subroutine reflection

  Subroutine refraction(ni,n,f, mrtemp, nt)
    Real(dp), Dimension(2), Intent(IN) :: ni, n, f 
    Real(dp) :: mrtemp !relative rafractive index, set the material of incident direction 1
    Real(dp), Dimension(2), Intent(Out) :: nt
    Real(dp) :: nin, norm, factor
    nin=n(1)*ni(1)+n(2)*ni(2)
    factor=mrtemp**2-1+nin**2
    If(factor>=0) then
       nt(1)=1/mrtemp*(ni(1)-nin*n(1)- &
            Sqrt(factor)*n(1))
       nt(2)=1/mrtemp*(ni(2)-nin*n(2)- &
            Sqrt(factor)*n(2))
       norm=nt(1)**2+nt(2)**2
       nt(1)=nt(1)/norm
       nt(2)=nt(2)/norm
    Else !there is no refraction light, we set nt=0.0
       nt(1)=0.0_dp
       nt(2)=0.0_dp
    End If
  End Subroutine refraction

  !> \details Fresnel formulas using adjusted refractive index
  Subroutine Fresnel(mrtemp,ni,n,Rl,Rr)
    Real(dp), Dimension(2), Intent(IN) :: ni, n
    Real(dp) :: mrtemp !relative rafractive index, set the material of incident direction 1
    Real(dp), Intent(out) :: Rl, Rr
    Real(dp) :: factor

    cosi=Abs(n(1)*ni(1)+n(2)*ni(2))
    factor=mrtemp**2-1+cosi**2
    If(factor>=0) then
       cost=Sqrt(factor)/mrtemp
       ! reflection for electric field
       Rl=(mrtemp*cosi-cost)/(mrtemp*cosi+cost)
       Rr=(cosi-mrtemp*cost)/(cosi+mrtemp*cost)
    Else !total reflection case, no trasmission
       Rl=1.0_dp
       Rr=1.0_dp
    End If
    !reflection for energy density
    Rl=Rl**2
    Rr=Rr**2
  End Subroutine Fresnel
  !------------------------------------------------------
  Subroutine PerpDrc(ni,nj)
    Real(dp), Dimension(2), Intent(IN) :: ni
    Real(dp), Dimension(2), Intent(OUT) :: nj
    Real(dp) :: norm !normalization factor for ni and nj
    !there are two opposite direction perpendicular to ni
    !here we choose the one nj x ni=nz
    nj(1)=ni(2)
    nj(2)=-ni(1)
    !normalize
    norm=ni(1)**2+ni(2)**2
    If(norm>1.0_dp) Then
       nj=nj/norm
    End If
  End Subroutine PerpDrc
  !---------------------------------------------------
  Subroutine ReachBoundary(ni, pti, n, f, pm, ns)
    Real(dp), Dimension(2), Intent(INOUT) :: ni, pti
    Real(dp), Dimension(2), Intent(IN) :: n, f, pm
    Real(dp), Dimension(2), Intent(Out) :: ns
    Real(dp), Dimension(2) :: wi, ws

    !for the transmitted ray, and scattered outside ray
    !(parallel, perpendicular)=(l,r)
    Real(dp), Dimension(2) :: nt, nr
    Real(dp) :: mrtemp
    Real(dp) :: Rl,Rr
    Real(dp) :: wabs, dis, mitemp
    Integer :: NrayPrint

    NrayPrint=1 ! the number of ray path printed
    wabs=1.0_dp
    If(l_inside) Then
       dis=Sqrt((pm(1)-pti(1))**2+(pm(2)-pti(2))**2)
       wabs=Exp(-mi*dis*300.0_dp) !size parameter pi d=lambda !check    
       !mi used here is updated from last time
       !       print *, dis, mi, wabs
    Else
       wabs=1.0_dp
    End If
    !------------------------------------
    pti=pm !for next ray 

    !--determine relative refractive index-----   
    If(.Not. l_inside) Then
       mrtemp=mr0
    Else If(l_inside) Then
       mrtemp=1/mr0
    Else
       Print *, "error for mrt"
    End If
    !------------------------------------------
    ! we adjust refractive index before each refraction, and reflection
    ! and fresnel
    Call AdjustIndex(l_inside, ni, n) 
    !----------------------------
    Call refraction(ni, n, f, mrtemp, nt)
    Call reflection(ni, n, f, nr)
    Call Fresnel(mrtemp, ni, n, Rl, Rr)

    If(.NOT. l_inside) Then
       ni=nt
       ns=nr
       wi(1)=1-Rl
       wi(2)=1-Rr
       ws(1)=Rl
       ws(2)=Rr
    Else If(l_inside) Then
       ni=nr
       ns=nt
       wi(1)=Rl
       wi(2)=Rr
       ws(1)=1.0_dp-Rl
       ws(2)=1.0_dp-Rr
    Else 
       Print *, "error l_inside"
    End If

    weights(1)=weighti(1)*ws(1)
    weights(2)=weighti(2)*ws(2)
    weighti(1)=weighti(1)*wi(1)
    weighti(2)=weighti(2)*wi(2)
    !------absorption-------------
    weighti=weighti*wabs
    !----------------------------
    ncol=ncol+1
    !--path and weight and direction------------------
    If(krayhit<=NrayPrint) Then
       Write(100,'((2x,I5), 2(2X, F6.2),4(2X, F10.6),4(2X, F10.8))') & 
            ncol, pm(1),pm(2),ni(1),ni(2),ns(1),ns(2),&
            weighti(1)+weighti(2),weights(1)+weights(2)
    End If
    l_inside=.True. 
  End Subroutine ReachBoundary

END MODULE m_Public
