Module m_Public
    Integer, Parameter :: dp = KIND(1.0D0)
    Real(dp), Parameter :: PI=3.141592653589793238462643383279502884197_dp
    Real(dp), Parameter :: Tiny1=0.0000000001_dp !10^(-10)
    TYPE Direction
       Real(dp) :: a,b,c, theta, phi
    END type Direction

    TYPE Position
       Real(dp) :: x,y,z
    END type Position

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

  !> \brief Update photon direction after Sampling scattering angle \f$ \theta_s \f$ for random orientation particles
  !> \details algorithm: 
  !> \author PengWang Zhai, Meng Gao
  !> \date 01/14/2011
  !> \todo put in the new form using rotation of frame
  !>\warning using the reference frame same as the limit of the fixed orientation,
  !> be careful when comparing with the original old code, it may be using different frame if not modified yet.
  SUBROUTINE UpdateDrcRan(ni, MuSample, PhiSample)
    implicit none
    TYPE(Direction), INTENT(INOUT) :: ni
    REAL(DP), INTENT(IN) :: MuSample, PhiSample
    REAL(DP) :: CANG, SCPHI, SANG, SA,SB,SC,SPHI,COSPHI, SINPHI, SINPC,&
         RK, SCS, SS, U2,V2
    CANG=MuSample
    !-- sampling phi------                                                                             
                                        
    SCPHI=PhiSample
    !PhiSample=SCPHI
    !---update a,b,c------------
    SANG=SQRT(ABS(1.0_DP-CANG**2))!mu=cang obtained by sampling mu                           
    !------before collision----                                
    SA=ni%A
    SB=ni%B
    SC=ni%C
    SPHI=ni%PHI
    !-----after collision-------                                                                       
    COSPHI=COS(SCPHI)
    SINPHI=SIN(SCPHI)
    SINPC=SQRT(ABS(1.0_DP-SC**2))

    IF(ABS(SINPC)<TINY1) then !sinpc=0                                                                 
       ni%A=SANG*COSPHI
       ni%B=SANG*SINPHI
       ni%C=CANG*ni%C
    ELSE
       RK=1.0_DP/SINPC
       SCS=SANG*SINPHI*ni%C
       SS=-SANG*COSPHI

       U2=(SCS*ni%A-SS*ni%B)*RK+CANG*ni%A
       V2=(SCS*ni%B+SS*ni%A)*RK+CANG*ni%B
       ni%C=ni%C*CANG-SINPC*SANG*SINPHI
       ni%A=U2
       ni%B=V2
    End IF
  END SUBROUTINE UpdateDrcRan

  SUBROUTINE FindDrc(ntemp)
    implicit none
    TYPE(Direction), INTENT(INOUT) :: ntemp
        ntemp%theta=ACos(ntemp%c)
        If(ntemp%theta>2*Pi) Then
           ntemp%theta=ntemp%theta-2*Pi
        Else If(ntemp%theta<0) Then
           ntemp%theta=ntemp%theta+2*Pi
        End If
    ! get rid of singularity !no matter what z is
    IF(ABS(ntemp%A)<TINY1 .AND. ABS(ntemp%B)<TINY1) THEN
       ntemp%phi=0.0_DP ! (0,0,+-1) along +-z directioni, phi=0
       RETURN
    ENDIF
    IF(ABS(ntemp%A)<TINY1)THEN
       IF(ntemp%B>0.0D0) THEN
          ntemp%phi=PI/2.0D0 ! (0,1,1) in the y-z plan, phi=pi/2
       ELSE IF(ntemp%B<0.0D0)THEN
          ntemp%phi=3.0D0*PI/2.0D0 !( 0, -1,1)
       ENDIF
       RETURN
    ENDIF
    ! using trigonometric function
    IF (ntemp%A<0.0_DP) THEN
       ntemp%phi=PI+ATAN(ntemp%B/ntemp%A)
    ELSE IF (ntemp%B>=0.0_DP) THEN
       ntemp%phi=ATAN(ntemp%B/ntemp%A)
    ELSE IF (ntemp%B<0.0_DP) THEN
       ntemp%phi=2*PI+ATAN(ntemp%B/ntemp%A)
    ENDIF
    IF(ABS(ntemp%phi-2.0_dp*Pi)<TINY1)ntemp%phi=0.0_DP ! large than 2pi get back to 0
  END SUBROUTINE FindDrc

END MODULE m_Public
