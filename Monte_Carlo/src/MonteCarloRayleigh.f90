!>\author Meng Gao
!>\date 04/07/2010
!>\details 3D Monte Carlo code for Rayleigh Scattering in a plane parallel system

Program MC
  use m_Public
  implicit none
  Integer :: Kphoton,Nphoton, ncol
  Integer :: IX !seed for ran2
  Real(dp) :: Weight, Weight_CutOff, AfterEscape, Ref, Tran, F0
  Real(dp) :: t1,t2,ttol
  Real(dp) :: xi, tau,tau1, tau_Max, albedo
  Real(dp) :: theta0, mu0,phi0, muSample, phiSample
  Real(dp) :: RN1, RN2, Ptest,PRayleigh
  Type(Direction) :: ni
  Type(Position) :: r
  Real(dp) :: theta_d,phi_d, dtheta, dOmega
  Integer :: Ntheta,i,j
  Real(dp), Dimension(:), Allocatable :: theta_Set, Radiance
  Real(dp) :: Estimation
  Type(Direction), Dimension(:),Allocatable :: nd
  Real(dp) :: mu !ni.nd for estimation
  Integer :: NphotonPrint 

  !the number of photon used to print trajectory
  NphotonPrint=1

  Call cpu_time(t1)
  !---input parameters:--------------------------------------------------------------
  ! total optical depth
  !  tau_Max=1.0_dp
  !  tau_Max=0.01_dp
  tau_Max=10.0_dp
  ! incident direction
  theta0=180_dp*Pi/180_dp
  phi0=0_dp*Pi/180_dp
  ! incident irradiance
  !  F0=1367
  F0=1.0
  ! total photon numbers
  Nphoton=10000
  ! albedo
  albedo=1.0_dp !check energy conservation
  !-- detector directions--in the system frame
  phi_d=0.0_dp*Pi/180_dp !for principle plane, sinc incident direction phi_i=0.0_dp
  Ntheta=180 !number of angle for upward refl

  !--initialization-------------------------------
  Allocate(theta_Set(Ntheta), nd(Ntheta), Radiance(Ntheta))
  Radiance=0.0_dp
  Ref=0.0_dp
  Tran=0.0_dp  
  ! photon weight cutoff
  Weight_CutOff=0.00000001_dp !10^(-8)
  ! random number seed
  IX=-439
  !we assign the angle for reflection and transmission seperately, take out 90
  Do i=1, 90
     theta_d=(i-1)*Pi/180_dp
     theta_Set(i)=theta_d
     nd(i)%a=Sin(theta_d)*Cos(phi_d)
     nd(i)%b=Sin(theta_d)*Sin(phi_d)
     nd(i)%c=Cos(theta_d)
  End Do
  Do i=91, 180
     theta_d=i*Pi/180_dp
     theta_Set(i)=theta_d
     nd(i)%a=Sin(theta_d)*Cos(phi_d)
     nd(i)%b=Sin(theta_d)*Sin(phi_d)
     nd(i)%c=Cos(theta_d)
  End Do
  !-------------------------------------------------------------------------
  !--open files----------
  open(unit=1200,file='trajectory.dat',status='unknown')
  open(unit=1201,file='radiance.dat',status='unknown')
  !-----------------------
  Print *, "---starting---"
  !new photon
  New_Photon : DO Kphoton=1,Nphoton  ! generate new photon
     !initial photon direction
     mu0=Cos(theta0)
     ni%a=Sin(theta0)*Cos(phi0)
     ni%b=Sin(theta0)*Sin(phi0)
     ni%c=mu0
     ni%phi=phi0
     !ni%theta=theta0
     ! initialize photon weight
     Weight=1.0_dp
     ! collision numbers
     ncol=0
     ! position in the units of tau
     r%x=0.0_dp
     r%y=0.0_dp
     r%z=tau_Max
     !------------
     If(Kphoton>100) Then
        If (Mod(Kphoton,Nphoton/100).EQ.0) &
             Write(*,*) Int(DBLE(Kphoton)/Nphoton*100), '/100' ! output percentage
     End If

     If(Kphoton<NphotonPrint) Then
        Write(1200,'(2(2X,I5),(2X, F10.8),3(2X,F15.8))') kphoton, ncol, Weight, r%x,r%y,r%z
     End If

     Next_Scattering : DO WHILE(Weight>Weight_CutOff)
        !---Hitting, sample tau from current position to the boundary-------
        ! tau is the total length, photon will traveled after this collision
        ! tau1 is the toal distance from current photon position along the current direction 
        ! to the boundary
        xi=RAN2(IX) !random number [0,1]       
        If(Abs(ni%c)<tiny1) Then
           tau=-Log(xi)
        Else If(ni%c>0)Then
           tau1=(tau_Max-r%z)/ni%c
           AfterEscape=1.0_DP-EXP(-tau1) ! force move some part of photon package out
           Ref=Ref+Weight*Exp(-tau1) 
           weight=weight*AfterEscape
           tau=-Log(1.0_dp-AfterEscape*xi)
        Else If(ni%c<0) Then
           tau1=(0.0_dp-r%z)/ni%c
           AfterEscape=1.0_DP-EXP(-tau1) ! force move some part of photon package out
           Tran=Tran+Weight*Exp(-tau1)
           weight=weight*AfterEscape
           tau=-Log(1.0_dp-AfterEscape*xi)

        End If
        Ncol=Ncol+1
        Weight=Weight*albedo
        !----------------------------------------------------------------------
        !---- update the position of photon in term of optical depth---------
        r%x=r%x+tau*ni%a
        r%y=r%y+tau*ni%b
        r%z=r%z+tau*ni%c
        !----check trajectory------------------------------------------------------------
        If(kphoton<=NphotonPrint) Then
           Write(1200,'(2(2X,I5),(2X, F10.8),3(2X,F15.8))') kphoton, ncol, Weight, r%x,r%y,r%z
        End If
        !--ESTIMATION--------
        Do i=1, Ntheta
           mu=ni%a*nd(i)%a+ni%b*nd(i)%b+ni%c*nd(i)%c         ! n1*n2=cos theta12
           If(Abs(mu)>1.0_DP) mu=mu/Abs(mu)
           PRayleigh=3.0_dp/(16.0_dp*Pi)*(1.0_dp+mu**2)
           If(ABS(nd(i)%c)<tiny1) then
              Print *, "detector angle is too small"
              STOP
           Else If(nd(i)%c>0)Then
              tau1=(tau_Max-r%z)/nd(i)%c
           Else If(nd(i)%c<0) Then
              tau1=(0-r%z)/nd(i)%c
           End If

           Estimation=Weight*PRayleigh*Exp(-tau1)/Abs(Cos(theta_Set(i)))
           Radiance(i)=Radiance(i)+Estimation
        End Do

        !-------------------------------------------------
        !        Call FindDrc(ni)

        !----SCATTERING----------------------------
        !sample scattering angles
100     RN1=RAN2(IX)
        RN2=RAN2(IX)
        muSample=2.0_DP*RN2-1.0_DP ! 
        PTEST=0.5_DP*(1.0_DP+muSample**2) ! make the maximum=1
        IF(PTEST<RN1) GOTO 100

        !     PRayleigh=0.5_DP*(1.0_DP+muSample**2)
        phiSample=RAN2(IX)*2*PI
        Call UpdateDrcRan(ni, muSample, phiSample)

     End Do Next_Scattering

  End Do New_Photon

  !--EXPORT--------------
  Ref=Ref/Nphoton
  Tran=Tran/Nphoton
  !get radiance---------
  Do i=1, Ntheta
     Radiance(i)=F0*Abs(Cos(theta0))*Radiance(i)/Nphoton !/dOmega*
     Write(1201,'((2X, F6.2),(2X,F12.6))') theta_Set(i)*180_dp/Pi, Radiance(i)
  End Do

  !timing------------
  Call cpu_time(t2)
  ttol=t2-t1
  Print *, "running time Total",ttol
  Print *, "Reflection", Ref, "Transmission", Tran, "Loss", 1-(Ref+Tran)

  Close(1200)
  Close(1201)
  Close(1202)
  Deallocate(Radiance, theta_Set)

End Program MC
