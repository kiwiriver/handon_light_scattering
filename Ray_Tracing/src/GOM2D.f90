!>\author Meng Gao
!>\date 05/03/2011
!>\details Geometric Optics Methods for 2D sphere

Program GOM2D
  use m_Public
  use m_Sphere
  implicit none
  Real(dp), Dimension(181):: P11
  Real(dp) :: l,h
  Real(dp),Dimension(2) :: n,f ! Drc n and f are for surface
  Real(dp),Dimension(2) :: ni0,pti0, ni, pti !incident direction, and one point on the ray
  Real(dp),Dimension(2) :: nit, nr, nt,ns ! nit temporary save ni
  !reflected and refracted light nr, nt, scattered light ns
  Real(dp),Dimension(2) :: pm, pm1,pm2 !ray intersect a surface at pm
  Integer :: IX, I,ti
  Real(dp) :: norm, thetai, x0, y0
  Real(dp) :: t1,t2,ttol
  Logical :: l_miss
  Real(dp) :: nis ! ni.ns
  !-----------------
  Real(dp) :: R ! radius of sphere
  !----initialize variables--------
  Call cpu_time(t1)

  Nray=1000000 ! the number of tested rays
  !--------------------
  !refractive index: real mr0, imaginary mi0
  !  mr0=1.33_dp 
  mr0=1.31
  mi0=0.0_dp
  !  mi0=0.02_dp 

  !----Do Not Modify After This Line----
  IX=233
  !  IX=178 !seed for ran2
  l_miss=.True.
  P11=0.0_dp 
  mr=mr0
  mi=mi0
  !--open files-----------
  open(unit=99,file='log.dat',status='unknown')
  open(unit=100,file='path.dat',status='unknown')
  open(unit=102,file='p11.dat',status='unknown')

  krayhit=0 ! record the number of rays hitting the particle
  New_Ray: Do kray=1, Nray
     If(Kray>100) Then
        If (Mod(Kray,Nray/100).EQ.0) &
             Write(*,*) Int(DBLE(Kray)/Nray*100), '/100' ! output percentage
     End If

     thetai=0.0_dp
     R=1.0_dp
     x0=-2.0_dp*R
     y0=2.0_dp*R*RAN2(IX)-R
     
     !-----------------
     ncol=0
     weighti=0.5_dp
     weights=0.5_dp
     !--------------
     ni0(1)=Cos(thetai)
     ni0(2)=Sin(thetai)
     pti0(1)=x0
     pti0(2)=y0
     ni=ni0
     pti=pti0

     ! whether the ray can hit the sphere
     Call HitNoSphere(ni, pti, R, l_miss)

     If(l_miss) Then
        Cycle New_Ray
     End If

     krayhit=krayhit+1
     !--initially the ray is outside-------------------------
     l_inside=.False. !outside
     Tracing_Next_Ray: DO WHILE(weighti(1)+weighti(2)>Weight_CutOff)

        Call HitWhichSphere(l_inside, ni, pti, R, n, f, pm)
        
        Call ReachBoundary(ni, pti, n,f, pm, ns)
        nis=ni0(1)*ns(1)+ni0(2)*ns(2)
        If(Abs(nis)>1.0_dp) Then 
           nis=nis/Abs(nis)
        End If
        ti=ANINT(ACos(nis)*180/pi)+1
        P11(ti)=P11(ti)+weights(1)+weights(2)
     End Do Tracing_Next_Ray

  End Do New_Ray
  !--Export------------------
  norm=0.0_dp
  Do I=1,180
     norm=norm+(P11(I)+P11(I+1))/2*(Cos(Real(I-1)/180.0_dp*Pi)-Cos(Real(I)/180.0_dp*Pi))
  End Do

  P11=2*P11/norm
  print *,"normalization factor", norm

  Do I=1,181
     Write(102,'((2X, I5),(2X,F12.6))') I-1, P11(I)
  End Do

  Close(100)
  Close(101)

  !---timing------------
  Call cpu_time(t2)
  ttol=t2-t1

  Write(99,'(2X, "total number of rays used:", 2X, I10)') Nray
  Write(99,'(2X, "total number of rays hitting the particle:", 2X, I10)') krayhit
  Write(99,'(2X, "refractive index real:", 2X, F10.5, 2X, "imaginary", F10.5)') mr0, mi0
  Write(99,'(2X, "running time in sec:", 2X, F10.5)') ttol
  Close(99)
  !--------------
  Print *, "total ray numbers used", Nray
  Print *, "total ray numbers hitted", krayhit
  Print *, "total running time (sec)",ttol
  !--------------------
End Program GOM2D
