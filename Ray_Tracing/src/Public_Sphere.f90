!>\author Meng Gao
!>\date 01/11/2011
!>\details public functions and subroutines for sphere
MODULE m_Sphere
  use m_Public

CONTAINS
  Subroutine HitNoSphere(ni, pti, R, l_miss)
    Logical, Intent(Out) :: l_miss
    Real(dp), Dimension(2), Intent(In) :: ni, pti
    Real(dp) :: R !radius of the sphere located at (0,0)
    Real(dp) :: factor
    factor=ni(1)**2*(1-(ni(2)*pti(1)-ni(1)*pti(2))**2)
    l_miss=.False.
    If(factor<=0) Then
       l_miss=.True. 
       !when factor=0, tangent to the boundary
       !when factor<0, no intersecting points
       !when factor>0 two solutions
    End If
  End Subroutine HitNoSphere

  Subroutine WhereTheyMeetSphere(ni, pti, R, pm1, pm2)
    Real(dp), Dimension(2), Intent(In) :: ni, pti
    Real(dp), Intent(In) :: R !radius of the sphere located at (0,0)
    Real(dp), Dimension(2), Intent(Out) :: pm1, pm2
    Real(dp) :: factor0,factor
    factor0=R**2-pm1(1)**2
    If(Abs(ni(1))<Tiny1) then
       If(factor0<0) Then
          Print *, "error to find y on the sphere"
          Return
       Else
          pm1(1)=pti(1)
          pm2(1)=pti(1)
          pm1(2)=factor0**(1/2)
          pm2(2)=-factor0**(1/2)
       End If
       Return
    End If

    factor=ni(1)**2*(1-(ni(2)*pti(1)-ni(1)*pti(2))**2)
    If(factor<=0) Then
       Print *, "error, fail to select the right ray hit on sphere" 
       Return
    End If
    pm1(1)=(ni(2)*(ni(2)*pti(1)-ni(1)*pti(2))-Sqrt(factor))
    pm2(1)=(ni(2)*(ni(2)*pti(1)-ni(1)*pti(2))+Sqrt(factor))
    pm1(2)=(-ni(2)*pti(1)+ni(1)*pti(2)+ni(2)*pm1(1))/(ni(1))
    pm2(2)=(-ni(2)*pti(1)+ni(1)*pti(2)+ni(2)*pm2(1))/(ni(1))

  End Subroutine WhereTheyMeetSphere

  Subroutine DetermineDrcSphere(pm,ni,n,f)
    Real(dp), Dimension(2), Intent(In) :: pm,ni
    Real(dp), Dimension(2), Intent(Out) :: n,f
    Real(dp) :: norm1
    n(1)=pm(1)
    n(2)=pm(2)
    norm1=Sqrt(n(1)**2+n(2)**2)
    n(1)=n(1)/norm1
    n(2)=n(2)/norm1
    If(n(1)*ni(1)+n(2)*ni(2)>0) Then
       n=-n
    End If !ortherwise it's fine
    Call PerpDrc(n,f) !get new f
  End Subroutine DetermineDrcSphere

  Subroutine HitWhichSphere(l_inside, ni, pti, R, n, f, pm)
    Logical, Intent(IN):: l_inside
    Real(dp), Intent(In) :: R ! radius
    Real(dp), Dimension(2), Intent(IN):: ni, pti
    Real(dp), Dimension(2), Intent(Out):: n,f
    Real(dp), Dimension(2), Intent(Out):: pm !the intersecting point
    Real(dp), Dimension(2) :: pm1, pm2 !the intersecting points
    Integer :: flag
    Integer :: I, J
    Real(dp) :: d, dtemp

    Call WhereTheyMeetSphere(ni, pti, R, pm1, pm2)
    If((pm1(1)-pm2(1))*ni(1)+(pm1(2)-pm2(2))*ni(2)>0) Then
       If(l_inside) Then
          pm=pm1
       Else
          pm=pm2
       End If
    Else
       If(l_inside) Then
          pm=pm2
       Else
          pm=pm1
       End If
    End If
    Call DetermineDrcSphere(pm, ni, n, f)
  End Subroutine HitWhichSphere

END MODULE m_Sphere
