      program ttt
       implicit real*8 (a-h,o-z)
C compile ifort ua.f -r8  -I./fftw-3.3.4/include -L./fftw-3.3.4/lib -lfftw3
       real*8 u(129),uu(129),x(129),wait(129),v(129),tt(129)
       real*8 us(129),vs(129),ws(129)
       real*8 w(129)

       km = 129
       CPI = DBLE(4.)*DATAN(DBLE(1))
       do i=1,km
        wait(i)=DBLE(0)
        x(i)=DCOS( (DBLE(i-1)*CPI)/DBLE(km-1) )
       enddo
       do i=1,km,2
         r=i-1
         wait(i)=-DBLE(1)/DBLE(r*r-1)
       enddo
       wait(1)=DBLE(1)/DBLE(2)

       do k=1,km
         us(k)=DBLE(0)
         vs(k)=DBLE(0)
         ws(k)=DBLE(0)
       enddo

C PROCESS MEAN STREAM VELOCITY FIRST
       do i=1,200000
         read(2111,'(512E24.16)') (u(K),K=1,KM)
C make symmetric
        do j=1,km/2+1
          tmp=(u(km+1-j)+u(j))/DBLE(2)
          u(j)=tmp
          u(km+1-j)=tmp
        enddo
C 1st derivative in spectral space
         call ctran(u)
         call chebdz(u,uu,km)
         call citran(uu)
         do k=1,km
           u(k)=uu(k)/DBLE(180)
         enddo
C integral of ^2 over wall-normal direction
         do k=1,km
           us(k)=us(k)+(u(k)-us(k))/DBLE(i)
           tt(k)=us(k)
         enddo
         call ctran(tt)
         su3=DBLE(0)
         do j=1,km
           su3=su3+wait(j)*tt(j)
         enddo
         sumu=su3
         do k=1,km
           tt(k)=(us(k)-sumu)**2
         enddo
         call ctran(tt)
         su3=DBLE(0)
         do j=1,km
           su3=su3+wait(j)*tt(j)
         enddo
         rmsu=su3

C PROCESS STREAMWISE REYNOLDS STRESS 
         read(2123,'(512E24.16)') (v(K),K=1,KM)
c make symmetric
        do j=1,km/2+1
          tmp=(v(km+1-j)-v(j))/DBLE(2)
          v(j)=-tmp
          v(km+1-j)=tmp
        enddo
C integral of ^2 over wall-normal direction
         do k=1,km
           vs(k)=vs(k)+(v(k)-vs(k))/DBLE(i)
           tt(k)=vs(k)
         enddo
         call ctran(tt)
         su3=DBLE(0)
         do j=1,km
           su3=su3+wait(j)*tt(j)
         enddo
         sumv=su3
         do k=1,km
           tt(k)=(vs(k)-sumv)**2
         enddo
         call ctran(tt)
         su3=DBLE(0)
         do j=1,km
           su3=su3+wait(j)*tt(j)
         enddo
         rmsv=su3

C add into time averaging sum
         do k=1,km
          w(k)=u(k)-v(k)+x(k)
          ws(k)=ws(k)+(w(k)-ws(k))/DBLE(i)
          tt(k)=ws(k)
         enddo
         call ctran(tt)
         su3=DBLE(0)
         do j=1,km
           su3=su3+wait(j)*tt(j)
         enddo
         sumw=su3
         do k=1,km
          tt(k)=(ws(k)-sumw)**2
         enddo
         call ctran(tt)
         su3=DBLE(0)
         do j=1,km
           su3=su3+wait(j)*tt(j)
         enddo
         rmsw=su3
         write(3117,'(i8,512E24.16)') i,
     & sumu**2+rmsu,
     & sumv**2+rmsv,
     & sumw**2+rmsw
       enddo

      end

      SUBROUTINE CHEBDZ( a, b, km  )
C     ------------------------------
C - THIS SUBR. PUTS INTO B THE CHEB EXPANSION
C - OF THE DERIVATIVE OF THE CHEB FIELD A.
C
      DIMENSION a(km), b(km)
        b(km)=DBLE(0)
        b(km-1)=a(km)*DBLE(2*(km-1))
        do k=km-2,1,-1
         b(k)=b(k+2)+DBLE(2*k)*a(k+1)
        enddo
      RETURN
      END

      SUBROUTINE CHEBDZINV( b, a, km  )
C     ------------------------------
C - INVERSE OF  CHEBDZ
C
      DIMENSION a(km), b(km)
        a(km)=b(km-1)/DBLE(2*(km-1))
        a(1)=DBLE(0)
        do k=2,km-1
         a(k)=(b(k-1)-b(k+1))/DBLE(2*(k-1))
        enddo
      RETURN
      END

      SUBROUTINE CITRAN(W)
C     --------------------
C - THIS SUBR. EVALUATES A CHEB EXPANSION AT COLLOC POINTS :
C - ZK = COS((K-1)*PI/(KM-1)). THE FIRST TERM IN THE EXPANSION
C - IS WILL BE MULTIPLED BY HALF ( BUT NOT THE LAST TERM).
C
C - THE DIMENSION OF W SHOULD BE AT LEAST 2*KM
C - KM SHOULD BE OF THE FORM 2**M + 1
C
      DIMENSION W(2*129)
      KM=129
C
      W(KM) = DBLE(2)*W(KM)
      CALL DCT(W,KM)
C
      RETURN
      END
      SUBROUTINE CTRAN(W)
C     -------------------
C - THIS SUBR. APPROXIMATES THE THE CHEB EXPANSION COEFFS.
C - OF A FUNCTION GIVEN ITS VALUE AT KM COLLOCATION POINTS.
C - THE POINTS ARE ZK = COS((K-1)*PI/(KM-1)).
C - W MUST BE DIMENSIONED WITH SIZE AT LEAST 2*KM
C - KM MUST BE OF THE FORM 2**M + 1.
C
      DIMENSION W(2*129)
      KM=129
C
      CALL DCT(W,KM)
C
      N = KM - 1
      FACTOR = DBLE(2)/DBLE(N)
      DO 10 K = 1, KM
         W(K) = FACTOR*W(K)
   10 CONTINUE
      W(KM) = W(KM)/DBLE(2)
C
      RETURN
      END

      SUBROUTINE DCT(X,KM)
C **********************************************************************
C *   THIS ROUTINE PERFORMS A DISCRETE CHEBYSHEV TRANSFORM      *
C *   ROUTINE IS NOT THREAD SAFE!                               *
C *   X      -    INPUT VECTOR DIM X(KM)                        *
C *   Y      -    OUTPUT VECTOR DIM Y(KM)                       *
C **********************************************************************
      INCLUDE 'fftw3.f'
      DIMENSION X(KM)
      DIMENSION XX(KM-1)
      DIMENSION YY(KM+1)
      DIMENSION WSAVE(KM)

       DT = DBLE(4)*DATAN(DBLE(1))/DBLE(KM-1)
       DO K=2,KM/2
         WSAVE(K) = DBLE(2)*DSIN(DBLE(K-1)*DT)
         WSAVE(KM+1-K) = DBLE(2)*DCOS(DBLE(K-1)*DT)
       ENDDO
       FAC=DBLE(1)/DBLE(2)
       MODN = MOD(KM,2)
       CALL DFFTW_PLAN_R2R_1D(PLAN,KM-1,XX,YY,FFTW_R2HC,
     *    FFTW_ESTIMATE+FFTW_UNALIGNED)
        C1 = X(1)-X(KM)
        XX(1) = X(1)+X(KM)
        DO K=2,KM/2
          T1 = X(K)+X(KM+1-K)
          T2 = X(K)-X(KM+1-K)
          C1 = C1+WSAVE(KM+1-K)*T2
          T2 = WSAVE(K)*T2
          XX(K) = T1-T2
          XX(KM+1-K) = T1+T2
        ENDDO
        IF (MODN .NE. 0) XX(KM/2+1)=X(KM/2+1)+X(KM/2+1)
        CALL DFFTW_EXECUTE_R2R(PLAN,XX,YY)
        X(1)=YY(1)*FAC
        X(2) = C1*FAC
        DO K=2,KM/2
          X(2*K) = X(2*K-2)-YY(KM+1-K)*FAC
          X(2*K-1) = YY(K)*FAC
        ENDDO
        IF (MODN .NE. 0) X(KM) = YY(KM/2+1)*FAC
      CALL DFFTW_DESTROY_PLAN(PLAN)
      RETURN
      END

