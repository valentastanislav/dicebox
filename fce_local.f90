module lokalni_fce
contains
!***********************************************************************
SUBROUTINE progress_1(i,j,celek,switch)
integer(kind=4) ::  i,j,celek,switch
      open(6)
      if (j.NE.celek) then
        if (switch.EQ.1) then
          call print_bar(i,j,celek)
        else
          call delete_bar(j,celek)
        endif
      else
        if (switch.EQ.1) then
          call print_bar(i,j,celek)
          write(6,*)
        endif
      endif
      close(6)
      return
END SUBROUTINE progress_1
!***********************************************************************
SUBROUTINE print_bar(i,j,celek)
integer(kind=4) ::  j,celek,k,i
character(len=1) :: bar, back, mezera
bar = '='
mezera = '-'
      ! print the percentage and the bar
      write(6,'(1x,1a11,1i4,2x,1i3,1a1,2x,1a1,256a1,1a1)', advance='no') 'Realizace #',i,100*j/celek,'%','[', &
      (bar, k =1,50*j/celek),(mezera, k=(50*j/celek+1),50),']'
END SUBROUTINE print_bar
!***********************************************************************
SUBROUTINE delete_bar(j,celek)
integer(kind=4) ::  j,celek,k
character(len=1) :: bar, back
back = char(8)
      ! delete the bar and the percentage
      write(6,'(256a1)', advance='no') (back, k =1,(50)+26)
END SUBROUTINE delete_bar
!***********************************************************************
SUBROUTINE WRITEMTRX(OUTNAME,MATRIX,SLOUPCE,NSREAL)        !should be OK
INTEGER::                             ILOCAL,JLOCAL
character(12)::                       OUTNAME
integer::                             SLOUPCE,NSREAL
integer,dimension(:,:),allocatable::  MATRIX
!
      if (.not.allocated(MATRIX)) then
       write(*,*) 'matice neni alokovana, nelze tedy zapsat'
      endif 
      OPEN(UNIT=4,FILE=OUTNAME,ACCESS='SEQUENTIAL',STATUS='OLD')
  102 FORMAT(I12,I12,I12,I12)
      DO ILOCAL=1,NSREAL
        WRITE(4,*) (MATRIX(JLOCAL,ILOCAL),JLOCAL=1,SLOUPCE)
      END DO
      CLOSE(4)
      RETURN
END SUBROUTINE WRITEMTRX
!***********************************************************************
SUBROUTINE CNTRLMTRX(MATRIX,SLOUPCE,NSREAL)        !should be OK
!   Let's generate the 4xN matrix of random seeds
INTEGER::                             NEXTSTEP,J,I,K
REAL::                                dummy
integer::                             SLOUPCE,NSREAL
integer,dimension(:,:),allocatable::  MATRIX
!
!      if (.not.allocated(MATRIX)) then
!       allocate(MATRIX(1:SLOUPCE,1:NSREAL))
!      endif 
      DO I=1,SLOUPCE
        NEXTSTEP=MATRIX(I,1)
        DO J=2,NSREAL
          DO K=1,100  
            dummy=RAN0(NEXTSTEP)
          ENDDO
          MATRIX(I,J)=NEXTSTEP-J
          IF (MATRIX(I,J).LE.0) MATRIX(I,J) = MATRIX(I,J) + 2 * J 
        ENDDO
      ENDDO
      RETURN
END SUBROUTINE CNTRLMTRX
!***********************************************************************
REAL FUNCTION RAN0(ir)
!   "Minimal" RND generator of Park and Miller - from Numerical Recipes
!   seed is restored after 2147483645 steps,
!   for correct function take seed [1;2147483646]
!   TODO proc by mela vadit hodnota 1?
      INTEGER*4 ir,IA,IM,IQ,IS,mask,K
      REAL*4    AM
      PARAMETER (IA=16807,IM=2147483647,IQ=127773,IS=2836,mask=123459876)
!  other posible values of constants are:
!      parameter (IA=48271,IM=2147483647,IQ=44488,IS=3399,mask=123459876)
!      parameter (IA=69621,IM=2147483647,IQ=30845,IS=23902,mask=123459876)
!
!      ir=ieor(ir,mask)
      AM=1./float(IM)
      K=ir/IQ
      ir=IA*(ir-K*IQ)-IS*K
      if (ir.LT.0) ir=ir+IM
      RAN0=AM*float(ir)
!      ir=ieor(ir,mask)
      return
      END FUNCTION RAN0
!***********************************************************************
INTEGER FUNCTION ISUBSC(A)
!     The functio ISUBSC is needed in such situations
!     when values of spin are to be transformed to values
!     that can serve as subscripts. The functio converts
!     values of a real variable A to integer values
!     according to the following rule:
!
!                A = 0.0  ISUBSC = 0
!                    0.5           0
!                    1.0           1
!                    1.5           1
!                    2.0           2
!                    2.5           2
!                     .            .
!                     .            .
!                    etc.         etc.
!
!     This functio is to be immune against a finite precission in
!     manipulations with real values.
!
      INTEGER*4 I
      REAL*4    X,A
!      
      I=INT(A)
      X=A-FLOAT(I)
      IF (ABS(X-.5)-0.25) 1,2,2
!
!     The case of half-integer spin:
!
    1 ISUBSC=I
      RETURN
!
!     The case of integer spin:
!
    2 ISUBSC=NINT(A)
      RETURN
      END FUNCTION ISUBSC
!***********************************************************************
INTEGER FUNCTION ITYPE(SP_IN,P_IN,SP_FI,P_FI)
!     The functio ITYPE determines the type of gamma-transition.
!     The meaning of the variables is evident.
!      ITYPE=1: Pure E1-transition
!            2: Mixed (M1+E2)-transition
!            3: Pure M1-transition
!            4: Pure E2-transition
!     This functio also checks wheather the spin of the final state
!     falls within the limits 0. to 50. or -- in case of odd product
!     nuclei -- within the limits 0.5 to 50.5.
!
      INTEGER*4 P_IN,P_FI,JL,JU
      REAL*4    SP_IN,SP_FI
!      
      JL=NINT(ABS(2.*SP_IN-2.*SP_FI))
      JU=NINT(2.*SP_IN+2.*SP_FI)
      IF((SP_FI.LT.(-.25)).OR.(SP_FI.GT.(50.75))) GOTO 3
      IF ((JL.LE.2).AND.(JU.GE.2)) GOTO 1
      IF ((JL.LE.4).AND.(JU.GE.4)) GOTO 2
    3 ITYPE=0
      RETURN
    2 IF (P_IN.NE.P_FI) GOTO 3
      ITYPE=4
      RETURN
    1 IF ((JL.LE.4).AND.(JU.GE.4)) GOTO 4
      IF (P_IN.EQ.P_FI) GOTO 5
    6 ITYPE=1
      RETURN
    5 ITYPE=3
      RETURN
    4 IF (P_IN.NE.P_FI) GOTO 6
      ITYPE=2
      RETURN
      END FUNCTION ITYPE
!***********************************************************************
REAL FUNCTION ALPH_TOT (EI,SPI,IPI,EF,SPF,IPF,DMISQ,NEN,ELEN,CONV)
!
!           DMISQ   - SQUARED (!) mixing amplitude
!
!           MAEL0   - the type of the lowest order of gamma radiation.
!                     for MAgnetic radiation MAEL0=1, for ELectric
!                     MAEL0=0
!
!           MUL0    - the lowest multipolarity contributing to the
!                     transition (say, 2 in case of E2+M3)
!
!           MUL1    - the next contributing multipolarity
!                     (3 in this case)
!
real,dimension(1:100):: ELEN
real,dimension(0:1,1:5,1:100):: CONV
      INTEGER*4 MUL0,MUL1,IPI,IPF,NEN,IAUX,MAEL0,MAEL1
      REAL*4    EI,SPI,EF,SPF,DMISQ,EG,CT0,CT1
!
      MUL0=NINT(ABS(SPI-SPF))
      IF (MUL0.EQ.0) MUL0=1
      MUL1=MUL0+1
      IAUX=(-1)**(IPI+IPF+MUL0+1)
      IF (IAUX.EQ.-1) THEN
        MAEL0=0    ! The dominating type is electric
      ELSE
        MAEL0=1    ! ... magnetic
      ENDIF
!
!     MAEL1 - the type  of the next-order contributing radiation
!
      MAEL1=1-MAEL0
      EG=EI-EF
      CT0=AICC(EG,ELEN,CONV,MAEL0,MUL0,NEN)
      CT1=AICC(EG,ELEN,CONV,MAEL1,MUL1,NEN)
!
      ALPH_TOT=(CT0+CT1*DMISQ)/(1.+DMISQ)
      RETURN
      END FUNCTION ALPH_TOT
!***********************************************************************
REAL FUNCTION AICC(ETRA,TABEN,TABICC,MAEL,MUL,N)
!***********************************************************************
!
!     This subroutine provides cubic interpolation of values
!     Y=TABICC(MAEL,MUL,I) for fixed MAEL and MUL and variable I.
!     These values Y for various I are assumed to represent internal
!     conversion coefficients (ICC's) for generally non-equdistant
!     transition energies that are specified by TABEN(I). MAEL stands for
!     the type of a transition (0 for M1 and 1 for E1) while MUL means
!     multipolarity. N is the length of the table for a fixed MAEL and MUL.
!     Be careful, if ETRA (the energy at which it is desirable to get
!     the ICC-coefficient) is lower than TABEN(1) or higher than TABEN(N),
!     then the interpolation changes to extrapolation and difficulties may
!     start ...
!                                                    Version from 6-OCT-95
!     updated on 21-MAY-2025 to avoid problems with "zeroes" (subthreshold for K- and pair-conversion)
!
      real,dimension(1:100):: TABEN
      real,dimension(0:1,1:5,1:100):: TABICC
      integer:: MAEL,MUL,N
      INTEGER:: I,J,K
      REAL*4 XX(4),YY(4),A(4),ETRA
      logical :: LOG_INTERPOLATE
!
      IF (ETRA.LE.TABEN(1)) THEN
        AICC=TABICC(MAEL,MUL,1)
        RETURN
      ELSEIF (ETRA.GE.TABEN(N)) THEN
        AICC=TABICC(MAEL,MUL,N)
        RETURN
      ELSEIF (ETRA.LE.TABEN(2)) THEN
        AICC=TABICC(MAEL,MUL,1)+(ETRA-TABEN(1))*(TABICC(MAEL,MUL,2)-TABICC(MAEL,MUL,1))/(TABEN(2)-TABEN(1))
        RETURN
      ELSE
        IF (ETRA.LE.TABEN(N-1)) THEN
          DO I=3,N-1
            IF (ETRA.LE.TABEN(I)) THEN
              K=I-3
              GO TO 1   !TODO ODSTRANIT GO TO
            ENDIF
          ENDDO !I
        ELSE
          K=N-4
        ENDIF
      ENDIF
!
    1 LOG_INTERPOLATE = .TRUE.
      do J = 1, 4
        if (TABICC(MAEL, MUL, K + J) .LE. 1e7*tiny(TABICC)) then
          LOG_INTERPOLATE = .FALSE.
          exit
        endif
      enddo
      IF(LOG_INTERPOLATE) THEN
        DO J=1,4
           XX(J)=log(TABEN(K+J))
           YY(J)=log(TABICC(MAEL,MUL,K+J))
        ENDDO !J
        AICC=0.
        DO I=1,4
          A(I)=YY(I)
          DO J=1,4
            IF (I.NE.J) A(I)=A(I)/(XX(I)-XX(J))
          ENDDO !J
          DO J=1,4
            IF (I.NE.J) A(I)=A(I)*(log(ETRA)-XX(J))
          ENDDO !J
          AICC=AICC+A(I)
        ENDDO !I
        AICC=exp(AICC)
      ELSE
      ! Fallback: linear interpolation in linear space using closest bracket
        do I = K+1, K+3
          if (ETRA .LE. TABEN(I+1)) then
            AICC=TABICC(MAEL,MUL,I)+(ETRA-TABEN(I))*(TABICC(MAEL,MUL,I+1)-TABICC(MAEL,MUL,I))/(TABEN(I+1)-TABEN(I))
            exit
          endif
        enddo
      ENDIF
      RETURN
      END FUNCTION AICC
!***********************************************************************
      SUBROUTINE rsm1(nm,n,a,w)
!
      integer n,nm,ierr
      integer k1,k2
      double precision a(nm,n),w(n),fwork2(n),fwork1(n)
!
!     this subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find all of the eigenvalues and some of the eigenvectors
!     of a real symmetric matrix.
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrix  a.
!
!        a  contains the real symmetric matrix.
!
!        m  the eigenvectors corresponding to the first m eigenvalues
!           are to be computed.
!           if m = 0 then no eigenvectors are computed.
!           if m = n then all of the eigenvectors are computed.
!
!     on output
!
!        w  contains all n eigenvalues in ascending order.
!
!        z  contains the orthonormal eigenvectors associated with
!           the first m eigenvalues.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for tqlrat,
!           imtqlv and tinvit.  the normal completion code is zero.
!
!        fwork  is a temporary storage array of dimension 8*n.
!
!        iwork  is an integer temporary storage array of dimension n.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
      ierr = 10 * n
      if (n .le. nm) then
        k1 = 1
        k2 = k1 + n
!     .......... find eigenvalues only ..........
!      call  tred1(nm,n,a,w,fwork(k1),fwork(k2))
!      call  tqlrat(n,w,fwork(k2),ierr)
        CALL  tred1(nm,n,a,w,fwork1,fwork2)
        CALL  tqlrat(n,w,fwork2,ierr)
      endif
      return
      END SUBROUTINE rsm1
!***********************************************************************
      SUBROUTINE tred1(nm,n,a,d,e,e2)
!
      integer i,j,k,l,n,ii,nm,jp1
      double precision a(nm,n),d(n),e(n),e2(n)
      double precision f,g,h,scale
!
!     this subroutine is a translation of the algol procedure tred1,
!     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
!     this subroutine reduces a real symmetric matrix
!     to a symmetric tridiagonal matrix using
!     orthogonal similarity transformations.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        a contains the real symmetric input matrix.  only the
!          lower triangle of the matrix need be supplied.
!
!     on output
!
!        a contains information about the orthogonal trans-
!          formations used in the reduction in its strict lower
!          triangle.  the full upper triangle of a is unaltered.
!
!        d contains the diagonal elements of the tridiagonal matrix.
!
!        e contains the subdiagonal elements of the tridiagonal
!          matrix in its last n-1 positions.  e(1) is set to zero.
!
!        e2 contains the squares of the corresponding elements of e.
!          e2 may coincide with e if the squares are not needed.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      do i = 1, n
         d(i) = a(n,i)
         a(n,i) = a(i,i)
      enddo
!     .......... for i=n step -1 until 1 do -- ..........
      do 300 ii = 1, n
        i = n + 1 - ii
        l = i - 1
        h = 0.0d0
        scale = 0.0d0
        if (l .ge. 1) then
!     .......... scale row (algol tol then not needed) ..........
          do k = 1, l
            scale = scale + dabs(d(k))
          enddo
!
         if (scale .ne. 0.0d0) go to 140
!
           do j = 1, l
              d(j) = a(l,j)
              a(l,j) = a(i,j)
              a(i,j) = 0.0d0
           enddo
!
        endif
        e(i) = 0.0d0
        e2(i) = 0.0d0
        go to 300
!
  140   do k = 1, l
          d(k) = d(k) / scale
          h = h + d(k) * d(k)
        enddo
!
        e2(i) = scale * scale * h
        f = d(l)
        g = -dsign(dsqrt(h),f)
        e(i) = scale * g
        h = h - f * g
        d(l) = f - g
        if (l .ne. 1) then
!     .......... form a*u ..........
          do j = 1, l
            e(j) = 0.0d0
          enddo
          do j = 1, l
            f = d(j)
            g = e(j) + a(j,j) * f
            jp1 = j + 1
            if (l .ge. jp1) then
              do k = jp1, l
                g = g + a(k,j) * d(k)
                e(k) = e(k) + a(k,j) * f
              enddo
            endif
            e(j) = g
          enddo
!     .......... form p ..........
          f = 0.0d0
          do j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
          enddo
          h = f / (h + h)
!     .......... form q ..........
          do j = 1, l
            e(j) = e(j) - h * d(j)
          enddo
!     .......... form reduced a ..........
          do j = 1, l
            f = d(j)
            g = e(j)
            do k = j, l
              a(k,j) = a(k,j) - f * e(k) - g * d(k)
            enddo
          enddo
        endif
        do j = 1, l
          f = d(j)
          d(j) = a(l,j)
          a(l,j) = a(i,j)
          a(i,j) = f * scale
        enddo
!
  300 continue
!
      return
      END SUBROUTINE tred1
!***********************************************************************

      SUBROUTINE SORT (N,X,XN)
!
!        ALGORITHM AS 304.8 APPL.STATIST. (1996), VOL.45, NO.3
!
!        Sorts the N values stored in array X in ascending order
!
      INTEGER N
      DOUBLE PRECISION X(N),XN(N),XTEMP(N)
!
      INTEGER I, J, INCR
      DOUBLE PRECISION TEMP
!
    	DO I = 1, N
       XTEMP(I) = X(I) 
	    ENDDO
!
      INCR = 1
!
!        Loop : calculate the increment
!
   10 INCR = 3 * INCR + 1
      IF (INCR .LE. N) GOTO 10

!
!        Loop : Shell-Metzner sort
!
   20 INCR = INCR / 3
      I = INCR + 1
   30 IF (I .GT. N) GOTO 60
      TEMP = X(I)
      J = I
   40 IF (X(J - INCR) .LT. TEMP) GOTO 50
      X(J) = X(J - INCR)
      J = J - INCR
      IF (J .GT. INCR) GOTO 40
   50 X(J) = TEMP
      I = I + 1
      GOTO 30
   60 IF (INCR .GT. 1) GOTO 20
!
	    DO I = 1, N
        XN(I) = X(I)
        X(I) = XTEMP(I) 
	    ENDDO
!
      RETURN
      END SUBROUTINE SORT
!***********************************************************************
double precision FUNCTION epslon (x)
double precision:: x
!
!     estimate unit roundoff in quantities of size x.
!
double precision:: a,b,c,eps
!
!     this program should functio properly on all systems
!     satisfying the following two assumptions,
!        1.  the base used in representing floating point
!            numbers is not a power of three.
!        2.  the quantity  a  in statement 10 is represented to
!            the accuracy used in floating point variables
!            that are stored in memory.
!     the statement number 10 and the go to 10 are intended to
!     force optimizing compilers to generate code satisfying
!     assumption 2.
!     under these assumptions, it should be true that,
!            a  is not exactly equal to four-thirds,
!            b  has a zero for its last bit or digit,
!            c  is not exactly equal to one,
!            eps  measures the separation of 1.0 from
!                 the next larger floating point number.
!     the developers of eispack would appreciate being informed
!     about any systems where these assumptions do not hold.
!
!     this version dated 4/6/83.
!
      a = 4.0d0/3.0d0
   10 b = a - 1.0d0
      c = b + b + b
      eps = dabs(c-1.0d0)
      if (eps .eq. 0.0d0) go to 10
      epslon = eps*dabs(x)
      return
      END FUNCTION epslon
!***********************************************************************
      double precision FUNCTION pythag(a,b)
      double precision a,b
!
!     finds dsqrt(a**2+b**2) without overflow or destructive underflow
!
      double precision p,r,s,t,u
!
      p = dmax1(dabs(a),dabs(b))
      if (p .eq. 0.0d0) go to 20
      r = (dmin1(dabs(a),dabs(b))/p)**2
   10 continue
         t = 4.0d0 + r
         if (t .eq. 4.0d0) go to 20
         s = r/t
         u = 1.0d0 + 2.0d0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      END FUNCTION pythag
!***********************************************************************
! not used in res_all
      SUBROUTINE CLEBSCH(AJ,BJ,CJ,AM,BM,CM,CG)
!      to calculate Clebsch-Gordan coefficients
!      You need to add a "NED(AJ,BJ,CJ,AM,BM,CM,CG)" in your main routine
!      Input:
!        AJ,BJ,CJ,AM,BM,CM (the usual Clebsch-Gordan indices)
!      Output:
!        CG=C-G(AJ,BJ,CJ,AM,BM,CM)
!
      DIMENSION Q(100,100)
      INTEGER ZZ,I,K
      REAL*4 AJ,BJ,CJ,AM,BM,CM,CG
!      
      ZZ=MAX(2*AJ+1,2*BJ+1,2*CJ+1,AJ+BJ+CJ,AJ+AM,BJ+BM,CJ+CM)+2
      DO I=1,ZZ
        Q(I,1)=1.D0
        Q(I,I)=1.D0
      ENDDO
      DO I=2,ZZ-1
        DO K=2,I
          Q(I+1,K)=Q(I,K-1)+Q(I,K)
        ENDDO
      ENDDO
      CG=0.D0
      JA=AJ+AM+1.01D0
      MA=AJ-AM+1.01D0
      JB=BJ+BM+1.01D0
      MB=BJ-BM+1.01D0
      JC=CJ+CM+1.01D0
      MC=CJ-CM+1.01D0
      LA=BJ+CJ-AJ+1.01D0
      LB=CJ+AJ-BJ+1.01D0
      LC=AJ+BJ-CJ+1.01D0
      LT=AJ+BJ+CJ+1.01D0
      D=ABS(AM+BM-CM)-0.01D0
      IF (D) 10,10,20
10    LD=MIN0(JA,JB,JC,MA,MB,MC,LA,LB,LC)
      IF (LD) 20,20,30
30    JA2=AJ+AJ+AM+AM
      JB2=BJ+BJ+BM+BM
      JC2=CJ+CJ-CM-CM
      I2=JA2+JB2+JC2-JA2/2*2-JB2/2*2-JC2/2*2
      IF (I2) 20,40,20
40    FN=Q(JA+MA-1,LC)/Q(LT,JC+MC-1)
      FN=FN*Q(JB+MB-1,LC)/Q(LT+1,2)
      FN=FN/Q(JA+MA-1,JA)
      FN=FN/Q(JB+MB-1,JB)
      FN=FN/Q(JC+MC-1,JC)
      K0=MAX(0,LC-JA,LC-MB)+1
      K1=MIN(LC,MA,JB)
      X=0.D0
      DO K=K0,K1
        X=-X-Q(LC,K)*Q(LB,MA-K+1)*Q(LA,JB-K+1)
      ENDDO
      IP=K1+LB+JC
      P=1-2*(IP-IP/2*2)
      CG=P*X*SQRT(FN)
! What we've calculated is a Wigner 3-j coefficient
! Next, we'll turn it into a Clebsch-Gordan coefficient
      CG=CG*SQRT(2*CJ+1)*(-1)**NINT(AJ-BJ-CM)
20    CONTINUE
      RETURN
      END SUBROUTINE CLEBSCH
!***********************************************************************

      SUBROUTINE TQLRAT(N,D,E2,IERR)
!
      INTEGER I,J,L,M,N,II,L1,MML,IERR
      DOUBLE PRECISION D(N),E2(N)
      DOUBLE PRECISION B,C,F,G,H,P,R,S,T!,EPSLON,PYTHAG
!
!     This subroutine is a translation of the Algol procedure tqlrat,
!     Algorithm 464, Comm. ACM 16, 689(1973) by Reinsch.
!
!     This subroutine finds the eigenvalues of a symmetric
!     tridiagonal matrix by the rational QL method.
!
!     On input
!
!        N is the order of the matrix.
!
!        D contains the diagonal elements of the input matrix.
!
!        E2 contains the squares of the subdiagonal elements of the
!          input matrix in its last N-1 positions.  E2(1) is arbitrary.
!
!      On output
!
!        D contains the eigenvalues in ascending order.  If an
!          error exit is made, the eigenvalues are correct and
!          ordered for indices 1,2,...IERR-1, but may not be
!          the smallest eigenvalues.
!
!        E2 has been destroyed.
!
!        IERR is set to
!          zero       for normal return,
!          J          if the J-th eigenvalue has not been
!                     determined after 30 iterations.
!
!     Calls PYTHAG for  DSQRT(A*A + B*B) .
!
!     Questions and comments should be directed to Burton S. Garbow,
!     Mathematics and Computer Science Div, Argonne National Laboratory
!
!     This version dated August 1987.
!     Modified by C. Moler to fix underflow/overflow difficulties,
!     especially on the VAX and other machines where epslon(1.0d0)**2
!     nearly underflows.  See the loop involving statement 102 and
!     the two statements just before statement 200.
!
!     ------------------------------------------------------------------
!

      IERR = 0
      IF (N .EQ. 1) GO TO 1001
!
      DO I = 2, N
        E2(I-1) = E2(I)
      ENDDO
!
      F = 0.0D0
      T = 0.0D0
      E2(N) = 0.0D0
!
      DO L = 1, N
         J = 0
         H = DABS(D(L)) + DSQRT(E2(L))
         IF (T .GT. H) GO TO 105
         T = H
         B = EPSLON(T)
         C = B * B
         if (c .ne. 0.0d0) go to 105
!        Spliting tolerance underflowed.  Look for larger value.
         do i = l, n
            h = dabs(d(i)) + dsqrt(e2(i))
            if (h .gt. t) t = h
         enddo
         b = epslon(t)
         c = b * b
!     .......... LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT ..........
  105    DO M = L, N
            IF (E2(M) .LE. C) GO TO 120
!     .......... E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
!                THROUGH THE BOTTOM OF THE LOOP ..........
         ENDDO
!
  120    IF (M .EQ. L) GO TO 210
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
!     .......... FORM SHIFT ..........
         L1 = L + 1
         S = DSQRT(E2(L))
         G = D(L)
         P = (D(L1) - G) / (2.0D0 * S)
         R = PYTHAG(P,1.0D0)
         D(L) = S / (P + DSIGN(R,P))
         H = G - D(L)
!
         DO I = L1, N
           D(I) = D(I) - H
         ENDDO
!
         F = F + H
!     .......... RATIONAL QL TRANSFORMATION ..........
         G = D(M)
         IF (G .EQ. 0.0D0) G = B
         H = G
         S = 0.0D0
         MML = M - L
!     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO II = 1, MML
            I = M - II
            P = G * H
            R = P + E2(I)
            E2(I+1) = S * R
            S = E2(I) / R
            D(I+1) = H + S * (H + D(I))
            G = D(I) - E2(I) / G
!           Avoid division by zero on next pass
            if (g .eq. 0.0d0) g = epslon(d(i))
            h = g * (p / r)
         ENDDO
!
         E2(L) = S * G
         D(L) = H
!     .......... GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST ..........
         IF (H .EQ. 0.0D0) GO TO 210
         IF (DABS(E2(L)) .LE. DABS(C/H)) GO TO 210
         E2(L) = H * E2(L)
         IF (E2(L) .NE. 0.0D0) GO TO 130
  210    P = D(L) + F
!     .......... ORDER EIGENVALUES ..........
         IF (L .EQ. 1) GO TO 250
!     .......... FOR I=L STEP -1 UNTIL 2 DO -- ..........
         DO II = 2, L
            I = L + 2 - II
            IF (P .GE. D(I-1)) GO TO 270
            D(I) = D(I-1)
         ENDDO
!
  250    I = 1
  270    D(I) = P
      ENDDO
!
      GO TO 1001
!     .......... SET ERROR -- NO CONVERGENCE TO AN
!                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END SUBROUTINE TQLRAT
end module lokalni_fce