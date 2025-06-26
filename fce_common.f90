!TODO gilbert cameron and spin distribution of initial states email jutta escher
!TODO NOPTE1.EQ.65 -> 75, added scaling with PAR_E1(1)
!TODO ETR deleted, replaced by PAR_E1(1) or PAR_M1(1) in relevant models
!TODO added TCONST via PAR_E1(2) parameter, now in models 39,46,56,74,75,76,79
!TODO deleted 78 as it was the 74
!TODO check about the models 48,49
module spolecne
use lokalni_fce
integer::                             nbin,LMODE,LDENP,LDSTAG,NLD,NGIGE,NLOWLOR,NGIGM,NGIGE2,NOPTE1,NOPTM1,NOPTE2
integer::                             max_decays,numlev,NOPTDE

real::                                factnrm,BN,AMASS,DELTA,PAIRING,FJ
real::                                ASHELL09,DEL09,TEMPER09,EZERO09,PAIRING09,SIG_CUSTOM,EZERO,TEMPER,DEL,ASHELL
real::                                TKpair, TKeshell, TKematch ! Kawano's stuff
real::                                DENLO,DENHI,DENPA,DENPB,DENPC,DENPD,ZNUM,DENPPC,DENPA0,DENPA1,DENPA2
real::                                FERMC,PAIR_PSF,DEG,DMG,QEL,EK0
real::                                EGZERO,DIPSLP,DIPZER,EFEC_E
integer, dimension(1:199)::           denum,ilowip,LVL_CLASS !TODO maybe make these allocatable
integer, dimension(0:49,0:1)::        NDIS,NEIGENVAL
integer, dimension(1:199,1:20)::      delev,deparity
integer,dimension(:,:),allocatable::  ityp

real,    dimension(1:2)::             spinc
real,    dimension(1:4)::             PAR_E1,PAR_M1
real,    dimension(1:5)::             ER,SIG,W0,ERM,SIGM,WM0,ERE,SIGE,WE0
real,    dimension(1:199)::           LVL_ENERGY,prim,errprim,elowlev,elowsp !TODO maybe make these allocatable
real,    dimension(0:330)::           TABENLD !TODO make this allocatable
real,    dimension(1:199,0:20)::      sal,errsal,alpha !TODO somehow smart determine the maximum number of decays in DIS and make these allocatable
real,    dimension(0:24,0:20)::       F4
real,    dimension(1:199,1:20)::      despin
real,    dimension(1:100,0:49,0:1)::  ENDIS
real,    dimension(0:330,0:49,0:1)::  TABLD !TODO make this allocatable
integer, dimension(1:3)::             NPSF
real,    dimension(1:3,0:400)::       TABENPSF,TABPSF
real,    dimension(0:1000,0:49,0:1):: EIGENVAL
real,    dimension(1:4,1:4,0:24,0:20):: Fk

contains
!***********************************************************************
SUBROUTINE ADJUST_NBIN(SPC,NBIN)
integer::                             NBIN
real::                                SPC
INTEGER,PARAMETER::                   MAXJC = 49
INTEGER::                             n_adjust,I,n_nul,J,IP
REAL::                                SP,E
!globalni BN,DELTA
        n_adjust=-1
        DO I=NBIN,1,-1
          SP=SPC-INT(SPC+.25)-1.
          n_nul=0
          E=BN+DELTA/2.-FLOAT(I)*DELTA
          DO J=0,MAXJC
            SP=SP+1.
              DO IP=0,1
                IF (DENSITY(E,SP,IP).EQ.0.) n_nul=n_nul+1
              ENDDO
          ENDDO
          IF (n_nul.EQ.((MAXJC+1)*2)) THEN
  235       FORMAT('in bin #',I3,' that is under energy',F11.6,' total density is zero')
            WRITE(*,235) I,(BN+(FLOAT(1-I))*DELTA)
            WRITE(*,*) I,(BN+(FLOAT(1-I))*DELTA),DELTA,NBIN
            n_adjust=I
          ENDIF
        ENDDO
        IF(n_adjust.NE.(-1)) THEN
          NBIN=n_adjust-1
        ENDIF
    RETURN
END SUBROUTINE ADJUST_NBIN
!***********************************************************************
SUBROUTINE READ_INT(sall,STDISa,IFLAG,U,IR4)   !should be OK
!***********************************************************************
integer::                             IFLAG,IR4
real::                                U
real,dimension(:,:),allocatable::     sall
real, dimension(0:2,0:20,-2:2,0:1)::  STDISa
INTEGER::                             I,J,K, spin_diff
REAL::                                sal_tp
!global despin,deparity,delev,denum,prim,errprim,alpha,ndis,endis,max_decays,sal,errsal
      if(.not.allocated(sall)) then
        allocate(sall(1:numlev,0:max_decays))
      endif
      DO I=1,numlev
        DO K=0,max_decays
          sall(I,K)=0.
        ENDDO
      ENDDO
      DO I=0,2
        DO J=0,20
          DO spin_diff=-2,2
            DO K=0,1
              STDISa(I,J,spin_diff,K)=0.0
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      sal_tp=0.
      IFLAG=0
      DO I=1,numlev
!   SECONDARIES part
        DO K=1,denum(I)
!         write(*,*) I,K
         DO J=1,ndis(ISUBSC(despin(I,K)),deparity(I,K))
           IF (endis(J,ISUBSC(despin(I,K)),deparity(I,K)).EQ.endis(delev(I,K),ISUBSC(despin(I,K)),deparity(I,K))) THEN
   51       sal_tp=sal(I,K)+errsal(I,K)*GAUSS(IR4,U,IFLAG) !TODO we might want to change this to better accomodate E0 transitions <-alpha here, not later
            if (sal_tp.LE.0.) goto 51  !!!TODO ask FB ??/MK .LT.
            sall(I,K)=sall(I,K-1)+sal_tp*(1+alpha(I,K))
           ENDIF
         ENDDO
        ENDDO
!   PRIMARIES part; TODO introduce fluctuations according to errprim(I)
        DO J=1,ndis(ISUBSC(elowsp(I)),ilowip(I)) ! we go over all levels of the same spin and parity to find THE level
          spin_diff = NINT(elowsp(I)+.25)-NINT(SPINc(1)+.25)
          IF ((LVL_ENERGY(I).EQ.endis(J,ISUBSC(elowsp(I)),ilowip(I))).AND.(abs(spin_diff).LE.2)) THEN
            ! write(*,*) 'in READ_INT ',I,J,spin_diff,ilowip(I),prim(I)
            STDISa(1,J,spin_diff,ilowip(I)) = &
            STDISa(1,J-1,spin_diff,ilowip(I)) + prim(I)*factnrm
          ENDIF
        ENDDO
      ENDDO ! I=1,numlev
      ! write(*,*) 'in READ_INT ',sum(prim)*factnrm,(maxval(STDISa(1,:,-2,0))+maxval(STDISa(1,:,-1,0)) &
      ! +maxval(STDISa(1,:,0,0))+maxval(STDISa(1,:,1,0))+maxval(STDISa(1,:,2,0))+maxval(STDISa(1,:,-2,1)) &
      ! +maxval(STDISa(1,:,-1,1))+maxval(STDISa(1,:,0,1))+maxval(STDISa(1,:,1,1))+maxval(STDISa(1,:,2,1)))
      RETURN
END SUBROUTINE READ_INT
!***********************************************************************
SUBROUTINE LEVELSCH(IR,IREAL,IFLAG,NTOTAL,ITID,SPC,U,LEVCON)         !should be OK
!***********************************************************************
!    - the Poisson distribution of neighbourhood level spacing is
!      assumed for LMODE=0 
!    - the Wigner distribution (with long-range correlations) is
!      assumed for LMODE=1 
!    - the "restricted Wigner distribution" - no long-range correlations
INTEGER,PARAMETER::                   MAXJC = 49                      
INTEGER::                             IP,J,I,IAUX,K,KLO,KHI
REAL::                                SP,E,AVNL,X,OMEGA,ALPHA,RA     
integer,dimension(:,:,:),allocatable::LEVCON
integer::                             IR,IREAL,IFLAG,NTOTAL,ITID
real::                                SPC,U
! vstup(BN,DELTA,LMODE,NBIN) vystup(NTOTAL,LEVCON,IFLAG)
      if (.not.allocated(LEVCON)) then
       allocate(LEVCON(0:NBIN,0:MAXJC,0:1))
      endif 
      NTOTAL=0
      DO IP=0,1
        DO J=0,MAXJC
          DO I=0,NBIN
            LEVCON(I,J,IP)=0
          ENDDO
        ENDDO
      ENDDO
      IF (LMODE.EQ.0) THEN              !Poisson distribution
       DO IP=0,1
        SP=SPC-INT(SPC+.25)-1.
        DO J=0,MAXJC
         SP=SP+1.
         DO I=1,NBIN
          E=BN+DELTA/2.-FLOAT(I)*DELTA
          AVNL=DELTA*DENSITY(E,SP,IP)    !neni nic nahodneho, pri paralelnim behu ma byt pocet pruchodu porad stejny
          LEVCON(I,J,IP)=NPOISS(IR,AVNL,IFLAG,U)
          NTOTAL=NTOTAL+LEVCON(I,J,IP)
         ENDDO !I
        ENDDO !J
       ENDDO !IP
!       WRITE(*,*) 'poissona zavolam',2*(MAXJC+1)*NBIN,' ale gausse jen',uziti_u  !should be OK
      ELSEIF (LMODE.EQ.1) THEN      !Wigner distribution with long-range correlations
       DO IP=0,1
        SP=SPC-FLOAT(INT(SPC+.25))-1.
        DO J=0,MAXJC
          SP=SP+1.
          X=10.*RAN0(IR)
          K=1
    8    IF (GOE_EIGEN_VAL(K,J,IP).LT.X) THEN
           K=K+1
           GOTO 8
         ENDIF
         KLO=K
         DO I=NBIN,1,-1
           E=BN+DELTA/2.-FLOAT(I)*DELTA
           AVNL=DELTA*DENSITY(E,SP,IP)
           X=X+AVNL
    9      IF (GOE_EIGEN_VAL(K,J,IP).LT.X) THEN 
             K=K+1
             GOTO 9
           ENDIF   
           KHI=K
           LEVCON(I,J,IP)=KHI-KLO        
           NTOTAL=NTOTAL+KHI-KLO
           KLO=KHI
         ENDDO !I
        ENDDO !J
       ENDDO !IP
      ELSE                         !Wigner distribution - no long-range correlations
       DO IP=0,1                   !pouzivam i pozici LEVCON(0,...) kvuli kumulativnimu ukladani
        SP=SPC-INT(SPC+.25)-1.
         DO J=0,MAXJC
          DO I=0,NBIN
           LEVCON(I,J,IP)=0
          ENDDO  !I
          SP=SP+1.
          OMEGA=0.
          ALPHA=0.
          I=0
!         OMEGA is a random sample drawn from the Wigner distribution;
!         the expectation value of average distance between
!         neighbouring levels of a given spin and parity is assumed
!         to be equal to 1.
!           The constant 1.1283791671 is equal to "two divided by
!           square root of pi"
    1     RA=RAN0(IR)
          IF (RA.LE.0.) GO TO 1
          OMEGA=OMEGA+1.1283791671*SQRT(-ALOG(RA))
    3     IF (OMEGA.LT.ALPHA) THEN
            LEVCON(I,J,IP)=LEVCON(I,J,IP)+1
            NTOTAL=NTOTAL+1
            GO TO 1
          ELSE
            I=I+1
            IF (I.GT.NBIN) GO TO 2
            E=BN+DELTA/2.-FLOAT(I)*DELTA
            ALPHA=ALPHA+DENSITY(E,SP,IP)*DELTA
            GO TO 3
          ENDIF
    2     CONTINUE
         ENDDO !J
       ENDDO !IP
      ENDIF
!debug line
!      CALL WRITELVLSCH(IREAL,NBIN,ITID,IFLAG,NTOTAL,LEVCON)
!     NTOTAL is the total number of generated levels.
!     At this moment for each value of I the variable LEVCON(...) contains
!     the number of those levels of a particular spin and parity that fall
!     within the corresponding energy bin (whose width is DELTA).
!     The following DO-loops, however, convert this differential distribution
!     of level energies into a CUMULATIVE (i.e. integral) form. This simple
!     conversion leads to a significant increase of the speed of
!     functio SEED.
      IAUX=0
      DO IP=0,1
       DO J=0,MAXJC
        LEVCON(0,J,IP)=IAUX
        DO I=1,NBIN
          LEVCON(I,J,IP)=LEVCON(I,J,IP)+LEVCON(I-1,J,IP)
        ENDDO 
        IAUX=LEVCON(NBIN,J,IP)
       ENDDO !J
      ENDDO !IP
      RETURN
END SUBROUTINE LEVELSCH
!***********************************************************************
REAL FUNCTION DENSITY(EEXC,SPIN,IPAR)         !should be OK
!
!     Explicit expressions for level density
!     NOPTDE= 0: CTF-model
!           = 1: Bethe's level-density formula following formulation
!               of T.von Egidy et al., Nucl.Phys. A (1988)
!           = 2: modified BSFG
!           = 3: modified CTF
!          =4,5: BSFG with modified spin cut-off parameter
!           = 6: BSFG from T. von Egidy 2005
!           =11: Goriely - tabulated level density
!
!     Changed a factor 1/2 in the row DENSITY=DENSITY*FJ*.5
!          => the parity-dependent level density allowed (PRC67,015803)
!
!***********************************************************************
REAL::                                EFEC_E,SIGSQ,FJ,PARDEP,PAIRS
real::                                EEXC,SPIN
integer::                             IPAR
!uses 'global' variables EZERO,TEMPER,AMASS,LDENP,DEL,ASHELL,various09,...
      DENSITY=0.
      IF (NOPTDE.EQ.0) THEN                                 ! CTF
        EFEC_E=EEXC-EZERO
        IF (EFEC_E.LE.0.) RETURN
        DENSITY=EXP(EFEC_E/TEMPER)/TEMPER
        SIGSQ=(.98*AMASS**.29)**2.
!        SIGSQ=(2.*AMASS**.29)**2.
      ELSEIF (NOPTDE.EQ.1) THEN                             ! BSFG
        EFEC_E=EEXC-DEL
        IF (EFEC_E.LE.0.) RETURN
        SIGSQ=.0888*SQRT(ASHELL*EFEC_E)*AMASS**.666667
        DENSITY=EXP(2.*SQRT(ASHELL*EFEC_E))/(16.9706*SQRT(SIGSQ)*ASHELL**.25*EFEC_E**1.25)
      ELSEIF (NOPTDE.EQ.2) THEN                             ! modified BSFG
        EFEC_E=EEXC-DEL
        IF (EFEC_E.LE.0.) RETURN
        SIGSQ=.0888*SQRT(ASHELL*EFEC_E)*AMASS**.666667
        DENSITY=EXP(2.*SQRT(ASHELL*EFEC_E))/(16.9706*SQRT(SIGSQ)*ASHELL**.25*EFEC_E**1.25)
        IF ((EEXC.GE.DENLO).AND.(EEXC.LE.DENHI)) THEN
         FCTDEN=DENPA+DENPB*EEXC+DENPC*EEXC**2+DENPD*EEXC**3
         DENSITY=DENSITY*FCTDEN
        ENDIF
      ELSEIF (NOPTDE.EQ.3) THEN                              ! modified CTF
        EFEC_E=EEXC-EZERO
        IF (EFEC_E.LE.0.) RETURN
        DENSITY=EXP(EFEC_E/TEMPER)/TEMPER
        SIGSQ=(.98*AMASS**.29)**2.
        IF ((EEXC.GE.DENLO).AND.(EEXC.LE.DENHI)) THEN
         FCTDEN=DENPA+DENPB*EEXC+DENPC*EEXC**2+DENPD*EEXC**3
         DENSITY=DENSITY*FCTDEN
        ENDIF
      ELSEIF (NOPTDE.EQ.4) THEN                             ! BSFG, s=0.1446 (Paar,Al-Quraishi)
        EFEC_E=EEXC-DEL
        IF (EFEC_E.LE.0.) RETURN
        SIGSQ=.1446*SQRT(ASHELL*EFEC_E)*AMASS**.666667
        DENSITY=EXP(2.*SQRT(ASHELL*EFEC_E))/(16.9706*SQRT(SIGSQ)*ASHELL**.25*EFEC_E**1.25)
      ELSEIF (NOPTDE.EQ.5) THEN                             ! BSFG, another s (Al-Quraishi)
        EFEC_E=EEXC-DEL
        IF (EFEC_E.LE.0.) RETURN
        SIGSQ=.0145*0.8*SQRT(EFEC_E/ASHELL)*AMASS**1.66667
        DENSITY=EXP(2.*SQRT(ASHELL*EFEC_E))/(16.9706*SQRT(SIGSQ)*ASHELL**.25*EFEC_E**1.25)
      ELSEIF (NOPTDE.EQ.6) THEN                             ! BSFG - Von Egidy (2006) cut-off
        EFEC_E=EEXC-DEL
        IF (EFEC_E.LE.0.) RETURN
        SIGSQ=.0146*(1+SQRT(1+4*ASHELL*EFEC_E))/2./ASHELL*AMASS**1.666667
        DENSITY=EXP(2.*SQRT(ASHELL*EFEC_E))/(16.9706*SQRT(SIGSQ)*ASHELL**.25*EFEC_E**1.25)
      ELSEIF (NOPTDE.EQ.66) THEN                             ! BSFG - Von Egidy (2006) cut-off
        EFEC_E=EEXC-DEL
        IF (EFEC_E.LE.0.) RETURN
        IF ((SPIN.EQ.2.5).AND.(EEXC.GE.DENPC).AND.(EEXC.LE.DENPD)) RETURN
        SIGSQ=.0146*(1+SQRT(1+4*ASHELL*EFEC_E))/2./ASHELL*AMASS**1.666667
        DENSITY=EXP(2.*SQRT(ASHELL*EFEC_E))/(16.9706*SQRT(SIGSQ)*ASHELL**.25*EFEC_E**1.25)
      ELSEIF (NOPTDE.EQ.67) THEN                             ! BSFG - Von Egidy (2006) cut-off
        EFEC_E=EEXC-DEL
        IF (EFEC_E.LE.0.) RETURN
        IF (((SPIN.EQ.2.5).OR.(SPIN.EQ.3.5)).AND.(EEXC.GE.DENPC).AND.(EEXC.LE.DENPD)) RETURN
        SIGSQ=.0146*(1+SQRT(1+4*ASHELL*EFEC_E))/2./ASHELL*AMASS**1.666667
        DENSITY=EXP(2.*SQRT(ASHELL*EFEC_E))/(16.9706*SQRT(SIGSQ)*ASHELL**.25*EFEC_E**1.25)
      ELSEIF (NOPTDE.EQ.68) THEN                             ! modified BSFG
        EFEC_E=EEXC-DEL
        IF (EFEC_E.LE.0.) RETURN
        SIGSQ=.0146*(1+SQRT(1+4*ASHELL*EFEC_E))/2./ASHELL*AMASS**1.666667
        DENSITY=EXP(2.*SQRT(ASHELL*EFEC_E))/(16.9706*SQRT(SIGSQ)*ASHELL**.25*EFEC_E**1.25)
        IF ((EEXC.GE.DENLO).AND.(EEXC.LE.DENHI)) THEN
          FCTDEN=DENPA+DENPB*EEXC+DENPC*EEXC**2+DENPD*EEXC**3
          DENSITY=DENSITY*FCTDEN
        ENDIF
      ELSEIF (NOPTDE.EQ.11) THEN                            ! Goriely
        IF (EEXC.LE.0.) RETURN
        DENSITY = ALD(EEXC,SPIN,IPAR)
        RETURN
      ELSEIF (NOPTDE.EQ.12) THEN                            ! Kawano
        IF (EEXC.LE.0.) RETURN
        DENSITY = ALD(EEXC,SPIN,IPAR)
        RETURN
      ELSEIF (NOPTDE.EQ.7) THEN                             ! Voinov Mo BSFG
        EFEC_E=EEXC-DEL
        IF (EFEC_E.LE.0.) RETURN
        SIGSQ=.0146*SQRT(EFEC_E/ASHELL)*AMASS**1.666667
        DENSITY=EXP(2.*SQRT(ASHELL*EFEC_E))/(16.9706*SQRT(SIGSQ)*ASHELL**.25*EFEC_E**1.25)
      ELSEIF (NOPTDE.EQ.8) THEN                             ! CTF von Egidy 09
        EFEC_E=EEXC-EZERO09
        IF (EFEC_E.LE.0.) RETURN
        IF (EEXC.LE..5*PAIRING09) RETURN
        DENSITY=EXP(EFEC_E/TEMPER09)/TEMPER09
        SIGSQ=.391*AMASS**.675*(EEXC-.5*PAIRING09)**.312
      ELSEIF (NOPTDE.EQ.9) THEN                             ! BSFG von Egidy 09
        EFEC_E=EEXC-DEL09
        IF (EFEC_E.LE.0.) RETURN
        IF (EEXC.LE..5*PAIRING09) RETURN
        SIGSQ=.391*AMASS**.675*(EEXC-.5*PAIRING09)**.312
        DENSITY=EXP(2.*SQRT(ASHELL09*EFEC_E))/(16.9706*SQRT(SIGSQ)*ASHELL09**.25*EFEC_E**1.25)
      ELSEIF (NOPTDE.EQ.13) THEN                             ! Oslo BS - CTF with BSFG von Egidy cut-off
        EEFF=EEXC-DEL
        IF (EEFF.LE.0.) RETURN
        SIGSQ=.0146*(1+SQRT(1+4*ASHELL*EEFF))/2./ASHELL*AMASS**1.666667
        EEFF=EEXC-EZERO
        DENSITY=EXP(EEFF/TEMPER)/TEMPER
      ELSEIF (NOPTDE.EQ.18) THEN                             ! CTF with custom SIGSQ
        EFEC_E=EEXC-EZERO09
        IF (EFEC_E.LE.0.) RETURN
        DENSITY=EXP(EFEC_E/TEMPER09)/TEMPER09
        SIGSQ=SIG_CUSTOM**2
      ELSEIF (NOPTDE.EQ.19) THEN                             ! BSFG with custom SIGSQ
        EFEC_E=EEXC-DEL09
        IF (EFEC_E.LE.0.) RETURN
        SIGSQ=SIG_CUSTOM**2
        DENSITY=EXP(2.*SQRT(ASHELL09*EFEC_E))/(16.9706*SQRT(SIGSQ)*ASHELL09**.25*EFEC_E**1.25)
      ELSEIF (NOPTDE.EQ.29) THEN                             ! BSFG custom for 168Er, inspired by von Egidy 09
        EFEC_E=EEXC !E1_ours=0.00; NOTE needs staggering - maximum effect up to at least 2.4~MeV, probably higher
        IF (EFEC_E.LE.0.) RETURN
        IF (EEXC.LE..5*PAIRING09) RETURN !pairing_ours=pairing_09
        SIGSQ=.28*AMASS**.695*(EEXC-.5*PAIRING09)**.14 !sigma dependence with modified coefficients
        DENSITY=EXP(2.*SQRT(15.00*EFEC_E))/(16.9706*SQRT(SIGSQ)*15.00**.25*EFEC_E**1.25) !a_ours=15.00
      ELSEIF (NOPTDE.EQ.39) THEN                             ! BSFG custom for 168Er, inspired by von Egidy 09
        EFEC_E=EEXC-DEL09                                    ! f(J,-) significantly wider at low energies
        IF (EFEC_E.LE.0.) RETURN                             ! same width as f(J,+) at neutron sep. energy
        IF (EEXC.LE..5*PAIRING09) RETURN
        IF (IPAR.EQ.0) THEN
          SIGSQ=.391*AMASS**.675*(EEXC-.5*PAIRING09)**.312
        ELSE
          SIGSQ=SIG_CUSTOM**2.
        ENDIF
        DENSITY=EXP(2.*SQRT(ASHELL09*EFEC_E))/(16.9706*SQRT(SIGSQ)*ASHELL09**.25*EFEC_E**1.25)
      ENDIF
!
      FJ=(SPIN+.5)*EXP(-(SPIN+.5)**2/(2.*SIGSQ))/SIGSQ
!     "staggering" as originally proposed by von Egidy (2009) for models 8 and 9
      IF (LDSTAG.EQ.1) THEN !staggering in + parity only
        IF (IPAR.EQ.0) THEN
          STAG = (EEXC - DENLO) / (DENHI - DENLO)
          IF (STAG.LE.0.0) STAG = 0.0
          IF (STAG.GE.1.0) STAG = 1.0
          IF (SPIN.LT.0.25) THEN
            FJ = FJ * (1.0 + 1.02 * (1.0 - STAG) )
          ELSEIF (MOD(INT(SPIN+0.25),2).EQ.0) THEN
            FJ = FJ * (1.0 + 0.227 * (1.0 - STAG) )
          ELSE
            FJ = FJ * (1.0 - 0.227 * (1.0 - STAG) )
          ENDIF
        ENDIF
      ELSEIF (LDSTAG.EQ.2) THEN !staggering in both parities
        STAG = (EEXC - DENLO) / (DENHI - DENLO)
        IF (STAG.LE.0.0) STAG = 0.0
        IF (STAG.GE.1.0) STAG = 1.0
        IF (SPIN.LT.0.25) THEN
          FJ = FJ * (1.0 + 1.02 * (1.0 - STAG) )
        ELSEIF (MOD(INT(SPIN+0.25),2).EQ.0) THEN
          FJ = FJ * (1.0 + 0.227 * (1.0 - STAG) )
        ELSE
          FJ = FJ * (1.0 - 0.227 * (1.0 - STAG) )
        ENDIF
      ENDIF

!
!     "Parity-dependence" term
!
      IF (LDENP.EQ.0) THEN
        PARDEP=0.5
      ELSEIF (LDENP.GE.2) THEN
        PAIRS=DENPA0+DENPA1/AMASS**DENPA2             !see PRC67, 015803  
      ELSEIF (LDENP.EQ.1) THEN                        !see PRC67, 015803
        IF (MOD(INT(AMASS+0.25),2).EQ.0) THEN
          IF (MOD(INT(ZNUM+0.25),2).EQ.0) THEN
            PAIRS= 1.34+75.22/AMASS**0.89             !E-E nucleus
          ELSE
            PAIRS=-0.90+75.22/AMASS**0.89             !O-O nucleus
          ENDIF
        ELSE
          IF (MOD(INT(ZNUM+0.25),2).EQ.0) THEN
            PAIRS=-0.08+75.22/AMASS**0.89             !E-O nucleus
          ELSE
            PAIRS=-0.42+75.22/AMASS**0.89             !O-E nucleus
          ENDIF
        ENDIF
      ENDIF
      IF ((LDENP.EQ.1).OR.(LDENP.EQ.2)) THEN
        IF (IPAR.EQ.0) THEN 
          PARDEP=0.5*(1+1/(1+EXP(DENPPC*(EEXC-PAIRS))))
        ELSE
          PARDEP=0.5*(1-1/(1+EXP(DENPPC*(EEXC-PAIRS))))
        ENDIF
      ELSEIF (LDENP.EQ.3) THEN
        IF (IPAR.EQ.0) THEN 
          PARDEP=0.5*(1-1/(1+EXP(DENPPC*(EEXC-PAIRS))))
        ELSE
          PARDEP=0.5*(1+1/(1+EXP(DENPPC*(EEXC-PAIRS))))
        ENDIF
      ENDIF           
!
      DENSITY=DENSITY*FJ*PARDEP
      RETURN
END FUNCTION DENSITY
!***********************************************************************
REAL FUNCTION ALD(EX,SPIN,IP)
!***********************************************************************
!
!     This subroutine provides cubic interpolation of values
!     of level density - based on the functio AICC.
!     At the lowest excitations (very low density) the linear interpolation 
!     is used while at higher energies a cubic interpolation is adopted 
!                                                    Version from 16-MAY-09
!
real,DIMENSION(1:4):: XX,YY,A
real::                EX,SPIN
integer::             IP
INTEGER::             I,J,K,NLMIN
!
      ALD = 0.0
!    Low level densities treated in a special way
      NLMIN=NLD
      DO WHILE ((TABLD(NLMIN,ISUBSC(SPIN),IP).GT.0.0).AND.(NLMIN.GT.0))
        NLMIN=NLMIN-1 !this can be anything from NLD down to 0
      ENDDO
      ! IF ((ISUBSC(SPIN).EQ.0).AND.(IP.EQ.0)) THEN
      !   write(*,*) 'E_exc ',EX,' NLMIN ',NLMIN,' E_tab ',TABENLD(NLMIN),TABENLD(NLMIN+1),' LD_tab ',TABLD(NLMIN,ISUBSC(SPIN),IP)&
      !   ,TABLD(NLMIN+1,ISUBSC(SPIN),IP)
      ! ENDIF
      IF (NLMIN.GE.(NLD-1)) THEN !this happens for highest spins which are non existent -> density is plain zero
        ! IF ((ISUBSC(SPIN).EQ.0).AND.(IP.EQ.0)) THEN
        !   write(*,*) 'lvl density comes out at #0 as ',ALD
        ! ENDIF
        RETURN
      ELSEIF ((NLMIN.EQ.0).OR.(EX.GT.TABENLD(NLMIN+2))) THEN ! NLMIN=0 means all densities are bigger than zero and any interpolation should be fine, otherwise we go two bins above the last zero density bin where the exp-log interpolation should be safe
       IF (EX.LE.TABENLD(2)) THEN
        K=0
       ELSE
        IF (EX.LE.TABENLD(NLD-1)) THEN
          DO I=3,NLD-1
            IF (EX.LE.TABENLD(I)) THEN
              K=I-3
              GO TO 1
            ENDIF
          ENDDO !I
        ELSE
        K=NLD-4
        ENDIF
       ENDIF
    1  EXLOG=log(EX)
       DO J=1,4
         XX(J)=log(TABENLD(K+J))
         YY(J)=log(TABLD(K+J,ISUBSC(SPIN),IP))
       ENDDO !J
       DO I=1,4
         A(I)=YY(I)
         DO J=1,4
            IF (I.NE.J) A(I)=A(I)/(XX(J)-XX(I))
         ENDDO !J
         DO J=1,4
            IF (I.NE.J) A(I)=A(I)*(XX(J)-EXLOG)
         ENDDO !J
         ALD=ALD+A(I)
       ENDDO !I
       ALD=exp(ALD)
      !  IF ((ISUBSC(SPIN).EQ.0).AND.(IP.EQ.0)) THEN
      !   write(*,*) 'lvl density comes out at #1 as ',ALD
      !  ENDIF
       RETURN
      ELSE !in regime near zero densities we use linear interpolation of two nearest bins
       IF (EX.LE.TABENLD(1)) THEN !this should never happen, the tables should start at ~0 MeV while the E_crit should be at least few states higher
         ALD = TABLD(1,ISUBSC(SPIN),IP)
        !  IF ((ISUBSC(SPIN).EQ.0).AND.(IP.EQ.0)) THEN
        !   write(*,*) 'lvl density comes out at #2 as ',ALD
        !  ENDIF
         RETURN
       ELSEIF (EX.GE.TABENLD(NLD)) THEN !this should never happen, the tables should go above the initial cascading energy
         ALD = TABLD(NLD,ISUBSC(SPIN),IP)
        !  IF ((ISUBSC(SPIN).EQ.0).AND.(IP.EQ.0)) THEN
        !   write(*,*) 'lvl density comes out at #3 as ',ALD
        !  ENDIF
         RETURN
       ELSE
         I=1
         DO WHILE (EX.GT.TABENLD(I))
           I=I+1
         ENDDO
         ALD = TABLD(I-1,ISUBSC(SPIN),IP) + (EX-TABENLD(I-1))*((TABLD(I,ISUBSC(SPIN),IP)-TABLD(I-1,ISUBSC(SPIN),IP))/&
         (TABENLD(I)-TABENLD(I-1)))
        !  IF ((ISUBSC(SPIN).EQ.0).AND.(IP.EQ.0)) THEN
        !   write(*,*) 'lvl density comes out at #4 as ',ALD,' with I =',I
        !  ENDIF
         RETURN
       ENDIF
      ENDIF
      ! IF ((ISUBSC(SPIN).EQ.0).AND.(IP.EQ.0)) THEN
      !   write(*,*) 'lvl density comes out at #5 as ',ALD
      ! ENDIF

      RETURN
      END FUNCTION ALD
!***********************************************************************
REAL FUNCTION APSF(EGX,MTYP)
!***********************************************************************
real,DIMENSION(1:4):: XX,YY,A
real::                EGX,EXLOG
INTEGER::             I,J,K,NLMIN
!
      APSF = 0.0
!    Low PSF treated in a special way
      NLMIN=NPSF(MTYP)
      DO WHILE ((TABPSF(MTYP,NLMIN).GT.0.0).AND.(NLMIN.GT.0))
        NLMIN=NLMIN-1 !this is the highest bin where PSF is zero, can be anything from NPSF(MTYP) down to 0
      ENDDO
      IF (NLMIN.GE.(NPSF(MTYP)-1)) THEN !if all but last two are zero, return zero
       RETURN
      ELSEIF ((NLMIN.EQ.0).OR.(EGX.GT.TABENPSF(MTYP,NLMIN+2))) THEN ! NLMIN=0 means all PSFs are bigger than zero and any interpolation should be fine, otherwise we go two bins above the last zero PSF bin where the exp-log interpolation should be safe
       IF (EGX.LE.TABENPSF(MTYP,2)) THEN
        K=0
       ELSE
        IF (EGX.LE.TABENPSF(MTYP,NPSF(MTYP)-1)) THEN
          DO I=3,NPSF(MTYP)-1
            IF (EGX.LE.TABENPSF(MTYP,I)) THEN
              K=I-3
              GO TO 1
            ENDIF
          ENDDO !I
        ELSE
        K=NPSF(MTYP)-4
        ENDIF
       ENDIF
    1  EXLOG=log(EGX)
       DO J=1,4
         XX(J)=log(TABENPSF(MTYP,K+J))
         YY(J)=log(TABPSF(MTYP,K+J))
       ENDDO !J
       DO I=1,4
         A(I)=YY(I)
         DO J=1,4
            IF (I.NE.J) A(I)=A(I)/(XX(J)-XX(I))
         ENDDO !J
         DO J=1,4
            IF (I.NE.J) A(I)=A(I)*(XX(J)-EXLOG)
         ENDDO !J
         APSF=APSF+A(I)
       ENDDO !I
       APSF=exp(APSF)
       RETURN
      ELSE !in regime near zero PSF we use linear interpolation of two nearest bins
       IF (EGX.LE.TABENPSF(MTYP,1)) THEN !this should never happen, the tables should start at ~0 MeV while the E_crit should be at least few states higher
         APSF = TABPSF(MTYP,1)
         RETURN
       ELSEIF (EGX.GE.TABENPSF(MTYP,NPSF(MTYP))) THEN !this should never happen, the tables should go above the initial cascading energy
         APSF = TABPSF(MTYP,NPSF(MTYP))
         RETURN
       ELSE
         I=1
         DO WHILE (EGX.GT.TABENPSF(MTYP,I))
           I=I+1
         ENDDO
         APSF = TABPSF(MTYP,I-1) +(EGX-TABENPSF(MTYP,I-1))*((TABPSF(MTYP,I)-TABPSF(MTYP,I-1))/(TABENPSF(MTYP,I)-TABENPSF(MTYP,I-1)))
         RETURN
       ENDIF
      ENDIF

      RETURN
      END FUNCTION APSF
!***********************************************************************
SUBROUTINE GENERATE_GOE_EIGEN_VAL(IR,N,IFLAG,U)         !should be OK
!***********************************************************************
!     Generates eigenvalues of random matrices - they are stored in
!     the EIGENVAL,NEIGENVAL and are used for generating level in the 
!     subroutine LEVELSCH via calling GOE_EIGEN_VAL
!
INTEGER,PARAMETER::                   MAXJC=49
DOUBLE PRECISION,PARAMETER::          DPI=3.141592653589793d0
INTEGER::                             IP,IS,I,J,M,L,MM
REAL::                                EIG0
REAL*8                                A(N,N),RA(N),RAORD(N),X
integer::                             IR,N,IFLAG
real::                                U
!globalni EIGENVAL,NEIGENVAL

!
      DO IP=0,1
       DO IS=0,MAXJC
        DO I=1,N
         A(I,I)=DBLE(GAUSS(IR,U,IFLAG)/SQRT(2.*FLOAT(N-1)))
         DO J=I+1,N
          A(I,J)=DBLE(GAUSS(IR,U,IFLAG)/SQRT(4.*FLOAT(N-1)))
          A(J,I)=A(I,J)
         ENDDO !J
        ENDDO !I
        CALL RSM1(N,N,A,RA)
        CALL SORT(N,RA,RAORD)
        DO I = 1,N
          RA(I) = RAORD(I)
        ENDDO

        M=-1
        DO L=1,N
         X=RA(L)
         IF (DABS(X).lt.0.9d0) THEN     !Only eigenval. in the middle are treated
          M=M+1
          EIGENVAL(M,IS,IP)=FLOAT(N-1)*(SNGL((X*DSQRT(1.d0-X**2)+DASIN(X))/DPI+.5d0))
         ENDIF
        ENDDO !L
        NEIGENVAL(IS,IP)=M
        EIG0=EIGENVAL(0,IS,IP)
        DO MM=M,0,-1
         EIGENVAL(MM,IS,IP)=EIGENVAL(MM,IS,IP)-EIG0
        ENDDO !MM
       ENDDO !IS
      ENDDO !IP
      RETURN
END SUBROUTINE GENERATE_GOE_EIGEN_VAL
!***********************************************************************
REAL FUNCTION GOE_EIGEN_VAL(N,IS,IP)  !integer aritmetics?
!***********************************************************************
!
!     Eigenvalues were generated at the beginning 
!     by SUBROUTINE GENERATE_GOE_EIGEN_VAL(IR,N,IFLAG,U)
!
INTEGER::                             NN,K
integer::                             N,IS,IP
!global NEIGENVAL,EIGENVAL
      NN=NEIGENVAL(IS,IP)
      K=N-(N/NN)*NN
      GOE_EIGEN_VAL=EIGENVAL(K,IS,IP)+FLOAT((N/NN)*NN)
      RETURN
END FUNCTION GOE_EIGEN_VAL
!***********************************************************************
REAL FUNCTION SGAMMA(EGAM,EINI,ITYP)
!
!   Photon strengths for E1, M1+E2, M1 and E2 transitions. "Photon
!   strength" does not mean "photon strength functio" here, but
!   photon strength functio multiplied by EGAM**(2*L+1) and, in
!   the case of M1+E2 transitions, summed over both XL-components.
!
!     ITYP is equal to 1, 2, 3 or 4 (see ITYPE functio)
!     EINI is the initial state energy in MeV
!     EGAM is gamma-ray energy in MeV
!
!   E1:  NOPTE1= 0: Single-particle approximation
!                1: Classical Lorentzian GDER
!                2: GDER with an energy and temperature dependent
!                   damping width (J.Kopecky, R.Chrien, Nucl.Phys.
!                   A468,p.285)
!                3: correct EGLO model - see 6
!                4: Kadmenskij-Markushev-Furman original Strength functio
!                   (no high energy approximation according to Chrien)
!                5: The Chrien's Strength functio (Nucl. Phys. A468, 285
!                   (1987)) only. In 3: is this model used only for
!                   high energy region
!                6: The strength functio according to Chrien with
!                   phenomenological temperature dependent damping
!                   proposed by Kopecky (Distribution of Radiative Strength
!                   in Gd-156, 157 and 158 Nuclei)
!                   - I have found an error, correct EGLO is 3:
!                7: GDER with phenomenological temperature dependent
!                   damping width proposed by Kopecky
!                8: 4: with the first resonance of Lorentz type
!                9: 6: with the first resonance of Lorentz type
!               10: KMF (4:) for EG<4 MeV; lin. combination of KMF and BA
!                   for 4 MeV<EG<8 MeV; BA (1:) for EG>8 MeV
!               31-40: correspond to 1-10 for high EGAM; for low EGAM original
!                      values are multiplied by a factor (given in input data)
!                      in between the PSF is a linear combination ... 
!                      motivation comes from Au
!               51: KMF (4:) without temperature-dependent term in damping
!                   width
!               52: EGLO (6:) without temperature-dependent term in damping
!                   width; temperature is taken into account only in the
!                   "second term" - FK*...
!               53: EELO (7:) without temperature-dependent term in damping
!                   width - i.e. no termperature dependence assumed
!               41: KMF according to Oslo group 
!
!   M1:  NOPTM1= 0: Single-particle approximation
!                1: Classical lorentzian GDMR
!                3: Scissors (first) resonance is build up only on states
!                   with excitation energy lower than PAR_M1(1)
!                4: Classical lorentzian build on the "background"
!                   that is described by the SP (constant functio)
!
!                ?: Enery of scissors resonance depends linearly on
!                   the energy of final state (and is build up only on states
!                   below certain excitation energy)
!                ?: Scissors resonance is considered only for primary transitions
!
!                ?: power dependence
!
!   E2:  NOPTE2= 0: Single-particle approximation
!                1: Classical Lorentzian GQER
!
!***********************************************************************
!
REAL,PARAMETER::  PIH=  8.673592583E-08,& ! 1/(3*(pi*hbar*c)**2)
                  PIHQ= 5.204155555E-08,& ! 1/(5*(pi*hbar*c)**2)
                  PI42=39.4784176       ! 4*pi**2
REAL::            SFCEE1,SFCEM1,SFCEE2,Q,QQ,TFIN,W,WPHEN,SLIM,x,FACTOR,ALPPL,ER0PL,WD,WDR,WR,FKs0,FNS,EFERMI,WWALL,Efinal
INTEGER::         I
integer::         ITYP
real::            EGAM,EINI

      SGAMMA=0.
      SFCEM1=0
      SFCEE2=0
      IF ((ITYP.GT.4).OR.(ITYP.LT.1).OR.(EGAM.LE.0.)) RETURN
!
!*****                       E1 component
!
      IF (ITYP.EQ.1) THEN
!
        IF     (NOPTE1.EQ.0) THEN   ! The single-particle approximation
          SGAMMA=DEG*EGAM**3
          RETURN
        ELSEIF (NOPTE1.EQ.1) THEN   ! Classical Lorentzian
          Q=0.
          DO I=1,NGIGE              ! loop over both GDR peaks
            QQ=SIG(I)*(EGAM*W0(I)**2/((EGAM**2-ER(I)**2)**2+(EGAM*W0(I))**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          RETURN
        ELSEIF (NOPTE1.EQ.11) THEN  ! Goriely tables
          SGAMMA = APSF(EGAM,1)
          SGAMMA = SGAMMA + EINI * PAR_E1(1) / (1.0+EXP(EGAM-PAR_E1(2))) 
          SGAMMA = SGAMMA * EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.50) THEN  ! Goriely tables
          SGAMMA = APSF(EGAM,1)*EGAM**3
          Q=0.
          DO I=1,NGIGE              ! loop over both GDR peaks
            QQ=SIG(I)*(EGAM*W0(I)**2/((EGAM**2-ER(I)**2)**2+(EGAM*W0(I))**2))
            Q=Q+QQ
          ENDDO
          SGAMMA = SGAMMA + PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.91) THEN  !SLO with exponential low energy enhancement as requested by Artemis
          Q=SIG(1)*EXP(-(EGAM-ER(1))*W0(1))
          DO I=2,NGIGE
            QQ=SIG(I)*(EGAM*W0(I)**2/((EGAM**2-ER(I)**2)**2+(EGAM*W0(I))**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          RETURN
        ELSEIF (NOPTE1.EQ.2) THEN   !ELO = GDER with E,T-dependent damping
          TFIN=TERM(EINI-EGAM)
          Q=0.
          DO I=1,NGIGE
            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2   !energy and temperature dependent width
            QQ=SIG(I)*W0(I)*(EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          RETURN
        ELSEIF (NOPTE1.EQ.3) THEN  !Empirical generalization of temperature
!          dependent damping according to Kopecky in Chrien model (EGLO)
          TFIN=TERM(EINI-EGAM)
          Q=0.
          DO I=1,NGIGE
            WPHEN=EK0+(1.-EK0)*(EGAM-EGZERO)/(ER(I)-EGZERO)
            W=WPHEN*W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2 !energy and temperature dependent width
            WPHENZ=EK0-(1.-EK0)*EGZERO/(ER(I)-EGZERO)
            SLIM=WPHENZ*FERMC*PI42*TFIN**2*W0(I)/ER(I)**5 !the non-zero limit at Egam-->0
            QQ=SIG(I)*W0(I)*(SLIM+EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          RETURN
        ELSEIF (NOPTE1.EQ.4) THEN   ! Pure Fermi liquid theory (Kadmenskij)
          TFIN=TERM(EINI-EGAM)          ! (no high energy approximation)
          Q=0.
          DO I=1,NGIGE
            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
            QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          RETURN
        ELSEIF (NOPTE1.EQ.5) THEN    !GLO as called nowdays, Pure Chrien model (similar to !OPT=3,but no Kadmenskij for low EGAM) viz Kopecky
          TFIN=TERM(EINI-EGAM)
          Q=0.
          DO I=1,NGIGE
            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2 !energy and temperature dependent width
            SLIM=FERMC*PI42*TFIN**2*W0(I)/ER(I)**5 !the non-zero limit at Egam-->0
            QQ=SIG(I)*W0(I)*(SLIM+EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          RETURN
        ELSEIF (NOPTE1.EQ.56) THEN  !GLO with constant temperature and low-energy enhancement (for Artemis)
          TFIN=PAR_E1(2)
          Q=SIG(1)*EXP(-W0(1)*(EGAM-ER(1)))
          DO I=2,NGIGE
            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
            SLIM=FERMC*PI42*TFIN**2*W0(I)/ER(I)**5
            QQ=SIG(I)*W0(I)*(SLIM+EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          RETURN
        ELSEIF (NOPTE1.EQ.58) THEN  ! GLO with constant temperature and pygmy resonance(s)
          TFIN=PAR_E1(2)
          Q=0.
          DO I=1,NLOWLOR
            QQ=SIG(I)*(EGAM*W0(I)**2/((EGAM**2-ER(I)**2)**2+(EGAM*W0(I))**2))
            Q=Q+QQ
          ENDDO
          DO I=NLOWLOR+1,NGIGE+NLOWLOR
            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
            SLIM=FERMC*PI42*TFIN**2*W0(I)/ER(I)**5
            QQ=SIG(I)*W0(I)*(SLIM+EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          RETURN
        ELSEIF (NOPTE1.EQ.6) THEN  !MGLO <-Empirical generalization of temperature dependent damping from EGLO(3)
          TFIN=TERM(EINI-EGAM)
          Q=0.
          DO I=1,NGIGE
            WPHEN=EK0+(1.-EK0)*(EGAM-EGZERO)/(ER(I)-EGZERO)
            W=WPHEN*W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2 !energy and temperature dependent width
            SLIM=FERMC*PI42*TFIN**2*W0(I)/ER(I)**5 !the non-zero limit at Egam-->0
            QQ=SIG(I)*W0(I)*(SLIM+EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          RETURN
        ELSEIF (NOPTE1.EQ.7) THEN   ! SMLO
          IF (EGAM.LE.EINI) THEN
            TFIN=SQRT((EINI-EGAM)/AMASS*10.0)
          ELSE
            TFIN=0.0
          ENDIF
          Q=0.
          DO I=1,NGIGE
            W=W0(I)*(EGAM*ER(I)+PI42*TFIN**2)/ER(I)**2
            QQ=SIG(I)*W*EGAM / (1.0-EXP(-EGAM/TFIN))/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2)
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.207) THEN  !Empirical generelization of temperature dependent damping according to Kopecky aplied to TD model (EELO)
          TFIN=TERM(EINI-EGAM)
          Q=0.
          DO I=1,NGIGE
            WPHEN=EK0+(1.-EK0)*(EGAM-EGZERO)/(ER(I)-EGZERO)
            W=WPHEN*W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2 !energy and temperature dependent width
            QQ=SIG(I)*W0(I)*(EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          RETURN
        ELSEIF (NOPTE1.EQ.8) THEN   ! Fermi liquid theory (Kadmenskij)
          TFIN=TERM(EINI-EGAM)      ! with 1st resonance of Lorentz. shape
          Q=0.
          QQ=SIG(1)*(EGAM*W0(1)**2/((EGAM**2-ER(1)**2)**2+(EGAM*W0(1))**2))
          Q=Q+QQ
          DO I=2,NGIGE
            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
            QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          RETURN
        ELSEIF (NOPTE1.EQ.9) THEN  !Our empirical generalization of temperature
!          dependent damping according to Kopecky in Chrien model (MGLO)
!          with the first resonance of Lorentzian type = pygmy
          TFIN=TERM(EINI-EGAM)
          Q=0.
          QQ=SIG(1)*(EGAM*W0(1)**2/((EGAM**2-ER(1)**2)**2+(EGAM*W0(1))**2))
          Q=Q+QQ
          DO I=2,NGIGE
            WPHEN=EK0+(1.-EK0)*(EGAM-EGZERO)/(ER(I)-EGZERO)
            W=WPHEN*W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2 !energy and temperature dependent width
            SLIM=FERMC*PI42*TFIN**2*W0(I)/ER(I)**5 !the non-zero limit at Egam-->0
            QQ=SIG(I)*W0(I)*(SLIM+EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          RETURN
        ELSEIF (NOPTE1.EQ.10) THEN   ! KMF for low energies
          TFIN=TERM(EINI-EGAM)       ! Mix KMF and BA for higher energies
          Q=0.
          IF (EGAM.LE.PAR_E1(1)) THEN
           DO I=1,NGIGE
            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
            QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
            Q=Q+QQ
           ENDDO
          ELSE
           x=(EGAM-PAR_E1(1))/(PAR_E1(2)-PAR_E1(1))         ! Admixture of BA to KMF
           IF (x.GT.1.) x=1.
           DO I=1,NGIGE
            QQ=SIG(I)*(EGAM*W0(I)**2/((EGAM**2-ER(I)**2)**2+(EGAM*W0(I))**2))
            Q=Q+x*QQ
            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
            QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
            Q=Q+(1.-x)*QQ
           ENDDO
          ENDIF
          SGAMMA=PIH*Q*EGAM**3
          RETURN
!
        ELSEIF (NOPTE1.EQ.31) THEN  ! Classical Lor.; suppressed for small EGAM
          Q=0.
          DO I=1,NGIGE              ! loop over both GDR peaks
            QQ=SIG(I)*(EGAM*W0(I)**2/((EGAM**2-ER(I)**2)**2+(EGAM*W0(I))**2))
            Q=Q+QQ
          ENDDO
          x=DIPSLP*EGAM+DIPZER 
          if (x.LT.PAR_E1(3)) x=PAR_E1(3)
          if (x.GT.1.0)    x=1.0
          SGAMMA=PIH*Q*EGAM**3*x
          RETURN
        ELSEIF (NOPTE1.EQ.34) THEN   ! Pure Fermi liquid theory (Kadmenskij)
!                                      suppressed for small EGAM
          TFIN=TERM(EINI-EGAM)          
          Q=0.
          DO I=1,NGIGE
            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
            QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
            Q=Q+QQ
          ENDDO
          x=DIPSLP*EGAM+DIPZER 
          if (x.LT.PAR_E1(3)) x=PAR_E1(3)
          if (x.GT.1.0)    x=1.0
          SGAMMA=PIH*Q*EGAM**3*x
          RETURN

        ELSEIF (NOPTE1.EQ.35) THEN   ! Pure Fermi liquid theory (Kadmenskij)
!                                      suppressed for small EGAM (34); no restriction for large Eg
          TFIN=TERM(EINI-EGAM)          
          Q=0.
          DO I=1,NGIGE
            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
            QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
            Q=Q+QQ
          ENDDO
          x=DIPSLP*EGAM+DIPZER 
          if (x.LT.PAR_E1(3)) x=PAR_E1(3)
!          if (x.GT.1.0)    x=1.0
          SGAMMA=PIH*Q*EGAM**3*x
          RETURN

        ELSEIF (NOPTE1.EQ.36) THEN  !TODO from which is this derived: EGLO (6 - incorrect); suppressed for small EGAM
          TFIN=TERM(EINI-EGAM)
          Q=0.
          DO I=1,NGIGE
            WPHEN=EK0+(1.-EK0)*(EGAM-EGZERO)/(ER(I)-EGZERO)
            W=WPHEN*W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
            SLIM=FERMC*PI42*TFIN**2*W0(I)/ER(I)**5
            QQ=SIG(I)*W0(I)*(SLIM+EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          x=DIPSLP*EGAM+DIPZER 
          if (x.LT.PAR_E1(3)) x=PAR_E1(3)
          if (x.GT.1.0)    x=1.0
          SGAMMA=PIH*Q*EGAM**3*x
          RETURN
        ELSEIF (NOPTE1.EQ.37) THEN  !EELO; suppressed for small EGAM
          TFIN=TERM(EINI-EGAM)
          Q=0.
          DO I=1,NGIGE
            WPHEN=EK0+(1.-EK0)*(EGAM-EGZERO)/(ER(I)-EGZERO)
            W=WPHEN*W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
            QQ=SIG(I)*W0(I)*(EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          x=DIPSLP*EGAM+DIPZER 
          if (x.LT.PAR_E1(3)) x=PAR_E1(3)
          if (x.GT.1.0)    x=1.0
          SGAMMA=PIH*Q*EGAM**3*x
          RETURN
! TODO change the PAR_E1(1) - talk to MK about the logic of this
        ELSEIF (NOPTE1.EQ.73) THEN  ! EGLO(3) with constant T
          TFIN=PAR_E1(2)
          Q=0.
          DO I=1,NGIGE
            WPHEN=EK0+(1.-EK0)*(EGAM-EGZERO)/(ER(I)-EGZERO)
            W=WPHEN*W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
            WPHENZ=EK0-(1.-EK0)*EGZERO/(ER(I)-EGZERO)
            SLIM=WPHENZ*FERMC*PI42*TFIN**2*W0(I)/ER(I)**5
            QQ=SIG(I)*W0(I)*(SLIM+EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.74) THEN   ! KMF with constant T
          TFIN=PAR_E1(2)
          Q=0.
          DO I=1,NGIGE
           W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
           QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
           Q=Q+QQ
          ENDDO
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3    
        RETURN
        ELSEIF (NOPTE1.EQ.75) THEN  !GLO with constant T
          TFIN=PAR_E1(2)
          Q=0.
          DO I=1,NGIGE
            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
            SLIM=FERMC*PI42*TFIN**2*W0(I)/ER(I)**5
            QQ=SIG(I)*W0(I)*(SLIM+EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3
          RETURN
        ELSEIF (NOPTE1.EQ.76) THEN  ! MGLO(6) with constant T
          TFIN=PAR_E1(2)
          Q=0.
          DO I=1,NGIGE
            WPHEN=EK0+(1.-EK0)*(EGAM-EGZERO)/(ER(I)-EGZERO)
            W=WPHEN*W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
            SLIM=FERMC*PI42*TFIN**2*W0(I)/ER(I)**5
            QQ=SIG(I)*W0(I)*(SLIM+EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3
          RETURN
        ELSEIF (NOPTE1.EQ.77) THEN  ! MGLO(6) with constant T and Lorentzian LLR pygmy
          TFIN=PAR_E1(2)
          Q=0.
          DO I=1,NLOWLOR
            QQ=SIG(I)*(EGAM*W0(I)**2/((EGAM**2-ER(I)**2)**2+(EGAM*W0(I))**2))
            Q=Q+QQ
          ENDDO
          DO I=NLOWLOR+1,NGIGE+NLOWLOR
            WPHEN=EK0+(1.-EK0)*(EGAM-EGZERO)/(ER(I)-EGZERO)
            W=WPHEN*W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
            SLIM=FERMC*PI42*TFIN**2*W0(I)/ER(I)**5
            QQ=SIG(I)*W0(I)*(SLIM+EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3
          RETURN
        ELSEIF (NOPTE1.EQ.66) THEN  !MGLO <-Empirical generalization of temperature dependent damping from EGLO(3)
! NLOWLOR resonances of Lorentzian shape pygmy, then NGIGE resonances of MGLO shape
          TFIN=TERM(EINI-EGAM)
          Q=0.
          DO I=1,NLOWLOR
            QQ=SIG(I)*(EGAM*W0(I)**2/((EGAM**2-ER(I)**2)**2+(EGAM*W0(I))**2))
            Q=Q+QQ
          ENDDO
          DO I=NLOWLOR+1,NGIGE+NLOWLOR
            WPHEN=EK0+(1.-EK0)*(EGAM-EGZERO)/(ER(I)-EGZERO)
            W=WPHEN*W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2 !energy and temperature dependent width
            SLIM=FERMC*PI42*TFIN**2*W0(I)/ER(I)**5 !the non-zero limit at Egam-->0
            QQ=SIG(I)*W0(I)*(SLIM+EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          RETURN
        ELSEIF (NOPTE1.EQ.39) THEN   ! Pure Fermi liquid theory (Kadmenskij)
!                                      suppressed for small EGAM - const. T !!!
          TFIN=PAR_E1(2)
          Q=0.
          DO I=1,NGIGE
            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
            QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
            Q=Q+QQ
          ENDDO
          x=DIPSLP*EGAM+DIPZER 
          if (x.LT.PAR_E1(3)) x=PAR_E1(3)
          if (x.GT.1.0)    x=1.0
          SGAMMA=PIH*Q*EGAM**3*x
          RETURN
        ELSEIF (NOPTE1.EQ.64) THEN   ! Pure Fermi liquid theory (Kadmenskij)
!                                      suppressed for small EGAM
!                                    ! 1st resonance of Lorentzian shape
          TFIN=TERM(EINI-EGAM)          
          Q=0.
          QQ=SIG(1)*(EGAM*W0(1)**2/((EGAM**2-ER(1)**2)**2+(EGAM*W0(1))**2))
          Q=Q+QQ
          DO I=2,NGIGE
            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
            QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
            Q=Q+QQ
          ENDDO
          x=DIPSLP*EGAM+DIPZER 
          if (x.LT.PAR_E1(3)) x=PAR_E1(3)
          if (x.GT.1.0)    x=1.0
          SGAMMA=PIH*Q*EGAM**3*x
          RETURN
        ELSEIF (NOPTE1.EQ.12) THEN   ! MLO1 (Original Plujko)
!         !tato procedura nepostihuje mozna uplne vsechny pripady (hodnoty parametru),
!         !ktere mohou podle 'Plujkovy teorie' nastat 
!         !urcite neni zakomponovana jeho moznost KEYSET=2 (je jine ALPPL) 
!         !hodnoty nekterych 'promennych' nastaveny zde
!
          TFIN=TERM(EINI-EGAM)
          FACTOR=1.0                        ! Recommended by Plujko
          ALPPL=185.659                     ! 4*pi^2*alpha_free - from Plujko
          ALPPL=ALPPL/FACTOR
          ER0PL=SMFREQ()
!
          Q=0.
          DO I=1,NGIGE
            WD=EINI*ER(I)/ALPPL
            WDR=ER(I)**2/ALPPL
            WR=2*WDR*(ER(I)**2+ER0PL**2)/((ER(I)**2-ER0PL**2)**2+4*(WDR*ER(I))**2)
            W=2*WD*W0(i)/WR*(ER(I)**2+ER0PL**2)/((ER(I)**2-ER0PL**2)**2+4*(WD*EGAM)**2)

            IF (TFIN.NE.0) THEN
             QQ=SIG(I)*W0(I)/(1-exp(-EGAM/TFIN))*(EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            ELSE
             QQ=SIG(I)*W0(I)*(EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            ENDIF
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          RETURN
        ELSEIF (NOPTE1.EQ.13) THEN   ! MLO1 (Original Plujko); Jina T!!!!
!         !tato procedura nepostihuje mozna uplne vsechny pripady (hodnoty parametru),
!         !ktere mohou podle 'Plujkovy teorie' nastat 
!         !urcite neni zakomponovana jeho moznost KEYSET=2 (je jine ALPPL) 
!         !hodnoty nekterych 'promennych' nastaveny zde
!
          TFIN=TERMDILG(EINI-EGAM)
          FACTOR=1.0                        ! Recommended by Plujko
          ALPPL=185.659                     ! 4*pi^2*alpha_free - from Plujko
          ALPPL=ALPPL/FACTOR
          ER0PL=SMFREQ()
!
          Q=0.
          DO I=1,NGIGE
            WD=EINI*ER(I)/ALPPL
            WDR=ER(I)**2/ALPPL
            WR=2*WDR*(ER(I)**2+ER0PL**2)/((ER(I)**2-ER0PL**2)**2+4*(WDR*ER(I))**2)
            W=2*WD*W0(i)/WR*(ER(I)**2+ER0PL**2)/((ER(I)**2-ER0PL**2)**2+4*(WD*EGAM)**2)

            IF (TFIN.NE.0) THEN
             QQ=SIG(I)*W0(I)/(1-exp(-EGAM/TFIN))*(EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            ELSE
             QQ=SIG(I)*W0(I)*(EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            ENDIF
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          RETURN
!
        ELSEIF (NOPTE1.EQ.14) THEN   ! MLO - RIPL2 (some errors)
!         !it is not clear if there is a factor EINI or (EGAM+a*TFIN**2) 
!          or (EGAM+a*TINI**2) in the width WD 
!
          TFIN=TERM(EINI-EGAM)
          TINI=TERM(EINI)
!          TINI=TERM1(EINI)
          ALPPL=0.00542              !in RIPL2 error in the order (5.42e-4)
          FACTOR=1.0                 !All parameters recommended by Plujko
          ALPPL=ALPPL/FACTOR
          FKs0=0.3                          
          FNS=1.0                           
!
          Q=0.
          DO I=1,NGIGE
            WD=EINI*ER(I)*ALPPL
!            WD=(EGAM+UINI)*ER(I)*ALPPL
            FKR=(W0(I)-ALPPL*ER(I)**2)/WWALL1()
            IF (EGAM.LE.(2*ER(I))) THEN
             FKS1=FKR+(FKs0-FKR)*ABS((EGAM-W0(I))/W0(I))**FNS             
            ELSE
             FKS1=FKs0
            ENDIF
            W=WD+FKS1*WWALL1()
            IF (TFIN.NE.0) THEN
             QQ=SIG(I)*W0(I)/(1-exp(-EGAM/TFIN))*(EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            ELSE
             QQ=SIG(I)*W0(I)*(EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            ENDIF
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          RETURN
!
        ELSEIF (NOPTE1.EQ.61) THEN   ! MLO - RIPL2 
!         !it is not clear if there is a factor EINI or (EGAM+a*TFIN**2) 
!          or (EGAM+a*TINI**2) in the width WD - here TFIN
!
          TFIN=TERM(EINI-EGAM)
!          TINI=TERM(EINI)
!          ALPPL=0.00542              !in RIPL2 error in the order (5.42e-4)
          ALPPL=0.000542              !in RIPL2 error in the order (5.42e-4)
          FACTOR=1.0                 !All parameters recommended by Plujko
          ALPPL=ALPPL/FACTOR
          FKs0=0.3                          
          FNS=1.0                           
!
          Q=0.
          DO I=1,NGIGE
            WD=(EGAM+EFEXC(EINI))*ER(I)*ALPPL
            FKR=(W0(I)-ALPPL*ER(I)**2)/WWALL1()
            IF (EGAM.LE.(2*ER(I))) THEN
             FKS1=FKR+(FKs0-FKR)*ABS((EGAM-ER(I))/ER(I))**FNS             
            ELSE
             FKS1=FKs0
            ENDIF
            W=WD+FKS1*WWALL1()
            IF (TFIN.NE.0) THEN
             QQ=SIG(I)*W0(I)/(1-exp(-EGAM/TFIN))*(EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            ELSE
             QQ=SIG(I)*W0(I)*(EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            ENDIF
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          RETURN
        ELSEIF (NOPTE1.EQ.62) THEN   ! MLO - RIPL2 
!         !it is not clear if there is a factor EINI or (EGAM+a*TFIN**2) 
!          or (EGAM+a*TINI**2) in the width WD - here TINI
!
          TFIN=TERM(EINI-EGAM)
!          ALPPL=0.00542              !in RIPL2 error in the order (5.42e-4)
          ALPPL=0.000542              !in RIPL2 error in the order (5.42e-4)
          FACTOR=1.0                 !All parameters recommended by Plujko
          ALPPL=ALPPL/FACTOR
          FKs0=0.3                          
          FNS=1.0                           
!
          Q=0.
          DO I=1,NGIGE
            WD=(EGAM+EINI)*ER(I)*ALPPL
            FKR=(W0(I)-ALPPL*ER(I)**2)/WWALL1()
            IF (EGAM.LE.(2*ER(I))) THEN
             FKS1=FKR+(FKs0-FKR)*ABS((EGAM-ER(I))/ER(I))**FNS             
            ELSE
             FKS1=FKs0
            ENDIF
            W=WD+FKS1*WWALL1()
            IF (TFIN.NE.0) THEN
             QQ=SIG(I)*W0(I)/(1-exp(-EGAM/TFIN))*(EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            ELSE
             QQ=SIG(I)*W0(I)*(EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            ENDIF
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          RETURN
        ELSEIF (NOPTE1.EQ.67) THEN   !GFL - Mughabghab+Dunford - formulation from RIPL2 
!
          S2=217.16/AMASS/AMASS      !global parametrization, see RIPL2 
          CDQ=1.05
          TFIN=TERM(EINI-EGAM)
          BETA2 = PAR_E1(3)*PAR_E1(3)    !square of deformation
!
          Q=0.
          DO I=1,NGIGE
            WDQ0= CDQ*SQRT(BETA2*ER(I)*ER(I)+ER(I)*S2)
            WDQ = CDQ*SQRT(BETA2*EGAM*EGAM+EGAM*S2)
            WC  = (W0(I)-WDQ0)/ER(I)/ER(I)*(EGAM**2+PI42*TFIN**2)
            W   = WDQ + WC 
            QQ=SIG(I)*W0(I)*FERMC*ER(I)*W/((EGAM**2-ER(I)**2)**2+FERMC*(EGAM*W)**2)
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          RETURN
        ELSEIF (NOPTE1.EQ.68) THEN   ! Mughabghab+Dunford - original formulation (PLB487,155)  
!
          S2=217.16/AMASS/AMASS      !global parametrization, see RIPL2 
          CDQ=1.05
          TFIN=TERM(EINI-EGAM)
          BETA2 = PAR_E1(3)*PAR_E1(3)    !square of deformation
!
          Q=0.
          DO I=1,NGIGE
            WDQ0= CDQ*SQRT(BETA2*ER(I)*ER(I)+ER(I)*S2)
            WDQ = CDQ*SQRT(BETA2*EGAM*EGAM+EGAM*S2)
            WM  = (W0(I)-WDQ0)/ER(I)/ER(I)*(EGAM**2+PI42*TFIN**2)
            W   = WDQ + WM 
            QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          RETURN
        ELSEIF (NOPTE1.EQ.69) THEN   ! Goriely (PLB436,10) - in RIPL2
          TFIN=TERM(EINI-EGAM)
          Q=0.
          DO I=1,NGIGE
            W=W0(I)*FERMC*(EGAM**2+PI42*TFIN**2)/ER(I)/EGAM !energy and temperature dependent width
            QQ=SIG(I)*W0(I)*(EGAM*W/((EGAM**2-ER(I)**2)**2+EGAM**2*W0(I)*W))
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          RETURN
!
!
        ELSEIF (NOPTE1.EQ.16) THEN   ! MLO2 (Original Plujko)
!         !tato procedura nepostihuje mozna uplne vsechny pripady (hodnoty parametru),
!         !ktere mohou podle 'Plujkovy teorie' nastat 
!         !urcite neni zakomponovana jeho moznost KEYSET=2 (je jine ALPPL) 
!         !hodnoty nekterych 'promennych' nastaveny zde
!
          TFIN=TERM(EINI-EGAM)
          BC=1.0                            ! Recommended by Plujko
          FACTOR=1.0                        ! Recommended by Plujko
          ALPPL=185.659                     ! 4*pi^2*alpha_free - from Plujko
          ALPPL=ALPPL/FACTOR
          EFERMI=37.                        ! Recommended by Plujko 
          WWALL=6.857764*sqrt(EFERMI)/R0()  ! origin of the const. is a puzzle for me
          FKs0=0.3                          ! Recommended by Plujko
          FNS=1.0                           ! Recommended by Plujko
!
          Q=0.
          DO I=1,NGIGE
            WD=EINI*ER(I)/ALPPL
            IF (NGIGE.EQ.1) THEN
              W=WD+FKS(EGAM,FNS,FKs0,WWALL,ALPPL)*WWALL
            ELSE
              W=WD+(W0(I)-ER(I)**2/ALPPL)*FKS(EGAM,FNS,FKs0,WWALL,ALPPL)/FKS(ER(I),FNS,FKs0,WWALL,ALPPL)
            ENDIF
            IF (TFIN.NE.0) THEN
             QQ=SIG(I)*W0(I)/(1-exp(-EGAM/TFIN))*(EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            ELSE
             QQ=SIG(I)*W0(I)*(EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            ENDIF
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          RETURN
!
        ELSEIF (NOPTE1.EQ.17) THEN   ! MLO3 (Original Plujko)
!         !tato procedura nepostihuje mozna uplne vsechny pripady (hodnoty parametru),
!         !ktere mohou podle 'Plujkovy teorie' nastat 
!         !urcite neni zakomponovana jeho moznost KEYSET=2 (je jine ALPPL) 
!         !hodnoty nekterych 'promennych' nastaveny zde
!
          TFIN=TERM(EINI-EGAM)
          BC=1.0                            ! Recommended by Plujko
          FACTOR=1.0                        ! Recommended by Plujko
          ALPPL=185.659                     ! 4*pi^2*alpha_free - from Plujko
          ALPPL=ALPPL/FACTOR
          EFERMI=37.                        ! Recommended by Plujko 
          WWALL=6.857764*sqrt(EFERMI)/R0()  ! origin of the const. is a puzzle for me
          FKs0=0.3                          ! Recommended by Plujko
          FNS=1.0                           ! Recommended by Plujko
!
          Q=0.
          DO I=1,NGIGE
            WD=(EGAM**2+PI42*TFIN**2)/ALPPL
            IF (NGIGE.EQ.1) THEN
              W=WD+FKS(EGAM,FNS,FKs0,WWALL,ALPPL)*WWALL
            ELSE
              W=WD+(W0(I)-ER(I)**2/ALPPL)*FKS(EGAM,FNS,FKs0,WWALL,ALPPL)/FKS(ER(I),FNS,FKs0,WWALL,ALPPL)
            ENDIF
            IF (TFIN.NE.0) THEN
             QQ=SIG(I)*W0(I)/(1-exp(-EGAM/TFIN))*(EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            ELSE
             QQ=SIG(I)*W0(I)*(EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            ENDIF
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          RETURN
!
        ELSEIF (NOPTE1.EQ.18) THEN   ! Sirotkin (1) - with KMF damping
          TFIN=TERM(EINI-EGAM)
          Q=0.
          DO I=1,NGIGE
            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
            O=FERMC*ER(I)
            IF (TFIN.NE.0) THEN
            QQ=SIG(I)*W0(I)*W*EGAM*(EGAM**2+O**2)/(EGAM**2-ER(I)**2)**2/O**2/(1-exp(-EGAM/TFIN))
            ELSE
            QQ=SIG(I)*W0(I)*W*EGAM*(EGAM**2+O**2)/(EGAM**2-ER(I)**2)**2/O**2
            ENDIF
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          RETURN
        ELSEIF (NOPTE1.EQ.19) THEN   ! Sirotkin (2) - with damping from Sir
          TFIN=TERM(EINI-EGAM)
          Q=0.
          DO I=1,NGIGE
!            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
            W=0.666*4./35.*(EGAM**2+PI42*TFIN**2)
            O=FERMC*ER(I)
            IF (TFIN.NE.0) THEN
            QQ=SIG(I)*W0(I)*W*EGAM*(EGAM**2+O**2)/(EGAM**2-ER(I)**2)**2/O**2/(1-exp(-EGAM/TFIN))
            ELSE
            QQ=SIG(I)*W0(I)*W*EGAM*(EGAM**2+O**2)/(EGAM**2-ER(I)**2)**2/O**2
            ENDIF
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          RETURN
        ELSEIF (NOPTE1.EQ.20) THEN   ! Pure Fermi liquid theory (Kadmenskij)
          TFIN=TERM(EINI-EGAM)       ! with theoretical damping from Sir
          Q=0.
          DO I=1,NGIGE
!            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
            W=0.666*4./35.*(EGAM**2+PI42*TFIN**2)
            QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          RETURN
!
        ELSEIF (NOPTE1.EQ.51) THEN    ! Fermi liquid theory (Kadmenskij)
                                      ! without dependence on T !!!
          Q=0.
          DO I=1,NGIGE
            W=W0(I)*(EGAM**2)/ER(I)**2
            QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          RETURN
        ELSEIF (NOPTE1.EQ.52) THEN  !Modified empirical generalization
!          of temperature dependent damping according to Kopecky
!          in Chrien model (EGLO); Temperature not assumed in damping
!          width - see above
          TFIN=TERM(EINI-EGAM)
          Q=0.
          DO I=1,NGIGE
            WPHEN=EK0+(1.-EK0)*(EGAM-EGZERO)/(ER(I)-EGZERO)
            W=WPHEN*W0(I)*(EGAM**2)/ER(I)**2
!          !     ^ .......... energy dependent width
            SLIM=FERMC*PI42*TFIN**2*W0(I)/ER(I)**5
!          !     ^ ...... the non-zero limit at Egam-->0
            QQ=SIG(I)*W0(I)*(SLIM+EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          RETURN
        ELSEIF (NOPTE1.EQ.53) THEN  !Modified empirical generelization
!        of temperature dependent damping according to Kopecky aplied
!        to TD model (EELO); temperature dependence removed - see above
          Q=0.
          DO I=1,NGIGE
            WPHEN=EK0+(1.-EK0)*(EGAM-EGZERO)/(ER(I)-EGZERO)
            W=WPHEN*W0(I)*(EGAM**2)/ER(I)**2 !energy and temperature dependent width
            QQ=SIG(I)*W0(I)*(EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          RETURN
!
        ELSEIF (NOPTE1.EQ.41) THEN   !!!!!! Oslo KMF for Yb
          TFIN=TERMOSLO(EINI-EGAM)        
          Q=0.
          DO I=1,NGIGE
           W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
           QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
           Q=Q+QQ
          ENDDO
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3    
        RETURN
        ELSEIF (NOPTE1.EQ.42) THEN   !!!!!! Oslo KMF - our temperature
          TFIN=TERM(EINI-EGAM)        
          Q=0.
          DO I=1,NGIGE
           W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
           QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
           Q=Q+QQ
          ENDDO
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3    
        RETURN
        ELSEIF (NOPTE1.EQ.57) THEN   ! Fermi liquid theory (Kadmenskij)
          TFIN=TERM(EINI-EGAM)      ! with 1st resonance of Lorentz. shape
          Q=0.
          DO I=1,2
            QQ=SIG(I)*(EGAM*W0(I)**2/((EGAM**2-ER(I)**2)**2+(EGAM*W0(I))**2))
            Q=Q+QQ
          ENDDO
          DO I=3,NGIGE
            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
            QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
            Q=Q+QQ
          ENDDO
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3
        RETURN
        ELSEIF (NOPTE1.EQ.43) THEN   !!!!!! Oslo KMF; including softpole
          TFIN=TERMOSLO(EINI-EGAM)        
          Q=0.
          DO I=1,NGIGE
           W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
           QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
           Q=Q+QQ
          ENDDO
          Q=Q+PAR_E1(2)/(EGAM**PAR_E1(3))
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3    
        RETURN
        ELSEIF (NOPTE1.EQ.44) THEN   !!!!!! Oslo KMF - our temperature; including softpole
          TFIN=TERM(EINI-EGAM)        
          Q=0.
          DO I=1,NGIGE
           W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
           QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
           Q=Q+QQ
          ENDDO
          Q=Q+PAR_E1(2)/(EGAM**PAR_E1(3))
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3    
        RETURN
        ELSEIF (NOPTE1.EQ.45) THEN   !!!!!! Oslo KMF - const T 
          TFIN=TERMOSLO(BN-EGAM)        
          Q=0.
          DO I=1,NGIGE
           W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
           QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
           Q=Q+QQ
          ENDDO
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3    
        RETURN
        ELSEIF (NOPTE1.EQ.46) THEN   !!!!!! Oslo KMF; including softpole; TODO OUTDATED
          TFIN=PAR_E1(2)       !Temperature independent PSF  
          Q=0.
          DO I=1,NGIGE
           W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
           QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
           Q=Q+QQ
          ENDDO
          Q=Q+PAR_E1(3)/(EGAM**PAR_E1(4))
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3    
        RETURN
        ELSEIF (NOPTE1.EQ.47) THEN   !!!!!! Oslo KMF; including softpole
          TFIN=TERMOSLO(EINI-EGAM)     !!!! const below PAR_E1(1)   
          Q=0.
          EGAMOLD=EGAM
          IF (EGAM.LT.PAR_E1(1)) EGAM=PAR_E1(1)
          DO I=1,NGIGE
           W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
           QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
           Q=Q+QQ
          ENDDO
          Q=Q+PAR_E1(2)/(EGAM**PAR_E1(3))
          EGAM=EGAMOLD
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3    
        RETURN
        ELSEIF (NOPTE1.EQ.48) THEN   !!!!!! Oslo KMF - T independent; including softpole
          TFIN=TERMOSLO(BN-EGAM)     !!!! const below PAR_E1(1)   
          Q=0.
          EGAMOLD=EGAM
          IF (EGAM.LT.PAR_E1(1)) EGAM=PAR_E1(1)
          DO I=1,NGIGE
           W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
           QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
           Q=Q+QQ
          ENDDO
          Q=Q+PAR_E1(2)/(EGAM**PAR_E1(3))
          EGAM=EGAMOLD
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3    
        RETURN
        ELSEIF (NOPTE1.EQ.49) THEN   !!!!!! Oslo KMF - T independent; including softpole
          TFIN=TERMOSLO(BN-EGAM)     !!!! zero below PAR_E1(1)   
          Q=0.
          DO I=1,NGIGE
           W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
           QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
           Q=Q+QQ
          ENDDO
          Q=Q+PAR_E1(2)/(EGAM**PAR_E1(3))
          IF (EGAM.LT.PAR_E1(1)) Q=0.0
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3    
        RETURN
        ELSEIF (NOPTE1.EQ.54) THEN   !!!!!! Oslo KMF - T constant, f_py in M1
          Q=0.
          DO I=1,NGIGE
           W=(EGAM**2+PI42*EK0**2)
           QQ=FERMC*SIG(I)*W0(I)**2*W/ER(I)/(EGAM**2-ER(I)**2)**2
           Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3    
        RETURN
        ELSEIF (NOPTE1.EQ.55) THEN   !!!!!! Oslo KMF - T constant, f_py in E1
          Q=0.
          Q=SIG(1)*(EGAM*W0(1)**2/((EGAM**2-ER(1)**2)**2+(EGAM*W0(1))**2))
          DO I=2,NGIGE
           W=(EGAM**2+PI42*EK0**2)
           QQ=FERMC*SIG(I)*W0(I)**2*W/ER(I)/(EGAM**2-ER(I)**2)**2
           Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3    
        RETURN
        ELSEIF (NOPTE1.EQ.79) THEN   !!!!!! Voinov Mo LLR; KMF with constant T and 1st lorentzian resonance, equal to 55 with diferently defined parameters on input (uses kappa)
          TFIN=PAR_E1(2)        
          Q=0.
          Q=SIG(1)*(EGAM*W0(1)**2/((EGAM**2-ER(1)**2)**2+(EGAM*W0(1))**2))/PAR_E1(1)  ! Lorentzian LLR
          DO I=2,NGIGE
           W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
           QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
           Q=Q+QQ
          ENDDO
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3    
        RETURN
        ELSEIF (NOPTE1.EQ.83) THEN  ! EGLO(3) with constant T + Lorentzian pygmy (Oslo - actinides)
          TFIN=PAR_E1(2)
          Q=SIG(1)*(EGAM*W0(1)**2/((EGAM**2-ER(1)**2)**2+(EGAM*W0(1))**2)) ! pygmy
          DO I=2,NGIGE
            WPHEN=EK0+(1.-EK0)*(EGAM-EGZERO)/(ER(I)-EGZERO)
            W=WPHEN*W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
            WPHENZ=EK0-(1.-EK0)*EGZERO/(ER(I)-EGZERO)
            SLIM=WPHENZ*FERMC*PI42*TFIN**2*W0(I)/ER(I)**5
            QQ=SIG(I)*W0(I)*(SLIM+EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.25) THEN    !Pure Chrien+Kopecky model but Toshihiko's interpretation of temperature, NOT exactly his temperature!!!
          TFIN=TERM(BN-EGAM)
          Q=0.
          DO I=1,NGIGE
            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
          !     ^ .......... energy and temperature dependent width
            SLIM=FERMC*PI42*TFIN**2*W0(I)/ER(I)**5
          !     ^ ...... the non-zero limit at Egam-->0
            QQ=SIG(I)*W0(I)*(SLIM+EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.28) THEN  ! CoH3 GL model (= Kawano's implementation of GLO model)
          TFIN=TERM_TK(BN-EGAM)
          Q=0.
          DO I=1,NGIGE
            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
          !     ^ .......... energy and temperature dependent width
            SLIM=FERMC*PI42*TFIN**2*W0(I)/ER(I)**5
          !     ^ ...... the non-zero limit at Egam-->0
            QQ=SIG(I)*W0(I)*(SLIM+EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ENDIF
      ENDIF
!
!*****                       M1 component
      IF ((ITYP.EQ.2).OR.( ITYP.EQ.3)) THEN
!
        IF (NOPTM1.EQ.0) THEN       ! The single-particle approximation
          SGAMMA=DMG*EGAM*EGAM*EGAM
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN
        ELSEIF (NOPTM1.EQ.1) THEN   ! Classical Lorentzian
          Q=0.
          DO I=1,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN
        ELSEIF (NOPTM1.EQ.3) THEN   ! Scissors mode is build up only on states
          Q=0.                      ! with Efinal lower than PAR_M1(1)
          DO I=2,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          Efinal=Eini-Egam
          IF (Efinal.LE.PAR_M1(1)) THEN
            QQ=SIGM(1)*EGAM*WM0(1)**2/((EGAM**2-ERM(1)**2)**2+(EGAM*WM0(1))**2)
            Q=Q+QQ
          ENDIF
          SGAMMA=PIH*Q*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN
        ELSEIF (NOPTM1.EQ.4) THEN   ! Classical Lorentzian on backgroung
          Q=0.                      ! which is described by SP up to PAR_M1(1)
          DO I=1,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          IF (EGAM.LE.PAR_M1(1)) THEN
            SGAMMA=(PIH*Q+DMG)*EGAM**3
          ELSE
            SGAMMA=PIH*Q*EGAM**3
          ENDIF
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN
        ELSEIF (NOPTM1.EQ.5) THEN  ! SP on states
          Q=0.                     ! with Efinal lower than PAR_M1(1)
          DO I=1,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          Efinal=EINI-EGAM
          IF (Efinal.LE.PAR_M1(1)) THEN
            SGAMMA=(PIH*Q+DMG)*EGAM**3
          ELSE
            SGAMMA=PIH*Q*EGAM**3
          ENDIF
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN
        ELSEIF (NOPTM1.EQ.7) THEN   ! SMLO = Classical Lorentzian on backgroung
          Q=0.                      ! which is described by exponencial fcion
          DO I=1,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          SGAMMA=(PIH*Q+PAR_M1(1)*EXP(-PAR_M1(2)*EGAM)*(1.0+PAR_M1(3)*EGAM**3))*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN
        ELSEIF (NOPTM1.EQ.6) THEN   ! Classical Lorentzian ...
          Q=0.
          DO I=1,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          Efinal=EINI-EGAM
          SGAMMA=(PIH*Q+DMG*exp(-Efinal/PAR_M1(1)))*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN
!
        ELSEIF (NOPTM1.EQ.8) THEN   ! Scissors mode is NOT build on
          Q=0.                      ! some of the terminal states
          DO I=1,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          Efinal=Eini-Egam
          DO i=1,numlev
            IF((Efinal.EQ.LVL_ENERGY(i)).AND.(1.EQ.LVL_CLASS(i))) then
              QQ=SIGM(1)*EGAM*WM0(1)**2/((EGAM**2-ERM(1)**2)**2+(EGAM*WM0(1))**2)
              Q=Q-QQ  
            ENDIF
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN
        ELSEIF (NOPTM1.EQ.9) THEN   ! Scissors mode is build up only on
          Q=0.                      ! some of the terminal states
          DO I=2,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          Efinal=Eini-Egam
          IF ((Efinal.EQ.LVL_ENERGY(1)).OR.(Efinal.EQ.LVL_ENERGY(2)).OR.(Efinal.EQ.LVL_ENERGY(3))) THEN
            QQ=SIGM(1)*EGAM*WM0(1)**2/((EGAM**2-ERM(1)**2)**2+(EGAM*WM0(1))**2)
            Q=Q+QQ
          ENDIF
          SGAMMA=PIH*Q*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN
        ELSEIF (NOPTM1.EQ.11) THEN   ! Goriely tables
          SGAMMA = APSF(EGAM,2)
          SGAMMA = SGAMMA + PAR_M1(1) * EXP(-PAR_M1(2)*EGAM) *(1.0+PAR_M1(3)*EGAM**PAR_M1(4))
          SGAMMA = SGAMMA * EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN
        ELSEIF (NOPTM1.EQ.91) THEN   ! Scissors mode is build up only on states
          Q=0.                      ! with Efinal lower than PAR_M1(1)
          DO I=2,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          Efinal=Eini-Egam
          IF (Efinal.LE.PAR_M1(1)) THEN
            QQ=SIGM(1)*EGAM*WM0(1)**2/((EGAM**2-ERM(1)**2)**2+(EGAM*WM0(1))**2)
            Q=Q+QQ
          ENDIF
          SGAMMA=(PIH*Q+DMG)*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN
        ELSEIF (NOPTM1.EQ.13) THEN   ! Lorentzian - energy of scissors
          Q=0.                      ! resonance (1) depends on final state
          ERORIG=ERM(1)
          DO I=1,NGIGM
            IF (I.EQ.1) ERM(1)=-0.1*(EINI-EGAM)+3.0
            QQ=SIGM(I)*EGAM*WM0(I)**2/((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          ERM(1)=ERORIG
          SGAMMA=PIH*Q*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN

        ELSEIF (NOPTM1.EQ.14) THEN   ! Width of Scissors resonance depends on Efin
          Efinal=Eini-Egam           ! => total strength increases
          WM0SC=WM0(1)
          WM0(1)=WM0(1)+PAR_M1(1)*Efinal
          Q=0.
          DO I=1,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEM1=SGAMMA
          WM0(1)=WM0SC               !restoring the original value
          IF (ITYP.EQ.3) RETURN
        ELSEIF (NOPTM1.EQ.15) THEN   ! Width of Scissors resonance depends on Efin
          Efinal=Eini-Egam           ! Total strength should remain 
          WM0SC=WM0(1)
          WM0(1)=WM0(1)+PAR_M1(1)*Efinal
          SIGMSC=SIGM(1)
          SIGM(1)=WM0SC*SIGMSC/WM0(1)
          Q=0.
          DO I=1,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEM1=SGAMMA
          WM0(1)=WM0SC               !restoring the original values
          SIGM(1)=SIGMSC
          IF (ITYP.EQ.3) RETURN
        ELSEIF (NOPTM1.EQ.16) THEN   ! cross section of Scissors resonance depends on Efin
          Efinal=Eini-Egam           
          SIGMSC=SIGM(1)
          SIGM(1)=SIGM(1)+PAR_M1(1)*Efinal
          Q=0.
          DO I=1,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEM1=SGAMMA
          WM0(1)=WM0SC               !restoring the original values
          IF (ITYP.EQ.3) RETURN
!
        ELSEIF (NOPTM1.EQ.21) THEN   ! Oslo - factor k
          Q=0.
          DO I=1,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          SGAMMA=PAR_M1(1)*PIH*Q*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN
        ELSEIF (NOPTM1.EQ.22) THEN   ! Oslo - factor k + softpole
          Q=0.
          DO I=1,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          Q=Q+PAR_M1(2)/(EGAM**PAR_M1(3))
          SGAMMA=PAR_M1(1)*PIH*Q*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN
        ELSEIF (NOPTM1.EQ.23) THEN   ! Oslo - factor k + softpole with a const below PAR_M1(1)
          Q=0.
          EGAMOLD=EGAM
          IF (EGAM.LT.PAR_M1(1)) EGAM=PAR_M1(1)
          Q=SIGM(1)*EXP(-WM0(1)*(EGAM-ERM(1)))
          DO I=2,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          EGAM=EGAMOLD
          SGAMMA=PIH*Q*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN
        ELSEIF (NOPTM1.EQ.24) THEN   ! SP - factor k + softpole 
          Q=0.
          EGAMOLD=EGAM
          IF (EGAM.LT.PAR_M1(1)) EGAM=PAR_M1(1)
          Q=Q+PAR_M1(2)/(EGAM**PAR_M1(3))
          EGAM=EGAMOLD
          SGAMMA=PAR_M1(1)*(PIH*Q+DMG)*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN
        ELSEIF (NOPTM1.EQ.25) THEN   ! Oslo - factor k + softpole with zero below PAR_M1(1)
          Q=0.
          DO I=1,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          Q=Q+PAR_M1(2)/(EGAM**PAR_M1(3))
          IF (EGAM.LT.PAR_M1(1)) Q=0
          SGAMMA=PAR_M1(1)*PIH*Q*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN

        ELSEIF (NOPTM1.EQ.27) THEN   ! Oslo - factor k not for LLR
          Q=0.
          DO I=1,1
            QQ=SIGM(I)*EGAM*WM0(I)**2/((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ/PAR_M1(1)
          ENDDO
          DO I=2,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          SGAMMA=PAR_M1(1)*PIH*Q*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN

        ENDIF
      ENDIF
!
!*****                       E2 component
      IF ((ITYP.EQ.2).OR.(ITYP.EQ.4)) THEN
        IF (NOPTE2.EQ.0) THEN
          SGAMMA=SGAMMA+QEL*EGAM**5
          SFCEE2=SGAMMA-SFCEM1
        ELSEIF (NOPTE2.EQ.11) THEN  ! Goriely tables
          SGAMMA = SGAMMA + APSF(EGAM,3)*EGAM**5
          SFCEE2=SGAMMA-SFCEM1
          RETURN
        ELSEIF(NOPTE2.EQ.1) THEN
          Q=0.
          DO I=1,NGIGE2
            QQ=(SIGE(I)*WE0(I)**2)/(EGAM*((EGAM**2-ERE(I)**2)**2+(EGAM*WE0(I))**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=SGAMMA+PIHQ*Q*EGAM**5
          SFCEE2=SGAMMA-SFCEM1
        ELSEIF(NOPTE2.EQ.4) THEN
          Q=0.
          DO I=1,NGIGE2
            QQ=(SIGE(I)*WE0(I)**2)/(EGAM*((EGAM**2-ERE(I)**2)**2+(EGAM*WE0(I))**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=SGAMMA+(PIHQ*Q+QEL)*EGAM**5
          SFCEE2=SGAMMA-SFCEM1
        ENDIF
        RETURN
      ENDIF
!
      RETURN
      END FUNCTION SGAMMA
!***********************************************************************
REAL FUNCTION TERM(EEXC)
!
!     The effective temperature corresponding to the exc.energy EEXC
!
!***********************************************************************
! attention!! : here corrected formula for EEF  / DEL ==> PAIRING !!! /
!
real::            EEXC
      EFEC_E=EEXC-PAIR_PSF
      IF (EFEC_E.LT.0.) EFEC_E=0.
      TERM=SQRT((EFEC_E)/ASHELL) ! TERM=SQRT((MAX(0.,EEXC-PAIR_PSF))/ASHELL)        
      RETURN
      END FUNCTION TERM
!***********************************************************************
REAL FUNCTION TERM1(EEXC)
!
!     The effective temperature corresponding to the exc.energy EEXC
!
!***********************************************************************
! attention!! : here corrected formula for EEF  / DEL ==> PAIR_PSF !!! /
!
real::            EEXC
      EFEC_E=EEXC
      IF (EFEC_E.LT.0.) EFEC_E=0.
      TERM1=SQRT((EFEC_E)/ASHELL)
      RETURN
      END FUNCTION TERM1
!***********************************************************************
REAL FUNCTION TERMOSLO(EEXC)
!
!     The effective temperature corresponding to the exc.energy EEXC
!
!***********************************************************************
! attention!! : here corrected formula for EEF  / DEL ==> PAIR_PSF !!! /
!
real::            EEXC
      EFEC_E=EEXC-PAIR_PSF+6.6/(AMASS**0.32)
      IF (EFEC_E.LT.0.) EFEC_E=0.
      TERMOSLO=SQRT((EFEC_E)/ASHELL)
      RETURN
      END FUNCTION TERMOSLO
!***********************************************************************
REAL FUNCTION EFEXC(EEXC)
!
!     The effective temperature corresponding to the exc.energy EEXC
!
!***********************************************************************
! attention!! : here corrected formula for EEF  / DEL ==> PAIR_PSF !!! /
!
real::            EEXC
      EFEXC=EEXC-PAIR_PSF
      RETURN
      END FUNCTION EFEXC
!***********************************************************************
REAL FUNCTION TERMDILG(EEXC)
!
!     The effective temperature corresponding to the exc.energy EEXC
!
!***********************************************************************
! attention!! : here corrected formula for EEF  / DEL ==> PAIR_PSF !!! /
!
real::            EEXC
      EFEC_E=EEXC-PAIR_PSF
      IF (EFEC_E.LT.0.) EFEC_E=0.
      TERMDILG=(1.0 + SQRT(1.0+4.0*ASHELL*EFEC_E))/2.0/ASHELL
      RETURN
      END FUNCTION TERMDILG
!***********************************************************************
REAL FUNCTION TERM_TK(EEXC)
      real:: EEXC
      REAL:: astar
      astar = 0.126181*AMASS + 7.52191e-05*AMASS*AMASS
      !let's determine astar
      if(TKematch.LE.0.0) then ! (not only) Kawano's Ematch is always 0.0 or higher, but let's be safe here
        astar = astar
      else
        if(EEXC.LT.TKematch) then
          if((TKematch-TKpair).GT.0.0) then
            astar = astar*(1.0+(1.0-exp(-0.31*(AMASS)**(-0.33333)*(TKematch-TKpair)))*TKeshell/(TKematch-TKpair)) !corresponds to ldShellCorrection(TKematch-TKpair,astar,TKeshell,AMASS)
          elseif((TKematch-TKpair).EQ.0.0) then
            astar = astar*(1.0+0.31*(AMASS)**(-0.33333)*TKeshell)
          else
            astar = astar
          endif
        else
          if((EEXC-TKpair).GT.0.0) then
            astar = astar*(1.0+(1.0-exp(-0.31*(AMASS)**(-0.33333)*(EEXC-TKpair)))*TKeshell/(EEXC-TKpair)) !corresponds to ldShellCorrection(EEXC-TKpair ,astar,TKeshell,AMASS)
          elseif((EEXC-TKpair).EQ.0.0) then
            astar = astar*(1.0+0.31*(AMASS)**(-0.33333)*TKeshell)
          else
            astar = astar
          endif
        endif
      endif !determine astar
      if (EEXC.LT.0.) EEXC=0.
      TERM_TK=SQRT(EEXC/astar) ! .cpp says t = (ex < 0.0) ? 0.0 : sqrt(ex/a); // but how can ex be <0? so this should be fine
      return
END FUNCTION TERM_TK
!***********************************************************************
REAL FUNCTION SMFREQ()
!
      SMFREQ=41./AMASS**(0.33333)
      RETURN
      END FUNCTION SMFREQ
!***********************************************************************
REAL FUNCTION EREL1()
!
      EREL1=31.2/(AMASS**(0.33333))+20.6/sqrt((AMASS**(0.33333)))
      RETURN
      END FUNCTION EREL1
!***********************************************************************
REAL FUNCTION WREL1()
!
      WREL1=0.026*EREL1()**1.91
      RETURN
      END FUNCTION WREL1
!***********************************************************************
REAL FUNCTION R0()
!
      R0=1.27*AMASS**(0.33333)
      RETURN
      END FUNCTION R0
!***********************************************************************
REAL FUNCTION FKS(EX,FKs0,FNS,WWALL,ALPPL)
real:: EX,FKs0,FNS,WWALL,ALPPL
      if (EX.GT.(2*EREL1())) then
        FKS=FKs0     
      else
        FKS=(WREL1()-EREL1()**2/ALPPL)/WWALL+(FKs0-(WREL1()-EREL1()**2/ALPPL)/WWALL)*ABS((EX-EREL1())/EREL1())**FNS
      endif
      RETURN
      END FUNCTION FKS
!***********************************************************************
REAL FUNCTION WWALL1()
!
          WWALL1=36.43*AMASS**(-0.3333333) 
      RETURN
      END FUNCTION WWALL1
!***********************************************************************
SUBROUTINE FactFk()
!   Returns values of standart functio Fk which are present in the
!   equation for angular correlation.
!   Data are from some russian book.
!   Standard factors F2 are stored in the form Fk(L,L',2*J1,2*J)
!   (reason is to include half-integer spins).
!   At this moment only combinations (L=1,L'=1), (L=1,L'=2) and
!   (L=2,L'=2) are assumed (no octupole mixing).
!   F4 factors are stored only for L=2, L'=2 in the form F4(2*J1,2*J).
!*********************************************************************
INTEGER::                             k,l,i,j
      DO k=0,24
       DO l=0,20
           F4(k,l)=0.
         DO i=1,4
          DO j=1,4
            Fk(i,j,k,l)=0.
          ENDDO
         ENDDO
       ENDDO
      ENDDO

      Fk(1,1,0,2)=0.707
      Fk(1,1,2,2)=-0.354
      Fk(1,1,4,2)=0.071
      Fk(1,1,2,4)=0.418
      Fk(1,1,4,4)=-0.418
      Fk(1,1,6,4)=0.120
      Fk(1,1,4,6)=0.346
      Fk(1,1,6,6)=-0.433
      Fk(1,1,8,6)=0.144
      Fk(1,1,6,8)=0.313
      Fk(1,1,8,8)=-0.439
      Fk(1,1,10,8)=0.160
      Fk(1,1,8,10)=0.294
      Fk(1,1,10,10)=-0.442
      Fk(1,1,12,10)=0.170
      Fk(1,1,10,12)=0.282
      Fk(1,1,12,12)=-0.443
      Fk(1,1,14,12)=0.177
      Fk(1,1,12,14)=0.273
      Fk(1,1,14,14)=-0.444
      Fk(1,1,16,14)=0.183
      Fk(1,1,14,16)=0.267
      Fk(1,1,16,16)=-0.445
      Fk(1,1,18,16)=0.187
      Fk(1,1,16,18)=0.262
      Fk(1,1,18,18)=-0.445
      Fk(1,1,20,18)=0.191
      Fk(1,1,18,20)=0.258
      Fk(1,1,20,20)=-0.446
      Fk(1,1,22,20)=0.194

      Fk(1,2,2,2)=-1.061
      Fk(1,2,4,2)=0.474
      Fk(1,2,2,4)=-0.935
      Fk(1,2,4,4)=-0.612
      Fk(1,2,6,4)=0.655
      Fk(1,2,4,6)=-0.949
      Fk(1,2,6,6)=-0.433
      Fk(1,2,8,6)=0.722
      Fk(1,2,6,8)=-0.940
      Fk(1,2,8,8)=-0.335
      Fk(1,2,10,8)=0.757
      Fk(1,2,8,10)=-0.931
      Fk(1,2,10,10)=-0.274
      Fk(1,2,12,10)=0.778
      Fk(1,2,10,12)=-0.923
      Fk(1,2,12,12)=-0.231
      Fk(1,2,14,12)=0.793
      Fk(1,2,12,14)=-0.917
      Fk(1,2,14,14)=-0.200
      Fk(1,2,16,14)=0.803
      Fk(1,2,14,16)=-0.912
      Fk(1,2,16,16)=-0.177
      Fk(1,2,18,16)=0.811
      Fk(1,2,16,18)=-0.907
      Fk(1,2,18,18)=-0.158
      Fk(1,2,20,18)=0.817
      Fk(1,2,18,20)=-0.904
      Fk(1,2,20,20)=-0.143
      Fk(1,2,22,20)=0.822

      Fk(2,2,2,2)=-0.354
      Fk(2,2,4,2)=0.354
      Fk(2,2,6,2)=-0.101
      Fk(2,2,0,4)=-0.598
      Fk(2,2,2,4)=-0.299
      Fk(2,2,4,4)=0.128
      Fk(2,2,6,4)=0.341
      Fk(2,2,8,4)=-0.171
      Fk(2,2,2,6)=-0.495
      Fk(2,2,4,6)=-0.124
      Fk(2,2,6,6)=0.227
      Fk(2,2,8,6)=0.309
      Fk(2,2,10,6)=-0.206
      Fk(2,2,4,8)=-0.448
      Fk(2,2,6,8)=-0.045
      Fk(2,2,8,8)=0.265
      Fk(2,2,10,8)=0.285
      Fk(2,2,12,8)=-0.228
      Fk(2,2,6,10)=-0.421
      Fk(2,2,8,10)=0.0
      Fk(2,2,10,10)=0.283
      Fk(2,2,12,10)=0.267
      Fk(2,2,14,10)=-0.243
      Fk(2,2,8,12)=-0.403
      Fk(2,2,10,12)=0.029
      Fk(2,2,12,12)=0.294
      Fk(2,2,14,12)=0.253
      Fk(2,2,16,12)=-0.253
      Fk(2,2,10,14)=-0.391
      Fk(2,2,12,14)=0.049
      Fk(2,2,14,14)=0.300
      Fk(2,2,16,14)=0.243
      Fk(2,2,18,14)=-0.261
      Fk(2,2,12,16)=-0.381
      Fk(2,2,14,16)=0.064
      Fk(2,2,16,16)=0.304
      Fk(2,2,18,16)=0.234
      Fk(2,2,20,16)=-0.268
      Fk(2,2,14,18)=-0.374
      Fk(2,2,16,18)=0.075
      Fk(2,2,18,18)=0.307
      Fk(2,2,20,18)=0.227
      Fk(2,2,22,18)=-0.273
      Fk(2,2,16,20)=-0.369
      Fk(2,2,18,20)=0.084
      Fk(2,2,20,20)=0.310
      Fk(2,2,22,20)=0.221
      Fk(2,2,24,20)=-0.277

      F4(0,4)=-1.069
      F4(2,4)=0.713
      F4(4,4)=-0.305
      F4(6,4)=0.076
      F4(8,4)=-0.008
      F4(2,6)=-0.447
      F4(4,6)=0.670
      F4(6,6)=-0.447
      F4(8,6)=0.149
      F4(10,6)=-0.020
      F4(4,8)=-0.304
      F4(6,8)=0.609
      F4(8,8)=-0.498
      F4(10,8)=0.194
      F4(12,8)=-0.030
      F4(6,10)=-0.243
      F4(8,10)=0.567
      F4(10,10)=-0.523
      F4(12,10)=0.224
      F4(14,10)=-0.037
      F4(8,12)=-0.209
      F4(10,12)=0.537
      F4(12,12)=-0.537
      F4(14,12)=0.246
      F4(16,12)=-0.043
      F4(10,14)=-0.187
      F4(12,14)=0.515
      F4(14,14)=-0.546
      F4(16,14)=0.263
      F4(18,14)=-0.048
      F4(12,16)=-0.173
      F4(14,16)=0.499
      F4(16,16)=-0.551
      F4(18,16)=0.276
      F4(20,16)=-0.059
      F4(14,18)=-0.162
      F4(16,18)=0.486
      F4(18,18)=-0.555
      F4(20,18)=0.286
      F4(22,18)=-0.056
      F4(16,20)=-0.154
      F4(18,20)=0.476
      F4(20,20)=-0.558
      F4(22,20)=0.295
      F4(24,20)=-0.059

!     Odd spins

      Fk(1,1,1,3)=0.500
      Fk(1,1,3,3)=-0.400
      Fk(1,1,5,3)=0.100
      Fk(1,1,3,5)=0.374
      Fk(1,1,5,5)=-0.428
      Fk(1,1,7,5)=0.134
      Fk(1,1,5,7)=0.327
      Fk(1,1,7,7)=-0.436
      Fk(1,1,9,7)=0.153
      Fk(1,1,7,9)=0.303
      Fk(1,1,9,9)=-0.440
      Fk(1,1,11,9)=0.165
      Fk(1,1,9,11)=0.288
      Fk(1,1,11,11)=-0.442
      Fk(1,1,13,11)=0.174
      Fk(1,1,11,13)=0.277
      Fk(1,1,13,13)=-0.444
      Fk(1,1,15,13)=0.180
      Fk(1,1,13,15)=0.270
      Fk(1,1,15,15)=-0.445
      Fk(1,1,17,15)=0.185
      Fk(1,1,15,17)=0.264
      Fk(1,1,17,17)=-0.445
      Fk(1,1,19,17)=0.189
      Fk(1,1,17,19)=0.260
      Fk(1,1,19,19)=-0.446
      Fk(1,1,21,19)=0.192

      Fk(1,2,1,3)=-0.866
      Fk(1,2,3,3)=-0.755
      Fk(1,2,5,3)=0.592
      Fk(1,2,3,5)=-0.949
      Fk(1,2,5,5)=-0.507
      Fk(1,2,7,5)=0.694
      Fk(1,2,5,7)=-0.945
      Fk(1,2,7,7)=-0.378
      Fk(1,2,9,7)=0.742
      Fk(1,2,7,9)=-0.935
      Fk(1,2,9,9)=-0.302
      Fk(1,2,11,9)=0.769
      Fk(1,2,9,11)=-0.927
      Fk(1,2,11,11)=-0.251
      Fk(1,2,13,11)=0.786
      Fk(1,2,11,13)=-0.920
      Fk(1,2,13,13)=-0.215
      Fk(1,2,15,13)=0.798
      Fk(1,2,13,15)=-0.914
      Fk(1,2,15,15)=-0.188
      Fk(1,2,17,15)=0.807
      Fk(1,2,15,17)=-0.910
      Fk(1,2,17,17)=-0.167
      Fk(1,2,19,17)=0.814
      Fk(1,2,17,19)=-0.906
      Fk(1,2,19,19)=-0.150
      Fk(1,2,21,19)=0.820

      Fk(2,2,1,3)=-0.500
      Fk(2,2,5,3)=0.357
      Fk(2,2,7,3)=-0.143
      Fk(2,2,1,5)=-0.535
      Fk(2,2,3,5)=-0.191
      Fk(2,2,5,5)=0.191
      Fk(2,2,7,5)=0.325
      Fk(2,2,9,5)=-0.191
      Fk(2,2,3,7)=-0.468
      Fk(2,2,5,7)=-0.078
      Fk(2,2,7,7)=0.249
      Fk(2,2,9,7)=0.296
      Fk(2,2,11,7)=-0.218
      Fk(2,2,5,9)=-0.433
      Fk(2,2,7,9)=-0.020
      Fk(2,2,9,9)=0.275
      Fk(2,2,11,9)=0.275
      Fk(2,2,13,9)=-0.236
      Fk(2,2,7,11)=-0.411
      Fk(2,2,9,11)=0.016
      Fk(2,2,11,11)=0.289
      Fk(2,2,13,11)=0.260
      Fk(2,2,15,11)=-0.248
      Fk(2,2,9,13)=-0.396
      Fk(2,2,11,13)=0.040
      Fk(2,2,13,13)=0.297
      Fk(2,2,15,13)=0.248
      Fk(2,2,17,13)=-0.258
      Fk(2,2,11,15)=-0.386
      Fk(2,2,13,15)=0.057
      Fk(2,2,15,15)=0.302
      Fk(2,2,17,15)=0.238
      Fk(2,2,19,15)=-0.265
      Fk(2,2,13,17)=-0.378
      Fk(2,2,15,17)=0.070
      Fk(2,2,17,17)=0.306
      Fk(2,2,19,17)=0.231
      Fk(2,2,21,17)=-0.270
      Fk(2,2,15,19)=-0.371
      Fk(2,2,17,19)=0.080
      Fk(2,2,19,19)=0.309
      Fk(2,2,21,19)=0.224
      Fk(2,2,23,19)=-0.275

      F4(1,5)=-0.617
      F4(3,5)=0.705
      F4(5,5)=-0.397
      F4(7,5)=0.118
      F4(9,5)=-0.015
      F4(3,7)=-0.358
      F4(5,7)=0.637
      F4(7,7)=-0.478
      F4(9,7)=0.174
      F4(11,7)=-0.025
      F4(5,9)=-0.268
      F4(7,9)=0.586
      F4(9,9)=-0.512
      F4(11,9)=0.210
      F4(13,9)=-0.034
      F4(7,11)=-0.224
      F4(9,11)=0.551
      F4(11,11)=-0.531
      F4(13,11)=0.236
      F4(15,11)=-0.041
      F4(9,13)=-0.197
      F4(11,13)=0.525
      F4(13,13)=-0.542
      F4(15,13)=0.255
      F4(17,13)=-0.046
      F4(11,15)=-0.179
      F4(13,15)=0.507
      F4(15,15)=-0.549
      F4(17,15)=0.270
      F4(19,15)=-0.051
      F4(13,17)=-0.167
      F4(15,17)=0.492
      F4(17,17)=-0.554
      F4(19,17)=0.281
      F4(21,17)=-0.054
      F4(15,19)=-0.158
      F4(17,19)=0.481
      F4(19,19)=-0.557
      F4(21,19)=0.291
      F4(23,19)=-0.058

      return
      END SUBROUTINE FactFk
!***********************************************************************
SUBROUTINE WRITELVLSCH(IGLOB,NBIN,ITID,IFLAG,NTOTAL,LEVCON)
!   This subroutin writes the level scheme, see code
character(12)::                       OUTNAME
integer::                             ILOCAL,IGLOB,NBIN,ITID,IFLAG,NTOTAL
integer,dimension(:,:,:),allocatable:: LEVCON
!
      OUTNAME='LVL_??.DAT'
      WRITE(OUTNAME(5:6),'(I2)') IGLOB
      IF (IGLOB.LE.9) WRITE(OUTNAME(5:5),'(A1)') '0'
      OPEN(UNIT=(IGLOB+10),FILE=OUTNAME)
      WRITE((IGLOB+10),*) 'NTOTAL =',NTOTAL,' IFLAG =',IFLAG
      DO ILOCAL=1,NBIN
        WRITE((IGLOB+10),*) ILOCAL,(LEVCON(ILOCAL,KLOCAL,0),KLOCAL=0,49)
      ENDDO
      RETURN
END SUBROUTINE WRITELVLSCH
!***********************************************************************
INTEGER FUNCTION NPOISS(IR,AM,IFLAG,U)
!     Poisson's random numbers
!
      INTEGER*4 IR,IFLAG
      REAL*4    U,AM,Q
!
      IF (AM.LE.50.) THEN
        Q=1.
        NPOISS=0
      2 Q=Q*RAN0(IR)
        IF (Q.LT.EXP(-AM)) RETURN
        NPOISS=NPOISS+1
        GOTO 2
      ENDIF
      NPOISS=INT(AM+SQRT(AM)*GAUSS(IR,U,IFLAG))
      RETURN
      END FUNCTION NPOISS
!***********************************************************************
REAL FUNCTION GAUSS(IR,U,IFLAG)
!     Normally distributed random numbers with
!     a zero mean and a unit variance
!     the IFLAG is 0-GAUSS->1-GAUSS->0
      INTEGER*4 IR,IFLAG
      REAL*4    U,X1,X2,SUM
!
      IF(IFLAG) 2,1,2
    1 X1=2.*RAN0(IR)-1.
      X2=2.*RAN0(IR)-1.
      SUM=X1*X1+X2*X2
      IF (SUM.GE.1.) GOTO 1
      SUM=SQRT(-2.*ALOG(SUM)/SUM)
      U=X1*SUM
      GAUSS=X2*SUM
      IFLAG=1
      RETURN
    2 IFLAG=0
      GAUSS=U
      RETURN
      END FUNCTION GAUSS
!***********************************************************************
REAL FUNCTION CHISQR(NOPTFL,IR,U,IFLAG)
!     chi^2 random numbers with DOF = NOPTFL
      INTEGER*4 IR,IFLAG,NOPTFL,I
      REAL*4    U,RND
!
      CHISQR=0.0
      DO I=1,NOPTFL
        RND=GAUSS(IR,U,IFLAG)
        CHISQR=CHISQR+RND**2
      ENDDO
      RETURN
      END FUNCTION CHISQR
!***********************************************************************
!TODO pridat cast SUBROUTINE WRITE_PARAMS(NAME)
end module spolecne
!***********************************************************************
