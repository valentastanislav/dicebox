!TODO gilbert cameron and spin distribution of initial states email jutta escher
!gamma_multiplicita make 

module vsechno
use lokalni_fce
use spolecne
implicit none

integer::                             NEVENTS,NREAL,NSUB,NDDD,NOPTFL,NLINc,ipinc,nfilev,NSTEPS,NENT,NENK,NEN_IPF,&
                                      KpopGS,ISWWR,ISWBN,ISWEL,ISWSP,ISWPA,ISWIC,ISWMX,ISWWI,ISWLS,NOPTCS,&
                                      N_MSC_FS,MIN_MULTIPLICITA,MAX_MULTIPLICITA,TRGT_PI

real::                                eall,EIN,EFI,ecrit,xrayk,xrayl,max_spin,factnrm,SUMNO,BIN_WIDTH,&
                                      D0,TRGT_SPIN

integer, dimension(1:199)::            ilowip
integer, dimension(0:49,0:1)::        levdis
integer, dimension(1:20,0:49,0:1)::   dekod
integer,dimension(:,:),allocatable::  isbspin

real,    dimension(:),allocatable::   MSC_FS
real,    dimension(1:2)::             CAPFR,spinc
real,    dimension(1:100)::           elent,elenk,ELEN_IPF
real,    dimension(1:199)::            elowlev,elowsp
real,    dimension(1:199,0:20)::       p_conv,p_conv_K,p_conv_IPF !TODO somehow smart determine the maximum number of decays in DIS and make these allocatable
real,    dimension(1:199,1:20)::       deltx
real,    dimension(0:1,1:5,1:100)::    CONVt,CONVk,CONV_IPF

contains
!***********************************************************************
SUBROUTINE READ_EV(NAME,lopopgs,KONTROLMATRIX)
!     reading of the input file
INTEGER,PARAMETER::     MAXJC  = 49
CHARACTER(80)::         TITLE1,TITLE2,TITLE3
character(80)::         NAME
logical::               lopopgs
integer,dimension(:,:),allocatable:: KONTROLMATRIX
INTEGER::               I,J,K,NMU,ipfi,ipar,control
REAL::                  enrg,spfi,enrgf,desp,dlt,alphak,alphaIPF,SPACRES,dummy,FSPAC,corrAlpha,corrDelta
      OPEN (UNIT=5,FILE=NAME,STATUS='OLD')
!     User's alphanumeric titles:
      READ (5,100) TITLE1
      READ (5,100) TITLE2
      READ (5,100) TITLE3
  100 FORMAT (A80)
!     The regime of run:
      READ (5,*)
      READ (5,*) ISWWR,ISWBN,ISWEL,ISWSP,ISWPA,ISWIC,ISWMX,ISWWI,ISWLS
      READ (5,*)
      Read (5,*) Nreal, NEVENTS, NSUB
      if (.not.allocated(KONTROLMATRIX)) then
       allocate(KONTROLMATRIX(1:4,1:NREAL))
      endif 
      READ (5,*)
      READ (5,*) NBIN,(KONTROLMATRIX(I,1),I=1,4)
      READ (5,*)
      READ (5,*) NOPTFL,NOPTE1,NOPTM1,NOPTE2,NOPTDE,LMODE,LDENP,LDSTAG
      READ (5,*)
      READ (5,*) N_MSC_FS, BIN_WIDTH
      if (N_MSC_FS.GE.1) then
      allocate (MSC_FS(1:N_MSC_FS))
      DO I=1,N_MSC_FS
        READ (5,*) MSC_FS(I)
      ENDDO
      endif
      MIN_MULTIPLICITA=1
      MAX_MULTIPLICITA=7
!
!     Giant Resonaces:
!
      NLOWLOR=0
      READ (5,*)
      IF((NOPTE1.EQ.66).OR.(NOPTE1.EQ.77)) THEN
        READ (5,*) NGIGE, NLOWLOR
        DO I=1, NGIGE+NLOWLOR
          READ (5,*) ER(I),W0(I),SIG(I)
        ENDDO
      ELSE
        READ (5,*) NGIGE
        DO I=1, NGIGE
          READ (5,*) ER(I),W0(I),SIG(I)
        ENDDO
      ENDIF
      READ (5,*)
      READ (5,*) NGIGM
      DO I=1, NGIGM
        READ (5,*) ERM(I),WM0(I),SIGM(I)
      ENDDO
      READ (5,*)
      READ (5,*) NGIGE2
      DO I=1, NGIGE2
        READ (5,*) ERE(I),WE0(I),SIGE(I)
      ENDDO
!
!     Other data needed for photon strength:
      READ (5,*)
      READ (5,*) DEG,DMG,QEL
      READ (5,*)
      READ (5,*) FERMC, TCONST, PAIR_PSF
      READ (5,*)
      READ (5,*) EK0,EGZERO
      READ (5,*)
      READ (5,*) (PAR_E1(I),I=1,3) ! DIPELO,DIPEHI,DIPSUP
        DIPSLP=(1.0-PAR_E1(3))/(PAR_E1(2)-PAR_E1(1))
        DIPZER=1.0-DIPSLP*PAR_E1(2)
      READ (5,*)
      READ (5,*) (PAR_M1(I),I=1,4)
!
!     Data needed for level density formulas:
!
      READ (5,*)
      READ (5,*) EZERO,DEL,TEMPER,ASHELL,AMASS,ZNUM,PAIRING
      READ (5,*)
      READ (5,*) ASHELL09,DEL09,TEMPER09,EZERO09,PAIRING09,SIG_CUSTOM
      READ (5,*)
      READ (5,*) DENPPC,DENPA0,DENPA1,DENPA2        !see PRC67, 015803
      IF ((LDENP.GT.0).AND.(DENPPC.LT.0.)) THEN
        WRITE(*,*) 'Negative DENPPC means parity asymmetry is appearing with excitation energy',&
         ' while oposite is usually the case. Are You sure about the DENPPC = ',DENPPC,' ?'
      ENDIF
      READ (5,*)
      READ (5,*) DENLO,DENHI,DENPA,DENPB,DENPC,DENPD
      IF (LDSTAG.NE.0) THEN
        IF ((MOD(INT(AMASS+0.25),2).NE.0).OR.(MOD(INT(ZNUM+0.25),2).NE.0)) THEN
            WRITE(*,*) 'Are You sure You want to use staggering in isotope which is not even-even?!'
        ENDIF
        WRITE(*,*) 'Effect of staggering linearly decrasing from energy of ',DENLO,' MeV to ',DENHI,' MeV'
      ENDIF
!
!     Data relating to the (thermal) neutron capturing state:
!
      READ (5,*)
      READ (5,*) BN,SPINc(1),IPINC,TRGT_SPIN,TRGT_PI
      NOPTCS=1
      NLINc=1
      CAPFR(1)=1.
!
!     Now read the tables of ICC-coefficients and other information
!     on electron conversion. The used notation:
!
!             XRAYK,XRAYL         -- electron binding energy for K-shell
!                                    (and L-shell) expressed in MeV
!             NENT,NENK           -- the number of electron energy points
!             NMU                 -- the highest multipolarity, e.g. 3 for
!                                    octupole (electric and magnetic)
!             CONVK(I,MU,K)       -- K-shell ICC for I-th type of radiation
!                                    (0 for electric, 1 for magnetic),
!                                    multipolarity MU and K-th electron
!                                    energy value
!             CONVT(I,M,K)        -- total ICC ... (as the CONVK)
!             ELENT(K),ELENK(K)   -- K-th value of electron energy
!
      READ (5,*)
      READ (5,*) xrayk,xrayl
      READ (5,*) NMU
      READ (5,*) NENT
      READ (5,*) (ELENT(K),K=1,NENT)
      READ (5,*) (((CONVT(I,J,K),I=0,1),J=1,NMU),K=1,NENT)
!     Let's handle alpha_K now
!     number of energies
      READ (5,*) NENK
!     add 3 energies to the beginning (corresponding to below threshold) and put zeroes to the coefficients
      DO K=1,3
        ELENK(K)=ELENT(NENT-NENK-3+K)
        DO I=0,1
          DO J=1,NMU
            CONVK(I,J,K)=0.0E+00
          ENDDO
        ENDDO
      ENDDO
!     now read the rest - input itself
      READ (5,*) (ELENK(K),K=1+3,NENK+3)
      READ (5,*) (((CONVK(I,J,K),I=0,1),J=1,NMU),K=1+3,NENK+3)
!     now tell the code we've added 3 energies
      NENK=NENK+3
!     Let's handle alpha_IPF now
!     number of energies
      READ (5,*) NEN_IPF
!     add 3 energies to the beginning (corresponding to below threshold) and put zeroes to the coefficients
      DO K=1,3
        ELEN_IPF(K)=ELENT(NENT-NEN_IPF-3+K)
        DO I=0,1
          DO J=1,NMU
            CONV_IPF(I,J,K)=0.0E+00
          ENDDO
        ENDDO
      ENDDO
!     now read the rest - input itself
      READ (5,*) (ELEN_IPF(K),K=1+3,NEN_IPF+3)
      READ (5,*) (((CONV_IPF(I,J,K),I=0,1),J=1,3),K=1+3,NEN_IPF+3)
!     now tell the code we've added 3 energies
      NEN_IPF=NEN_IPF+3
!
!     Data related to the discrete levels (J, pi, Eexc, primary intensities
!     and branchings):
!
      READ (5,*)
      READ (5,*) EALL,ECRIT ! od ground state do EALL zname "uplne" vsechno o hladinach, mezi EALL a ECRIT "jen" energie, spin a paritu, a statisticky (do)generujeme jejich rozpad
      read (5,*)
      LOpopGS=.false.
      DO J=0,MAXJC
        DO K=0,1
          NDIS(J,K)=0
        ENDDO
      ENDDO
      max_decays=0
      max_spin=0.
      READ (5,*) numlev !TODO zde bude druha promenna udavajici pocet hladin mezi EALL a ECRIT
      if (.not.allocated(ityp)) then
        allocate (ityp(1:numlev,1:2))
      endif  
      if (.not.allocated(isbspin)) then
        allocate (isbspin(1:numlev,1:2))
      endif  
      DO i=1,numlev
       READ (5,*) enrg,spfi,ipfi,denum(i),dummy,dummy,LVL_CLASS(i)
       write(*,*) 'reading lvl # ',i,' at energy ',enrg
       IF((ipfi.NE.0).AND.(ipfi.NE.1)) THEN
         WRITE(*,*) 'parity can be either 0 (+) or 1 (-)'
         STOP
       ENDIF
       LVL_ENERGY(i)=enrg
       ndis(ISUBSC(spfi),ipfi)=ndis(ISUBSC(spfi),ipfi)+1
       endis(ndis(ISUBSC(spfi),ipfi),ISUBSC(spfi),ipfi)=enrg
       dekod(ndis(ISUBSC(spfi),ipfi),ISUBSC(spfi),ipfi)=i
       IF (denum(i).GT.0) then
        DO k=1,denum(i)
!         write(*,*) i,k,enrg
         READ (5,*) enrgf,sal(i,k),errsal(i,k),desp,ipar,dlt,alpha(i,k)
         IF((ipar.NE.0).AND.(ipar.NE.1)) THEN
           WRITE(*,*) 'parity can be either 0 (+) or 1 (-)'
           STOP
         ENDIF
         control=0
         DO j=1,ndis(ISUBSC(desp),ipar)
           IF(elowlev(dekod(j,ISUBSC(desp),ipar)).NE.enrgf) control=control+1
         ENDDO
         IF ((enrgf.GT.enrg).OR.(control.EQ.ndis(ISUBSC(desp),ipar))) THEN
          WRITE(*,*) 'wrong decay pattern at level',enrg,'to level',enrgf,desp,ipar
          STOP
         ENDIF
         DO j=1,ndis(ISUBSC(desp),ipar)
           IF (endis(j,ISUBSC(desp),ipar).EQ.enrgf) THEN
            delev(i,k)=j
            despin(i,k)=desp
            deparity(i,k)=ipar
            deltx(i,k)=dlt
            if (alpha(i,k).LE.1e-6) alpha(i,k)=ALPH_TOT(enrg,spfi,ipfi,enrgf,desp,ipar,dlt**2,nent,elent,convt)
            alphak=ALPH_TOT(enrg,spfi,ipfi,enrgf,desp,ipar,dlt**2,nenk,elenk,convk)
            alphaIPF=ALPH_TOT(enrg,spfi,ipfi,enrgf,desp,ipar,dlt**2,NEN_IPF,ELEN_IPF,CONV_IPF)
            if (alphak.GT.alpha(i,k)) alphak=alpha(i,k)
            p_conv(i,k)=alpha(i,k)/(1+alpha(i,k))
            p_conv_K(i,k)=alphak/(1+alpha(i,k))
            p_conv_IPF(i,k)=(alphaIPF+alphak)/(1+alpha(i,k))
           ENDIF
         ENDDO
        ENDDO
       ENDIF
       elowlev(i)=enrg
       IF(elowlev(i).EQ.0.) THEN
        LOpopGS=.TRUE.
        KpopGS=i
       ENDIF
       elowsp(i)=spfi
       ilowip(i)=ipfi
       IF(denum(i).GT.max_decays) max_decays=denum(i)
       IF(elowsp(i).GT.max_spin) max_spin=elowsp(i)
!
       DO j=1,NLINc
         isbspin(i,j)=NINT(spfi+.25)-NINT(spinc(j)+.25)
         ityp(i,j)=ITYPE(spinc(j),ipinc,spfi,ipfi)
       ENDDO
       IF (NLINc.EQ.1) ityp(i,2)=0
      ENDDO !i over numlev
!
      DO ipfi=0,1
       DO j = 0,ISUBSC(max_spin)
        spfi = spinc(1)-int(spinc(1)+.25) + float(j)
        levdis(ISUBSC(spfi),ipfi)=nddd+ndis(ISUBSC(spfi),ipfi)
        nddd=levdis(ISUBSC(spfi),ipfi)
       ENDDO
      ENDDO
!
      CLOSE (UNIT=5,STATUS='KEEP')
      WRITE(*,*) 'nddd ',nddd,' numlev ',numlev

      IF ((NOPTDE.EQ.8).OR.(NOPTDE.EQ.9)) THEN
        IF (MOD(INT(AMASS+0.25),2).EQ.0) THEN
          IF (MOD(INT(ZNUM+0.25),2).EQ.0) THEN
            WRITE(*,*) 'beware, LD model is going to use staggering!'
          ENDIF
        ENDIF
      ENDIF
!
!     Tabulated level density (Kawano)
!
      IF (NOPTDE.EQ.12) THEN
        NAME = "LDTAB_K.DAT"
        OPEN (UNIT=5,FILE=NAME,STATUS='OLD')
        READ(5,*) 
        READ(5,*) NLD
        READ(5,*) 
        DO I = NLD, 1, -1
          READ(5,*) TABENLD(I),DUMMY,(TABLD(I,J,0),J=0,9) 
        ENDDO
        DO I = 1, NLD
          DO J = 0, 9
            TABLD(I,J,0) = TABLD(I,J,0) / 2.0 ! Kawano gives total NLD (p-independent)
            TABLD(I,J,1) = TABLD(I,J,0) 
          ENDDO
          DO J = 10, MAXJC
            TABLD(I,J,0) = 0.0
            TABLD(I,J,1) = 0.0            
          ENDDO
        ENDDO
        CLOSE (UNIT=5,STATUS='KEEP')
        IF (TABENLD(NLD).LT.BN) THEN
          WRITE(*,*) 'tabulated LD does not span up to the initial state: ',TABENLD(NLD),' vs ',BN
          STOP
        ENDIF
      ENDIF
!
!     Tabulated level density (Goriely)
!
      IF (NOPTDE.EQ.11) THEN
        NAME = "LDTAB.DAT"
        OPEN (UNIT=5,FILE=NAME,STATUS='OLD')
        READ(5,*) 
        READ(5,*) NLD, SPACRES,spfi,ipfi,corrAlpha,corrDelta
        READ(5,*) 
        DO I = 1, NLD
          READ(5,*) TABENLD(I),DUMMY,DUMMY,DUMMY,DUMMY,(TABLD(I,J,0),J=0,MAXJC)
        ENDDO
        READ(5,*) 
        READ(5,*) 
        READ(5,*) 
        DO I = 1, NLD
          READ(5,*) TABENLD(I),DUMMY,DUMMY,DUMMY,DUMMY,(TABLD(I,J,1),J=0,MAXJC) 
        ENDDO
        IF (TABENLD(NLD).LT.BN) THEN
          WRITE(*,*) 'tabulated LD does not span up to the initial state: ',TABENLD(NLD),' vs ',BN
          STOP
        ENDIF
!
!     "Normalization to experimental level spacing" - spfi/ipfi needed
!
        IF (spfi.EQ.0.0) THEN !matching the s-wave spacing for even-even target with gs spin 0
         FSPAC = (1.0 / SPACRES * 1e6) / DENSITY(BN,spfi+0.5,ipfi)
         corrAlpha=0.0
         corrDelta=0.0
         write(*,*) 'using tabulated NLD matched to s-wave resonance spacing of ',SPACRES
        ELSEIF (spfi.EQ.-1.0) THEN !constant renormalization by one common factor
         write(*,*) 'using tabulated NLD with renormalization by one given factor of ',SPACRES
         FSPAC = SPACRES
         corrAlpha=0.0
         corrDelta=0.0
        ELSEIF (spfi.EQ.-2.0) THEN !renormalization by Goriely PRC 78 064307 (2008): \rho(u,J,p)=\exp(alpha x \sqrt(U-delta)) * \rho(U-delta,J,p)
         write(*,*) 'using tabulated NLD with renormalization by Goriely PRC 78 064307 (2008) from ',TABENLD(1),&
         ' to ', TABENLD(NLD),' with alpha ',corrAlpha,' and delta ',corrDelta
         FSPAC = 1
        ELSE !matching the s-wave spacing for target with nonzero spin
         FSPAC = (1.0 / SPACRES * 1e6) /(DENSITY(BN,spfi-0.5,ipfi) + DENSITY(BN,spfi+0.5,ipfi))
         corrAlpha=0.0
         corrDelta=0.0
         write(*,*) 'tabulated NLD gives s-wave resonance spacing of ',1e6/(DENSITY(BN,spfi-0.5,ipfi) + DENSITY(BN,spfi+0.5,ipfi))
         write(*,*) 'using tabulated NLD matched to s-wave resonance spacing of ',SPACRES,' by a factor ',FSPAC
        ENDIF
        DO I = 1, NLD !Standa's version
         TABENLD(I)=TABENLD(I)+corrDelta
         DO J = 0, MAXJC
          DO K = 0, 1
            if (corrAlpha.eq.0.0) then
              TABLD(I,J,K) = TABLD(I,J,K) * FSPAC
            else
              TABLD(I,J,K) = TABLD(I,J,K) * FSPAC * exp(corrAlpha*sqrt(TABENLD(I)-corrDelta))
            endif
          ENDDO
         ENDDO
        ENDDO
        write(*,*) 'normalization of tabulated NLD done'
        CLOSE (UNIT=5,STATUS='KEEP')
      ENDIF

      IF (TRGT_SPIN.EQ.0.0) THEN
        D0=1e6/DENSITY(BN,TRGT_SPIN+0.5,TRGT_PI)
      ELSE
        D0=1e6/(DENSITY(BN,TRGT_SPIN-0.5,TRGT_PI) + DENSITY(BN,TRGT_SPIN+0.5,TRGT_PI))
      ENDIF
!
!     Tabulated PSF
!
      IF ((NOPTE1.EQ.11).OR.(NOPTE1.EQ.50)) THEN
        NAME = "PSFE1.DAT"
        OPEN (UNIT=5,FILE=NAME,STATUS='OLD')
        READ(5,*) 
        READ(5,*) NPSF(1), dummy
        READ(5,*)
        TABPSF(1,0)=0.
        DO I = 1, NPSF(1)
          READ(5,*) TABENPSF(1,I),TABPSF(1,I)
          TABPSF(1,I)=dummy*TABPSF(1,I)
        ENDDO
        CLOSE(5)
        IF (TABENPSF(1,NPSF(1)).LT.BN) THEN
          WRITE(*,*) 'tabulated E1 does not span up to the initial state: ',TABENPSF(1,NPSF(1)),' vs ',BN
          STOP
        ENDIF
      ENDIF

      IF (NOPTM1.EQ.11) THEN
        NAME = "PSFM1.DAT"
        OPEN (UNIT=5,FILE=NAME,STATUS='OLD')
        READ(5,*) 
        READ(5,*) NPSF(2), dummy
        READ(5,*) 
        TABPSF(2,0)=0.
        DO I = 1, NPSF(2)
          READ(5,*) TABENPSF(2,I),TABPSF(2,I)
          TABPSF(2,I)=dummy*TABPSF(2,I)
        ENDDO
        CLOSE(5)
        IF (TABENPSF(2,NPSF(2)).LT.BN) THEN
          WRITE(*,*) 'tabulated M1 does not span up to the initial state: ',TABENPSF(2,NPSF(2)),' vs ',BN
          STOP
        ENDIF
      ENDIF

      IF (NOPTE2.EQ.11) THEN
        NAME = "PSFE2.DAT"
        OPEN (UNIT=5,FILE=NAME,STATUS='OLD')
        READ(5,*) 
        READ(5,*) NPSF(3), dummy
        READ(5,*) 
        TABPSF(1,0)=0.
        DO I = 1, NPSF(3)
          READ(5,*) TABENPSF(3,I),TABPSF(3,I)
          TABPSF(3,I)=dummy*TABPSF(3,I)
        ENDDO
        CLOSE(5)
        IF (TABENPSF(3,NPSF(3)).LT.BN) THEN
          WRITE(*,*) 'tabulated E2 does not span up to the initial state: ',TABENPSF(3,NPSF(3)),' vs ',BN
          STOP
        ENDIF
      ENDIF
!
      DELTA=(BN-ECRIT)/FLOAT(NBIN)
      WRITE(*,*) ' Inputs have been successfully loaded.'

      RETURN
END SUBROUTINE READ_EV
!***********************************************************************
SUBROUTINE GERMS(IR,NTOTAL,NDDD,IRCONc,IRCON)     !should be OK
!
!     This subroutine generates a random-generator random seed for
!     each individual level. The seeds obtained are stored in IRCON(K).
!     For a level of interest, characterized by {IP,SP,IB,IL}, the
!     corresponding seed can be fetched in a simple way. This is evident
!     from the body of the functio SEED that performs such an operation.
!
!***********************************************************************
! vstup(NTOTAL_tp,NDDD_c,NOPTCS_c) vystup(IRCON,IRCONc)
INTEGER::                             I,IG,K,N,NOLD,NDIF,KAUX,JAUX
INTEGER,dimension(:),allocatable::    IAUX
REAL::                                dummy
integer::                             IR,NTOTAL,NDDD
integer,dimension(:),allocatable::    IRCONc,IRCON
!globalni NLINc
      if (allocated(IRCON)) then
        deallocate(IRCON)
      endif  
      allocate(IRCON(1:NTOTAL+NDDD))
      if (allocated(IRCONc)) then
        deallocate(IRCONc)
      endif  
      allocate(IRCONc(1:NLINc))
      
!     Seed for 1st capturing state
      NOLD = 0
      NDIF = 0
      IRCONc(1)=IR
!     Initialization of IRCON
      DO K=1,NTOTAL+NDDD
        IRCON(K)=0
      ENDDO
!     IRCON(.) is seeded in a random fashion by seeds IR produced
!     consecutively by repeated calls of RAN0(IR), see the body of the
!     routine. If two or more seeds are falling to the same word of
!     IRCON(.), then only the last of them is kept there.
    4 DO IG=1,NTOTAL+NDDD
        dummy=RAN0(IR)
        K=CEILING(DBLE(NTOTAL+NDDD)*dummy)  !obcas dela problem predelano by SV
        IF (K.GT.(NTOTAL+NDDD)) write(*,*) dummy, NTOTAL+NDDD, DBLE((NTOTAL+NDDD)*dummy)
        dummy=RAN0(IR)
        IRCON(K)=IR
      ENDDO
!     Determine the number of not seeded IRCON(K)
      N=0
      DO K=1,NTOTAL+NDDD
        IF (IRCON(K).EQ.0) N=N+1
      ENDDO
      NDIF = ABS(NOLD - N)
      NOLD = N
      IF ((NOLD.GT.2000).AND.(NDIF.GT.1000)) GOTO 4
      IF (NOLD.EQ.0) GOTO 7
      allocate(IAUX(1:NOLD))
!     If every site of the IRCON(.) is seeded, the possible seeding of 2nd capture state 
!     and RETURN follows. If not, ONLY EMPTY positions are additionally
!     seeded. Hopefully, seeds are again being distributed randomly.
!     This new seeding to a restricted area ensures that ALL
!     positions are seeded. After that the (2nd CS &) RETURN follows.
      N=0
      DO K=1,NTOTAL+NDDD
        IF (IRCON(K).EQ.0) THEN
          N=N+1
          IAUX(N)=K
        ENDIF
      ENDDO
      IF (N.NE.NOLD) WRITE(*,*) 'prusvih velikost'  !debug line
      DO KAUX=1,NOLD
        I=INT(FLOAT(N)*RAN0(IR))+1
        dummy=RAN0(IR)
        IRCON(IAUX(I))=IR
        DO JAUX=I+1,N
          IAUX(JAUX-1)=IAUX(JAUX)
        ENDDO
        N=N-1
      ENDDO
      IF (N.NE.0) WRITE(*,*) 'prusvih nenulovost'  !debug line
    7 IF (NLINc.NE.1) THEN
        dummy=RAN0(IR)
        IRCONc(2)=IR         ! 2nd capture state is now seeded
      ENDIF
      if (allocated(IAUX)) then
        deallocate(IAUX)
      endif  
      RETURN
END SUBROUTINE GERMS
!***********************************************************************
FUNCTION SEEDS(MODE,SP,IP,IL,IB,LEVCON,IRCON,IRCONc)     !should be OK
!     This functio yields the seed that is needed for the intrinsic
!     functio RAN inside the subroutine WIDTHS immediately before
!     this subroutine starts to generate the values of the subscripted
!     variable STCON.
!***********************************************************************
INTEGER::                             SEEDS
integer::                             MODE,IP,IL,IB
real::                                SP
integer,dimension(:),allocatable::    IRCONc,IRCON
integer,dimension(:,:,:),allocatable::LEVCON
!no globals, pro-ly strictly threadprivate
      IF (MODE.GT.0) THEN
        SEEDS=IRCONc(MODE)
      ELSE
        SEEDS=IRCON(LEVCON(IB-1,ISUBSC(SP),IP)+IL)
      ENDIF
      RETURN
END FUNCTION SEEDS
!***********************************************************************
!should be OK, just TODOs for chi2 random number generation
SUBROUTINE WIDTHS_R (MODE,IPIN,SPIN,IBIN,ILIN,TOTCON,STCON,GACON,ISCON,TOTDIS,STDIS,GADIS,ISDIS,LEVCON,&
IRCON,IRCONc,IFLAG,U,IREGI,EIN,EFI)
!***********************************************************************
double precision,dimension(0:2,-2:2,0:1)::     GACON
double precision,dimension(:,:,:,:),allocatable::  STCON
double precision,dimension(0:2)::     TOTCON
integer,dimension(:,:,:,:),allocatable:: ISCON
real, dimension(0:2,0:20,-2:2,0:1)::  STDIS
real, dimension(0:2,-2:2,0:1)::       GADIS
real, dimension(0:2)::                TOTDIS
real::                                U,SPIN,EIN,EFI
integer::                             MODE,IPIN,IBIN,ILIN,IFLAG,IREGI
integer,dimension(:),allocatable::    IRCONc,IRCON
integer,dimension(:,:,:),allocatable::LEVCON
integer,dimension(:,:,:,:),allocatable::ISDIS
INTEGER::                             IPFI,ISPFI,ISP,IP,IS,I,NLEV,ISEED,ISEEDGS,ISEEDORIG,ISBS,IT,IT1,IT2,ITT,IL,NL,ICOR
REAL::                                SPFI,Q,SPAC,SP,EG,Z,GSQ,G,CLEB,alpha,DMIX2
REAL,dimension(1:2)::                 SIMPL,GG
!globalni NOPTFL,DELTA,corri,re2res,im2res,NDIS,ENDIS,BN,NBIN
!
!     If you use MODE=1 or 2 you have to set IBIN=1 and ILIN=1 !
!
      IF (IREGI.EQ.0) THEN
        IF (MODE.NE.0) THEN
          EIN=BN
          Q=0.0
        ELSE
          NLEV=LEVCON(IBIN,ISUBSC(SPIN),IPIN)-LEVCON(IBIN-1,ISUBSC(SPIN),IPIN)
          Q=1.-FLOAT(2*ILIN-1)/FLOAT(2*NLEV)
          EIN=BN-(FLOAT(IBIN)-Q)*DELTA
        ENDIF
        IF (NOPTFL.GE.1) ISEED=SEEDS(MODE,SPIN,IPIN,ILIN,IBIN,LEVCON,IRCON,IRCONc) !originaly .EQ.
      ELSE
        EIN=EFI
      ENDIF !end IF (IREGI.EQ.0) THEN
      SPAC=1./DENSITY(EIN,SPIN,IPIN)
      IF (SPAC.LE.0.0) SPAC=0.0000001  !nikdy byt nemusi
      TOTCON(MODE)=0.D+0
      TOTDIS(MODE)=0.
!
      DO IS=-2,2
       DO IP=0,1
        DO I=0,NBIN
         STCON(MODE,I,IS,IP)=0.D+0
        ENDDO
!        IF (MODE.EQ.0) THEN !what the hell
        DO I=0,20
         STDIS(MODE,I,IS,IP)=0.
        ENDDO
!        ENDIF
        GACON(MODE,IS,IP)=0.D+0
        GADIS(MODE,IS,IP)=0.
       ENDDO !IP
      ENDDO !IS
!
      SP  = SPIN - INT(SPIN + .25)
      ISP = INT(SPIN + .25)
      DO IPFI=0,1
        DO ISPFI= ISP-2, ISP+2
          SPFI = SP + FLOAT(ISPFI)
          ISBS=NINT(SPFI+.25)-NINT(SPIN+.25)
          IT=ITYPE(SPIN,IPIN,SPFI,IPFI)
          IF (IT.EQ.0) GOTO 1
          IF (IT.EQ.2) THEN
            IT1=3
            IT2=4
          ELSE
           IT1=IT
           IT2=IT
          ENDIF
      STCON(MODE,IBIN,ISBS,IPFI)=0.D+0
      IF (iregi.eq.1) then
       ISEED=SEEDS(mode,spin,ipin,ilin,ibin,LEVCON,IRCON,IRCONc)
       GOTO 24
      ENDIF
!
      DO I=IBIN+1,NBIN
       NL=LEVCON(I,ISUBSC(SPFI),IPFI)-LEVCON(I-1,ISUBSC(SPFI),IPFI)
       EG=(FLOAT(I-IBIN-1)+Q+.5)*DELTA    ! Q=0.5 ... fixed; trans. between mid-bins
       Z=0.
       DO ITT=IT1,IT2
        SIMPL(ITT-IT1+1)=SGAMMA(EG,EIN,ITT)
       ENDDO
       IF (NOPTFL.LT.1) THEN !originaly .NE.
        DO ITT=IT1,IT2
         Z=Z+FLOAT(NL)*SIMPL(ITT-IT1+1)
         GG(ITT-IT1+1)=FLOAT(NL)*SIMPL(ITT-IT1+1)
        ENDDO
        !TODO can I get a decent speed back when I introduce some IFs here?
        IF (GG(1).GT.0.0) THEN
         DMIX2 = GG(2) / GG(1)
        ELSE
         DMIX2 = 0.0
        ENDIF
        alpha = ALPH_TOT(EG,spfi,ipfi,0.0,spin,ipin,DMIX2,nent,elent,convt)
        Z = Z + (1. + alpha) * (GG(1)+GG(2))
       ELSE  !TODO modify for chi2 with more DOF
        ISCON(MODE,I,ISBS,IPFI)=ISEED
        IFLAG=0
        IF (MODE.EQ.0) THEN
         DO IL=1,NL  !neprobehne kdyz NL je 0, coz ma za nasledek nulovou intenzitu do prazdnych binu
          DO ITT=IT1,IT2
           G=GAUSS(ISEED,U,IFLAG)
          !  Z=Z+G*G*SIMPL(ITT-IT1+1)
           GG(ITT-IT1+1)=G*G*SIMPL(ITT-IT1+1)
          ENDDO
          !TODO can I get a decent speed back when I introduce some IFs here?
          IF ((GG(2).GT.0.0).OR.(GG(1).GT.0.0)) THEN
            IF (GG(1).GT.0.0) THEN
             DMIX2 = GG(2) / GG(1)
            ELSE
             DMIX2 = 0.0
            ENDIF
            alpha = ALPH_TOT(EG,spfi,ipfi,0.0,spin,ipin,DMIX2,nent,elent,convt)
            Z = Z + (1. + alpha) * (GG(1)+GG(2))
          ENDIF
         ENDDO  !IL
        ELSE                    ! Primary transitions (the same fluctuation)
         DO IL=1,NL
          DO ITT=IT1,IT2
           GSQ=CHISQR(NOPTFL,ISEED,U,IFLAG) !originaly GSQ=GAUSS(ISEED,U,IFLAG)
          !  Z=Z+GSQ*SIMPL(ITT-IT1+1)
           GG(ITT-IT1+1)=GSQ*SIMPL(ITT-IT1+1)
          ENDDO
          IF ((GG(2).GT.0.0).OR.(GG(1).GT.0.0)) THEN
            IF (GG(1).GT.0.0) THEN
              DMIX2 = GG(2) / GG(1)
            ELSE
              DMIX2 = 0.0
            ENDIF
            alpha = ALPH_TOT(EG,spfi,ipfi,0.0,spin,ipin,DMIX2,nent,elent,convt)
            Z = Z + (1. + alpha) * (GG(1)+GG(2))
          ENDIF
         ENDDO  !IL
        ENDIF !IF (MODE.EQ.0)
       ENDIF !IF (NOPTFL.LT.1)
!
!  fix equivalent to NL.EQ.0 cause of construction of Z
       IF (DENSITY(EIN,SPIN,IPIN).GT.0.0) THEN
         STCON(MODE,I,ISBS,IPFI)=STCON(MODE,I-1,ISBS,IPFI)+DBLE(Z*SPAC)
       ELSE
         STCON(MODE,I,ISBS,IPFI)=STCON(MODE,I-1,ISBS,IPFI)
       ENDIF
!       IF (IT.EQ.1) THEN
!         GE1=GE1+Z
!       ELSEIF (IT.EQ.3) THEN
!         GM1=GM1+Z
!       ELSEIF (IT.EQ.4) THEN
!         GE2=GE2+Z
!       ELSEIF (IT.EQ.2) THEN          ! This is not fully correct (strength fluctuate)
!         GM1=GM1+Z*SIMPL(1)/(SIMPL(1)+SIMPL(2))   ! But it is acceptable approximation
!         GE2=GE2+Z*SIMPL(2)/(SIMPL(1)+SIMPL(2))
!       ENDIF
      ENDDO !I=IBIN+1,NBIN
!
!     For primary transitions it is assumed that MODE=0. In such
!     a case it is understood that the values of the ACTUAL
!     subscripted variable that replaces the FICTIVE variable STDIS
!     are simply derived from input data without the need of Monte
!     Carlo simulation.
!
!      IF (MODE.NE.0) GOTO 22    !!!!!!!! Compare to WIDTHS()
   24 CONTINUE
      DO I=1,NDIS(ISUBSC(SPFI),IPFI)
       EG=EIN-ENDIS(I,ISUBSC(SPFI),IPFI)
       IF (NOPTFL.GE.1) ISDIS(MODE,I,ISUBSC(SPFI),IPFI)=ISEED !originaly .EQ.
       IFLAG=0
       Z=0.
       DO ITT=IT1,IT2
        SIMPL(ITT-IT1+1)=SGAMMA(EG,EIN,ITT)
       ENDDO
       IF (NOPTFL.GE.1) THEN !originaly .EQ.
        IF (MODE.EQ.0) THEN
         DO ITT=IT1,IT2
          G=GAUSS(ISEED,U,IFLAG)
          ! Z=Z+G*G*SIMPL(ITT-IT1+1)
          GG(ITT-IT1+1)=G*G*SIMPL(ITT-IT1+1)
         ENDDO
         IF (GG(1).GT.0.0) THEN
          DMIX2 = GG(2) / GG(1)
         ELSE
          DMIX2 = 0.0
         ENDIF
         alpha = ALPH_TOT(EG,spfi,ipfi,0.0,spin,ipin,DMIX2,nent,elent,convt)
         Z = Z + (1. + alpha) * (GG(1)+GG(2))
        ELSE !MODE.NE.0             ! Primary transitions (the same fluctuation)
         DO ITT=IT1,IT2
!          G1=GAUSS(ISEED)
!          G2=GAUSS(ISEED)+CORRI(MODE)*G1
!          GSQ=(Re2Res(mode)*G2*G2/(1+CORRI(mode)**2)+
!     *         Im2Res(mode)*G1*G1)/(Re2Res(mode)+Im2Res(mode))
          GSQ=CHISQR(NOPTFL,ISEED,U,IFLAG) !originaly GSQ=GAUSS(ISEED,U,IFLAG)
          Z=Z+GSQ*SIMPL(ITT-IT1+1)
          GG(ITT-IT1+1)=GSQ*SIMPL(ITT-IT1+1)
         ENDDO
         IF (GG(1).GT.0.0) THEN
          DMIX2 = GG(2) / GG(1)
         ELSE
          DMIX2 = 0.0
         ENDIF
         alpha = ALPH_TOT(EG,spfi,ipfi,0.0,spin,ipin,DMIX2,nent,elent,convt)
         Z = Z + (1. + alpha) * (GG(1)+GG(2))
        ENDIF !end of MODE.EQ.0
       ELSE ! in IF (NOPTFL.GE.1) THEN
        DO ITT=IT1,IT2
         Z=Z+SIMPL(ITT-IT1+1)
         GG(ITT-IT1+1)=SIMPL(ITT-IT1+1)
        ENDDO
        IF (GG(1).GT.0.0) THEN
         DMIX2 = GG(2) / GG(1)
        ELSE
         DMIX2 = 0.0
        ENDIF
        alpha = ALPH_TOT(EG,spfi,ipfi,0.0,spin,ipin,DMIX2,nent,elent,convt)
        Z = Z + (1. + alpha) * (GG(1)+GG(2))
       ENDIF !end of IF (NOPTFL.GE.1) THEN
!
!       IF (IT.EQ.1) THEN
!         GE1=GE1+Z
!       ELSEIF (IT.EQ.3) THEN
!         GM1=GM1+Z
!       ELSEIF (IT.EQ.4) THEN
!         GE2=GE2+Z
!       ELSEIF (IT.EQ.2) THEN           ! This is not fully correct (strength fluctuate)
!         GM1=GM1+Z*SIMPL(1)/(SIMPL(1)+SIMPL(2))    ! But it is acceptable approximation
!         GE2=GE2+Z*SIMPL(2)/(SIMPL(1)+SIMPL(2))
!       ENDIF
!
       IF (DENSITY(EIN,SPIN,IPIN).GT.0.0) THEN
!         WRITE(*,*) MODE,I,ISBS,IPFI,STDIS(MODE,I-1,ISBS,IPFI)
!         write(*,*) 'stdis',stdis(mode,0,isbs,ipfi)
         STDIS(MODE,I,ISBS,IPFI)=STDIS(MODE,I-1,ISBS,IPFI)+Z*SPAC
       ELSE
         STDIS(MODE,I,ISBS,IPFI)=STDIS(MODE,I-1,ISBS,IPFI)
       ENDIF
      ENDDO !I=1,NDIS
      
   22 TOTCON(MODE)=TOTCON(MODE)+STCON(MODE,NBIN,ISBS,IPFI)
      GACON(MODE,ISBS,IPFI)=TOTCON(MODE)
      TOTDIS(MODE)=TOTDIS(MODE)+STDIS(MODE,NDIS(ISUBSC(SPFI),IPFI),ISBS,IPFI)
      GADIS(MODE,ISBS,IPFI)=TOTDIS(MODE)
!      WRITE(*,*) TOTDIS(MODE),STDIS(MODE,NDIS(ISUBSC(SPFI),IPFI),ISBS,IPFI),Z,SPAC
!
    1     CONTINUE
        ENDDO
      ENDDO
      RETURN
END SUBROUTINE WIDTHS_R
!***********************************************************************      
SUBROUTINE ONESTEP(MODE,IPIN,SPIN,IBIN,ILIN,TOTCON,STCON,GACON,ISCON,TOTDIS,STDIS,GADIS,ISDIS,IPFI,&
SPFI,IBFI,ILFI,DMIX2,SIGN,IR,IRX,LEVCON,sall,U,IFLAG,EIN,EFI,IREGI,IC_type,IRCON,IRCONc)
!should be OK, changes EIN,EFI, outputs DMIX2,IC_type,sign
!***********************************************************************
REAL::                                SP,SPF,ALPHA,ALPHAK,ALPHAIPF,SPAC,EG,Z,G1,G2,GSQ,Q1,dummy,RN,G
DOUBLE PRECISION::                    AUX0,DRN,GAC
INTEGER::                             IPF,ISPF,ISEED,ISEEDEN,deaux,ISP,ISBS,IT,I,NL,IL,ITT,IT1,IT2,NLEV1
REAL,dimension(1:2)::                 SIMPL,GG
real::                                U,SPIN,EIN,EFI,SPFI,DMIX2,SIGN
integer::                             MODE,IPIN,IBIN,ILIN,IFLAG,IREGI,IPFI,IBFI,ILFI,IR,IRX,IC_type
double precision,dimension(0:2,-2:2,0:1)::     GACON
double precision,dimension(:,:,:,:),allocatable::  STCON
double precision,dimension(0:2)::     TOTCON
integer,dimension(:,:,:,:),allocatable:: ISCON
real, dimension(0:2,0:20,-2:2,0:1)::  STDIS
real, dimension(0:2,-2:2,0:1)::       GADIS
real, dimension(0:2)::                TOTDIS
integer,dimension(:),allocatable::    IRCONc,IRCON
integer,dimension(:,:,:),allocatable::LEVCON
real,dimension(:,:),allocatable::     sall
integer,dimension(:,:,:,:),allocatable::ISDIS
!globalni NBIN,BN,DELTA,NOPTFL,CORRI,Re2Res,Im2Res,ecrit,NENT,ELENT,CONVT,NENK,ELENK,CONVK,NEN_IPF,ELEN_IPF,CONV_IPF,NDIS,dekod,delev,despin,deparity,deltx
      IF (iregi.eq.2) GOTO 23
      SPAC=1./DENSITY(EIN,SPIN,IPIN)
      IF (SPAC.LE.0.0) SPAC=0.00001
!
    5 DRN=(DBLE(INT(DBLE(RAN0(IR))/1.D-4))+DBLE(RAN0(IR)))*1.D-4*(TOTCON(MODE)+DBLE(TOTDIS(MODE)))
      IF (DRN-TOTCON(MODE)) 1,2,2
!
!     Label #1 means that the transition ends in the continuum; this
!     can be learnt from IREGI=0.
!
    1 IREGI=0
      SP  = SPIN - INT(SPIN + .25)
      ISP = INT(SPIN + .25)
      DO IPF=0,1
       DO ISPF = ISP-2, ISP+2
        SPF = SP + FLOAT(ISPF)
        ISBS=NINT(SPF+.25)-NINT(SPIN+.25)
        IT=ITYPE(SPIN,IPIN,SPF,IPF)
        IF (IT.NE.0) THEN
          GAC=GACON(MODE,ISBS,IPF)
          IF(DRN.LT.GAC) GOTO 4
        ENDIF
       ENDDO !ISPF
      ENDDO !IPF
!
!     GOTO 5 statements are here for security reasons
!
      GOTO 5
!
!     Now IPFI,SPFI (& ISBS) are already determined:
!
    4 IPFI=IPF
      SPFI=SPF
      DRN=DRN-GAC+STCON(MODE,NBIN,ISBS,IPFI)
      IF (IT.EQ.2) THEN
        IT1=3
        IT2=4
      ELSE
        IT1=IT
        IT2=IT
      ENDIF
!
!     A DO-loop that follows should be later replaced by a faster
!     algorithm
!
      DO I=IBIN+1,NBIN
        IF (STCON(MODE,I,ISBS,IPFI).GT.DRN) GOTO 7
      ENDDO
      GOTO 5
    7 IBFI=I
!
!     Now even IBFI is known. It remains to determine ILFI.
!
      ISEED=ISCON(MODE,IBFI,ISBS,IPFI)
      IFLAG=0
!
!     The appropriate ISEED is now fetched and IFLAG reset.
!     The algorithm for finding ILFI can start.
!
      EG=EIN-BN+(FLOAT(IBFI)-0.5)*DELTA
!      EG=(FLOAT(IBFI-IBIN))*DELTA	       !Energies between mid-bins
      NL=LEVCON(IBFI,ISUBSC(SPFI),IPFI)-LEVCON(IBFI-1,ISUBSC(SPFI),IPFI)
      AUX0=STCON(MODE,IBFI-1,ISBS,IPFI)
      Z=0.
      DO ITT=IT1,IT2
        SIMPL(ITT-IT1+1)=SGAMMA(EG,EIN,ITT)
      ENDDO
      DO IL=1,NL
       DO ITT=IT1,IT2
        IF (NOPTFL.GE.1) THEN !originaly .EQ.
         IF (MODE.EQ.0) THEN
          G=GAUSS(ISEED,U,IFLAG)
          Z=Z+G*G*SIMPL(ITT-IT1+1)
          GG(ITT-IT1+1)=G*G*SIMPL(ITT-IT1+1)
         ELSE                                 !Primary transitions
!          G1=GAUSS(ISEED)
!          G2=GAUSS(ISEED)+CORRI(MODE)*G1
!          GSQ=(Re2Res(mode)*G2*G2/(1+CORRI(mode)**2)+
!     *         Im2Res(mode)*G1*G1)/(Re2Res(mode)+Im2Res(mode))
          GSQ=CHISQR(NOPTFL,ISEED,U,IFLAG) !originaly GSQ=GAUSS(ISEED,U,IFLAG)
          Z=Z+GSQ*SIMPL(ITT-IT1+1)
          GG(ITT-IT1+1)=GSQ*SIMPL(ITT-IT1+1)
         ENDIF
        ELSE
          Z=Z+SIMPL(ITT-IT1+1)
          GG(ITT-IT1+1)=SIMPL(ITT-IT1+1)
        ENDIF
       ENDDO !ITT
       IF ((GG(2).GT.0.0).OR.(GG(1).GT.0.0)) THEN
         IF (GG(1).GT.0.0) THEN
           DMIX2 = GG(2) / GG(1)
         ELSE
           DMIX2 = 0.0
         ENDIF
         alpha = ALPH_TOT(EG,spfi,ipfi,0.0,spin,ipin,DMIX2,nent,elent,convt)
         Z = Z + (1. + alpha) * (GG(1)+GG(2))
       ENDIF
       IF ((AUX0+DBLE(Z*SPAC)).GT.DRN) GOTO 9
      ENDDO !IL
      GO TO 5
!
!     The resulting ILFI:
!
    9 ILFI=IL
!
      NLEV1=LEVCON(IBFI,ISUBSC(SPFI),IPFI)-LEVCON(IBFI-1,ISUBSC(SPFI),IPFI)
      Q1=1.-FLOAT(2*ILFI-1)/FLOAT(2*NLEV1)
!      ISEEDEN=SEEDS(0,SPFI,IPFI,ILFI,IBFI,LEVCON,IRCON,IRCONc)
!      Q1=RAN0(ISEEDEN)   ! Energy of the final level random in bin
      EFI=BN-(FLOAT(IBFI)-Q1)*DELTA
      IF (efi.le.ecrit) IREGI=2
!
!       delta**2 and conversion
!
      IF (NOPTFL.LT.1) THEN !NOTE originaly .NE.
        GG(1)=1.
        GG(2)=1.
      ENDIF
      IF (IT.EQ.2) THEN
        IF (GG(1).NE.0.) THEN
          DMIX2=GG(2)**2*SIMPL(2)/GG(1)**2/SIMPL(1)
          IF (GG(1).GT.0) THEN
            SIGN=1.
          ELSE
            SIGN=-1.
          ENDIF
        ELSE
          DMIX2=1.E+6
          SIGN=1.
        ENDIF
      ELSE
        DMIX2=0.
        SIGN=1.
      ENDIF
      dummy=RAN0(IRX)
      ALPHA=ALPH_TOT(EIN,SPIN,IPIN,EFI,SPFI,IPFI,DMIX2,NENT,ELENT,CONVT)
      IF (dummy.GT.(ALPHA/(1+ALPHA))) THEN
        IC_type=0
      ELSE
        ALPHAK=ALPH_TOT(EIN,SPIN,IPIN,EFI,SPFI,IPFI,DMIX2,NENK,ELENK,CONVK)
        ALPHAIPF=ALPH_TOT(EIN,SPIN,IPIN,EFI,SPFI,IPFI,DMIX2,NEN_IPF,ELEN_IPF,CONV_IPF)+ALPHAK
        IF (dummy.LE.(ALPHAK/(1+ALPHA))) THEN
          IC_type=1
        ELSEIF (dummy.LE.(ALPHAIPF/(1+ALPHA))) THEN
          IC_type=3
        ELSE
          IC_type=2
        ENDIF
      ENDIF
      RETURN
!
    2 IREGI=1
      RN=SNGL(DRN-TOTCON(MODE))
      SP  = SPIN - INT(SPIN + .25)
      ISP = INT(SPIN + .25)
      DO IPF=0,1
       DO ISPF = ISP-2, ISP+2
        SPF = SP + FLOAT(ISPF)
        ISBS=NINT(SPF+.25)-NINT(SPIN+.25)
        IT=ITYPE(SPIN,IPIN,SPF,IPF)
        IF (IT.NE.0) THEN
         IF (RN.LT.GADIS(MODE,ISBS,IPF)) GOTO 11
        ENDIF
       ENDDO !SPF
      ENDDO !IPF
      GOTO 5
   11 SPFI=SPF
      IPFI=IPF
!
      RN=RN-GADIS(MODE,ISBS,IPFI)+STDIS(MODE,NDIS(ISUBSC(SPFI),IPFI),ISBS,IPFI)
      DO I=1,NDIS(ISUBSC(SPFI),IPFI)
        IF (RN.LT.STDIS(MODE,I,ISBS,IPFI)) GOTO 15
      ENDDO
      GOTO 5
   15 ILFI=I
      ibin=0                                        !!!!Why
      EFI=ENDIS(ILFI,ISUBSC(SPFI),IPFI)
      IF (EFI.LE.ecrit) IREGI=2
!
!     delta**2 and conversion
!
      ISEED=ISDIS(MODE,ILFI,ISUBSC(SPFI),IPFI)
      IFLAG=0
      IF (IT.EQ.2) THEN
         IT1=3
         IT2=4
         ELSE
         IT1=IT
         IT2=IT
      ENDIF
      EG=EIN-EFI
      DO ITT=IT1,IT2
        SIMPL(ITT-IT1+1)=SGAMMA(EG,EIN,ITT)
      ENDDO
      IF (NOPTFL.GE.1) THEN  !originaly .EQ. ,TODO modify for chi2 with more DOF
       IF (MODE.EQ.0) THEN
        DO ITT=IT1,IT2
         G=GAUSS(ISEED,U,IFLAG)
         GG(ITT-IT1+1)=G
        ENDDO
       ELSE                             !Primary transitions
        DO ITT=IT1,IT2
!         G1=GAUSS(ISEED)
!         G2=GAUSS(ISEED)+CORRI(MODE)*G1
!         GSQ=(Re2Res(mode)*G2*G2/(1+CORRI(mode)**2)+
!     *        Im2Res(mode)*G1*G1)/(Re2Res(mode)+Im2Res(mode))
!         GG(ITT-IT1+1)=sqrt(GSQ)
!         IF (G1.LT.0) GG(ITT-IT1+1)=-1.*GG(ITT-IT1+1)
         G=CHISQR(NOPTFL,ISEED,U,IFLAG) !originaly G=GAUSS(ISEED,U,IFLAG)
         GG(ITT-IT1+1)=SQRT(G)
        ENDDO
       ENDIF
      ENDIF
      IF (NOPTFL.LT.1) THEN !originaly .NE.
        GG(1)=1.
        GG(2)=1.
      ENDIF
      IF (IT.EQ.2) THEN
        IF (GG(1).NE.0.) THEN
          DMIX2=GG(2)**2*SIMPL(2)/GG(1)**2/SIMPL(1)
          IF (GG(1).GT.0) THEN
            SIGN=1.
          ELSE
            SIGN=-1.
          ENDIF
        ELSE
          DMIX2=1.E+6
          SIGN=1.
        ENDIF
      ELSE
        DMIX2=0.
        SIGN=1.
      ENDIF
      dummy=RAN0(IRX)
      ALPHA=ALPH_TOT(EIN,SPIN,IPIN,EFI,SPFI,IPFI,DMIX2,NENT,ELENT,CONVT)
      IF (dummy.GT.(ALPHA/(1+ALPHA))) THEN
        IC_type=0
      ELSE
        ALPHAK=ALPH_TOT(EIN,SPIN,IPIN,EFI,SPFI,IPFI,DMIX2,NENK,ELENK,CONVK)
        ALPHAIPF=ALPH_TOT(EIN,SPIN,IPIN,EFI,SPFI,IPFI,DMIX2,NEN_IPF,ELEN_IPF,CONV_IPF)+ALPHAK
        IF (dummy.LE.(ALPHAK/(1+ALPHA))) THEN
          IC_type=1
        ELSEIF (dummy.LE.(ALPHAIPF/(1+ALPHA))) THEN
          IC_type=3
        ELSE
          IC_type=2
        ENDIF
      ENDIF
      RETURN
!
!
!     Note that this time ILFI means order number of an actual
!     level, not the order number of a level in an energy bin of
!     the energy continuum.
!
!     End in a part, where one knows every branching ratios
!
  23  EIN=EFI
      deaux=dekod(ilfi,isubsc(spfi),ipfi)
      rn=ran0(ir)*sall(deaux,denum(deaux))
      DO i=1,denum(deaux)
        IF (rn.lt.sall(deaux,i)) GOTO 24
      ENDDO
      GOTO 23
  24  ilfi=delev(deaux,i)
      spfi=despin(deaux,i)
      ipfi=deparity(deaux,i)
      ibin=0
      efi=endis(ilfi,ISUBSC(spfi),ipfi)
      dmix2=deltx(deaux,i)**2
      IF (deltx(deaux,i).LT.0) THEN
        sign=-1.
      ELSE
        sign=1.
      ENDIF
!
!     conversion
!      IC_type !0 = gamma, 1 = K-shell, 2 = higher-shell, 3 = pair
      if (p_conv(deaux,i).gt.0.0) then
         dummy=ran0(irx)
         if (dummy.gt.p_conv(deaux,i)) then
           IC_type=0
         else
           if (dummy.le.p_conv_K(deaux,i)) then
             IC_type=1
           elseif (dummy.le.p_conv_IPF(deaux,i)) then
             IC_type=3
           else
             IC_type=2
           endif  
         endif
      endif
      RETURN
END SUBROUTINE ONESTEP
!***********************************************************************
SUBROUTINE SPECTRA (INR,ELQQ,SPQQ,DMQQ,IPQQ,ICQQ,NR_STEPS,intermediate2,intermediate3,XSC_work,gamma_multiplicita)
!***********************************************************************
integer:: INR
real,dimension(:,:),allocatable:: ELQQ,SPQQ,DMQQ
integer,dimension(:),allocatable::   NR_STEPS
integer,dimension(:,:),allocatable:: IPQQ,ICQQ,gamma_multiplicita
integer,dimension(:,:,:),allocatable:: intermediate2
integer,dimension(:,:,:,:),allocatable:: intermediate3,XSC_work
INTEGER::              IEVENTS,J,K,I_MSC_FS,MULTIPLICITA,REAL_MULT,SUM_CONV
REAL,DIMENSION(:),allocatable::    E_G
!
      if (.not.allocated(E_G)) then
        allocate(E_G(1:126))
      endif
      DO I_MSC_FS=1,N_MSC_FS
        DO IEVENTS=1,NEVENTS
          REAL_MULT=NR_STEPS(IEVENTS)
          DO J=1,NR_STEPS(IEVENTS)
            E_G(J)=ELQQ(IEVENTS,J-1)-ELQQ(IEVENTS,J)
            IF (ICQQ(IEVENTS,J).NE.0) THEN
              E_G(J)=0.0 !asi nepotrebny radek
              REAL_MULT=REAL_MULT-1
            ENDIF
          ENDDO
          IF (I_MSC_FS.EQ.1) THEN
            gamma_multiplicita(INR,MIN(REAL_MULT,MAX_MULTIPLICITA))=gamma_multiplicita(INR,MIN(REAL_MULT,MAX_MULTIPLICITA))+1
          ENDIF

          DO MULTIPLICITA=MIN_MULTIPLICITA,MAX_MULTIPLICITA
            IF (NR_STEPS(IEVENTS).GE.MULTIPLICITA) THEN
              SUM_CONV=0
              DO J=1,NR_STEPS(IEVENTS)
                SUM_CONV=SUM_CONV+ICQQ(IEVENTS,J)
                IF ((ELQQ(IEVENTS,J).EQ.MSC_FS(I_MSC_FS)).AND.(J.EQ.MULTIPLICITA).AND.(SUM_CONV.EQ.0)) THEN !dostal jsem se na finalni hladinu v J krocich a nemel jsem konverzni elektron
!cut                IF ((ELQQ(IEVENTS,J).EQ.MSC_FS(I_MSC_FS)).AND.(J.EQ.MULTIPLICITA).AND.(SUM_CONV.EQ.0).AND.(E_G(J).GE.0.2)) THEN !dostal jsem se na finalni hladinu ve 3 krocich a nemel jsem konverzni elektron
                  IF (MULTIPLICITA.EQ.2) THEN
                    intermediate2(I_MSC_FS,INR,INT(ELQQ(IEVENTS,1)/BIN_WIDTH))=&
intermediate2(I_MSC_FS,INR,INT(ELQQ(IEVENTS,1)/BIN_WIDTH))+1
                  ELSEIF (MULTIPLICITA.EQ.3) THEN
                    intermediate3(I_MSC_FS,INR,INT(ELQQ(IEVENTS,1)/BIN_WIDTH),&
INT(ELQQ(IEVENTS,2)/BIN_WIDTH))=intermediate3(I_MSC_FS,INR,INT(ELQQ(IEVENTS,1)/BIN_WIDTH),INT(ELQQ(IEVENTS,2)/BIN_WIDTH))+1
                  ENDIF
                  DO K=1,MULTIPLICITA
                    XSC_work(I_MSC_FS,INR,MULTIPLICITA,INT(E_G(K)/BIN_WIDTH))=&
XSC_work(I_MSC_FS,INR,MULTIPLICITA,INT(E_G(K)/BIN_WIDTH))+1
                  ENDDO
                ENDIF
              ENDDO
            ENDIF
          ENDDO !MULTIPLICITA
        ENDDO !NEVENTS
      ENDDO !N_MSC_FS
      RETURN
END SUBROUTINE SPECTRA
!************************************************************************
SUBROUTINE WR_SPECTRA (intermediate2,intermediate3,XSC_work,gamma_multiplicita)

integer,dimension(:,:),allocatable::  gamma_multiplicita
integer,dimension(:,:,:),allocatable:: intermediate2
integer,dimension(:,:,:,:),allocatable:: intermediate3,XSC_work
INTEGER,DIMENSION(:),allocatable:: POCET_BINU
REAL,DIMENSION(:),ALLOCATABLE::  G_MULT,G_MULT_ERR
REAL,DIMENSION(:,:),ALLOCATABLE::  XSC_MEAN,XSC_ERR  !todo pripadne DP
REAL,DIMENSION(:,:,:),ALLOCATABLE::  int_mean
CHARACTER(1)::   cislosouboru,cislostavu
CHARACTER(7)::   final_mult
INTEGER::        I_MSC_FS,MULTIPLICITA,I,J,INR

      if ((.not.allocated(POCET_BINU)).AND.(N_MSC_FS.GE.1)) then
        allocate(POCET_BINU(1:N_MSC_FS))
      endif
      if (.not.allocated(G_MULT)) then
        allocate(G_MULT(0:MAX_MULTIPLICITA))
      endif
      if (.not.allocated(G_MULT_ERR)) then
        allocate(G_MULT_ERR(0:MAX_MULTIPLICITA))
      endif
      DO MULTIPLICITA=0,MAX_MULTIPLICITA
        G_MULT(MULTIPLICITA)=0.0
        G_MULT_ERR(MULTIPLICITA)=0.0
      ENDDO
!      write(*,*) (G_MULT(MULTIPLICITA),MULTIPLICITA=0,MAX_MULTIPLICITA)
      write(final_mult,'(F7.5)') MSC_FS(1)
      OPEN (UNIT=12,FILE='multiplicity_'//final_mult//'.dat',ACCESS='SEQUENTIAL',STATUS='UNKNOWN')
      DO MULTIPLICITA=0,MAX_MULTIPLICITA
!        write(*,*) MULTIPLICITA,':'
        DO INR=1,(NREAL*NSUB)
          G_MULT(MULTIPLICITA)=G_MULT(MULTIPLICITA)+gamma_multiplicita(INR,MULTIPLICITA)
          G_MULT_ERR(MULTIPLICITA)=G_MULT_ERR(MULTIPLICITA)+gamma_multiplicita(INR,MULTIPLICITA)**2
!          write(*,*) gamma_multiplicita(INR,MULTIPLICITA)
        ENDDO
!        write(*,*) G_MULT(MULTIPLICITA)
        G_MULT(MULTIPLICITA)=G_MULT(MULTIPLICITA)/(NREAL*NSUB)
        G_MULT_ERR(MULTIPLICITA)=G_MULT_ERR(MULTIPLICITA)/(NREAL*NSUB)-G_MULT(MULTIPLICITA)**2
        WRITE(12,*) MULTIPLICITA,G_MULT(MULTIPLICITA),SQRT(G_MULT_ERR(MULTIPLICITA)*FLOAT((NREAL*NSUB))/FLOAT((NREAL*NSUB)-1))
      ENDDO
      CLOSE (12)

      DO I_MSC_FS=1,N_MSC_FS
        POCET_BINU(I_MSC_FS)=INT((BN-MSC_FS(I_MSC_FS))/BIN_WIDTH)+1
        if (.not.allocated(XSC_MEAN)) then
          allocate(XSC_MEAN(MIN_MULTIPLICITA:MAX_MULTIPLICITA,0:(POCET_BINU(I_MSC_FS)-1)))
        else 
          write(*,*) 'neco je spatne'
        endif
        if (.not.allocated(XSC_ERR)) then
          allocate(XSC_ERR(MIN_MULTIPLICITA:MAX_MULTIPLICITA,0:(POCET_BINU(I_MSC_FS)-1)))
        else 
          write(*,*) 'neco je spatne'
        endif
        if (.not.allocated(int_mean)) then
          allocate(int_mean(MIN_MULTIPLICITA:MAX_MULTIPLICITA,0:INT(BN/BIN_WIDTH),0:INT(BN/BIN_WIDTH)))
        else 
          write(*,*) 'neco je spatne'
        endif
          DO MULTIPLICITA=MIN_MULTIPLICITA,MAX_MULTIPLICITA
             !TODO vypsat int_mean
                  IF (MULTIPLICITA.EQ.2) THEN
                    DO I=0,INT(BN/BIN_WIDTH)
                      DO INR=1,(NREAL*NSUB)
                        int_mean(MULTIPLICITA,I,0)=int_mean(MULTIPLICITA,I,0)+intermediate2(I_MSC_FS,INR,I)
                      ENDDO
                    ENDDO
                  ELSEIF (MULTIPLICITA.EQ.3) THEN
                    DO I=0,INT(BN/BIN_WIDTH)
                    DO J=0,INT(BN/BIN_WIDTH)
                      DO INR=1,(NREAL*NSUB)
                        int_mean(MULTIPLICITA,I,J)=int_mean(MULTIPLICITA,I,J)+intermediate3(I_MSC_FS,INR,I,J)
                      ENDDO
                    ENDDO
                    ENDDO
                  ENDIF
            write(cislosouboru,'(I1.1)') MULTIPLICITA
            write(cislostavu,'(I1.1)') I_MSC_FS
            OPEN (UNIT=12,FILE=cislostavu//'.'//cislosouboru,ACCESS='SEQUENTIAL',STATUS='UNKNOWN')
            DO I=0,(POCET_BINU(I_MSC_FS)-1)
              XSC_MEAN(MULTIPLICITA,I)=0.0
              XSC_ERR(MULTIPLICITA,I)=0.0
              DO INR=1,(NREAL*NSUB)
                XSC_MEAN(MULTIPLICITA,I)=XSC_MEAN(MULTIPLICITA,I)+XSC_work(I_MSC_FS,INR,MULTIPLICITA,I)
                XSC_ERR(MULTIPLICITA,I)=XSC_ERR(MULTIPLICITA,I)+XSC_work(I_MSC_FS,INR,MULTIPLICITA,I)**2
              ENDDO
              XSC_MEAN(MULTIPLICITA,I)=XSC_MEAN(MULTIPLICITA,I)/(NREAL*NSUB)
              XSC_ERR(MULTIPLICITA,I)=XSC_ERR(MULTIPLICITA,I)/(NREAL*NSUB)-XSC_MEAN(MULTIPLICITA,I)**2
              XSC_ERR(MULTIPLICITA,I)=SQRT(XSC_ERR(MULTIPLICITA,I)*FLOAT((NREAL*NSUB))/FLOAT((NREAL*NSUB)-1))
              WRITE(12,*) I,XSC_MEAN(MULTIPLICITA,I),XSC_ERR(MULTIPLICITA,I)
            ENDDO !I
            CLOSE (12)
          ENDDO !MULTIPLICITA
        deallocate(XSC_MEAN)
        deallocate(XSC_ERR)
        deallocate(int_mean)
      ENDDO !I_MSC_FS

      if (allocated(POCET_BINU)) then
        deallocate (POCET_BINU)
      else
        write(*,*) 'sth funky goin on'
      endif
      RETURN
END SUBROUTINE WR_SPECTRA
!************************************************************************
SUBROUTINE DO_IT (INS,INR,ELQQ,SPQQ,DMQQ,IPQQ,ICQQ,WIQQ,NR_STEPS,ITID)
!***********************************************************************
integer:: INS,INR,ITID
real,dimension(:,:),allocatable:: ELQQ,SPQQ,DMQQ,WIQQ
integer,dimension(:,:),allocatable:: IPQQ,ICQQ
integer,dimension(:),allocatable::   NR_STEPS
INTEGER::              IEV,J,REAL_MULT
CHARACTER*3  EXTR,EXTS
CHARACTER*4  EXT
REAL,DIMENSION(:),allocatable::    E_G
!
      if (.not.allocated(E_G)) then
        allocate(E_G(1:126))
      endif
!
!      IERR=0
      IF (ISWBN.LT.0.OR.ISWBN.GT.4) THEN
!         IERR=1
         RETURN
      ENDIF
      WRITE (EXTS,101) INS
      WRITE (EXTR,101) INR
      WRITE (EXT,102) (INS+(INR-1)*NREAL)
  101 FORMAT (I3.3)
  102 FORMAT (I4.4)
      IF (ISWBN.EQ.1) THEN
         OPEN (UNIT=12+ITID,FILE='EVENTS.S'//EXTS//'.R'//EXTR,FORM='UNFORMATTED',ACCESS='SEQUENTIAL',STATUS='NEW')
      ELSEIF (ISWBN.EQ.2) THEN
         OPEN (UNIT=12+ITID,FILE='EVS.'//EXT,STATUS='UNKNOWN')
      ELSE
         OPEN (UNIT=12+ITID,FILE='EVENTS.S'//EXTS//'.R'//EXTR,STATUS='UNKNOWN')
      ENDIF
      DO IEV=1,NEVENTS
       IF (ISWBN.EQ.0) THEN !used for DANCE so far, needed conversion by db_g4.exe
!         WRITE (12,*)   ' '
         WRITE (12+ITID,100) SPQQ(IEV,0),IPQQ(IEV,0),NR_STEPS(IEV),ELQQ(IEV,0)
         IF (ISWEL.EQ.1) WRITE (12+ITID,200) (ELQQ(IEV,J),J=1,NR_STEPS(IEV))
         IF (ISWSP.EQ.1) WRITE (12+ITID,300) (SPQQ(IEV,J),J=1,NR_STEPS(IEV))
         IF (ISWPA.EQ.1) WRITE (12+ITID,400) (IPQQ(IEV,J),J=1,NR_STEPS(IEV))
         IF (ISWIC.EQ.1) WRITE (12+ITID,400) (ICQQ(IEV,J),J=1,NR_STEPS(IEV))
         IF (ISWMX.EQ.1) WRITE (12+ITID,500) (DMQQ(IEV,J),J=1,NR_STEPS(IEV))
         IF (ISWWI.EQ.1) WRITE (12+ITID,600) (WIQQ(IEV,J),J=1,NR_STEPS(IEV))
       ELSEIF (ISWBN.EQ.1) THEN !TODO ask MK about this
         WRITE (12+ITID) SPQQ(IEV,0),IPQQ(IEV,0),NR_STEPS(IEV),ELQQ(IEV,0)
         IF (ISWEL.EQ.1) WRITE (12+ITID) (ELQQ(IEV,J),J=1,NR_STEPS(IEV))
         IF (ISWSP.EQ.1) WRITE (12+ITID) (SPQQ(IEV,J),J=1,NR_STEPS(IEV))
         IF (ISWPA.EQ.1) WRITE (12+ITID) (IPQQ(IEV,J),J=1,NR_STEPS(IEV))
         IF (ISWIC.EQ.1) WRITE (12+ITID) (ICQQ(IEV,J),J=1,NR_STEPS(IEV))
         IF (ISWMX.EQ.1) WRITE (12+ITID) (DMQQ(IEV,J),J=1,NR_STEPS(IEV))
         IF (ISWWI.EQ.1) WRITE (12+ITID) (WIQQ(IEV,J),J=1,NR_STEPS(IEV))
       ELSEIF (ISWBN.EQ.2) THEN !DANCE GEANT4 input format
         REAL_MULT=0 ! IC_type 0 = gamma, 1 = K-shell, 2 = higher-shell, 3 = pair
         DO J=1,NR_STEPS(IEV)
           IF (ICQQ(IEV,J).EQ.0) THEN
             REAL_MULT=REAL_MULT+1
             E_G(REAL_MULT)=ELQQ(IEV,J-1)-ELQQ(IEV,J)
           ELSEIF (ICQQ(IEV,J).EQ.1) THEN
             REAL_MULT=REAL_MULT+1
             E_G(REAL_MULT)=xrayk-xrayl
           ELSEIF (ICQQ(IEV,J).EQ.3) THEN
             REAL_MULT=REAL_MULT+2
             E_G(REAL_MULT-1)=.511
             E_G(REAL_MULT)=.511
           ENDIF
         ENDDO
         IF (REAL_MULT.GT.0) WRITE (12+ITID,700) REAL_MULT, (E_G(J),J=1,REAL_MULT)
       ELSEIF (ISWBN.EQ.4) THEN !DANCE GEANT4 input format with no X-rays or anihilation gammas
         REAL_MULT=0 ! IC_type 0 = gamma, 1 = K-shell, 2 = higher-shell, 3 = pair
         DO J=1,NR_STEPS(IEV)
           IF (ICQQ(IEV,J).EQ.0) THEN
             REAL_MULT=REAL_MULT+1
             E_G(REAL_MULT)=ELQQ(IEV,J-1)-ELQQ(IEV,J)
           ENDIF
         ENDDO
         IF (REAL_MULT.GT.0) WRITE (12+ITID,700) REAL_MULT, (E_G(J),J=1,REAL_MULT)
       ELSEIF (ISWBN.EQ.3) THEN !TAC @ n_TOF tabuled GEANT4 input format
         REAL_MULT=0
         DO J=1,NR_STEPS(IEV)
           IF (ICQQ(IEV,J).EQ.0) THEN
             REAL_MULT=REAL_MULT+1
             E_G(REAL_MULT)=ELQQ(IEV,J-1)-ELQQ(IEV,J)
           ELSEIF (ICQQ(IEV,J).EQ.1) THEN
             REAL_MULT=REAL_MULT+1
             E_G(REAL_MULT)=xrayk-xrayl
           ELSEIF (ICQQ(IEV,J).EQ.3) THEN
             REAL_MULT=REAL_MULT+2
             E_G(REAL_MULT-1)=.511
             E_G(REAL_MULT)=.511
           ENDIF
         ENDDO
         IF (REAL_MULT.GT.0) THEN
           WRITE (12+ITID,800) REAL_MULT,0
           WRITE (12+ITID,900) (1000*E_G(J),J=1,REAL_MULT)
         ENDIF !TODO co se deje kdyz je kaskada tvorena jen elektrony
       ENDIF
      ENDDO
  100 FORMAT (F5.1,I2,I4,2X,F9.5)
  200 FORMAT (2X,126F9.5)
  300 FORMAT (2X,126F9.1)
  400 FORMAT (2X,126I9)
  500 FORMAT (2X,126E11.3)
  600 FORMAT (2X,126E9.2)
  700 FORMAT (1X,I2,126F9.5)
  800 FORMAT (2I3)
  900 FORMAT (126F10.3)

      CLOSE (UNIT=12+ITID,STATUS='KEEP')
      RETURN
END SUBROUTINE DO_IT
!
!**********************************************************************
SUBROUTINE LABELS(tdens,tsfe1,tsfm1,tsfe2)
!************************************************************************
CHARACTER*14 tdens,tsfe1,tsfm1,tsfe2
  !
  !    Conversion from the model number to model's string label
  !
        tdens='???'
        tsfe1='???'
        tsfm1='???'
        tsfe2='???'
  
        IF (NOPTDE.EQ.0) THEN
          tdens='CTF'
        ELSEIF (NOPTDE.EQ.8) THEN
          tdens='CTF-vE2009' 
        ELSEIF ((NOPTDE.EQ.1).OR.(NOPTDE.EQ.4).OR.(NOPTDE.EQ.5)) THEN
          tdens='BSFG'
        ELSEIF (NOPTDE.EQ.6) THEN
          tdens='BSFG-vE2005'
        ELSEIF (NOPTDE.EQ.9) THEN
          tdens='BSFG-vE2009'
        ELSEIF (NOPTDE.EQ.11) THEN
          tdens='(Goriely) tabs' 
        ENDIF
  
        IF (NOPTE1.EQ.0) THEN
          tsfe1='SP'
        ELSEIF (NOPTE1.EQ.1) THEN
          tsfe1='SLO'
        ELSEIF (NOPTE1.EQ.2) THEN
          tsfe1='ELO'
        ELSEIF (NOPTE1.EQ.3) THEN
          tsfe1='EGLO'
        ELSEIF (NOPTE1.EQ.4) THEN
          tsfe1='KMF'
        ELSEIF (NOPTE1.EQ.5) THEN
          tsfe1='GLO'
        ELSEIF (NOPTE1.EQ.6) THEN
          tsfe1='MGLO'
        ELSEIF (NOPTE1.EQ.7) THEN
          tsfe1='EELO'
        ELSEIF (NOPTE1.EQ.8) THEN
          tsfe1='1SLO+KMF'
        ELSEIF (NOPTE1.EQ.9) THEN
          tsfe1='1SLO+MGLO'
        ELSEIF (NOPTE1.EQ.10) THEN
          tsfe1='KMF->SLO'
        ELSEIF (NOPTE1.EQ.11) THEN
          tsfe1='(Goriely) tabs'
        ENDIF  !TODO continue updating these and decide if to be used
  
        IF (NOPTM1.EQ.0) THEN
          tsfm1='SP'
        ELSEIF (NOPTM1.EQ.1) THEN
          tsfm1='SLO'
        ELSEIF (NOPTM1.EQ.4) THEN
          tsfm1='SLO+SP'
        ELSEIF (NOPTM1.EQ.11) THEN
          tsfm1='(Goriely) tabs'
        ENDIF
  
        IF (NOPTE2.EQ.0) THEN
          tsfe2='SP'
        ELSEIF (NOPTE2.EQ.1) THEN
          tsfe2='SLO'
        ELSEIF (NOPTE2.EQ.11) THEN
          tsfe2='(Goriely) tabs'
        ENDIF
  
END SUBROUTINE LABELS 
!************************************************************************
SUBROUTINE WRITE_PARAMS
!***********************************************************************
CHARACTER*14                      tdens,tsfe1,tsfm1,tsfe2
CHARACTER*1                       CHPAR
INTEGER::                         i
!
!        Writing down some input parameters to output file
!
      CALL LABELS(tdens,tsfe1,tsfm1,tsfe2)
      OPEN(9,FILE='DICE.PRO',STATUS='UNKNOWN')
      WRITE(9,118) INT(ZNUM),INT(AMASS)
  118 FORMAT('Simulation for compound nucleus Z= ',I3,' A= ',I3)
      WRITE(9,*) 'Used models:'
      WRITE(9,*) ' LD             E1             M1             E2'
      WRITE(9,116) NOPTDE, NOPTE1, NOPTM1, NOPTE2
  116 FORMAT(1X,I3,12X,I3,12X,I3,12X,I3)
      WRITE(9,119) tdens,tsfe1,tsfm1,tsfe2
  119 FORMAT(1X,A14,1X,A14,1X,A14,1X,A14)
      WRITE(9,*) 'Parameters of PSFs:'
      WRITE(9,*) 'E1 (GDER): Erez[MeV]  Width[MeV]   Sigma[mb]'
      DO i=1,NGIGE+NLOWLOR
        WRITE(9,111) ER(i),W0(i),SIG(i)
      ENDDO
  111 FORMAT('      ',2F11.2,F13.1)
      WRITE(9,*) 'M1 (GDMR): Erez[MeV]  Width[MeV]   Sigma[mb]'
      DO i=1,NGIGM
        WRITE(9,111) ERM(i),WM0(i),SIGM(i)
      ENDDO
      WRITE(9,*) 'E2 (GQER): Erez[MeV]  Width[MeV]   Sigma[mb]'
      DO i=1,NGIGE2
       WRITE(9,111) ERE(i),WE0(i),SIGE(i)
      ENDDO
      WRITE(9,112) DEG,DMG,QEL
  112 FORMAT(' SP strength:   E1:',E9.2,'   M1:',E9.2,'   E2:',E9.2)
      WRITE(9,120) D0, 1e6/DENSITY(BN,TRGT_SPIN+0.5,TRGT_PI)
  120 FORMAT(' chosen LD yields D_0 =',F10.3,' eV and D_0^+ =',F10.3)    
      IF(NOPTDE.EQ.8) THEN
        WRITE(9,114) TEMPER09,EZERO09
      ELSE
        WRITE(9,114) TEMPER,EZERO
      ENDIF
  114 FORMAT(' Density parameters:  CTF:   Temper:',F6.3,'     E0:',F6.3)
      IF(NOPTDE.EQ.9) THEN
        WRITE(9,115) ASHELL09,DEL09,PAIRING09
      ELSE
        WRITE(9,115) ASHELL,DEL,PAIRING
      ENDIF
  115 FORMAT(' Density parameters: BSFG:   Ashell:',F7.3,'  Del/E1:',F6.3,'  Pairing:',F6.3)
      WRITE(9,113) fermc,ek0,EGZERO
  113 FORMAT(' Fermi liq. par.:',F4.1,'            k0:',F4.1,'     Eg0:',F4.1)
      IF (IPINc.EQ.0) THEN
       CHPAR='+'
      ELSE
       CHPAR='-'
      ENDIF
      WRITE(9,117) BN,SPINc(1),CHPAR,CAPFR(1)
  117 FORMAT(' Binding energy:',F7.4,'         Jpi:',F4.1,A1,'  Fract:',F6.3)
      WRITE(9,*)'**********************************************************'
      CLOSE(UNIT=9,STATUS='KEEP')
      RETURN
END SUBROUTINE WRITE_PARAMS
!***********************************************************************
SUBROUTINE WRITE_DICE_PRO(RADAVG,RADVAR,RADFLUCT,NDEAD,NISOM)
!***********************************************************************
integer::                         NDEAD,NISOM
real::                            RADAVG,RADVAR,RADFLUCT
      OPEN(9,FILE='DICE.PRO',ACCESS='APPEND',STATUS='OLD')
      WRITE(9,*)
      WRITE(9,*) 'Suprarealization #', NREAL,' with ',NSUB,' realizations'
      WRITE(9,*) 'Number of events',NEVENTS
      WRITE(9,*)
      WRITE(9,*) 'Capt.state tot.rad.width (MeV): '
      WRITE(9,*) 'global average +/- overall fluctuation (suprareal , realization)'
      WRITE(9,197) RADAVG,SQRT(RADVAR+RADFLUCT),SQRT(RADVAR),SQRT(RADFLUCT)
  197 format(E13.6,' +/-',E13.6,' (',E13.6,' , ',E13.6,') ')   
      WRITE(9,*)
      WRITE(9,*) 'Number of cascades terminating at isomeric state: ',NISOM
      WRITE(9,*)
      WRITE(9,*) 'Number of cascades terminating at a dead end: ',NDEAD
      WRITE(9,*)
      CLOSE(9)
      RETURN
END SUBROUTINE WRITE_DICE_PRO
!***********************************************************************
SUBROUTINE WRITE_DICE_POPS(POPULT,POPERT,POPVAR,COVAP,POPULS,POPERS,POPSVAR,COVAS)
!TODO navazat pripadne lepsi vypis pro popul_depopul grafy
!***********************************************************************
real,dimension(:,:),allocatable:: POPULT,POPERT,POPULS,POPERS
real,dimension(:),allocatable::   POPVAR,POPSVAR
real,dimension(:,:),allocatable:: COVAP,COVAS
REAL,dimension(:,:),allocatable:: CONO
REAL::                            poper
INTEGER::                         I,K,L,ILVL
      if (.not.allocated(CONO)) then
        allocate(CONO(1:NUMLEV,1:NUMLEV))
      endif

      OPEN(9,FILE='DICE.PRO',ACCESS='APPEND',STATUS='OLD')
      WRITE(9,*) 'Population of low-lying states'
      WRITE(9,*) ' #      E       pop       s_all    s_supr    s_real    J pi'
      DO ILVL=1,numlev
        poper=SQRT(POPVAR(ILVL)+POPERT(ILVL,0))
        WRITE(9,198) ILVL,elowlev(ILVL),POPULT(ILVL,0),poper,SQRT(POPVAR(ILVL)),SQRT(POPERT(ILVL,0)),elowsp(ILVL),ilowip(ILVL)
      ENDDO

      WRITE(9,*) 'Direct population of low-lying states from continuum'
      WRITE(9,*) ' #      E       pop       s_all    s_supr    s_real    J pi'
      DO ILVL=1,numlev
        poper=SQRT(POPSVAR(ILVL)+POPERS(ILVL,0))
        WRITE(9,198) ILVL,elowlev(ILVL),POPULS(ILVL,0),poper,SQRT(POPSVAR(ILVL)),SQRT(POPERS(ILVL,0)),elowsp(ILVL),ilowip(ILVL)
      ENDDO

  198 format(I3,F10.6,F11.7,' +- ',F9.7,' (',F9.7,' , ',F9.7,') ',F4.1,I2)

      ! !TODO vysledky s populacemi nizkolezicich hladin
      ! DO K=1,NUMLEV
      !   POPULT(K)=0.0
      !   POPERT(K)=0.0
      !   POPULS(K)=0.0
      !   POPERS(K)=0.0
      !   DO I=1,NUMLEV
      !     COVAP(K,I)=0.0
      !     COVAS(K,I)=0.0
      !   ENDDO
      ! ENDDO
      ! DO K=1,NUMLEV
      !   DO NUC=1,NREAL*NSUB
      !     POPULT(K)=POPULT(K)+POPTLEV(K,NUC)/FLOAT(NEVENTS)
      !     POPULS(K)=POPULS(K)+POPSLEV(K,NUC)/FLOAT(NEVENTS)
      !     POPERT(K)=POPERT(K)+(POPTLEV(K,NUC)/FLOAT(NEVENTS))**2
      !     POPERS(K)=POPERS(K)+(POPSLEV(K,NUC)/FLOAT(NEVENTS))**2
      !     DO I=1,K
      !       COVAP(K,I)=COVAP(K,I)+POPTLEV(K,NUC)/FLOAT(NEVENTS)*POPTLEV(I,NUC)/FLOAT(NEVENTS)
      !       COVAS(K,I)=COVAS(K,I)+POPSLEV(K,NUC)/FLOAT(NEVENTS)*POPSLEV(I,NUC)/FLOAT(NEVENTS)
      !     ENDDO
      !   ENDDO
      !   POPULT(K)=POPULT(K)/(NREAL*NSUB)
      !   POPULS(K)=POPULS(K)/(NREAL*NSUB)
      !   POPERT(K)=POPERT(K)/(NREAL*NSUB)-POPULT(K)**2
      !   POPERS(K)=POPERS(K)/(NREAL*NSUB)-POPULS(K)**2
      !   DO I=1,K
      !     COVAP(K,I)=COVAP(K,I)/(NREAL*NSUB)-POPULT(K)*POPULT(I)
      !     COVAS(K,I)=COVAS(K,I)/(NREAL*NSUB)-POPULS(K)*POPULS(I)
      !   ENDDO
      !   !write(*,*) elowlev(k), popult(k), (poptlev(k,i), i=1,NREAL)
      ! ENDDO
!         write(9,*)
!         WRITE(9,*) 'Direct population of low-lying states from continuum'
!         do i=1,numlev
!           if (popers(i).GT.0.0) then
!             poper=sqrt(popers(i)*rn1)
!           else
!             poper=0.0  
!           endif
!           WRITE(9,198) i,elowlev(i),populs(i),poper,elowsp(i),ilowip(i)
!         enddo
!         write(9,*)
! !
!         WRITE(9,*) ' Population - covariance matrix'
!         DO K=1,NUMLEV
!           DO L=1,K
!             IF ((COVAP(K,K)*COVAP(L,L)).GT.0.) THEN
!               CONO(K,L)=COVAP(K,L)/SQRT((COVAP(K,K)*COVAP(L,L)))
!             ELSE
!               CONO(K,L)=1.0
!             ENDIF
!           ENDDO
!         ENDDO
!         DO K=1,NUMLEV
! !          WRITE(9,105) (CONO(K,L),L=1,NUMLEV)
!           ! WRITE(9,*) (CONO(K,L),L=1,NUMLEV)
! !  105     FORMAT(<NUMLEV>F7.3)
!         ENDDO
!         WRITE(9,*)

!         WRITE(9,*) ' Sidefeeding - covariance matrix'
!         DO K=1,NUMLEV
!           DO L=1,K
!             IF ((COVAS(K,K)*COVAS(L,L)).GT.0.) THEN
!               CONO(K,L)=COVAS(K,L)/SQRT((COVAS(K,K)*COVAS(L,L)))
!             ELSE
!               CONO(K,L)=1.0
!             ENDIF
!           ENDDO
!         ENDDO
!         DO K=1,NUMLEV
! !          WRITE(9,105) (CONO(K,L),L=1,NUMLEV)
!           ! WRITE(9,*) (CONO(K,L),L=1,NUMLEV)
!         ENDDO
      CLOSE(9)
      RETURN
END SUBROUTINE WRITE_DICE_POPS
!***********************************************************************
SUBROUTINE INIT_POPS(POPULT,POPERT,POPVAR,COVAP,POPULS,POPERS,POPSVAR,COVAS)
!***********************************************************************
real,dimension(:,:),allocatable::    POPULT,POPERT,POPULS,POPERS
real,dimension(:),allocatable::      POPVAR,POPSVAR
real,dimension(:,:),allocatable::    COVAP,COVAS
INTEGER::                            IREAL,ILVL,JLVL

! total feeding variables
      allocate(POPULT(1:NUMLEV,0:NREAL))
      allocate(POPERT(1:NUMLEV,0:NREAL))
      allocate(POPVAR(1:NUMLEV))
      allocate(COVAP(1:NUMLEV,1:NUMLEV))

! side-feeding variables
      allocate(POPULS(1:NUMLEV,0:NREAL))
      allocate(POPERS(1:NUMLEV,0:NREAL))
      allocate(POPSVAR(1:NUMLEV))
      allocate(COVAS(1:NUMLEV,1:NUMLEV))

      DO ILVL=1,numlev
        DO IREAL=0,NREAL
          POPULT(ILVL,IREAL)=0.0
          POPERT(ILVL,IREAL)=0.0
          POPULS(ILVL,IREAL)=0.0
          POPERS(ILVL,IREAL)=0.0
        ENDDO
        POPVAR(ILVL)=0.0
        POPSVAR(ILVL)=0.0
        DO JLVL=1,numlev
          COVAP(ILVL,JLVL)=0.0
          COVAS(ILVL,JLVL)=0.0
        ENDDO
      ENDDO

      RETURN
END SUBROUTINE INIT_POPS
!***********************************************************************
SUBROUTINE CALC_POPS(POPTLEV,POPULT,POPERT,POPVAR,COVAP)
!***********************************************************************
real,dimension(:,:,:),allocatable::  POPTLEV
real,dimension(:,:),allocatable::    POPULT,POPERT
real,dimension(:),allocatable::      POPVAR
real,dimension(:,:),allocatable::    COVAP
INTEGER::                            IREAL,ISUB,ILVL
      
      DO ILVL=1,numlev
        DO IREAL=1,NREAL
          DO ISUB=1,NSUB
            POPULT(ILVL,IREAL)=POPULT(ILVL,IREAL)+poptlev(ILVL,IREAL,ISUB)/FLOAT(NEVENTS)
            POPERT(ILVL,IREAL)=POPERT(ILVL,IREAL)+(poptlev(ILVL,IREAL,ISUB)/FLOAT(NEVENTS))**2
          ENDDO
          POPULT(ILVL,IREAL)=POPULT(ILVL,IREAL)/NSUB
          POPERT(ILVL,IREAL)=POPERT(ILVL,IREAL)/NSUB-POPULT(ILVL,IREAL)**2
          POPULT(ILVL,0)=POPULT(ILVL,0)+POPULT(ILVL,IREAL)
          POPVAR(ILVL)=POPVAR(ILVL)+POPULT(ILVL,IREAL)**2
          POPERT(ILVL,0)=POPERT(ILVL,0)+POPERT(ILVL,IREAL)
        ENDDO
        POPULT(ILVL,0)=POPULT(ILVL,0)/NREAL !this holds the global average
        POPVAR(ILVL)=POPVAR(ILVL)/NREAL-POPULT(ILVL,0)**2 !this holds the variance of averages
        POPERT(ILVL,0)=POPERT(ILVL,0)/NREAL !this holds the average of variances within suprarealizations
      ENDDO
      RETURN
END SUBROUTINE CALC_POPS
!***********************************************************************
SUBROUTINE INICIALIZACE(NDEAD,NISOM,RADW,RADWID,RADWDI,POPTLEV,POPSLEV)
!TODO promenne dodelat vzhledem k jejich vypisu
integer::                            NDEAD,NISOM
real,dimension(:),allocatable::      RADWID,RADWDI
real,dimension(:,:),allocatable::    RADW
real,dimension(:,:,:),allocatable::  POPTLEV,POPSLEV
INTEGER::                            IREAL,ISUB,ILVL
!
      if (.not.allocated(RADWID)) then
        allocate(RADWID(0:NREAL))
      endif  
      if (.not.allocated(RADWDI)) then
        allocate(RADWDI(0:NREAL))
      endif  
      if (.not.allocated(RADW)) then
        allocate(RADW(1:NREAL,1:NSUB))
      endif  
      if (.not.allocated(POPTLEV)) then
        allocate(POPTLEV(1:numlev,1:NREAL,1:NSUB))
      endif 
      if (.not.allocated(POPSLEV)) then
        allocate(POPSLEV(1:numlev,1:NREAL,1:NSUB))
      endif 

      NDEAD=0
      NISOM=0
      DO IREAL=0,NREAL
        RADWID(IREAL)=0.0
        RADWDI(IREAL)=0.0
      ENDDO
      DO IREAL=1,NREAL
        DO ISUB=1,NSUB
          RADW(IREAL,ISUB)=0.0
          DO ILVL=1,numlev
            POPTLEV(ILVL,IREAL,ISUB) = 0.
            POPSLEV(ILVL,IREAL,ISUB) = 0.
          ENDDO
        ENDDO
      ENDDO
!
      RETURN
END SUBROUTINE INICIALIZACE
!***********************************************************************
SUBROUTINE WRITELEVELS(IGLOB,NBIN,LEVCON)
!   This subroutin writes the levels, see code
integer::                             IGLOB,NBIN
integer,dimension(:,:,:),allocatable:: LEVCON
CHARACTER(3)::                        EXT
REAL::                                S0,SPFI,ENRG,Q
INTEGER::                             IPFI,IS,NL,IL,IBIN
INTEGER,PARAMETER::                   MAXJC  = 49
!global shared SPINC(1), ndis, endis, BN, DELTA
      WRITE (EXT,211) IGLOB
      S0=SPINC(1)-FLOAT(INT(SPINC(1)))
      IF (S0.LE.1e-3) THEN  ! only for security reason
        S0 = 0.0
      ELSE
        S0 = 0.5
      ENDIF
      OPEN (UNIT=(IGLOB+10),FILE='LEVELS.'//EXT,STATUS='UNKNOWN')
      DO IS = 0, MAXJC
        DO IPFI = 0, 1
          SPFI = FLOAT(IS) + S0
          NL = ndis(ISUBSC(spfi),ipfi)
          DO IL = 1, NL
            ENRG = endis(IL,ISUBSC(spfi),ipfi)
            WRITE((IGLOB+10),212) ENRG, SPFI, IPFI
          ENDDO
          DO IBIN = NBIN, 1, -1
            NL = LEVCON(IBIN,IS,IPFI)-LEVCON(IBIN-1,IS,IPFI)
            DO IL = 1, NL
              Q=FLOAT(2*IL-1)/FLOAT(2*NL) !TODO swap from lowest to highest
              ENRG = BN - (FLOAT(IBIN)-Q)*DELTA
              WRITE((IGLOB+10),212) ENRG, SPFI, IPFI
            ENDDO
          ENDDO ! IBIN
          WRITE((IGLOB+10),*)
          WRITE((IGLOB+10),*)
        ENDDO !IPFI
      ENDDO !IS
      CLOSE((IGLOB+10))
  211 FORMAT (I3.3)
  212 FORMAT(F8.5,F5.1,I2)  ! End of the part writing generated levels
      RETURN
END SUBROUTINE WRITELEVELS
!***********************************************************************
end module vsechno
!***********************************************************************
