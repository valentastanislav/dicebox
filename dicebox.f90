!***********************************************************************
PROGRAM DICE_EVENT
use omp_lib
use lokalni_fce
use spolecne
use vsechno
!
implicit none
INTEGER,PARAMETER::    MAXJC  = 49
!PRIVATE
CHARACTER(80)::                       NAME
LOGICAL::        log_check,sidlev
INTEGER::        IC_type !0 = gamma, 1 = K-shell, 2 = higher-shell, 3 = pair
INTEGER::        IR0,IR1,IR2,IR3,IR4,ITID,IFLAG,NTOTAL,IREGI,IBIN,ILIN
INTEGER::        ilinc,ip,is,il,i,K,IEV,STEPS,IPFI,IBFI,ILFI,IPIN,ISUB
REAL::           U,dummy,SPFI,DMIX2,SIGN,SPIN, eg,EG_MAX,EG_STEP
real,dimension(:,:),allocatable::     sall
integer,dimension(:,:,:),allocatable::LEVCON
integer,dimension(:),allocatable::    IRCONc,IRCON
integer,dimension(:,:,:,:),allocatable::ISDIS

double precision,dimension(0:2,-2:2,0:1)::     GACON
double precision,dimension(:,:,:,:),allocatable::  STCON
double precision,dimension(0:2)::     TOTCON
integer,dimension(:,:,:,:),allocatable:: ISCON
real, dimension(0:2,0:20,-2:2,0:1)::  STDIS
real, dimension(0:2,-2:2,0:1)::       GADIS
real, dimension(0:2)::                TOTDIS

real,dimension(:,:),allocatable::    ELQQ,SPQQ,DMQQ,WIQQ
integer,dimension(:,:),allocatable:: IPQQ,ICQQ
integer,dimension(:),allocatable::   NR_STEPS

!SHARED
LOGICAL::        lopopgs
INTEGER::        NVLAKEN,NDEAD,NUC,NISOM
integer,dimension(:,:),allocatable:: KONTROLMATRIX,gamma_multiplicita
integer,dimension(:,:,:),allocatable:: intermediate2
integer,dimension(:,:,:,:),allocatable:: intermediate3,XSC_work

real,dimension(:,:),allocatable::    RADW
real,dimension(:,:,:),allocatable::  POPTLEV,POPSLEV

!FINAL SINGLE THREAD !TODO if NaN keep occuring, maybe try double precision?
real::                            RADVAR
real,dimension(:),allocatable::   RADWID,RADWDI
real,dimension(:,:),allocatable:: POPULT,POPERT,POPULS,POPERS
real,dimension(:),allocatable::   POPVAR,POPSVAR
real,dimension(:,:),allocatable:: COVAP,COVAS
real:: start, finish

!$OMP PARALLEL DEFAULT(PRIVATE) &
!$OMP SHARED(lopopgs,kpopgs,NVLAKEN,NDEAD,NISOM,KONTROLMATRIX,RADW,POPTLEV,POPSLEV,&
!$OMP ISWWR,ISWBN,ISWEL,ISWSP,ISWPA,ISWIC,ISWMX,ISWWI,ISWLS,&
!$OMP NOPTFL,NOPTE1,NOPTM1,NOPTE2,NOPTDE,LMODE,LDENP,LDSTAG,&
!$OMP NREAL, NEVENTS, NUMLEV, NSUB,&
!$OMP NGIGE,NLOWLOR,ER,W0,SIG,NGIGM,ERM,WM0,SIGM,NGIGE2,ERE,WE0,SIGE,&
!$OMP DEG,DMG,QEL,FERMC,TCONST,PAIR_PSF,EK0,EGZERO,&
!$OMP PAR_E1,PAR_M1,DIPSLP,DIPZER,&
!$OMP EZERO,DEL,TEMPER,ASHELL,AMASS,ZNUM,PAIRING,&
!$OMP ASHELL09,DEL09,TEMPER09,EZERO09,PAIRING09,&
!$OMP DENPPC,DENPA0,DENPA1,DENPA2,&
!$OMP DENLO,DENHI,DENPA,DENPB,DENPC,DENPD,&
!$OMP BN,SPINc,IPINC,NOPTCS,NLINC,CAPFR,&
!$OMP XRAYK,XRAYL,NENT,ELENT,CONVT,NENK,ELENK,CONVK,&
!$OMP ECRIT,EALL,max_decays,max_spin,&
!$OMP ndis,endis,dekod,denum,LVL_CLASS,LVL_ENERGY,delev,despin,deparity,deltx,&
!$OMP sal,errsal,alpha,p_conv,p_conv_K,p_conv_IPF,elowlev,elowsp,ilowip,isbspin,ityp,nddd,&
!$OMP TABENLD,TABLD,NLD,&
!$OMP TABENPSF,TABPSF,NPSF,&
!$OMP NBIN,DELTA,gamma_multiplicita,N_MSC_FS,MIN_MULTIPLICITA,MAX_MULTIPLICITA,MSC_FS,BIN_WIDTH,&
!$OMP RADWID,RADWDI,RADVAR,&
!$OMP intermediate2,intermediate3,XSC_work)
      ITID = OMP_GET_THREAD_NUM()
      IFLAG=0 !TODO ocesat od IFLAG
      U=0.
      call cpu_time(start)
!$OMP BARRIER
      IF (ITID.EQ.0) THEN
!******initial reading**************************************************
        CALL get_command_argument(1, NAME)
        IF (LEN_TRIM(NAME).EQ.0) THEN
          STOP 'Invalid input: need an input name as the first argument'
        ENDIF
        NAME=TRIM(NAME)
        INQUIRE(FILE=NAME, EXIST=log_check)
        IF (.not.log_check) THEN
          STOP 'Invalid input: file not found'
        ENDIF
        CALL READ_EV(NAME,lopopgs,KONTROLMATRIX)
        OPEN (UNIT=12,FILE='PSF_GS.DAT',STATUS='UNKNOWN')  !to GS
        OPEN (UNIT=13,FILE='PSF_INI.DAT',STATUS='UNKNOWN') !from initial level
        write(12,*) 'E_\gamma  PSF(E1)[MeV^-3]  PSF(M1)[MeV^-3]  PSF(E2)[MeV^-5]'
        write(13,*) 'E_\gamma  PSF(E1)[MeV^-3]  PSF(M1)[MeV^-3]  PSF(E2)[MeV^-5]'
        EG_MAX = 21.0
        EG_STEP = 0.05
        DO I = 1, INT(EG_MAX/EG_STEP)
          eg = FLOAT(I)*EG_STEP
          write(12,201) eg,sgamma(eg,eg,1)/eg**3,sgamma(eg,eg,3)/eg**3,sgamma(eg,eg,4)/eg**5
          write(13,201) eg,sgamma(eg,BN,1)/eg**3,sgamma(eg,BN,3)/eg**3,sgamma(eg,BN,4)/eg**5
    201   format(f7.3,4e13.5)
        ENDDO
        CLOSE (12)
        CLOSE (13)
!******adjustace*nbin***************************************************
        CALL ADJUST_NBIN(SPINC(1),NBIN)
        NVLAKEN=OMP_GET_NUM_THREADS()
        if (.not.allocated(gamma_multiplicita)) then
         allocate(gamma_multiplicita(1:NREAL*NSUB,0:MAX_MULTIPLICITA))
        endif
        DO NUC=1,NREAL*NSUB
          DO I=0,MAX_MULTIPLICITA
            gamma_multiplicita(NUC,I)=0
          ENDDO
        ENDDO
        if ((.not.allocated(intermediate2)).AND.(N_MSC_FS.GE.1)) then
         allocate(intermediate2(1:N_MSC_FS,1:NREAL*NSUB,0:INT(BN/BIN_WIDTH)))
        endif 
        if ((.not.allocated(intermediate3)).AND.(N_MSC_FS.GE.1)) then
         allocate(intermediate3(1:N_MSC_FS,1:NREAL*NSUB,0:INT(BN/BIN_WIDTH),0:INT(BN/BIN_WIDTH)))
        endif 
        if ((.not.allocated(XSC_work)).AND.(N_MSC_FS.GE.1)) then
         allocate(XSC_work(1:N_MSC_FS,1:NREAL*NSUB,MIN_MULTIPLICITA:MAX_MULTIPLICITA,0:INT(BN/BIN_WIDTH)))
        endif 
        CALL CNTRLMTRX(KONTROLMATRIX,4,NREAL)  !should be OK
        CALL INICIALIZACE(NDEAD,NISOM,RADW,RADWID,RADWDI,POPTLEV,POPSLEV)
!        IR0=KONTROLMATRIX(1,NREAL)
        IF (LMODE.EQ.1) THEN
          CALL GENERATE_GOE_EIGEN_VAL(IR1,700,IFLAG,U) !Maximum allowed dimension (2nd parameter) is 1000
        ENDIF
      ENDIF
!$OMP BARRIER
        if (.not.allocated(ELQQ)) then
          allocate(ELQQ(1:NEVENTS,0:126))
        endif
        if (.not.allocated(SPQQ)) then
          allocate(SPQQ(1:NEVENTS,0:126))
        endif
        if (.not.allocated(IPQQ)) then
          allocate(IPQQ(1:NEVENTS,0:126))
        endif
        if (.not.allocated(ICQQ)) then
          allocate(ICQQ(1:NEVENTS,0:126))
        endif
        if (.not.allocated(DMQQ)) then
          allocate(DMQQ(1:NEVENTS,0:126))
        endif
        if (.not.allocated(WIQQ)) then
          allocate(WIQQ(1:NEVENTS,0:126))
        endif
        if (.not.allocated(NR_STEPS)) then
          allocate(NR_STEPS(1:NEVENTS))
        endif
!$OMP DO SCHEDULE(DYNAMIC)
      DO NUC=1,NREAL
       IR1=KONTROLMATRIX(1,NUC) !level scheme
       IR2=KONTROLMATRIX(2,NUC) !radiation widths in form of precursors and low-lying intensity fluctuations
       IR3=KONTROLMATRIX(3,NUC) !MC of cascades - actual search for the final state
       IR4=KONTROLMATRIX(4,NUC) !MC of cascades - the seed for a) mixing of primaries b) "coin flip" of internal conversion
       write(*,*) 'starting NREAL with these seeds: ',NUC,IR1,IR2,IR3,IR4 !TODO delete after testing of reproducibility
!
       CALL LEVELSCH (IR1,NUC,IFLAG,NTOTAL,ITID,SPINC(1),U,LEVCON)
       IF (ISWLS.EQ.1) THEN   ! Writing generated levels if requested
         CALL WRITELEVELS(NUC,NBIN,LEVCON)
       ENDIF ! ISWLS
       CALL GERMS(IR2,NTOTAL,NDDD,IRCONc,IRCON)
!      Intensities of low-lying transitions can fluctuate
       CALL READ_INT(sall,IFLAG,U,IR2)
       write(*,*) 'processing NREAL with these seeds: ',NUC,IR1,IR2,IR3,IR4 !TODO delete after testing of reproducibility
       DO ISUB=1,NSUB !TODO alokace a pocitani pozorovatelnych
!
!      The following DO loop serves for computing of mixing ratios
!      delta for primary transitions (E2 admixture is probably not
!      very important, so this is done in a very simple way - and
!      probably not fully correctly in the case NLINc=2)
!
        if (.not.allocated(ISDIS)) then
         allocate(ISDIS(0:NLINc,1:20,0:ISUBSC(max_spin),0:1))
        endif 
        DO ilinc=1,nlinc
         DO ip=0,1
          DO is=0,8
           DO il=1,NDIS(is,ip)
            DO i=1,100
             dummy=ran0(IR4)
            ENDDO
            isdis(ilinc,il,is,ip)=IR4
           ENDDO
          ENDDO
         ENDDO
        ENDDO
!
!     Here, the procedure 'WIDTHS' is called in order to evaluate a total
!     radiative width and proper values of STDISp.
!
        if (.not.allocated(STCON)) then
         allocate(STCON(0:2,0:NBIN,-2:2,0:1))
        endif 
        if (.not.allocated(ISCON)) then
         allocate(ISCON(0:2,0:NBIN,-2:2,0:1))
        endif 
        IREGI=0
        IBIN=0
        ILIN=1
        IF (RADW(NUC,ISUB).NE.0.0) write(*,*) 'problem with routine inicializace' !TODO delete if working properly
!
        write(*,*) 'before WIDTHS_R for: ',ISUB,' of ',NUC
        DO ILINc=1,NLINc
          CALL WIDTHS_R(ILINc,IPINC,SPINC(ILINc),IBIN,ILIN,TOTCON,STCON,GACON,ISCON,TOTDIS,STDIS,GADIS,&
                        ISDIS,LEVCON,IRCON,IRCONc,IFLAG,U,IREGI,EIN,EFI)
        ENDDO
        DO ILINc=1,NLINc
          RADW(NUC,ISUB)=RADW(NUC,ISUB)+sngl(TOTCON(ILINc))+TOTDIS(ILINc)
        ENDDO
        IREGI=0
        IBIN=0
        ILIN=1
        write(*,*) 'cascading starting for: ',ISUB,' of ',NUC
!
!       The master DO-loop
        DO IEV=1,NEVENTS
          IF (MOD(log10(REAL(IEV-1)),1.).EQ.0) WRITE(*,5410) ITID,IEV
 5410     FORMAT('+',21X,I6,2X,I10)
          IC_type=0 !0 = gamma, 1 = K-shell, 2 = higher-shell, 3 = pair
          sidlev=.FALSE.
          EIN=BN
          STEPS=1
!
!         Determination of spin of capture state
          ILINc=1
          ELQQ(IEV,0)=BN
          SPQQ(IEV,0)=SPINc(ILINc)
          IPQQ(IEV,0)=IPINc
          WIQQ(IEV,0)=RADW(NUC,ISUB)
          IREGI=0
          IBIN=0
          ILIN=1
          CALL ONESTEP(ILINc,IPINC,SPINC(ILINc),IBIN,ILIN,TOTCON,STCON,GACON,ISCON,TOTDIS,STDIS,GADIS,&
                       ISDIS,IPFI,SPFI,IBFI,ILFI,DMIX2,sign,IR3,IR4,LEVCON,sall,U,IFLAG,EIN,EFI,IREGI,&
                       IC_type,IRCON,IRCONc)
          DO WHILE (EFI.GT.0.)
            ELQQ(IEV,steps)=EFI
            SPQQ(IEV,steps)=SPFI
            IPQQ(IEV,steps)=IPFI
            DMQQ(IEV,steps)=sign*sqrt(dmix2)
            IF (IREGI.GT.0) THEN
              do I=1,numlev
                if (efi.eq.elowlev(I)) then
                  poptlev(I,NUC,ISUB)=poptlev(I,NUC,ISUB)+1.
                  if (.NOT.sidlev) then
                    popslev(I,NUC,ISUB)=popslev(I,NUC,ISUB)+1.
                    sidlev=.TRUE.
                  endif
                endif
              enddo
            ENDIF
            ICQQ(IEV,steps)=IC_type
            IC_type=0
            IPIN=IPFI
            SPIN=SPFI
            IBIN=IBFI
            ILIN=ILFI
            IF (IREGI.LT.2) THEN
              CALL WIDTHS_R(0,IPIN,SPIN,IBIN,ILIN,TOTCON,STCON,GACON,ISCON,TOTDIS,STDIS,GADIS,ISDIS,LEVCON,&
                            IRCON,IRCONc,IFLAG,U,IREGI,EIN,EFI)
              IF ((SNGL(TOTCON(0))+TOTDIS(0)).LE.0.) THEN
                NDEAD=NDEAD+1
                GO TO 5
              ENDIF
            ENDIF
            if (iregi.eq.2) then 
              if (denum(dekod(ilfi,isubsc(spfi),ipfi)).eq.0) then
                NISOM=NISOM+1
                !WRITE(*,*) 'so you declared isomeric state that does not decay'
                GO TO 6
              endif
            endif
            WIQQ(IEV,steps)=SNGL(TOTCON(0))+TOTDIS(0)
            STEPS=STEPS+1
            CALL ONESTEP(0,IPIN,SPIN,IBIN,ILIN,TOTCON,STCON,GACON,ISCON,TOTDIS,STDIS,GADIS,ISDIS,IPFI,SPFI,&
                IBFI,ILFI,dmix2,sign,IR3,IR4,LEVCON,sall,U,IFLAG,EIN,EFI,IREGI,IC_type,IRCON,IRCONc)
          ENDDO  !WHILE EFI
   6      ELQQ(IEV,steps)=EFI
          SPQQ(IEV,steps)=SPFI
          IPQQ(IEV,steps)=IPFI
          WIQQ(IEV,steps)=0.0
          DMQQ(IEV,steps)=sign*sqrt(dmix2)
!   **** feeding of ground state ****
          if (lopopgs) then
            poptlev(kpopgs,NUC,ISUB)=poptlev(kpopgs,NUC,ISUB)+1.
            if (.NOT.sidlev) then
              popslev(kpopgs,NUC,ISUB)=popslev(kpopgs,NUC,ISUB)+1.
            endif
          endif
          ICQQ(IEV,steps)=IC_type
          IC_type=0
          NR_STEPS(IEV)=STEPS
   5      CONTINUE
        ENDDO
!       writin the cascades
!        write(*,*) 'cascading ended'
        CALL SPECTRA((ISUB+(NUC-1)*NSUB),ELQQ,SPQQ,DMQQ,IPQQ,ICQQ,NR_STEPS,intermediate2,intermediate3,XSC_work,gamma_multiplicita)
        IF (ISWWR.EQ.1) CALL DO_IT(NUC,ISUB,ELQQ,SPQQ,DMQQ,IPQQ,ICQQ,WIQQ,NR_STEPS,ITID)
        IRCONc(1)=IR3 !TODO staci takhle z hlediska reprodukovatelnosti?
       ENDDO !DO ISUB=1,NSUB
      ENDDO !DO NUC=1,NREAL
!$OMP END DO
!$OMP END PARALLEL
      IF(N_MSC_FS.GT.0) CALL WR_SPECTRA (intermediate2,intermediate3,XSC_work,gamma_multiplicita)

      RADVAR=0.0
      DO NUC=1,NREAL
        DO ISUB=1,NSUB
          RADWID(NUC)=RADWID(NUC)+RADW(NUC,ISUB)
          RADWDI(NUC)=RADWDI(NUC)+RADW(NUC,ISUB)**2
        ENDDO
        RADWID(NUC)=RADWID(NUC)/NSUB !this holds the average within given suprarealization
        RADWDI(NUC)=RADWDI(NUC)/NSUB-RADWID(NUC)**2 !this holds the variance within given suprarealization
        RADWID(0)=RADWID(0)+RADWID(NUC)
        RADVAR=RADVAR+RADWID(NUC)**2
        RADWDI(0)=RADWDI(0)+RADWDI(NUC)
      ENDDO
      RADWID(0)=RADWID(0)/NREAL !this holds the global average
      RADVAR=RADVAR/NREAL-RADWID(0)**2 !this holds the variance of averages
      RADWDI(0)=RADWDI(0)/NREAL !this holds the average of variances within suprarealizations
      CALL INIT_POPS(POPULT,POPERT,POPVAR,COVAP,POPULS,POPERS,POPSVAR,COVAS)
      CALL CALC_POPS(POPTLEV,POPULT,POPERT,POPVAR,COVAP)
      CALL CALC_POPS(POPSLEV,POPULS,POPERS,POPSVAR,COVAS)
      CALL WRITE_PARAMS
      CALL WRITE_DICE_PRO(RADWID(0),RADVAR,RADWDI(0),NDEAD,NISOM)
      CALL WRITE_DICE_POPS(POPULT,POPERT,POPVAR,COVAP,POPULS,POPERS,POPSVAR,COVAS)
      call cpu_time(finish)
      print '("Run Time = ",f6.3," seconds.")',finish-start
      END PROGRAM DICE_EVENT
      