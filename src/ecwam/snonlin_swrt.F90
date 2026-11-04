! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      MODULE SNONLIN_SWRT_MOD

! ----------------------------------------------------------------------

!**** *SNONLIN_SWRT_MOD* - NONLINEAR TRANSFER BY SPARSE WRT (SWRT).

!     PORT OF THE WAVEWATCH III NL6 MODULE W3SNL6MD (SWRT, DEFAULT
!     LOCUS-GROUPED KERNEL) TO ECWAM. SELECTED WITH INLSCHEME = 1.

!*    PURPOSE.
!     --------

!       SWRT: THE WRT PAIR SUM RESTRICTED TO K MIRROR-PAIR UNITS OF ITS
!       QUADRATURE TABLE, ONE GAIN PER UNIT, DEEP WATER. IT REPLACES
!       THE DEEP-WATER DIA TRANSFER OF SNONLIN; THE SHALLOW WATER
!       SCALING ENH OF SNONLIN (ISNONLIN) IS APPLIED ON TOP.

!**   INTERFACE.
!     ----------

!       *CALL* *INISNONLIN_SWRT*
!          READS THE OPERATOR FILE NAMED BY THE ENVIRONMENT VARIABLE
!          WAM_SWRT_OPERATOR (WW3-SWRT-OPERATOR FORMAT VERSION 1, SEE
!          export_operator.py OF THE SWRT PACKAGE) AND CHECKS IT AGAINST
!          THE MODEL GRID: NANG = na, NFRE = nk, FR = sigma/(2 pi) TO A
!          RELATIVE 1.E-5. WAM_SWRT_DIAG=0 SWITCHES THE DIAGONAL OFF.

!       *CALL* *SNONLIN_SWRT (KIJS, KIJL, FL1, FLD, SL, WAVNUM, DEPTH, AKMEAN)*
!          SAME ARGUMENTS AS SNONLIN.
!          *FL1*  - SPECTRUM F(f,theta) (m**2/(Hz rad)).
!          *FLD*  - DIAGONAL OF dS/dF (1/s), INCREMENTED.
!          *SL*   - TOTAL SOURCE dF/dt, INCREMENTED.

!     METHOD.
!     -------

!       THE KERNEL IS THAT OF W3SNL6 (LOCUSLOOP), IN SINGLE PRECISION AS
!       IN WAVEWATCH III. E(f,theta) IS FL1 ITSELF (W3SNL6 WORKS ON
!       E = A/CONX, THE SAME QUANTITY), SO S IS ADDED TO SL AS IT IS AND
!       THE DIAGONAL dS/dE = dS/dN * CONV IS ADDED TO FLD.

!       SNONLIN_SWRT (DEFAULT) RUNS ALL GRID POINTS KIJS:KIJL TOGETHER:
!       THE INNERMOST LOOPS RUN OVER THE POINT INDEX (CONTIGUOUS, AS IN
!       SNONLIN), WITH ANCHOR FREQUENCY, CURVE, ROW AND DIRECTION OUTSIDE.
!       IT DOES THE SAME OPERATIONS IN THE SAME ORDER ON EACH POINT AS
!       SNONLIN_SWRT_REF (ONE SPECTRUM AT A TIME, INNER LOOPS OVER
!       DIRECTIONS). WORK THAT THE DENSITY FILTER ZEROES IS SKIPPED PER
!       CURVE AND PER DIRECTION WHEN IT IS ZERO AT ALL POINTS. THE ONLY
!       DIFFERENCE: SNONLIN_SWRT FLUSHES SUBNORMAL SINGLE PRECISION
!       NUMBERS (< 1.2E-38) TO ZERO, FOR SPEED; WITH THE SAME UNDERFLOW
!       MODE THE TWO ARE BIT-IDENTICAL. WAM_SWRT_REF=1 SELECTS
!       SNONLIN_SWRT_REF INSTEAD (FOR CHECKS).

!       DIRECTIONS: THE OPERATOR'S DIRECTION OFFSETS ARE FOR AN INDEX
!       INCREASING ANTI-CLOCKWISE (WAVEWATCH III). ECWAM'S TH(K) INCREASES
!       CLOCKWISE (FROM NORTH, GOING TO), SO ECWAM BIN K IS OPERATOR BIN
!       J = NANG+1-K. THE ABSOLUTE ORIGIN (INCLUDING THE HALF-BIN OFFSET
!       OF TH) DOES NOT ENTER: THE KERNEL ONLY USES INDEX DIFFERENCES ON
!       A UNIFORM PERIODIC GRID.

!       DEPTH: CONTRIBUTIONS OF THE PAIRS ANCHORED AT FREQUENCY M ARE
!       SCALED BY ENH(IJ,M) OF SNONLIN (ISNONLIN = 0, 1, 2), AS THE DIA
!       SCALES EACH QUADRUPLET BY ENH AT ITS CENTRAL FREQUENCY.

!     REFERENCE.
!     ----------

!       WAVEWATCH III w3snl6md.F90 (NL6 = SWRT), keydoc/swrt/ww3/nl6.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB

      IMPLICIT NONE
      PRIVATE
      PUBLIC :: INISNONLIN_SWRT, SNONLIN_SWRT, SNONLIN_SWRT_REF

!     SINGLE PRECISION, AS THE WAVEWATCH III KERNEL
      INTEGER, PARAMETER :: JWRS = SELECTED_REAL_KIND(6,37)

!     OPERATOR (READ ONCE BY INISNONLIN_SWRT, READ-ONLY AFTERWARDS)
      LOGICAL :: LSWRTREAD = .FALSE.
      LOGICAL :: LSWRTDIAG = .TRUE.
      LOGICAL :: LSWRTREF = .FALSE.
      INTEGER(KIND=JWIM) :: NR, NG
      INTEGER(KIND=JWIM), ALLOCATABLE :: DI(:), DJ(:), J2(:), J4(:), GS(:)
      REAL(KIND=JWRS), ALLOCATABLE :: G2(:), G4(:), W(:,:), FA(:,:), FP(:,:)
      REAL(KIND=JWRS), ALLOCATABLE :: CONV(:), KPOW(:)
      REAL(KIND=JWRS) :: QFRAC
!     FREQUENCY CORNERS OF K2 AND K4 PER (ANCHOR, ROW), AS CORNER OF W3SNL6
      INTEGER(KIND=JWIM), ALLOCATABLE :: KA2(:,:), KB2(:,:), KA4(:,:), KB4(:,:)
      REAL(KIND=JWRS), ALLOCATABLE :: U2(:,:), U4(:,:)

      CONTAINS

! ----------------------------------------------------------------------

      SUBROUTINE INISNONLIN_SWRT

      USE YOWFRED  , ONLY : FR
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : ZPI
      USE YOWTEST  , ONLY : IU06
      USE YOWABORT , ONLY : WAM_ABORT
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

      IMPLICIT NONE

      CHARACTER(LEN=1024) :: CLFILE
      CHARACTER(LEN=256)  :: CLLINE, CLSRC, CLMSG
      CHARACTER(LEN=16)   :: CLDIAG
      INTEGER(KIND=JWIM) :: IU, IERR, ILEN, IVER, INK, INA, IKU, M, R, MW
      INTEGER(KIND=JWIM), ALLOCATABLE :: I2(:), I4(:)
      REAL(KIND=JWRS), ALLOCATABLE :: F2(:), F4(:)
      REAL(KIND=JWRB) :: RDIF, RMAX
      REAL(KIND=JWRB), ALLOCATABLE :: SIGOP(:)
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

      IF (LHOOK) CALL DR_HOOK('INISNONLIN_SWRT',0,ZHOOK_HANDLE)

      IF (LSWRTREAD) THEN
        IF (LHOOK) CALL DR_HOOK('INISNONLIN_SWRT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF

!*    1. FILE NAME AND HEADER.
!        ---------------------

      CLFILE = ' '
      CALL GET_ENVIRONMENT_VARIABLE('WAM_SWRT_OPERATOR', CLFILE, ILEN, IERR)
      IF (IERR /= 0 .OR. ILEN == 0) CALL SWRTERR('environment variable WAM_SWRT_OPERATOR is not set '// &
     &   '(it must name the SWRT operator file, .ww3op)')

      OPEN(NEWUNIT=IU, FILE=TRIM(CLFILE), STATUS='OLD', ACTION='READ', FORM='FORMATTED', IOSTAT=IERR)
      IF (IERR /= 0) CALL SWRTERR('cannot open the file')
      READ(IU,'(A)',IOSTAT=IERR) CLLINE
      IF (IERR /= 0 .OR. TRIM(ADJUSTL(CLLINE)) /= 'WW3-SWRT-OPERATOR') CALL SWRTERR('first line is not WW3-SWRT-OPERATOR')
      READ(IU,*,IOSTAT=IERR) IVER
      IF (IERR /= 0 .OR. IVER /= 1) CALL SWRTERR('unknown format version (this code reads 1)')
      READ(IU,*,IOSTAT=IERR) INK, INA, NR, IKU
      IF (IERR /= 0 .OR. MIN(INK,INA,NR) < 1) CALL SWRTERR('cannot read the line nk na nr K')
      READ(IU,'(A)',IOSTAT=IERR) CLSRC
      IF (IERR /= 0) CALL SWRTERR('cannot read the source name')

!*    2. SPECTRAL GRID: NANG, NFRE AND FR.
!        ---------------------------------

      IF (INK /= NFRE .OR. INA /= NANG) THEN
        WRITE(CLMSG,'(A,I0,A,I0,A,I0,A,I0)') 'spectral grid mismatch: model NFRE = ', NFRE, ', NANG = ', NANG, &
     &   '; operator nk = ', INK, ', na = ', INA
        CALL SWRTERR(CLMSG)
      ENDIF
      ALLOCATE(SIGOP(NFRE))
      CALL RDLAB('sigma')
      READ(IU,*,IOSTAT=IERR) SIGOP
      IF (IERR /= 0) CALL SWRTERR('cannot read block sigma')
      RMAX = 0.0_JWRB
      MW = 1
      DO M=1,NFRE
        RDIF = ABS(ZPI*FR(M)-SIGOP(M))/MAX(ABS(SIGOP(M)),TINY(1.0_JWRB))
        IF (RDIF > RMAX) THEN
          RMAX = RDIF
          MW = M
        ENDIF
      ENDDO
      IF (RMAX > 1.E-5_JWRB) THEN
        WRITE(CLMSG,'(A,I0,A,ES14.7,A,ES14.7,A,ES9.2,A)') 'frequency mismatch: FR(', MW, ') = ', FR(MW), &
     &   ' Hz, operator sigma/(2 pi) = ', SIGOP(MW)/ZPI, ' Hz (relative ', RMAX, ' > 1.E-5)'
        CALL SWRTERR(CLMSG)
      ENDIF

!*    3. ROW ARRAYS.
!        -----------

      ALLOCATE(DI(NR), DJ(NR), I2(NR), F2(NR), J2(NR), G2(NR), I4(NR), F4(NR), J4(NR), G4(NR))
      ALLOCATE(W(NFRE,NR), FA(NFRE,NR), FP(NFRE,NR), CONV(NFRE), KPOW(NFRE))
      CALL RDLAB('di')
      READ(IU,*,IOSTAT=IERR) DI
      CALL RDLAB('dj')
      READ(IU,*,IOSTAT=IERR) DJ
      CALL RDLAB('i2')
      READ(IU,*,IOSTAT=IERR) I2
      CALL RDLAB('f2')
      READ(IU,*,IOSTAT=IERR) F2
      CALL RDLAB('j2')
      READ(IU,*,IOSTAT=IERR) J2
      CALL RDLAB('g2')
      READ(IU,*,IOSTAT=IERR) G2
      CALL RDLAB('i4')
      READ(IU,*,IOSTAT=IERR) I4
      CALL RDLAB('f4')
      READ(IU,*,IOSTAT=IERR) F4
      CALL RDLAB('j4')
      READ(IU,*,IOSTAT=IERR) J4
      CALL RDLAB('g4')
      READ(IU,*,IOSTAT=IERR) G4
      CALL RDLAB('w')
      READ(IU,*,IOSTAT=IERR) W
      CALL RDLAB('fa')
      READ(IU,*,IOSTAT=IERR) FA
      CALL RDLAB('fp')
      READ(IU,*,IOSTAT=IERR) FP
      CALL RDLAB('conv')
      READ(IU,*,IOSTAT=IERR) CONV
      CALL RDLAB('qfrac')
      READ(IU,*,IOSTAT=IERR) QFRAC
      CALL RDLAB('kpow')
      READ(IU,*,IOSTAT=IERR) KPOW
      IF (IERR /= 0) CALL SWRTERR('cannot read block kpow')
      CLOSE(IU)

      IF (MINVAL(DI) < 0) CALL SWRTERR('negative offset di')
      IF (MINVAL(J2) < -NANG .OR. MAXVAL(J2) > NANG-1 .OR. MINVAL(J4) < -NANG .OR. MAXVAL(J4) > NANG-1 .OR. &
     &    MAXVAL(ABS(DJ)) > NANG) CALL SWRTERR('direction offsets too large for the kernel')

!     FREQUENCY CORNERS (CORNER OF W3SNL6): LOWER CORNER CLAMPED TO THE
!     GRID WITH FRACTION 0 WHEN IT IS OFF THE GRID, UPPER = LOWER+1.
      ALLOCATE(KA2(NFRE,NR), KB2(NFRE,NR), U2(NFRE,NR), KA4(NFRE,NR), KB4(NFRE,NR), U4(NFRE,NR))
      DO R=1,NR
        DO M=1,NFRE
          CALL CORNER(M+I2(R), F2(R), KA2(M,R), KB2(M,R), U2(M,R))
          CALL CORNER(M+I4(R), F4(R), KA4(M,R), KB4(M,R), U4(M,R))
        ENDDO
      ENDDO

!     ROWS GROUPED BY RESONANCE CURVE (SAME DI, DJ); THE FILE LISTS THEM TOGETHER
      ALLOCATE(GS(NR+1))
      NG = 1
      GS(1) = 1
      DO R=2,NR
        IF (DI(R) /= DI(R-1) .OR. DJ(R) /= DJ(R-1)) THEN
          NG = NG+1
          GS(NG) = R
        ENDIF
      ENDDO
      GS(NG+1) = NR+1

!*    4. DIAGONAL TERM ON (DEFAULT) OR OFF (WAM_SWRT_DIAG = 0).
!        ------------------------------------------------------

      CALL GET_ENVIRONMENT_VARIABLE('WAM_SWRT_DIAG', CLDIAG, ILEN, IERR)
      LSWRTDIAG = .NOT. (IERR == 0 .AND. TRIM(ADJUSTL(CLDIAG)) == '0')
!     WAM_SWRT_REF = 1: PER-SPECTRUM REFERENCE KERNEL (CHECKS ONLY)
      CALL GET_ENVIRONMENT_VARIABLE('WAM_SWRT_REF', CLDIAG, ILEN, IERR)
      LSWRTREF = (IERR == 0 .AND. TRIM(ADJUSTL(CLDIAG)) == '1')

      WRITE(IU06,*) ' '
      WRITE(IU06,*) ' SNONLIN_SWRT: SWRT OPERATOR ', TRIM(CLFILE)
      WRITE(IU06,'(3A,I0,A,I0,A,I0,A,I0,A,I0,3A)') '   (', TRIM(CLSRC), ', K = ', IKU, ' CONFIGURATIONS, ', NR, ' ROWS, ', &
     &   NG, ' OFFSETS, NFRE = ', NFRE, ', NANG = ', NANG, ', DIAGONAL TERM ', TRIM(MERGE('ON ','OFF',LSWRTDIAG)), ')'
      IF (LSWRTREF) WRITE(IU06,*) '  WAM_SWRT_REF = 1: PER-SPECTRUM REFERENCE KERNEL SNONLIN_SWRT_REF'
      WRITE(IU06,*) ' '

      DEALLOCATE(SIGOP, I2, F2, I4, F4)
      LSWRTREAD = .TRUE.

      IF (LHOOK) CALL DR_HOOK('INISNONLIN_SWRT',1,ZHOOK_HANDLE)

      CONTAINS

        SUBROUTINE SWRTERR(CDTXT)
        CHARACTER(LEN=*), INTENT(IN) :: CDTXT
        WRITE(IU06,*) ' *** ERROR IN INISNONLIN_SWRT (SWRT OPERATOR) ***'
        WRITE(IU06,*) '     OPERATOR FILE : ', TRIM(CLFILE)
        WRITE(IU06,*) '     ', TRIM(CDTXT)
        CALL WAM_ABORT('INISNONLIN_SWRT: operator file '//TRIM(CLFILE)//': '//TRIM(CDTXT))
        END SUBROUTINE SWRTERR

        SUBROUTINE RDLAB(CDLAB)
!       NAME LINE OF THE NEXT BLOCK; ALSO CATCHES A FAILED READ OF THE PREVIOUS BLOCK
        CHARACTER(LEN=*), INTENT(IN) :: CDLAB
        IF (IERR /= 0) CALL SWRTERR('cannot read the block before '//CDLAB)
        READ(IU,'(A)',IOSTAT=IERR) CLLINE
        IF (IERR /= 0 .OR. TRIM(ADJUSTL(CLLINE)) /= CDLAB) CALL SWRTERR('expected block '//CDLAB//', found: '//TRIM(CLLINE))
        END SUBROUTINE RDLAB

        SUBROUTINE CORNER(I0, FIN, IA, IB, FU)
        INTEGER(KIND=JWIM), INTENT(IN)  :: I0
        REAL(KIND=JWRS),    INTENT(IN)  :: FIN
        INTEGER(KIND=JWIM), INTENT(OUT) :: IA, IB
        REAL(KIND=JWRS),    INTENT(OUT) :: FU
        IA = I0
        FU = FIN
        IF (IA < 1 .OR. IA > NFRE) THEN
          FU = 0.0_JWRS
          IA = MIN(MAX(IA,1),NFRE)
        ENDIF
        IB = MIN(IA+1,NFRE)
        END SUBROUTINE CORNER

      END SUBROUTINE INISNONLIN_SWRT

! ----------------------------------------------------------------------

      SUBROUTINE SNONLIN_SWRT (KIJS, KIJL, FL1, FLD, SL, WAVNUM, DEPTH, AKMEAN)

!     SWRT OVER ALL POINTS KIJS:KIJL AT ONCE; THE INNERMOST LOOPS RUN OVER
!     THE POINT INDEX. PER POINT THE ARITHMETIC IS THAT OF SNONLIN_SWRT_REF
!     (SAME OPERATIONS, SAME ORDER OF SUMMATION).

      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWABORT , ONLY : WAM_ABORT
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY : IEEE_SUPPORT_UNDERFLOW_CONTROL, IEEE_GET_UNDERFLOW_MODE, &
     &                                          IEEE_SET_UNDERFLOW_MODE

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(INOUT):: FLD, SL
      REAL(KIND=JWRB), DIMENSION(KIJL,NFRE), INTENT(IN) :: WAVNUM
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: DEPTH, AKMEAN

      INTEGER(KIND=JWIM) :: IJ, K, M, J, G, RR, R0, K1, K3, JJ, IA, NJA, O3, JD3, J2A, J4A
      INTEGER(KIND=JWIM) :: KA2_, KB2_, KA4_, KB4_
      INTEGER(KIND=JWIM), DIMENSION(NANG) :: JACT
      LOGICAL :: LACT, LUFCTL, LGRADUAL

      REAL(KIND=JWRB), DIMENSION(KIJL,NFRE) :: ENH

!     SINGLE PRECISION WORK ARRAYS, POINT INDEX FIRST (LOCAL: THREAD-SAFE).
!     NX, SX, DX HAVE A PADDED PERIODIC DIRECTION AXIS 1-NANG:2*NANG IN THE
!     OPERATOR'S ORDER, AS IN SNONLIN_SWRT_REF; THEY ARE ALLOCATED (HEAP)
!     BECAUSE THEY ARE TOO LARGE FOR THE STACK OF AN OPENMP THREAD.
      REAL(KIND=JWRS), ALLOCATABLE :: NX(:,:,:), SX(:,:,:), DX(:,:,:)
      REAL(KIND=JWRS), DIMENSION(KIJL,NANG) :: MM, AA, BB, AS1, AD1
      REAL(KIND=JWRS), DIMENSION(KIJL,NFRE) :: RMX
      REAL(KIND=JWRS), DIMENSION(KIJL) :: NMAX, QM, E1, WA, WP
      REAL(KIND=JWRS) :: XN, X1, X3, X2, X4, A2, B2, A4, B4, WR, FA_, FP_, U2_, U4_, V2_, V4_, TS, MMX

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SNONLIN_SWRT',0,ZHOOK_HANDLE)

      IF (.NOT. LSWRTREAD) CALL WAM_ABORT('SNONLIN_SWRT: operator not read (INISNONLIN_SWRT not called)')

      IF (LSWRTREF) THEN
        CALL SNONLIN_SWRT_REF (KIJS, KIJL, FL1, FLD, SL, WAVNUM, DEPTH, AKMEAN)
        IF (LHOOK) CALL DR_HOOK('SNONLIN_SWRT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF

!*    1. SHALLOW WATER SCALING, AS IN SNONLIN (FREQUENCIES 1 TO NFRE).
!        -------------------------------------------------------------

      CALL SWRT_ENH(KIJS, KIJL, FL1, WAVNUM, DEPTH, AKMEAN, ENH)

!     SUBNORMAL NUMBERS FLUSHED TO ZERO IN SECTIONS 2 AND 3 (THE MODE OF
!     THE CALLER IS RESTORED BEFORE SECTION 4). PRODUCTS OF SMALL SPECTRAL
!     DENSITIES IN SINGLE PRECISION ARE OFTEN SUBNORMAL (< 1.2E-38), AND
!     GRADUAL UNDERFLOW MAKES THIS KERNEL ABOUT 5 TIMES SLOWER. ONLY TERMS
!     BELOW 1.2E-38 ARE AFFECTED.
      LUFCTL = IEEE_SUPPORT_UNDERFLOW_CONTROL(1.0_JWRS)
      IF (LUFCTL) THEN
        CALL IEEE_GET_UNDERFLOW_MODE(LGRADUAL)
        CALL IEEE_SET_UNDERFLOW_MODE(.FALSE.)
      ENDIF

!*    2. N = E * CONV, OPERATOR DIRECTION ORDER (BIN J = ECWAM BIN NANG+1-J),
!        PADDED PERIODIC; ROW AND SPECTRUM MAXIMA FOR THE DENSITY FILTER.
!        ------------------------------------------------------------------

      ALLOCATE(NX(KIJL,1-NANG:2*NANG,NFRE), SX(KIJL,1-NANG:2*NANG,NFRE), DX(KIJL,1-NANG:2*NANG,NFRE))

      NMAX(KIJS:KIJL) = -HUGE(1.0_JWRS)
      DO M=1,NFRE
        RMX(KIJS:KIJL,M) = -HUGE(1.0_JWRS)
        DO J=1,NANG
!$OMP SIMD PRIVATE(XN)
          DO IJ=KIJS,KIJL
            XN = REAL(FL1(IJ,NANG+1-J,M),JWRS)*CONV(M)
            NX(IJ,J,M) = XN
            NX(IJ,J-NANG,M) = XN
            NX(IJ,J+NANG,M) = XN
            RMX(IJ,M) = MAX(RMX(IJ,M), XN)
          ENDDO
        ENDDO
        DO IJ=KIJS,KIJL
          NMAX(IJ) = MAX(NMAX(IJ), RMX(IJ,M))
        ENDDO
      ENDDO
      SX(KIJS:KIJL,:,:) = 0.0_JWRS
      DX(KIJS:KIJL,:,:) = 0.0_JWRS

!*    3. SWRT (LOCUSLOOP OF W3SNL6), ALL POINTS TOGETHER.
!        ------------------------------------------------

      DO K1=1,NFRE
        DO IJ=KIJS,KIJL
          E1(IJ) = REAL(ENH(IJ,K1),JWRS)
        ENDDO
        AS1(KIJS:KIJL,:) = 0.0_JWRS
        AD1(KIJS:KIJL,:) = 0.0_JWRS

        DO G=1,NG
          R0 = GS(G)
          O3 = DJ(R0)
          K3 = K1 + DI(R0)
          IF (K3 > NFRE) CYCLE
!         NO ROW OF THIS CURVE HAS A WEIGHT AT K1
          IF (ALL(W(K1,GS(G):GS(G+1)-1) == 0.0_JWRS)) CYCLE

!         DENSITY FILTER. THE REFERENCE KERNEL SKIPS THE CURVE AT A POINT
!         WHEN BOTH FREQUENCY ROWS ARE BELOW QM THERE; THE MASK MM THEN IS
!         0 IN ALL DIRECTIONS, SO SKIPPING IT ONLY WHEN THAT HOLDS AT ALL
!         POINTS, AND SKIPPING DIRECTIONS WHERE MM = 0 AT ALL POINTS, DOES
!         NOT CHANGE THE RESULT.
          LACT = .FALSE.
          DO IJ=KIJS,KIJL
            QM(IJ) = QFRAC * NMAX(IJ) * KPOW(K3)
            LACT = LACT .OR. RMX(IJ,K1) >= QM(IJ) .OR. RMX(IJ,K3) >= QM(IJ)
          ENDDO
          IF (.NOT. LACT) CYCLE

          NJA = 0
          DO JJ=1,NANG
            JD3 = JJ+O3
            MMX = 0.0_JWRS
!$OMP SIMD REDUCTION(MAX:MMX)
            DO IJ=KIJS,KIJL
!             SAME TEST AS X1 >= QM .OR. X3 >= QM, WITHOUT A BRANCH
              MM(IJ,JJ) = MERGE(1.0_JWRS, 0.0_JWRS, MAX(NX(IJ,JJ,K1), NX(IJ,JD3,K3)) >= QM(IJ))
              MMX = MAX(MMX, MM(IJ,JJ))
            ENDDO
            IF (MMX > 0.0_JWRS) THEN
              NJA = NJA + 1
              JACT(NJA) = JJ
              DO IJ=KIJS,KIJL
                AA(IJ,JJ) = 0.0_JWRS
                BB(IJ,JJ) = 0.0_JWRS
              ENDDO
            ENDIF
          ENDDO
          IF (NJA == 0) CYCLE

!         ROWS OF THE CURVE: INTERPOLATED PARTNERS K2, K4 AND THE SUMS
!         AA = SUM W (N4 - N2), BB = SUM W N2 N4, IN ROW ORDER
          DO RR=GS(G),GS(G+1)-1
            WR = W(K1,RR)
            IF (WR == 0.0_JWRS) CYCLE
            B2 = G2(RR)
            A2 = 1.0_JWRS - B2
            B4 = G4(RR)
            A4 = 1.0_JWRS - B4
            KA2_ = KA2(K1,RR)
            KB2_ = KB2(K1,RR)
            U2_  = U2(K1,RR)
            V2_  = 1.0_JWRS-U2_
            KA4_ = KA4(K1,RR)
            KB4_ = KB4(K1,RR)
            U4_  = U4(K1,RR)
            V4_  = 1.0_JWRS-U4_
            DO IA=1,NJA
              JJ  = JACT(IA)
              J2A = JJ+J2(RR)
              J4A = JJ+J4(RR)
!$OMP SIMD PRIVATE(X2,X4)
              DO IJ=KIJS,KIJL
                X2 = V2_ * ( A2*NX(IJ,J2A,KA2_) + B2*NX(IJ,J2A+1,KA2_) )   &
     &             + U2_ * ( A2*NX(IJ,J2A,KB2_) + B2*NX(IJ,J2A+1,KB2_) )
                X4 = V4_ * ( A4*NX(IJ,J4A,KA4_) + B4*NX(IJ,J4A+1,KA4_) )   &
     &             + U4_ * ( A4*NX(IJ,J4A,KB4_) + B4*NX(IJ,J4A+1,KB4_) )
                AA(IJ,JJ) = AA(IJ,JJ) + WR * ( X4 - X2 )
                BB(IJ,JJ) = BB(IJ,JJ) + WR * X2 * X4
              ENDDO
            ENDDO
          ENDDO

!         DEPOSITS AT K1 (SUMMED OVER CURVES IN AS1, AD1) AND AT K3
          FA_ = FA(K1,R0)
          FP_ = FP(K1,R0)
          DO IJ=KIJS,KIJL
            WA(IJ) = FA_ * E1(IJ)
            WP(IJ) = FP_ * E1(IJ)
          ENDDO
          DO IA=1,NJA
            JJ  = JACT(IA)
            JD3 = JJ+O3
!$OMP SIMD PRIVATE(X1,X3,TS)
            DO IJ=KIJS,KIJL
              X1 = NX(IJ,JJ,K1)
              X3 = NX(IJ,JD3,K3)
              TS = MM(IJ,JJ) * ( X1*X3*AA(IJ,JJ) + (X3-X1)*BB(IJ,JJ) )
              AS1(IJ,JJ) = AS1(IJ,JJ) + TS * WA(IJ)
              SX(IJ,JD3,K3) = SX(IJ,JD3,K3) - TS * WP(IJ)
            ENDDO
            IF (LSWRTDIAG) THEN
!$OMP SIMD PRIVATE(X1,X3)
              DO IJ=KIJS,KIJL
                X1 = NX(IJ,JJ,K1)
                X3 = NX(IJ,JD3,K3)
                AD1(IJ,JJ) = AD1(IJ,JJ) + MM(IJ,JJ) * ( X3*AA(IJ,JJ) - BB(IJ,JJ) ) * WA(IJ)
                DX(IJ,JD3,K3) = DX(IJ,JD3,K3) - MM(IJ,JJ) * ( X1*AA(IJ,JJ) + BB(IJ,JJ) ) * WP(IJ)
              ENDDO
            ENDIF
          ENDDO
        ENDDO

        DO J=1,NANG
!$OMP SIMD
          DO IJ=KIJS,KIJL
            SX(IJ,J,K1) = SX(IJ,J,K1) + AS1(IJ,J)
            DX(IJ,J,K1) = DX(IJ,J,K1) + AD1(IJ,J)
          ENDDO
        ENDDO
      ENDDO

      IF (LUFCTL) CALL IEEE_SET_UNDERFLOW_MODE(LGRADUAL)

!*    4. FOLD THE PADDING BACK, BACK TO ECWAM DIRECTION ORDER, ADD TO SL
!        AND FLD (dS/dE = dS/dN * CONV).
!        ---------------------------------------------------------------

      DO M=1,NFRE
        DO J=1,NANG
          K = NANG+1-J
!$OMP SIMD
          DO IJ=KIJS,KIJL
            SL(IJ,K,M)  = SL(IJ,K,M)  + REAL(SX(IJ,J,M) + SX(IJ,J-NANG,M) + SX(IJ,J+NANG,M), JWRB)
            FLD(IJ,K,M) = FLD(IJ,K,M) + REAL((DX(IJ,J,M) + DX(IJ,J-NANG,M) + DX(IJ,J+NANG,M)) * CONV(M), JWRB)
          ENDDO
        ENDDO
      ENDDO

      DEALLOCATE(NX, SX, DX)

      IF (LHOOK) CALL DR_HOOK('SNONLIN_SWRT',1,ZHOOK_HANDLE)

      END SUBROUTINE SNONLIN_SWRT

! ----------------------------------------------------------------------

      SUBROUTINE SNONLIN_SWRT_REF (KIJS, KIJL, FL1, FLD, SL, WAVNUM, DEPTH, AKMEAN)

      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWABORT , ONLY : WAM_ABORT
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(INOUT):: FLD, SL
      REAL(KIND=JWRB), DIMENSION(KIJL,NFRE), INTENT(IN) :: WAVNUM
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: DEPTH, AKMEAN

      INTEGER(KIND=JWIM) :: IJ, K, M, J, G, RR, R0, K1, K3, JJ, O3, O2, O4, NU
      INTEGER(KIND=JWIM) :: KA2_, KB2_, KA4_, KB4_

      REAL(KIND=JWRB), DIMENSION(KIJL,NFRE) :: ENH

!     SINGLE PRECISION WORK ARRAYS OF ONE SPECTRUM (LOCAL: THREAD-SAFE)
      REAL(KIND=JWRS) :: NX(1-NANG:2*NANG,NFRE), SX(1-NANG:2*NANG,NFRE), DX(1-NANG:2*NANG,NFRE)
      REAL(KIND=JWRS) :: X1(NANG), X3(NANG), MM(NANG), AA(NANG), BB(NANG), AS1(NANG), AD1(NANG), RMX(NFRE)
      REAL(KIND=JWRS) :: XN, NMAX, QM, X2, X4, A2, B2, A4, B4, WR, WA, WP, U2_, U4_, TS, E1

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SNONLIN_SWRT_REF',0,ZHOOK_HANDLE)

      IF (.NOT. LSWRTREAD) CALL WAM_ABORT('SNONLIN_SWRT_REF: operator not read (INISNONLIN_SWRT not called)')

!*    1. SHALLOW WATER SCALING, AS IN SNONLIN (FREQUENCIES 1 TO NFRE).
!        -------------------------------------------------------------

      CALL SWRT_ENH(KIJS, KIJL, FL1, WAVNUM, DEPTH, AKMEAN, ENH)

!*    2. SWRT, ONE SPECTRUM AT A TIME (LOCUSLOOP OF W3SNL6).
!        ---------------------------------------------------

      DO IJ=KIJS,KIJL

!       N(kx,ky) = E * CONV ON A PADDED PERIODIC DIRECTION AXIS, IN THE
!       OPERATOR'S (ANTI-CLOCKWISE) ORDER: OPERATOR BIN J = ECWAM BIN NANG+1-J
        DO M=1,NFRE
          DO J=1,NANG
            XN = REAL(FL1(IJ,NANG+1-J,M),JWRS)*CONV(M)
            NX(J,M) = XN
            NX(J-NANG,M) = XN
            NX(J+NANG,M) = XN
          ENDDO
          RMX(M) = MAXVAL(NX(1:NANG,M))
        ENDDO
        NMAX = MAXVAL(RMX)
        SX(:,:) = 0.0_JWRS
        DX(:,:) = 0.0_JWRS

        DO K1=1,NFRE
          E1 = REAL(ENH(IJ,K1),JWRS)
          X1(1:NANG) = NX(1:NANG,K1)
          AS1(:) = 0.0_JWRS
          AD1(:) = 0.0_JWRS
          DO G=1,NG
            R0 = GS(G)
            O3 = DJ(R0)
            K3 = K1 + DI(R0)
            IF (K3 > NFRE) CYCLE
            QM = QFRAC * NMAX * KPOW(K3)
            IF (RMX(K1) < QM .AND. RMX(K3) < QM) CYCLE
            DO JJ=1,NANG
              X3(JJ) = NX(JJ+O3,K3)
              MM(JJ) = MERGE(1.0_JWRS, 0.0_JWRS, X1(JJ) >= QM .OR. X3(JJ) >= QM)
            ENDDO
            AA(:) = 0.0_JWRS
            BB(:) = 0.0_JWRS
            NU = 0
            DO RR=GS(G),GS(G+1)-1
              WR = W(K1,RR)
              IF (WR == 0.0_JWRS) CYCLE
              NU = NU + 1
              O2 = J2(RR)
              O4 = J4(RR)
              B2 = G2(RR)
              A2 = 1.0_JWRS - B2
              B4 = G4(RR)
              A4 = 1.0_JWRS - B4
              KA2_ = KA2(K1,RR)
              KB2_ = KB2(K1,RR)
              U2_  = U2(K1,RR)
              KA4_ = KA4(K1,RR)
              KB4_ = KB4(K1,RR)
              U4_  = U4(K1,RR)
              DO JJ=1,NANG
                X2 = (1.0_JWRS-U2_) * ( A2*NX(JJ+O2,KA2_) + B2*NX(JJ+O2+1,KA2_) )   &
     &             +           U2_  * ( A2*NX(JJ+O2,KB2_) + B2*NX(JJ+O2+1,KB2_) )
                X4 = (1.0_JWRS-U4_) * ( A4*NX(JJ+O4,KA4_) + B4*NX(JJ+O4+1,KA4_) )   &
     &             +           U4_  * ( A4*NX(JJ+O4,KB4_) + B4*NX(JJ+O4+1,KB4_) )
                AA(JJ) = AA(JJ) + WR * ( X4 - X2 )
                BB(JJ) = BB(JJ) + WR * X2 * X4
              ENDDO
            ENDDO
            IF (NU == 0) CYCLE
            WA = FA(K1,R0) * E1
            WP = FP(K1,R0) * E1
            DO JJ=1,NANG
              TS = MM(JJ) * ( X1(JJ)*X3(JJ)*AA(JJ) + (X3(JJ)-X1(JJ))*BB(JJ) )
              AS1(JJ) = AS1(JJ) + TS * WA
              SX(JJ+O3,K3) = SX(JJ+O3,K3) - TS * WP
            ENDDO
            IF (LSWRTDIAG) THEN
              DO JJ=1,NANG
                AD1(JJ) = AD1(JJ) + MM(JJ) * ( X3(JJ)*AA(JJ) - BB(JJ) ) * WA
                DX(JJ+O3,K3) = DX(JJ+O3,K3) - MM(JJ) * ( X1(JJ)*AA(JJ) + BB(JJ) ) * WP
              ENDDO
            ENDIF
          ENDDO
          SX(1:NANG,K1) = SX(1:NANG,K1) + AS1(:)
          DX(1:NANG,K1) = DX(1:NANG,K1) + AD1(:)
        ENDDO

!       FOLD THE PADDING BACK, BACK TO ECWAM DIRECTION ORDER, ADD TO SL
!       AND FLD (dS/dE = dS/dN * CONV)
        DO M=1,NFRE
          DO J=1,NANG
            K = NANG+1-J
            SL(IJ,K,M)  = SL(IJ,K,M)  + REAL(SX(J,M) + SX(J-NANG,M) + SX(J+NANG,M), JWRB)
            FLD(IJ,K,M) = FLD(IJ,K,M) + REAL((DX(J,M) + DX(J-NANG,M) + DX(J+NANG,M)) * CONV(M), JWRB)
          ENDDO
        ENDDO

      ENDDO

      IF (LHOOK) CALL DR_HOOK('SNONLIN_SWRT_REF',1,ZHOOK_HANDLE)

      END SUBROUTINE SNONLIN_SWRT_REF

! ----------------------------------------------------------------------

      SUBROUTINE SWRT_ENH (KIJS, KIJL, FL1, WAVNUM, DEPTH, AKMEAN, ENH)

!     SHALLOW WATER SCALING ENH(IJ,M), M = 1..NFRE, EXACTLY AS IN SNONLIN
!     FOR THE ACTIVE ISNONLIN (THE CODE OF ITS THREE CASES IS COPIED).

      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWSTAT  , ONLY : ISNONLIN

      IMPLICIT NONE
#include "peak_ang.intfb.h"
#include "transf.intfb.h"
#include "transf_snl.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJL,NFRE), INTENT(IN) :: WAVNUM
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: DEPTH, AKMEAN
      REAL(KIND=JWRB), DIMENSION(KIJL,NFRE), INTENT(OUT) :: ENH

      INTEGER(KIND=JWIM) :: IJ, M
      REAL(KIND=JWRB), PARAMETER :: ENH_MAX=10.0_JWRB
      REAL(KIND=JWRB), PARAMETER :: ENH_MIN=0.1_JWRB
      REAL(KIND=JWRB), DIMENSION(KIJL) :: ENHFR, XNU, SIG_TH

      SELECT CASE (ISNONLIN)
      CASE(0)
        DO IJ=KIJS,KIJL
          ENHFR(IJ) = MAX(0.75_JWRB*DEPTH(IJ)*AKMEAN(IJ), 0.5_JWRB)
          ENHFR(IJ) = 1.0_JWRB+(5.5_JWRB/ENHFR(IJ)) * (1.0_JWRB-.833_JWRB*ENHFR(IJ)) * EXP(-1.25_JWRB*ENHFR(IJ))
        ENDDO
        DO M=1,NFRE
          DO IJ=KIJS,KIJL
            ENH(IJ,M) = ENHFR(IJ)
          ENDDO
        ENDDO

      CASE(1)
        DO M=1,NFRE
          DO IJ = KIJS, KIJL
            ENH(IJ,M) = MAX(MIN(ENH_MAX,TRANSF(WAVNUM(IJ,M),DEPTH(IJ))),ENH_MIN)
          ENDDO
        ENDDO

      CASE(2)
        CALL PEAK_ANG(KIJS, KIJL, FL1, XNU, SIG_TH)
        DO M=1,NFRE
          DO IJ = KIJS, KIJL
            ENH(IJ,M) = TRANSF_SNL(WAVNUM(IJ,M), DEPTH(IJ), XNU(IJ), SIG_TH(IJ))
          ENDDO
        ENDDO
      END SELECT

      END SUBROUTINE SWRT_ENH

      END MODULE SNONLIN_SWRT_MOD
