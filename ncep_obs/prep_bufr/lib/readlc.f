      SUBROUTINE READLC(LUNIT,CHR,STR)

C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    READLC
C   PRGMMR: WOOLLEN          ORG: NP20       DATE: 2003-11-04
C
C ABSTRACT: THIS SUBROUTINE RETURNS A CHARACTER DATA ELEMENT ASSOCIATED
C   WITH A PARTICULAR SUBSET MNEMONIC FROM THE INTERNAL MESSAGE BUFFER
C   (ARRAY MBAY IN COMMON BLOCK /BITBUF/).  IT IS DESIGNED TO BE USED
C   TO RETURN CHARACTER ELEMENTS GREATER THAN THE USUAL LENGTH OF EIGHT
C   BYTES.  IT CURRENTLY WILL NOT WORK FOR COMPRESSED BUFR MESSAGES.
C
C PROGRAM HISTORY LOG:
C 2003-11-04  J. WOOLLEN -- ORIGINAL AUTHOR
C 2003-11-04  D. KEYSER  -- UNIFIED/PORTABLE FOR WRF; ADDED
C                           DOCUMENTATION; OUTPUTS MORE COMPLETE
C                           DIAGNOSTIC INFO WHEN ROUTINE TERMINATES
C                           ABNORMALLY OR UNUSUAL THINGS HAPPEN
C 2004-08-09  J. ATOR    -- MAXIMUM MESSAGE LENGTH INCREASED FROM
C                           20,000 TO 50,000 BYTES
C
C USAGE:    CALL READLC (LUNIT, CHR, STR)
C   INPUT ARGUMENT LIST:
C     LUNIT    - INTEGER: FORTRAN LOGICAL UNIT NUMBER FOR BUFR FILE
C     STR      - CHARACTER*(*): STRING (I.E., MNEMONIC)
C
C   OUTPUT ARGUMENT LIST:
C     CHR      - CHARACTER*(*): UNPACKED CHARACTER STRING (I.E.,
C                CHARACTER DATA ELEMENT GREATER THAN EIGHT BYTES)
C
C   OUTPUT FILES:
C     UNIT 06  - STANDARD OUTPUT PRINT
C
C REMARKS:
C    THIS ROUTINE CALLS:        BORT     PARSEQ   STATUS   UPC
C    THIS ROUTINE IS CALLED BY: UFDUMP
C                               Normally not called by any application
C                               programs.
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN 77
C   MACHINE:  PORTABLE TO ALL PLATFORMS
C
C$$$

      INCLUDE 'bufrlib.prm'

      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(MXMSGLD4),MBYT(NFILES),
     .                MBAY(MXMSGLD4,NFILES)
      COMMON /TABLES/ MAXTAB,NTAB,TAG(MAXJL),TYP(MAXJL),KNT(MAXJL),
     .                JUMP(MAXJL),LINK(MAXJL),JMPB(MAXJL),
     .                IBT(MAXJL),IRF(MAXJL),ISC(MAXJL),
     .                ITP(MAXJL),VALI(MAXJL),KNTI(MAXJL),
     .                ISEQ(MAXJL,2),JSEQ(MAXJL)
      COMMON /USRINT/ NVAL(NFILES),INV(MAXJL,NFILES),VAL(MAXJL,NFILES)
      COMMON /USRBIT/ NBIT(MAXJL),MBIT(MAXJL)
      COMMON /UNPTYP/ MSGUNP(NFILES)
      COMMON /QUIET / IPRT

      CHARACTER*(*) CHR,STR
      CHARACTER*128 BORT_STR
      CHARACTER*10  TAG,TGS(100)
      CHARACTER*8   CTAG
      CHARACTER*3   TYP
      REAL*8        VAL

      DATA MAXTG /100/

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C  CHECK THE FILE STATUS
C  ---------------------

      CALL STATUS(LUNIT,LUN,IL,IM)
      IF(IL.EQ.0) GOTO 900
      IF(IL.GT.0) GOTO 901
      IF(IM.EQ.0) GOTO 902

C  CHECK FOR TAGS (MNEMONICS) IN INPUT STRING (THERE CAN ONLY BE ONE)
C  ------------------------------------------------------------------

      CALL PARSEQ(STR,TGS,MAXTG,NTG)
      IF(NTG.GT.1) GOTO 903
      CTAG = TGS(1)
      CHR = ' '

C  FIND THE TAG IN THE SUBSET OR RETURN A BLANK STRING
C  ---------------------------------------------------

      DO N=1,NVAL(LUN)
      NOD = INV(N,LUN)
      IF(CTAG.EQ.TAG(NOD)) GOTO 1
      ENDDO
      IF(IPRT.GE.0) THEN
      PRINT*
      PRINT*,'+++++++++++++++++BUFR ARCHIVE LIBRARY++++++++++++++++++++'
      PRINT*, 'BUFRLIB: READLC - MNEMONIC ',CTAG,' NOT LOCATED IN ',
     .        'REPORT SUBSET - RETURN WITH BLANK STRING FOR CHARACTER ',
     .        'DATA ELEMENT'
      PRINT*,'+++++++++++++++++BUFR ARCHIVE LIBRARY++++++++++++++++++++'
      PRINT*
      ENDIF
      GOTO 100
1     IF(ITP(NOD).NE.3) GOTO 904

C  DECIPHER THE LONG CHARACTER
C  ---------------------------

      IF(MSGUNP(LUN).EQ.0.OR.MSGUNP(LUN).EQ.1) THEN
         NCHR = NBIT(N)/8
         KBIT = MBIT(N)
         CALL UPC(CHR,NCHR,MBAY(1,LUN),KBIT)
      ELSEIF(MSGUNP(LUN).EQ.2) THEN
c  .... compressed message
         GOTO 905
      ELSE
         GOTO 906
      ENDIF

C  EXITS
C  -----

100   RETURN
900   CALL BORT('BUFRLIB: READLC - INPUT BUFR FILE IS CLOSED, IT MUST'//
     . ' BE OPEN FOR INPUT')
901   CALL BORT('BUFRLIB: READLC - INPUT BUFR FILE IS OPEN FOR '//
     . 'OUTPUT, IT MUST BE OPEN FOR INPUT')
902   CALL BORT('BUFRLIB: READLC - A MESSAGE MUST BE OPEN IN INPUT '//
     . 'BUFR FILE, NONE ARE')
903   WRITE(BORT_STR,'("BUFRLIB: READLC - THERE CANNOT BE MORE THAN '//
     . 'ONE MNEMONIC IN THE INPUT STRING (",A,") (HERE THERE ARE ",'//
     . 'I3,")")') STR,NTG
      CALL BORT(BORT_STR)
904   WRITE(BORT_STR,'("BUFRLIB: READLC - MNEMONIC ",A," DOES NOT '//
     . 'REPRESENT A CHARACTER ELEMENT (ITP=",I2,")")') CTAG,ITP(NOD)
      CALL BORT(BORT_STR)
905   CALL BORT('BUFRLIB: READLC - NOT ENABLED FOR COMPRESSED BUFR '//
     . 'MESSAGES')
906   WRITE(BORT_STR,'("BUFRLIB: READLC - MESSAGE UNPACK TYPE",I3,'//
     . '" IS NOT RECOGNIZED")') MSGUNP
      CALL BORT(BORT_STR)
      END