      SUBROUTINE UFBINX(LUNIT,IMSG,ISUB,USR,I1,I2,IRET,STR)

C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    UFBINX
C   PRGMMR: WOOLLEN          ORG: NP20       DATE: 2003-11-04
C
C ABSTRACT: THIS SUBROUTINE EITHER OPENS A BUFR FILE CONNECTED TO
C   LOGICAL UNIT LUNIT FOR INPUT OPERATIONS (IF IT IS NOT ALREADY
C   OPENED AS SUCH), OR SAVES ITS POSITION AND REWINDS IT TO THE FIRST
C   DATA MESSAGE (IF BUFR FILE ALREADY OPENED), THEN (VIA A CALL TO
C   BUFR ARCHIVE LIBRARY SUBROUTINE UFBINT) READS SPECIFIED VALUES FROM
C   INTERNAL SUBSET ARRAYS ASSOCIATED WITH A PARTICULAR SUBSET FROM A
C   PARTICULAR BUFR MESSAGE IN A MESSAGE BUFFER.  THE PARTICULAR SUBSET
C   AND BUFR MESSAGE ARE BASED BASED ON THE SUBSET NUMBER IN THE
C   MESSAGE AND THE MESSAGE NUMBER IN THE BUFR FILE.  FINALLY, THIS
C   SUBROUTINE EITHER CLOSES THE BUFR FILE IN LUNIT (IF IS WAS OPENED
C   HERE) OR RESTORES IT TO ITS PREVIOUS READ/WRITE STATUS AND POSITION
C   (IF IT WAS NOT OPENED HERE).  SEE UFBINT FOR MORE INFORMATION ON
C   THE READING OF VALUES OUT OF A BUFR MESSAGE SUBSET.  NOTE: THE
C   MESSAGE NUMBER HERE DOES NOT INCLUDE THE DICTIONARY MESSAGES AT THE
C   BEGINNING OF THE FILE.
C
C PROGRAM HISTORY LOG:
C 2003-11-04  J. WOOLLEN -- ORIGINAL AUTHOR (WAS IN VERIFICATION
C                           VERSION BUT MAY HAVE BEEN IN THE PRODUCTION
C                           VERSION AT ONE TIME AND THEN REMOVED)
C 2003-11-04  D. KEYSER  -- UNIFIED/PORTABLE FOR WRF; ADDED
C                           DOCUMENTATION; OUTPUTS MORE COMPLETE
C                           DIAGNOSTIC INFO WHEN ROUTINE TERMINATES
C                           ABNORMALLY
C 2004-08-09  J. ATOR    -- MAXIMUM MESSAGE LENGTH INCREASED FROM
C                           20,000 TO 50,000 BYTES
C DART $Id$
C
C USAGE:    CALL UFBINX (LUNIT, IMSG, ISUB, USR, I1, I2, IRET, STR)
C   INPUT ARGUMENT LIST:
C     LUNIT    - INTEGER: FORTRAN LOGICAL UNIT NUMBER FOR BUFR FILE
C     IMSG     - INTEGER: POINTER TO BUFR MESSAGE NUMBER TO READ IN
C                BUFR FILE
C     ISUB     - INTEGER: POINTER TO SUBSET NUMBER TO READ IN BUFR
C                MESSAGE
C     I1       - INTEGER: LENGTH OF FIRST DIMENSION OF USR OR THE
C                NUMBER OF BLANK-SEPARATED MNEMONICS IN STR (FORMER
C                MUST BE AT LEAST AS LARGE AS LATTER)
C     I2       - INTEGER: LENGTH OF SECOND DIMENSION OF USR
C     STR      - CHARACTER*(*): STRING OF BLANK-SEPARATED TABLE B
C                MNEMONICS IN ONE-TO-ONE CORRESPONDENCE WITH FIRST
C                DIMENSION OF USR {THIS CAN ALSO BE A SINGLE TABLE D
C                (SEQUENCE) MNEMONIC WITH EITHER 8- OR 16-BIT DELAYED
C                REPLICATION (SEE REMARKS 1 IN UFBINT DOCBLOCK)}
C
C   OUTPUT ARGUMENT LIST:
C     USR      - REAL*8: (I1,I2) STARTING ADDRESS OF DATA VALUES READ
C                FROM DATA SUBSET
C     IRET     - INTEGER: NUMBER OF "LEVELS" OF DATA VALUES READ FROM
C                DATA SUBSET (MUST BE NO LARGER THAN I2)
C
C   INPUT FILES:
C     UNIT "LUNIT" - BUFR FILE
C
C REMARKS:
C    THIS ROUTINE CALLS:        BORT     CLOSBF   OPENBF   READMG
C                               READSB   REWNBF   STATUS   UFBINT
C                               UPB
C    THIS ROUTINE IS CALLED BY: None
C                               Normally called only by application
C                               programs.
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN 77
C   MACHINE:  PORTABLE TO ALL PLATFORMS
C
C$$$

      INCLUDE 'bufrlib.prm'

      COMMON /MSGCWD/ NMSG(NFILES),NSUB(NFILES),MSUB(NFILES),
     .                INODE(NFILES),IDATE(NFILES)
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(MXMSGLD4),MBYT(NFILES),
     .                MBAY(MXMSGLD4,NFILES)

      CHARACTER*(*) STR
      CHARACTER*128 BORT_STR
      CHARACTER*8   SUBSET
      LOGICAL       OPENIT
      REAL*8        USR(I1,I2)

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      CALL STATUS(LUNIT,LUN,IL,IM)
      OPENIT = IL.EQ.0

      IF(OPENIT) THEN

C  OPEN BUFR FILE CONNECTED TO UNIT LUNIT IF IT IS NOT ALREADY OPEN
C  ----------------------------------------------------------------

         CALL OPENBF(LUNIT,'IN',LUNIT)
      ELSE

C  IF BUFR FILE ALREADY OPENED, SAVE POSITION & REWIND TO FIRST DATA MSG
C  ---------------------------------------------------------------------

         CALL REWNBF(LUNIT,0)
      ENDIF

C  SKIP TO MESSAGE # IMSG
C  ----------------------

      DO I=1,IMSG-1
      READ(LUNIT,ERR=900,END=901)
      ENDDO

      CALL READMG(LUNIT,SUBSET,JDATE,JRET)
      IF(JRET.NE.0) GOTO 901

C  POSITION AT SUBSET # ISUB
C  -------------------------

      DO I=1,ISUB-1
      IF(NSUB(LUN).GT.MSUB(LUN)) GOTO 902
      IBIT = MBYT(LUN)*8
      CALL UPB(NBYT,16,MBAY(1,LUN),IBIT)
      MBYT(LUN) = MBYT(LUN) + NBYT
      NSUB(LUN) = NSUB(LUN) + 1
      ENDDO

      CALL READSB(LUNIT,JRET)
      IF(JRET.NE.0) GOTO 902

      CALL UFBINT(LUNIT,USR,I1,I2,IRET,STR)

      IF(OPENIT) THEN

C  CLOSE BUFR FILE IF IT WAS OPENED HERE
C  -------------------------------------

         CALL CLOSBF(LUNIT)
      ELSE


C  RESTORE BUFR FILE TO PREV. STATUS & POSITION IF NOT ORIG. OPENED HERE
C  ---------------------------------------------------------------------

         CALL REWNBF(LUNIT,1)
      ENDIF

C  EXITS
C  -----

      RETURN
900   WRITE(BORT_STR,'("BUFRLIB: UFBINX - ERROR READING MESSAGE '//
     . '(RECORD) NUMBER",I5," IN INPUT BUFR FILE CONNECTED TO UNIT",'//
     . 'I4)')  I,LUNIT
      CALL BORT(BORT_STR)
901   WRITE(BORT_STR,'("BUFRLIB: UFBINX - HIT END OF FILE BEFORE '//
     . 'READING REQUESTED MESSAGE NO.",I5," IN BUFR FILE CONNECTED TO'//
     . ' UNIT",I4)')  IMSG,LUNIT
      CALL BORT(BORT_STR)
902   WRITE(BORT_STR,'("BUFRLIB: UFBINX - ALL SUBSETS READ BEFORE '//
     . 'READING REQ. SUBSET NO.",I3," IN REQ. MSG NO.",I5," IN BUFR '//
     . 'FILE CONNECTED TO UNIT",I4)') ISUB,IMSG,LUNIT
      CALL BORT(BORT_STR)
      END