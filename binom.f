	

	subroutine binom(n,p,m)
	implicit double precision (a-h,o-z)

	real u(1)

c	xrr=0.99
c	if(p.ge.xrr) then  ! this should never be used
c	p=0.9
c	endif
	
	call random_number(r1)

	xlog=-dlog(1.0d0-p)
	
	xsum=0.0d0

	i=1
c  10	call uniran(1,iseed,u)

10      call random_number(r1)
	
c	xx1=u(1)

	xx1=r1

	xx2=dfloat(n-i+1)
	hx=-dlog(xx1)/xx2
	xsum=xsum+hx
	
	if(xsum.gt.xlog) goto 100

	i=i+1
	goto 10	
	
100	m=i-1

	
	return
	end 
c*****************************************************************************************

	subroutine expnum(xnav,nmax,nu,iseed)
	! exponential distributed integer from 0 to nmax with average nav
	implicit double precision (a-h,o-z)
	real x(1)
	

	xlambda=1.0/xnav
	
	call uniran(1,iseed,x)
		xx=x(1)

	u=-dlog(1.0d0-xx)/xlambda

	nu=int(u)

	if(nu.gt.nmax) then
	nu=nmax
	endif 

	return
	end 



c ******************************************************************************************

      SUBROUTINE UNIRAN(N,ISEED,X)
C
C     PURPOSE--THIS SUBROUTINE GENERATES A RANDOM SAMPLE OF SIZE N
C              FROM THE UNIFORM (RECTANGULAR)
C              DISTRIBUTION ON THE UNIT INTERVAL (0,1).
C              THIS DISTRIBUTION HAS MEAN = 0.5
C              AND STANDARD DEVIATION = SQRT(1/12) = 0.28867513.
C              THIS DISTRIBUTION HAS THE PROBABILITY
C              DENSITY FUNCTION F(X) = 1.
C
C     INPUT  ARGUMENTS--N      = THE DESIRED INTEGER NUMBER 
C                                OF RANDOM NUMBERS TO BE
C                                GENERATED.
C                     --ISEED  = AN INTEGER ISEED VALUE
C     OUTPUT ARGUMENTS--X      = A SINGLE PRECISION VECTOR
C                                (OF DIMENSION AT LEAST N)
C                                INTO WHICH THE GENERATED
C                                RANDOM SAMPLE WILL BE PLACED.
C     OUTPUT--A RANDOM SAMPLE OF SIZE N 
C             FROM THE RECTANGULAR DISTRIBUTION ON (0,1).
C     PRINTING--NONE UNLESS AN INPUT ARGUMENT ERROR CONDITION EXISTS. 
C     RESTRICTIONS--THERE IS NO RESTRICTION ON THE MAXIMUM VALUE
C                   OF N FOR THIS SUBROUTINE.
C     OTHER           SUBROUTINES NEEDED--NONE.
C     FORTRAN LIBRARY SUBROUTINES NEEDED--NONE.
C     MODE OF INTERNAL OPERATIONS--SINGLE PRECISION.
C     LANGUAGE--ANSI FORTRAN (1977)
C
C     ALGORITHM--FIBONACCI GENERATOR
C                AS DEFINED BY GEORGE MARSAGLIA.
C
C     NOTE--THIS GENERATOR IS TRANSPORTABLE.
C           IT IS NOT MACHINE-INDEPENDENT
C           IN THE SENSE THAT FOR A GIVEN VALUE
C           OF THE INPUT SEED ISEED AND FOR A GIVEN VALUE
C           OF MDIG (TO BE DEFINED BELOW),
C           THE SAME SEQUENCE OF UNIRFORM RANDOM
C           NUMBERS WILL RESULT ON DIFFERENT COMPUTERS
C           (VAX, PRIME, PERKIN-ELMER, IBM, UNIVAC, HONEYWELL, ETC.)
C
C     NOTE--IF MDIG = 32 AND IF ISEED = 305,
C           THEN THE OUTPUT FROM THIS GENERATOR SHOULD BE AS FOLLOWS--
C           THE FIRST      NUMBER TO RESULT IS .4771580...
C           THE SECOND     NUMBER TO RESULT IS .4219293...
C           THE THIRD      NUMBER TO RESULT IS .6646181...
C           ...
C           THE THOUSANDTH NUMBER TO RESULT IS .2036834...
C
C     NOTE--IF MDIG = 16 AND IF ISEED = 305,
C           THEN THE OUTPUT FROM THIS GENERATOR SHOULD BE AS FOLLOWS--
C           THE FIRST      NUMBER TO RESULT IS .027832881...
C           THE SECOND     NUMBER TO RESULT IS .56102176... 
C           THE THIRD      NUMBER TO RESULT IS .41456343... 
C           ...
C           THE THOUSANDTH NUMBER TO RESULT IS .19797357... 
C
C     NOTE--IT IS RECOMMENDED THAT UPON 
C           IMPLEMENTATION OF DATAPLOT, THE OUTPUT
C           FROM UNIRAN BE CHECKED FOR AGREEMENT
C           WITH THE ABOVE SAMPLE OUTPUT.
C           ALSO, THERE ARE MANY ANALYSIS AND DIAGNOSTIC
C           TOOLS IN DATAPLOT THAT WILL ALLOW THE 
C           TESTING OF THE RANDOMNESS AND UNIFORMITY
C           OF THIS GENERATOR.
C           SUCH CHECKING IS ESPECIALLY IMPORTANT 
C           IN LIGHT OF THE FACT THAT OTHER DATAPLOT RANDOM 
C           NUMBER GENERATOR SUBROUTINES (NORRAN--NORMAL,
C           LOGRAN--LOGISTIC, ETC.) ALL MAKE USE OF INTERMEDIATE
C           OUTPUT FROM UNIRAN.
C
C     NOTE--THE OUTPUT FROM THIS SUBROUTINE DEPENDS
C           ON THE INPUT SEED (ISEED) AND ON THE
C           VALUE OF MDIG.
C           MDIG MAY NOT BE SMALLER THAN 16.
C           MDIG MAY NOT BE LARGER THAN MAX INTEGER ON YOUR COMPUTER. 
C
C     NOTE--BECAUSE OF THE PREPONDERANCE OF MAINFRAMES
C           WHICH HAVE WORDS OF 32 BITS AND LARGER
C           (E.G, VAX (= 32 BITS), UNIVAC (= 36 BITS), CDC (= 60 BITS), ETC.)
C           MDIG HAS BEEN SET TO 32.
C           THUS THE SAME SEQUENCE OF RANDOM NUMBERS SHOULD RESULT
C           ON ALL OF THESE COMPUTERS.
C
C     NOTE--FOR SMALLER WORD SIZE COMPUTERS (E.G., 24-BIT AND 16-BIT),
C           THE VALUE OF MDIG SHOULD BE CHANGED TO 24 OR 16.
C           IN SUCH CASE, THE OUTPUT WILL NOT BE IDENTICAL TO
C           THE OUTPUT WHEN MDIG = 32.
C
C     NOTE--THE CYCLE OF THE RANDOM NUMBERS DEPENDS ON MDIG.
C           THE CYCLE FROM MDIG = 32 IS LONG ENOUGH FOR MOST
C           PRACTICAL APPLICATIONS.
C           IF A LONGER CYCLE IS DESIRED, THEN INCREASE MDIG.
C
C     NOTE--THE SEED MAY BE ANY POSITIVE INTEGER. 
C           NO APPRECIABLE DIFFERENCE IN THE QUALITY
C           OF THE RANDOM NUMBERS HAS BEEN NOTED
C           BY THE CHOICE OF THE SEED.  THERE IS NO
C           NEED TO USE PRIMES, NOR TO USE EXCEPTIONALLY
C           LARGE NUMBERS, ETC.
C
C     REFERENCES--MARSAGLIA G., "COMMENTS ON THE PERFECT UNIFORM RANDOM
C                 NUMBER GENERATOR", UNPUBLISHED NOTES, WASH S. U.
C               --JOHNSON AND KOTZ, CONTINUOUS UNIVARIATE
C                 DISTRIBUTIONS--2, 1970, PAGES 57-74.
C     WRITTEN BY--JAMES BLUE
C                 SCIENTIFIC COMPUTING DIVISION
C                 CENTER FOR APPLIED MATHEMATICS
C                 NATIONAL BUREAU OF STANDARDS
C                 WASHINGTON, D. C. 20234
C               --DAVID KAHANER
C                 SCIENTIFIC COMPUTING DIVISION
C                 CENTER FOR APPLIED MATHEMATICS
C                 NATIONAL BUREAU OF STANDARDS
C               --GEORGE MARSAGLIA
C                 COMPUTER SCIENCE DEPARTMENT
C                 WASHINGTON STATE UNIVERSITY
C               --JAMES J. FILLIBEN
C                 STATISTICAL ENGINEERING DIVISION
C                 CENTER FOR APPLIED MATHEMATICS
C                 NATIONAL BUREAU OF STANDARDS
C
C     LANGUAGE--ANSI FORTRAN (1977)
C     ORIGINAL VERSION--JUNE      1972. 
C     UPDATED         --AUGUST    1974. 
C     UPDATED         --SEPTEMBER 1975. 
C     UPDATED         --NOVEMBER  1975. 
C     UPDATED         --NOVEMBER  1981. 
C     UPDATED         --MAY       1982. 
C     UPDATED         --MARCH     1984. 
C
C-----CHARACTER STATEMENTS FOR NON-COMMON VARIABLES-------------------
C
C---------------------------------------------------------------------
C
      DIMENSION X(*)
C
      DIMENSION M(17)
C
C---------------------------------------------------------------------
C
CCCCC CHARACTER*4 IFEEDB
CCCCC CHARACTER*4 IPRINT
C
CCCCC COMMON /MACH/IRD,IPR,CPUMIN,CPUMAX,NUMBPC,NUMCPW,NUMBPW
CCCCC COMMON /PRINT/IFEEDB,IPRINT
C
C-----SAVE STATEMENTS-------------------------------------------------
C
      SAVE I,J,M,M1,M2
C
C-----DATA STATEMENTS-------------------------------------------------
C
      DATA M(1),M(2),M(3),M(4),M(5),M(6),M(7),M(8),M(9),M(10),M(11),
     1     M(12),M(13),M(14),M(15),M(16),M(17)
     1/    30788,23052,2053,19346,10646,19427,23975,
     1     19049,10949,19693,29746,26748,2796,23890,
     1     29168,31924,16499/ 
      DATA M1,M2,I,J / 32767,256,5,17 / 
C
      IPR=6
C
C-----START POINT-----------------------------------------------------
C
C               ********************************************
C               **  STEP 1--                              **
C               **  CHECK THE INPUT ARGUMENTS FOR ERRORS  **
C               ********************************************
C
      IF(N.GE.1)GOTO90
      WRITE(IPR,999)
  999 FORMAT(1H )
      WRITE(IPR,51) 
   51 FORMAT(1H ,'***** ERROR IN UNIRAN--')
      WRITE(IPR,52) 
   52 FORMAT(1H ,'      THE INPUT NUMBER OF OBSERVATIONS IS ',
     1'NON-POSITIVE.')
      WRITE(IPR,53)N
   53 FORMAT(1H ,'      N = ',I8)
      GOTO9000
   90 CONTINUE
C
C               *******************************************************
C               **  STEP 2--                                         **
C               **  IF A POSITIVE INPUT SEED HAS BEEN GIVEN,         **
C               **  THEN THIS INDICATES THAT THE GENERATOR           **
C               **  SHOULD HAVE ITS INTERNAL M(.) ARRAY REDEFINED--  **
C               **  DO SO IN THIS SECTION.                           **
C               **  IF A NON-POSITIVE INPUT SEED HAS BEEN GIVEN,     **
C               **  THEN THIS INDICATES THAT THE GENERATOR           **
C               **  SHOULD CONTINUE ON FROM WHERE IT LEFT OFF,       **
C               **  AND THEREFORE THIS SECTION IS SKIPPED.           **
C               *******************************************************
C
      IF(ISEED.LE.0)GOTO290
C
CCCCC MDIG=16
      MDIG=32
C
      M1=2**(MDIG-2)+(2**(MDIG-2)-1)
      M2=2**(MDIG/2)
CCCCC ISEED3=MIN0(IABS(ISEED),M1)
      ISEED3=IABS(ISEED)
      IF(M1.LT.IABS(ISEED))ISEED3=M1
      IF(MOD(ISEED3,2).EQ.0)ISEED3=ISEED3-1
      K0=MOD(9069,M2)
      K1=9069/M2
      J0=MOD(ISEED3,M2)
      J1=ISEED3/M2
C
      DO200I=1,17
      ISEED3=J0*K0
      J1=MOD(ISEED3/M2+J0*K1+J1*K0,M2/2)
      J0=MOD(ISEED3,M2)
      M(I)=J0+M2*J1 
  200 CONTINUE
C
      I=5 
      J=17
C
  290 CONTINUE
C
C               *************************************
C               **  STEP 3--                       **
C               **  GENERATE THE N RANDOM NUMBERS  **
C               *************************************
C
      DO300L=1,N
      K=M(I)-M(J)
      IF(K.LT.0)K=K+M1
      M(J)=K
      I=I-1
      IF(I.EQ.0)I=17
      J=J-1
      IF(J.EQ.0)J=17
      AK=K
      AM1=M1
      X(L)=AK/AM1
  300 CONTINUE
C
C               ***************************************************** 
C               **  STEP 4--                                       ** 
C               **  REGARDLESS OF THE VALUE OF THE INPUT SEED,     ** 
C               **  REDEFINE THE VALUE OF ISEED UPON EXIT HERE     ** 
C               **  TO -1 WITH THE NET EFFECT THAT                 ** 
C               **  IF THE USER DOES NOT REDEFINE THE SEED         ** 
C               **  VALUE BEFORE THE NEXT CALL TO THIS GENERATOR,  ** 
C               **  THEN THIS GENERATOR WILL PICK UP               ** 
C               **  WHERE IT LEFT OFF.                             ** 
C               ***************************************************** 
C
      ISEED=(-1)
C
C               *****************
C               **  STEP 90--  **
C               **  EXIT       **
C               *****************
C
 9000 CONTINUE
      RETURN
CCCCC DEBUG TRACE,INIT
CCCCC AT 90
CCCCC TRACE ON
      END 
