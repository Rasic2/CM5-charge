C     MAY, 2015
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (MAT=30000)
      DIMENSION Q(3,MAT),IZ(MAT),CHIR(MAT),CCM5(MAT),CRCM5(MAT),CN(MAT)
      CHARACTER SWITCH*10,W1*1,W10*10,W26*26,W50*50,W256*256
C
      CALL GETARG(1,SWITCH)
      IF (SWITCH=='1' .OR. SWITCH=='2') THEN
       CALL CM5GAUSS(SWITCH)
      ELSE IF (SWITCH=='3'.OR.SWITCH=='4') THEN
       CALL CM5PS(SWITCH)
      ELSE
        WRITE(*,'(A)')'CM5ERROR'
      ENDIF
      END
C
      SUBROUTINE CM5GAUSS(SWITCH)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (MAT=30000)
      DIMENSION Q(3,MAT),IZ(MAT),CHIR(MAT),CCM5(MAT),CRCM5(MAT),CN(MAT)
      CHARACTER SWITCH*10,W1*1,W10*10,W26*26,W50*50,W256*256
      NAT=0
5     CONTINUE
      READ (*,'(A)',END=100) W256
      WRITE (*,'(A)') W256
      W50=W256
      IF (W50.NE.'                         Standard orientation:')
     $ GOTO 5
      READ (*,'(A)') W256
      WRITE (*,'(A)') W256
      READ (*,'(A)') W256
      WRITE (*,'(A)') W256
      READ (*,'(A)') W256
      WRITE (*,'(A)') W256
      READ (*,'(A)') W256
      WRITE (*,'(A)') W256
      DO K=1,MAT 
      READ (*,'(A)') W256
      WRITE (*,'(A)') W256
      W10=W256
      IF (W10.EQ.' ---------') GOTO 10                                                         
      READ (W256,*) I00,IZ(K),I00,(Q(J,K),J=1,3)
      ENDDO
10    CONTINUE
      NAT=K-1
20    CONTINUE
      READ (*,'(A)',END=90) W256
      WRITE (*,'(A)') W256
      W50=W256
      IF (W50.NE.
     $' Hirshfeld spin densities, charges and dipoles usi') GOTO 20
      READ (*,'(A)') W256
      WRITE (*,'(A)') W256
      DO K=1,NAT
      READ (*,'(A)') W256
      WRITE (*,'(A)') W256
      READ (W256,*) I0,W1,A1,CHIR(K)
      ENDDO
21    CONTINUE
      READ (*,'(A)',END=90) W256
      WRITE (*,'(A)') W256
      W26=W256
      IF (W26.NE.
     $' Sum of Hirshfeld charges=') GOTO 21
C
      CALL CM5MOD(NAT,IZ,CHIR,Q,CCM5,CRCM5,DHIRX,DHIRY,DHIRZ,
     $ DCM5X,DCM5Y,DCM5Z,DRCM5X,DRCM5Y,DRCM5Z,CN)
      DRCM5=DSQRT(DRCM5X**2+DRCM5Y**2+DRCM5Z**2)
      DCM5=DSQRT(DCM5X**2+DCM5Y**2+DCM5Z**2)
      DHIR=DSQRT(DHIRX**2+DHIRY**2+DHIRZ**2)
C
      IF (SWITCH=='1') THEN
      WRITE (*,'(/,A,/,A)')
     $ ' Charges (in A.U.) from CM5PAC version 2015',
     $ ' -----------------------------------------------'
      WRITE (*,'(A,/,A)')
     $ ' Center     Atomic      CM5         Hirshfeld',            
     $ ' Number     Number      Charge      Charge'             
      WRITE (*,'(A)')
     $ ' -----------------------------------------------'
        DO I=1,NAT    
        WRITE (*,'(I5,6X,I5,5X,F11.6,X,F11.6)') 
     $  I,IZ(I),CCM5(I),CHIR(I)            
        ENDDO 
      WRITE (*,'(A)')
     $ ' -----------------------------------------------'
C
      WRITE (*,'(/,A,/,A)')
     $ ' Dipole moment (in Debye)',
     $ ' -----------------------------------------------'
      WRITE (*,'(A,/,A)') 
     $ '                 X        Y        Z     Total', 
     $ ' -----------------------------------------------'
      WRITE (*,'(A,5F9.4)') ' CM5       ',DCM5X,DCM5Y,DCM5Z,DCM5
      WRITE (*,'(A,5F9.4)') ' Hirshfeld ',DHIRX,DHIRY,DHIRZ,DHIR
      WRITE (*,'(A,/)')
     $ ' -----------------------------------------------'
C
      ELSE 
      WRITE (*,'(/,A,/,A)')
     $ ' Charges (in A.U.) from CM5PAC version 2015',
     $ ' -----------------------------------------------'
      WRITE (*,'(A,/,A)')
     $ ' Center     Atomic      CM5M        Hirshfeld',
     $ ' Number     Number      Charge      Charge'
      WRITE (*,'(A)')
     $ ' -----------------------------------------------'
        DO I=1,NAT
        WRITE (*,'(I5,6X,I5,5X,F11.6,X,F11.6)')
     $  I,IZ(I),CRCM5(I),CHIR(I)
        ENDDO
      WRITE (*,'(A)')
     $ ' -----------------------------------------------'
C
      WRITE (*,'(/,A,/,A)')
     $ ' Dipole moment (in Debye)',
     $ ' -----------------------------------------------'
      WRITE (*,'(A,/,A)')
     $ '                 X        Y        Z     Total',
     $ ' -----------------------------------------------'
      WRITE (*,'(A,5F9.4)') ' CM5M      '
     $                        ,DRCM5X,DRCM5Y,DRCM5Z,DRCM5
      WRITE (*,'(A,5F9.4)') ' Hirshfeld ',DHIRX,DHIRY,DHIRZ,DHIR
      WRITE (*,'(A,/)')
     $ ' -----------------------------------------------'
      ENDIF
C
      GOTO 5
      GOTO 100
90    CONTINUE
      WRITE (*,'(A,A)') 'CM5ERROR: Hirshfeld charges are not found.',
     $ ' Perhaps, you need to specify pop=Hirshfeld in the input.'
100   CONTINUE
      IF (NAT.EQ.0) THEN
      WRITE (*,'(A)') 'CM5ERROR: Cartesian coodinates are not found'
      ENDIF
      END
C
      SUBROUTINE CM5MOD(NAT,IZ,CHIR,Q,CCM5,CRCM5,DHIRX,DHIRY,DHIRZ,
     $ DCM5X,DCM5Y,DCM5Z,DRCM5X,DRCM5Y,DRCM5Z,CN)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (MZ=118)
      DIMENSION Q(3,*),IZ(*),CHIR(*),CCM5(*),CRCM5(*)
      DIMENSION RAD(MZ),A0(MZ),D(MZ,MZ),CN(*),ZLAMDA(MZ)
C
C COVALENT RADII
C
C based on "Atomic Radii of the Elements," M. Mantina, R. Valero, C. J. Cramer, and D. G. Truhlar,
C in CRC Handbook of Chemistry and Physics, 91st Edition (2010-2011),
C edited by W. M. Haynes (CRC Press, Boca Raton, FL, 2010), pages 9-49-9-50;
C corrected Nov. 17, 2010 for the 92nd edition.
C
      RAD(1)=0.32D0
      RAD(2)=0.37D0
      RAD(3)=1.30D0
      RAD(4)=0.99D0
      RAD(5)=0.84D0
      RAD(6)=0.75D0
      RAD(7)=0.71D0
      RAD(8)=0.64D0
      RAD(9)=0.60D0
      RAD(10)=0.62D0
      RAD(11)=1.60D0
      RAD(12)=1.40D0
      RAD(13)=1.24D0
      RAD(14)=1.14D0
      RAD(15)=1.09D0
      RAD(16)=1.04D0
      RAD(17)=1.00D0
      RAD(18)=1.01D0
      RAD(19)=2.00D0
      RAD(20)=1.74D0
      RAD(21)=1.59D0
      RAD(22)=1.48D0
      RAD(23)=1.44D0
      RAD(24)=1.30D0
      RAD(25)=1.29D0
      RAD(26)=1.24D0
      RAD(27)=1.18D0
      RAD(28)=1.17D0
      RAD(29)=1.22D0
      RAD(30)=1.20D0
      RAD(31)=1.23D0
      RAD(32)=1.20D0
      RAD(33)=1.20D0
      RAD(34)=1.18D0
      RAD(35)=1.17D0
      RAD(36)=1.16D0
      RAD(37)=2.15D0
      RAD(38)=1.90D0
      RAD(39)=1.76D0
      RAD(40)=1.64D0
      RAD(41)=1.56D0
      RAD(42)=1.46D0
      RAD(43)=1.38D0
      RAD(44)=1.36D0
      RAD(45)=1.34D0
      RAD(46)=1.30D0
      RAD(47)=1.36D0
      RAD(48)=1.40D0
      RAD(49)=1.42D0
      RAD(50)=1.40D0
      RAD(51)=1.40D0
      RAD(52)=1.37D0
      RAD(53)=1.36D0
      RAD(54)=1.36D0
      RAD(55)=2.38D0
      RAD(56)=2.06D0
      RAD(57)=1.94D0
      RAD(58)=1.84D0
      RAD(59)=1.90D0
      RAD(60)=1.88D0
      RAD(61)=1.86D0
      RAD(62)=1.85D0
      RAD(63)=1.83D0
      RAD(64)=1.82D0
      RAD(65)=1.81D0
      RAD(66)=1.80D0
      RAD(67)=1.79D0
      RAD(68)=1.77D0
      RAD(69)=1.77D0
      RAD(70)=1.78D0
      RAD(71)=1.74D0
      RAD(72)=1.64D0
      RAD(73)=1.58D0
      RAD(74)=1.50D0
      RAD(75)=1.41D0
      RAD(76)=1.36D0
      RAD(77)=1.32D0
      RAD(78)=1.30D0
      RAD(79)=1.30D0
      RAD(80)=1.32D0
      RAD(81)=1.44D0
      RAD(82)=1.45D0
      RAD(83)=1.50D0
      RAD(84)=1.42D0
      RAD(85)=1.48D0
      RAD(86)=1.46D0
      RAD(87)=2.42D0
      RAD(88)=2.11D0
      RAD(89)=2.01D0
      RAD(90)=1.90D0
      RAD(91)=1.84D0
      RAD(92)=1.83D0
      RAD(93)=1.80D0
      RAD(94)=1.80D0
      RAD(95)=1.73D0
      RAD(96)=1.68D0
      RAD(97)=1.68D0
      RAD(98)=1.68D0
      RAD(99)=1.65D0
      RAD(100)=1.67D0
      RAD(101)=1.73D0
      RAD(102)=1.76D0
      RAD(103)=1.61D0
      RAD(104)=1.57D0
      RAD(105)=1.49D0
      RAD(106)=1.43D0
      RAD(107)=1.41D0
      RAD(108)=1.34D0
      RAD(109)=1.29D0
      RAD(110)=1.28D0
      RAD(111)=1.21D0
      RAD(112)=1.22D0
      RAD(113)=1.36D0
      RAD(114)=1.43D0
      RAD(115)=1.62D0
      RAD(116)=1.75D0
      RAD(117)=1.65D0
      RAD(118)=1.57D0
C
C CM5 MODEL PARAMETERS
C
      DO I=1,MZ
      A0(I)=0.D0
      DO J=1,MZ
      D(I,J)=0.D0
      ENDDO
      ENDDO
C ATOMWISE PARAMETERS
      A0(  1)= 0.0056
      A0(  2)=-0.1543
      A0(  3)= 0.0000
      A0(  4)= 0.0333
      A0(  5)=-0.1030
      A0(  6)=-0.0446
      A0(  7)=-0.1072
      A0(  8)=-0.0802
      A0(  9)=-0.0629
      A0( 10)=-0.1088
      A0( 11)= 0.0184
      A0( 12)= 0.0000
      A0( 13)=-0.0726
      A0( 14)=-0.0790
      A0( 15)=-0.0756
      A0( 16)=-0.0565
      A0( 17)=-0.0444
      A0( 18)=-0.0767
      A0( 19)= 0.0130
      A0( 20)= 0.0000
      A0( 21)= 0.0000
      A0( 22)= 0.0000
      A0( 23)= 0.0000
      A0( 24)= 0.0000
      A0( 25)= 0.0000
      A0( 26)= 0.0000
      A0( 27)= 0.0000
      A0( 28)= 0.0000
      A0( 29)= 0.0000
      A0( 30)= 0.0000
      A0( 31)=-0.0512
      A0( 32)=-0.0557
      A0( 33)=-0.0533
      A0( 34)=-0.0399
      A0( 35)=-0.0313
      A0( 36)=-0.0541
      A0( 37)= 0.0092
      A0( 38)= 0.0000
      A0( 39)= 0.0000
      A0( 40)= 0.0000
      A0( 41)= 0.0000
      A0( 42)= 0.0000
      A0( 43)= 0.0000
      A0( 44)= 0.0000
      A0( 45)= 0.0000
      A0( 46)= 0.0000
      A0( 47)= 0.0000
      A0( 48)= 0.0000
      A0( 49)=-0.0361
      A0( 50)=-0.0393
      A0( 51)=-0.0376
      A0( 52)=-0.0281
      A0( 53)=-0.0220
      A0( 54)=-0.0381
      A0( 55)= 0.0065
      A0( 56)= 0.0000
      A0( 57)= 0.0000
      A0( 58)= 0.0000
      A0( 59)= 0.0000
      A0( 60)= 0.0000
      A0( 61)= 0.0000
      A0( 62)= 0.0000
      A0( 63)= 0.0000
      A0( 64)= 0.0000
      A0( 65)= 0.0000
      A0( 66)= 0.0000
      A0( 67)= 0.0000
      A0( 68)= 0.0000
      A0( 69)= 0.0000
      A0( 70)= 0.0000
      A0( 71)= 0.0000
      A0( 72)= 0.0000
      A0( 73)= 0.0000
      A0( 74)= 0.0000
      A0( 75)= 0.0000
      A0( 76)= 0.0000
      A0( 77)= 0.0000
      A0( 78)= 0.0000
      A0( 79)= 0.0000
      A0( 80)= 0.0000
      A0( 81)=-0.0255
      A0( 82)=-0.0277
      A0( 83)=-0.0265
      A0( 84)=-0.0198
      A0( 85)=-0.0155
      A0( 86)=-0.0269
      A0( 87)= 0.0046
      A0( 88)= 0.0000
      A0( 89)= 0.0000
      A0( 90)= 0.0000
      A0( 91)= 0.0000
      A0( 92)= 0.0000
      A0( 93)= 0.0000
      A0( 94)= 0.0000
      A0( 95)= 0.0000
      A0( 96)= 0.0000
      A0( 97)= 0.0000
      A0( 98)= 0.0000
      A0( 99)= 0.0000
      A0(100)= 0.0000
      A0(101)= 0.0000
      A0(102)= 0.0000
      A0(103)= 0.0000
      A0(104)= 0.0000
      A0(105)= 0.0000
      A0(106)= 0.0000
      A0(107)= 0.0000
      A0(108)= 0.0000
      A0(109)= 0.0000
      A0(110)= 0.0000
      A0(111)= 0.0000
      A0(112)= 0.0000
      A0(113)=-0.0179
      A0(114)=-0.0195
      A0(115)=-0.0187
      A0(116)=-0.0140
      A0(117)=-0.0110
      A0(118)=-0.0189
C
      DO K1=1,MZ
      DO K2=K1+1,MZ
      D(K1,K2)=A0(K1)-A0(K2)
      ENDDO
      ENDDO
C PAIRWISE PARAMETERS
      D( 1, 6)= 0.0502                   
      D( 1, 7)= 0.1747              
      D( 1, 8)= 0.1671             
      D( 6, 7)= 0.0556             
      D( 6, 8)= 0.0234             
      D( 7, 8)=-0.0346              
C
      DO I=1,MZ
      DO J=I+1,MZ
      D(J,I)=-D(I,J)
      ENDDO
      ENDDO
C ALPHA
      ALP=2.4740
C LAMDA in RCM5
      ZLAMDA(  1)= 0.0000
      ZLAMDA(  2)= 0.0000
      ZLAMDA(  3)= 0.0000
      ZLAMDA(  4)= 0.0000
      ZLAMDA(  5)= 0.0000
      ZLAMDA(  6)= 0.0000
      ZLAMDA(  7)= 0.0000
      ZLAMDA(  8)= 0.0000
      ZLAMDA(  9)= 0.0000
      ZLAMDA( 10)= 0.0000
      ZLAMDA( 11)= 0.0000
      ZLAMDA( 12)= 0.0000
      ZLAMDA( 13)= 0.0000
      ZLAMDA( 14)= 0.0000
      ZLAMDA( 15)= 0.0000
      ZLAMDA( 16)= 0.0000
      ZLAMDA( 17)= 0.0000
      ZLAMDA( 18)= 0.0000
      ZLAMDA( 19)= 0.0000
      ZLAMDA( 20)= 0.0000
      ZLAMDA( 21)= 0.0000
      ZLAMDA( 22)= 0.0000
      ZLAMDA( 23)= 0.0000
      ZLAMDA( 24)= 0.0000
      ZLAMDA( 25)= 0.0000
      ZLAMDA( 26)= 0.0000
      ZLAMDA( 27)= 0.0000
      ZLAMDA( 28)= 0.0000
      ZLAMDA( 29)= 0.0000
      ZLAMDA( 30)= 0.0000
      ZLAMDA( 31)= 0.0000
      ZLAMDA( 32)= 0.0000
      ZLAMDA( 33)= 0.0000
      ZLAMDA( 34)= 0.0000
      ZLAMDA( 35)= 0.0000
      ZLAMDA( 36)= 0.0000
      ZLAMDA( 37)= 0.0000
      ZLAMDA( 38)= 0.0000
      ZLAMDA( 39)= 0.0000
      ZLAMDA( 40)= 0.0000
      ZLAMDA( 41)= 0.0000
      ZLAMDA( 42)= 0.0000
      ZLAMDA( 43)= 0.0000
      ZLAMDA( 44)= 0.0000
      ZLAMDA( 45)= 0.0000
      ZLAMDA( 46)= 0.0000
      ZLAMDA( 47)=-0.0800
      ZLAMDA( 48)= 0.0000
      ZLAMDA( 49)= 0.0000
      ZLAMDA( 50)= 0.0000
      ZLAMDA( 51)= 0.0000
      ZLAMDA( 52)= 0.0000
      ZLAMDA( 53)= 0.0000
      ZLAMDA( 54)= 0.0000
      ZLAMDA( 55)= 0.0000
      ZLAMDA( 56)= 0.0000
      ZLAMDA( 57)= 0.0000
      ZLAMDA( 58)= 0.0000
      ZLAMDA( 59)= 0.0000
      ZLAMDA( 60)= 0.0000
      ZLAMDA( 61)= 0.0000
      ZLAMDA( 62)= 0.0000
      ZLAMDA( 63)= 0.0000
      ZLAMDA( 64)= 0.0000
      ZLAMDA( 65)= 0.0000
      ZLAMDA( 66)= 0.0000
      ZLAMDA( 67)= 0.0000
      ZLAMDA( 68)= 0.0000
      ZLAMDA( 69)= 0.0000
      ZLAMDA( 70)= 0.0000
      ZLAMDA( 71)= 0.0000
      ZLAMDA( 72)= 0.0000
      ZLAMDA( 73)= 0.0000
      ZLAMDA( 74)= 0.0000
      ZLAMDA( 75)= 0.0000
      ZLAMDA( 76)= 0.0000
      ZLAMDA( 77)= 0.0000
      ZLAMDA( 78)= 0.0000
      ZLAMDA( 79)= 0.0000
      ZLAMDA( 80)= 0.0000
      ZLAMDA( 81)= 0.0000
      ZLAMDA( 82)= 0.0000
      ZLAMDA( 83)= 0.0000
      ZLAMDA( 84)= 0.0000
      ZLAMDA( 85)= 0.0000
      ZLAMDA( 86)= 0.0000
      ZLAMDA( 87)= 0.0000
      ZLAMDA( 88)= 0.0000
      ZLAMDA( 89)= 0.0000
      ZLAMDA( 90)= 0.0000
      ZLAMDA( 91)= 0.0000
      ZLAMDA( 92)= 0.0000
      ZLAMDA( 93)= 0.0000
      ZLAMDA( 94)= 0.0000
      ZLAMDA( 95)= 0.0000
      ZLAMDA( 96)= 0.0000
      ZLAMDA( 97)= 0.0000
      ZLAMDA( 98)= 0.0000
      ZLAMDA( 99)= 0.0000
      ZLAMDA(100)= 0.0000
      ZLAMDA(101)= 0.0000
      ZLAMDA(102)= 0.0000
      ZLAMDA(103)= 0.0000
      ZLAMDA(104)= 0.0000
      ZLAMDA(105)= 0.0000
      ZLAMDA(106)= 0.0000
      ZLAMDA(107)= 0.0000
      ZLAMDA(108)= 0.0000
      ZLAMDA(109)= 0.0000
      ZLAMDA(110)= 0.0000
      ZLAMDA(111)= 0.0000
      ZLAMDA(112)= 0.0000
      ZLAMDA(113)= 0.0000
      ZLAMDA(114)= 0.0000
      ZLAMDA(115)= 0.0000
      ZLAMDA(116)= 0.0000
      ZLAMDA(117)= 0.0000
      ZLAMDA(118)= 0.0000
 
C C-COEFFICIENT: 0.7050   ! ALREADY INCLUDED IN A0
C
      A=0.36
C
      DO K=1,NAT
       CN(K)=0.0
       DO K1=1,NAT
        IF ((IZ(K).EQ.IZ(K1)).AND.(K.NE.K1)) THEN
         DIS=DSQRT((Q(1,K)-Q(1,K1))**2+(Q(2,K)-Q(2,K1))**2+
     $   (Q(3,K)-Q(3,K1))**2)
         BKK=DEXP(-ALP*(DIS-RAD(IZ(K))-RAD(IZ(K1))))
         CN(K)=CN(K)+BKK
        ENDIF
       ENDDO
      ENDDO
C
       DO K=1,NAT   
       CCM5(K)=CHIR(K)   
        DO K1=1,NAT   
         IF (IZ(K).NE.IZ(K1)) THEN
         DIS=DSQRT((Q(1,K)-Q(1,K1))**2+(Q(2,K)-Q(2,K1))**2+
     $   (Q(3,K)-Q(3,K1))**2)
          BKK=DEXP(-ALP*(DIS-RAD(IZ(K))-RAD(IZ(K1))))
          CCM5(K)=CCM5(K)+BKK*D(IZ(K),IZ(K1))
         ENDIF
        ENDDO
       ENDDO
C      
      DO K=1,NAT
      CRCM5(K)=CCM5(K)
       DO K1=1,NAT
        IF (IZ(K).EQ.IZ(K1)) THEN
         DIS=DSQRT((Q(1,K)-Q(1,K1))**2+(Q(2,K)-Q(2,K1))**2+
     $   (Q(3,K)-Q(3,K1))**2)
         BKK=DEXP(-ALP*(DIS-RAD(IZ(K))-RAD(IZ(K1))))
         IF ((CN(K)+CN(K1)).NE.0) THEN
         CRCM5(K)=CRCM5(K)+BKK*ZLAMDA(IZ(K))*
     $   (TANH(A*CN(K))-TANH(A*CN(K1)))
         ENDIF
        ENDIF
       ENDDO
      ENDDO

      DHIRX=0.D0
      DHIRY=0.D0
      DHIRZ=0.D0
      DCM5X=0.D0
      DCM5Y=0.D0
      DCM5Z=0.D0
      DRCM5X=0.D0
      DRCM5Y=0.D0
      DRCM5Z=0.D0
       DO J=1,NAT   
       DHIRX=DHIRX+Q(1,J)*CHIR(J)*4.803242D0
       DHIRY=DHIRY+Q(2,J)*CHIR(J)*4.803242D0
       DHIRZ=DHIRZ+Q(3,J)*CHIR(J)*4.803242D0
       DCM5X=DCM5X+Q(1,J)*CCM5(J)*4.803242D0
       DCM5Y=DCM5Y+Q(2,J)*CCM5(J)*4.803242D0
       DCM5Z=DCM5Z+Q(3,J)*CCM5(J)*4.803242D0
       DRCM5X=DRCM5X+Q(1,J)*CRCM5(J)*4.803242D0
       DRCM5Y=DRCM5Y+Q(2,J)*CRCM5(J)*4.803242D0
       DRCM5Z=DRCM5Z+Q(3,J)*CRCM5(J)*4.803242D0
       ENDDO
      RETURN
      END
C
      SUBROUTINE CM5PS(SWITCH)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (MAT=30000)
      DIMENSION Q(3,MAT),CELL(3,3),IZ(MAT),CHIR(MAT),
     $ CCM5(MAT),CRCM5(MAT),CN(MAT),QNEW(3,MAT),CHIRNEW(MAT)
      CHARACTER W1*1,W10*10,W50*50,W256*256,TYPE*1,SWITCH*10
      READ (*,'(A)') W256
      !WRITE (*,'(A)') W256
      READ (*,'(A)') W256
      !WRITE (*,'(A)') W256
      READ (W256,*) SCAL
      READ (*,'(A)') W256
      !WRITE (*,'(A)') W256
      READ (W256,*) (CELL(J,1),J=1,3)
      READ (*,'(A)') W256
      !WRITE (*,'(A)') W256
      READ (W256,*) (CELL(J,2),J=1,3)
      READ (*,'(A)') W256
      !WRITE (*,'(A)') W256
      READ (W256,*) (CELL(J,3),J=1,3)  
      DO J=1,3
        CELL(J,1)=CELL(J,1)*SCAL
        CELL(J,2)=CELL(J,2)*SCAL
        CELL(J,3)=CELL(J,3)*SCAL
      ENDDO
102   CONTINUE
      READ (*,'(A)') W256
      !WRITE (*,'(A)') W256
      READ (W256,*) W1 
      IF ((W1.NE.'D').AND.(W1.NE.'C').AND.(W1.NE.'d').AND.(W1.NE.'c'))
     $ GOTO 102
      TYPE=W1
      
      DO K=1,MAT
        READ (*,'(A)') W256
        !WRITE (*,'(A)') W256
        W1=W256 
        IF (W1.EQ.'-') GOTO 103
        READ (W256,*) IZ(K),(Q(J,K),J=1,3)
      ENDDO
103   CONTINUE
      NAT=K-1
      
      DO K=1,NAT
      READ (*,'(A)') W256
      !WRITE (*,'(A)') W256
      READ (W256, *) CHIR(K)
      ENDDO       
C
      KNEW=1
      IF (TYPE=='C'.OR.TYPE=='c') THEN
      DO N1=1,5
       DO N2=1,5
        DO N3=1,5
         DO K=1,NAT
          DO J=1,3
            QNEW(J,KNEW)=Q(J,K)+(N1-3)*CELL(J,1)+
     $      (N2-3)*CELL(J,2)+(N3-3)*CELL(J,3)
C            WRITE (*,*) Q(J,K),QNEW(J,KNEW)
          ENDDO
         CHIRNEW(KNEW)=CHIR(K)
         IZ(KNEW)=IZ(K)
         KNEW=KNEW+1
         ENDDO
        ENDDO
       ENDDO  
      ENDDO
      ELSE IF (TYPE=='D'.OR.TYPE=='d') THEN
      DO N1=1,5
       DO N2=1,5
        DO N3=1,5
         DO K=1,NAT
            QNEW(1,KNEW)=Q(1,K)*CELL(1,1)+Q(2,K)*CELL(1,2)+
     $                   Q(3,K)*CELL(1,3)+(N1-3)*CELL(1,1)+
     $                   (N2-3)*CELL(1,2)+(N3-3)*CELL(1,3)
            QNEW(2,KNEW)=Q(1,K)*CELL(2,1)+Q(2,K)*CELL(2,2)+
     $                   Q(3,K)*CELL(2,3)+(N1-3)*CELL(2,1)+
     $                   (N2-3)*CELL(2,2)+(N3-3)*CELL(2,3)
            QNEW(3,KNEW)=Q(1,K)*CELL(3,1)+Q(2,K)*CELL(3,2)+
     $                   Q(3,K)*CELL(3,3)+(N1-3)*CELL(3,1)+
     $                   (N2-3)*CELL(3,2)+(N3-3)*CELL(3,3)
            CHIRNEW(KNEW)=CHIR(K)
            IZ(KNEW)=IZ(K)
            KNEW=KNEW+1
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      ELSE 
      WRITE (*,*) 'ERROR' 
      ENDIF
      NATNEW=KNEW-1
      CALL CM5MOD(NATNEW,IZ,CHIRNEW,QNEW,CCM5,CRCM5,DHIRX,DHIRY,DHIRZ,
     $ DCM5X,DCM5Y,DCM5Z,DRCM5X,DRCM5Y,DRCM5Z,CN)
C
      IF (SWITCH=='3') THEN
      WRITE (*,'(A,/,A)')
     $ ' Charges (in A.U.) from CM5PAC version 2015',
     $ ' -----------------------------------------------'
      WRITE (*,'(A,/,A)')
     $ ' Center     Atomic      CM5         Hirshfeld',
     $ ' Number     Number      Charge      Charge'
      WRITE (*,'(A)')
     $ ' -----------------------------------------------'
        DO I=62*NAT+1,63*NAT
        WRITE (*,'(I5,6X,I5,5X,F11.6,X,F11.6)')
     $  I-62*NAT,IZ(I),CCM5(I),CHIRNEW(I)
        ENDDO
      WRITE (*,'(A)')
     $ ' -----------------------------------------------'
C
      ELSE
      WRITE (*,'(/,A,/,A)')
     $ ' Charges (in A.U.) from CM5PAC version 2015',
     $ ' -----------------------------------------------'
      WRITE (*,'(A,/,A)')
     $ ' Center     Atomic      CM5M        Hirshfeld',
     $ ' Number     Number      Charge      Charge'
      WRITE (*,'(A)')
     $ ' -----------------------------------------------'
        DO I=62*NAT+1,63*NAT 
        WRITE (*,'(I5,6X,I5,5X,F11.6,X,F11.6)')
     $  I-62*NAT,IZ(I),CRCM5(I),CHIRNEW(I)
        ENDDO
      WRITE (*,'(A)')
     $ ' -----------------------------------------------'
      ENDIF
      END


