		REAL*8 FUNCTION refl(isection,iside,theta,lambda,ipol)

c	calculates the supermirror guide reflectivity. 

		INCLUDE 'parameters.inc'
		INCLUDE 'common.inc'
		INCLUDE 'c_guide.inc'

		REAL*8	m, ki
		LOGICAL polarising
		CHARACTER*1 Answer

		is=MIN(iside,5)
		m=mnumber(isection,is)
		polarising=.FALSE.
		IF (m.LT.0.) THEN
			polarising=.TRUE.
			m=ABS(m)
		END IF
		R0=refl0(isection,is)
		Rm=reflm(isection,is)
		theta0=0.1*lambda*pi/180. ! critical angle of Ni in radians
		thetam=m*theta0
		IF (m.LT.1.) theta0=thetam

		IF (show) THEN
			PRINT*,'refl: lambda =',REAL(lambda),
     &'  theta =',REAL(theta*180./pi)
			PRINT*,'refl: theta0 =',REAL(theta0*180./pi),
     &'  thetam =',REAL(thetam*180./pi)
			PRINT*,'refl: polarising=',polarising,' ipol=',ipol
			PRINT*,'refl: R0=',R0
		END IF

		IF (R0.GT.0.) THEN	! use model reflectivity curve
			IF (theta.LT.theta0) THEN
				P1=R0
			ELSE IF (theta.LT.thetam) THEN
				P1=R0+(theta-theta0)*(Rm-R0)/(thetam-theta0)
			ELSE
				P1=0.
			END IF
			IF (ipol.EQ.1) THEN
				Prefl=P1
			ELSE
				IF (polarising) THEN
					Prefl=0.
				ELSE
					Prefl=P1
				END IF
			END IF
		ELSE				! use reflectivity curve(s) from data file(s)
			ki=2.*pi/lambda
			phi=2.*ABS(theta)
			Q=ki*SQRT(2.*(1.-COS(phi)))
			IF (show) PRINT*,'refl: Q=',Q
			IF (show) PRINT*,'isection=',isection,' is=',is,' nQ=',
     &nQ(isection,is)
			DO i=1,nQ(isection,is)
				IF (show) PRINT*,'Q=',Q,' Qarr=',Qarr(isection,is,i),
     &' R=',Parr(isection,is,ipol,i)
				IF (Q.LT.Qarr(isection,is,i)) THEN
					Prefl=Parr(isection,is,ipol,i-1)
     &+(Q-Qarr(isection,is,i-1))
     &*(Parr(isection,is,ipol,i)-Parr(isection,is,ipol,i-1))
     &/(Qarr(isection,is,i)-Qarr(isection,is,i-1))
					GOTO 300
				END IF
				Prefl=0.	! if Q is outside range of file
			END DO
  300		CONTINUE
		END IF

		IF (show) WRITE(*,310) theta*180./pi, ipol, Prefl
  310	FORMAT(' refl: theta=',F5.3,'deg: ipol=',I1,' Prefl=',F5.3)

		refl=Prefl
		RETURN
		END
