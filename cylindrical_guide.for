		SUBROUTINE cylindrical_guide(lambda,x,y,z,ay,by,az,bz,phiy,phiz,
     &in,through,L,nrefl,xloss,isection,UCN,intube,iseed,ipol)
		INCLUDE 'parameters.inc'
		INCLUDE 'common.inc'
		INCLUDE 'c_guide.inc'
		INCLUDE 'c_plot.inc'

		DIMENSION xloss(2,2)
		REAL*8	L, ki, kix, kiy, kiz, kf, kfx, kfy, kfz
		LOGICAL	through, UCN, intube, in
		CHARACTER*1 Answer

		IF (UCN) L=-1.
		xloss(1,1)=-1.
		xloss(2,1)=1.
		xloss(1,2)=-1.
		xloss(2,2)=1.

		IF (show) WRITE(*,10) isection, x, y, z, 
     &phiy*180./pi, phiz*180./pi
   10	FORMAT(' Section',I2,' : starting pt (x,y,z) = (',F8.2,',',
     &F7.2,',',F7.2,')  phiy=',F6.2,' deg   phiz=',F6.2,' deg')

		ag0=ag(isection,1)		! straight line of cylinder axis
		bg0=bg(isection,1)		! is along y=ag0*x+bg0
		R=ag(isection,2)		! radius of cylindrical guide
		ki=2.*pi/lambda
		ag02=ag0*ag0

		xstart0=xc1(isection)
		ystart0=yc1(isection)
		zstart0=0.
		xend0=xc2(isection)
		yend0=yc2(isection)
		zend0=0.

		IF (show) WRITE(*,20) xstart0, ystart0, zstart0, 
     &xend0, yend0, zend0
   20	FORMAT(' Cylinder axis : (',F8.2,',',F7.2,',',F7.2,')  =>  (',
     &F8.2,',',F7.2,',',F7.2,')')

		c=cw(isection,1)
		d=dw(isection,1)
		xentrance=(by*c+d)/(1.-ay*c)
		yentrance=ay*xentrance+by
		zentrance=az*xentrance+bz
		dr=SQRT((xentrance-xstart0)**2+(yentrance-ystart0)**2
     &+(zentrance-zstart0)**2)

		IF (show) WRITE(*,30) xentrance, yentrance, zentrance
   30	FORMAT(' Entrance point : (',F8.2,',',F7.2,',',F7.2,')')

		IF (dr.GE.R) THEN
			IF (show) WRITE(*,50) yentrance, zentrance
   50		FORMAT(' At guide entrance: (y,z) = (',F6.2,',',F7.2,')'/
     &' == neutron missed entrance window ==')
			IF (pn) THEN
				CALL pgpanl(1,1)
				CALL pgswin(xmin0,xmax0,ymin0,ymax0)
				CALL pgarro(REAL(x),REAL(y),REAL(xentrance),REAL(yentrance))
				CALL pgpanl(1,2)
				CALL pgswin(xmin0,xmax0,zmin0,zmax0)
				CALL pgarro(REAL(x),REAL(z),REAL(xentrance),REAL(zentrance))
			END IF
			DO ipol=1,2
				xloss(1,ipol)=xentrance
				xloss(2,ipol)=0.
			END DO
			through=.FALSE.
			intube=.FALSE.
			in=.FALSE.
			GOTO 1000
		END IF
		in=.TRUE.

  100	CONTINUE

		c=cw(isection,2)
		d=dw(isection,2)
		xexit=(by*c+d)/(1.-ay*c)
		yexit=ay*xexit+by
		zexit=az*xexit+bz

		IF (show) WRITE(*,110) xexit, yexit, zexit
  110	FORMAT(' Heading towards exit point : (',
     &F8.2,',',F7.2,',',F7.2,')')

		c1=1.-(1.+ag0*ay)/(1.+ag02)
		c2=(ag0*bg0-ag0*by)/(1.+ag02)
		c3=ay-(ag0+ag02*ay)/(1.+ag02)
		c4=by-bg0+(ag02*bg0-ag02*by)/(1.+ag02)
		c5=az
		c6=bz
		A=c1*c1+c3*c3+c5*c5
		B=2.*(c1*c2+c3*c4+c5*c6)
		C=c2*c2+c4*c4+c6*c6-R*R
		Det=B*B-4.*A*C
		IF (Det.LE.0.) THEN
			PRINT*,'Warning: cylindrical guide - Det < 0'
			GOTO 1000
		ELSE
			xcross=(-B+SQRT(Det))/(2.*A)
		END IF
	
		IF (xcross.GT.xexit) THEN	! no more reflections
			IF (show) PRINT*,
     &'No further reflections from section',isection
			IF (UCN) THEN
				IF (L.LT.0.) THEN
					L=SQRT((xexit-xentrance)**2+(yexit-yentrance)**2
     &+(zexit-zentrance)**2)
				ELSE
					L=L+SQRT((xexit-x)**2+(yexit-y)**2+(zexit-z)**2)
				END IF
			END IF
		ELSE
			xold=x
			yold=y
			zold=z
			x=xcross
			y=ay*x+by
			z=az*x+bz
			IF (show) WRITE(*,150) isection, x, y, z
  150		FORMAT(' Neutron is reflected from guide at section',I2,
     &' at (x,y,z)=(',F8.2,',',F7.2,',',F7.2,')')
			IF (pn.AND.through) THEN
				CALL pgpanl(1,1)
				CALL pgswin(xmin0,xmax0,ymin0,ymax0)
				CALL pgarro(REAL(xold),REAL(yold),REAL(x),REAL(y))
				CALL pgpanl(1,2)
				CALL pgswin(xmin0,xmax0,zmin0,zmax0)
				CALL pgarro(REAL(xold),REAL(zold),REAL(x),REAL(z))
			END IF

			IF (UCN) THEN
				IF (L.LT.0.) THEN
					L=SQRT((x-xentrance)**2+(y-yentrance)**2
     &+(z-zentrance)**2)
				ELSE
					L=L+SQRT((x-xold)**2+(y-yold)**2+(z-zold)**2)
				END IF
			ELSE
				L=L+SQRT((x-xold)**2+(y-yold)**2+(z-zold)**2)
			END IF

			IF (show) WRITE(*,210) 
     &SQRT((x-xold)**2+(y-yold)**2+(z-zold)**2)
  210		FORMAT(' Flight path from last reflection =',F8.2,' mm')

			x0=(x+ag0*y-ag0*bg0)/(1.+ag02)
			y0=ag0*x0+bg0
			z0=0.
			kix=x-xold
			kiy=y-yold
			kiz=z-zold
			rki=SQRT(kix*kix+kiy*kiy+kiz*kiz)
			kix=ki*kix/rki
			kiy=ki*kiy/rki
			kiz=ki*kiz/rki
			Qx=x-x0
			Qy=y-y0
			Qz=z-z0
			IF (show) PRINT*,'R =',REAL(R),
     &'   Rcalc =',REAL(SQRT(Qx**2+Qy**2+Qz**2))
			Qx=Qx/R
			Qy=Qy/R
			Qz=Qz/R
			theta=pi/2.-ABS(ACOS((kix*Qx+kiy*Qy+kiz*Qz)/ki))
			IF (show) PRINT*,' theta =',REAL(theta*180./pi),' degrees'

			IF (Wwav(isection,1).GT.0.) THEN
				sigma=Wwav(isection,1)/2.355
  215			dphig=10.*sigma*(random(iseed)-1.)
				IF (EXP(-0.5*(dphig/sigma)**2).LT.random(iseed)) GOTO 215
			ELSE
				dphig=0.
			END IF
			theta=theta+dphig
			IF (show) WRITE(*,220) dphig*180./pi, theta*180./pi
 220		FORMAT(' Waviness angle =',F6.2,'  => theta=',F6.2,'deg')

			Prefl=refl(isection,1,ABS(theta),lambda,ipol)
  235		ran0=random(iseed)
			IF (REAL(ran0).EQ.0.) GOTO 235

			IF (ran0.GT.Prefl) THEN
				IF (through) xloss(1,ipol)=x
				through=.FALSE.
			END IF

			IF (show) WRITE(*,240) theta*180./pi, Prefl, ran0
  240		FORMAT(' theta =',F5.2,' deg,  Prefl =',F6.3,
     &',  random(iseed) =',F6.3)
			IF (.NOT.through) THEN
				IF (show) PRINT*,' == neutron absorbed =='
				GOTO 1000
			END IF

			Qmod=2.*ki*SIN(theta)
			Qx=Qx*Qmod
			Qy=Qy*Qmod
			Qz=Qz*Qmod
			kfx=kix-Qx
			kfy=kiy-Qy
			kfz=kiz-Qz
			IF (show) PRINT*,' ki =',REAL(ki)
			IF (show) PRINT*,' kf =',REAL(SQRT(kfx**2+kfy**2+kfz**2))

			ay=kfy/kfx
			by=y-ay*x
			phiy=ATAN(ay)
			az=kfz/kfx
			bz=z-az*x
			phiz=ATAN(az)

			IF (show) WRITE(*,250) phiy*180./pi, phiz*180./pi
  250		FORMAT(' Neutron reflected. New angles: phiy =',F6.2,
     &' deg   phiz =',F6.2,' deg')
			nrefl=nrefl+1

			IF (pn.AND.show) THEN
				WRITE(*,255)
  255			FORMAT(' Hit return to continue : ',$)
				READ(*,'(A1)') Answer
			END IF

			GOTO 100
		END IF

 1000	RETURN
		END

