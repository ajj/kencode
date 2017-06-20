		SUBROUTINE guide_section(lambda,v,x,y,z,ay,by,gy,az,bz,gz,
     &in,through,L,nrefl,xloss,isection,iseed,
     &curved,hparabolic,vparabolic,icavity,ipol)
		INCLUDE 'parameters.inc'
		INCLUDE 'common.inc'
		INCLUDE 'c_guide.inc'
		INCLUDE 'c_plot.inc'

		DIMENSION xloss(2,2), x0(2), y0(2), R0(2)
		REAL*8	L
		LOGICAL	in, curved, hparabolic, vparabolic, outside, through, UCN
		CHARACTER*1 Answer
		REAL*4	xp(101), yp(101), zp(101), random

		xold=x
		yold=y
		zold=z
		ilast=0

		UCN=.FALSE.
		IF (L.LT.0.) UCN=.TRUE.
		xloss(1,1)=-1.
		xloss(1,2)=-1.
		ilast=0

		icase=0
		IF (curved) icase=1
		IF (hparabolic) icase=2
		IF (vparabolic) icase=3

		IF (curved) THEN
			x0(1)=xc1(isection)
			x0(2)=xc2(isection)
			y0(1)=yc1(isection)
			y0(2)=yc2(isection)
			R0(1)=Rguide(isection,1)
			R0(2)=Rguide(isection,2)
		END IF

		IF (hparabolic .OR. vparabolic) THEN
			IF (hparabolic) SINphi=SIN(phiguideH(isection))
			IF (vparabolic) SINphi=SIN(phiguideV(isection))
			COSphi=SQRT(1.-SINphi*SINphi)
			TANphi=SINphi/COSphi
			alpha=apara(isection)
			beta=TANphi
			x00=xfocus(isection)
		END IF

		IF (show) THEN
			slopey=-gy*x/v**2+ay
			slopez=-gz*x/v**2+az
			WRITE(*,10) isection, x, y, z, 
     &ATAN(slopey)*180./pi, ATAN(slopez)*180./pi
   10		FORMAT(' Section',I2,' : starting pt (x,y,z) = (',F8.2,',',
     &F7.2,',',F7.2,')  phiy=',F7.3,' deg   phiz=',F7.3,' deg')
			PRINT*,'gy=',REAL(gy),' gz=',REAL(gz)
		END IF

		c=cw(isection,1)
		d=dw(isection,1)
		IF (gy.EQ.0.) THEN
			xentrance=(by*c+d)/(1.-ay*c)
			yentrance=ay*xentrance+by
		ELSE
			F=-gy/(2.*v*v)
			IF (c.EQ.0.) THEN
				xentrance=d
				yentrance=F*xentrance**2+ay*xentrance+by
			ELSE
				AA=F*c*c
				BB=2.*F*c*d+ay*c-1.
				CC=F*d*d+ay*d+by
				y1=(-BB-SQRT(BB*BB-4.*AA*CC))/(2.*AA)
				y2=(-BB+SQRT(BB*BB-4.*AA*CC))/(2.*AA)
				IF (show) PRINT*,'Two solutions for guide entrance: y =',y1, y2
				yentrance=y1	! choose solution with smallest ABS(y)
				IF (ABS(y2).LT.ABS(y1)) yentrance=y2
				xentrance=c*yentrance+d
			END IF
		END IF
		zentrance=az*xentrance+bz-0.5*gz*(xentrance/v)**2

		IF (show) THEN
			slopey=-gy*xentrance/v**2+ay
			slopez=-gz*xentrance/v**2+az
			WRITE(*,20) xentrance, yentrance, zentrance, 
     &ATAN(slopey)*180./pi, ATAN(slopez)*180./pi
   20		FORMAT(' At entrance of section: (x,y,z) = (',F8.2,',',
     &F9.4,',',F9.4,')  phiy=',F7.3,' deg   phiz=',F7.3,' deg')
			PRINT*,'zentrance=',zentrance
			PRINT*,'zstart3=',zstart(isection,3)
			PRINT*,'zstart4=',zstart(isection,4)
		END IF

		offsetmulti=0.
		IF (2*(nmultiple(isection)/2).EQ.nmultiple(isection)) 
     &offsetmulti=-0.5

		imulti=-10000
		SINphi=SIN(phiguideH(isection))
		COSphi=COS(phiguideH(isection))
		DO i=1,nmultiple(isection)
			rm=REAL(i-(nmultiple(isection)+1)/2)+offsetmulti
			ym=rm*(Wguide(isection,1)+gmultiple(isection))*COSphi
			IF (yentrance.GT.ystart(isection,1)+ym .AND. 
     &yentrance.LT.ystart(isection,2)+ym .AND. 
     &zentrance.GT.zstart(isection,3) .AND. 
     &zentrance.LT.zstart(isection,4)) THEN
	 			imulti=i
				IF (show) PRINT*,'ymin,max=',ystart(isection,1)+ym, 
     &ystart(isection,2)+ym
				GOTO 40
			END IF
		END DO

   40	CONTINUE
		IF (imulti.EQ.-10000) THEN
			IF (show) WRITE(*,50) yentrance, zentrance
   50		FORMAT(' At guide entrance: (y,z) = (',F6.2,',',F7.2,')'/
     &' == neutron missed entrance window ==')
			IF (pn) THEN
				CALL pgpanl(1,1)
				CALL pgswin(xmin0,xmax0,ymin0,ymax0)
				IF (gy.EQ.0.) THEN
					CALL pgarro(REAL(xold),REAL(yold),
     &REAL(xentrance),REAL(yentrance))
				ELSE
					dx=(xentrance-xold)/100.
					DO j=1,101
						x=xold+REAL(j-1)*dx
						y=ay*x+by-0.5*gy*(x/v)**2
						xp(j)=x
						yp(j)=y
					END DO
					CALL pgline(101,xp,yp)
				END IF
				CALL pgpt1(REAL(xold),REAL(yold),17)
				CALL pgpt1(REAL(xentrance),REAL(yentrance),17)
				CALL pgpanl(1,2)
				CALL pgswin(xmin0,xmax0,zmin0,zmax0)
				IF (gz.EQ.0.) THEN
					CALL pgarro(REAL(xold),REAL(zold),
     &REAL(xentrance),REAL(zentrance))
				ELSE
					dx=(xentrance-xold)/100.
					DO j=1,101
						x=xold+REAL(j-1)*dx
						z=az*x+bz-0.5*gz*(x/v)**2
						xp(j)=x
						zp(j)=z
					END DO
					CALL pgline(101,xp,zp)
				END IF
				CALL pgpt1(REAL(xold),REAL(zold),17)
				CALL pgpt1(REAL(xentrance),REAL(zentrance),17)
			END IF
			xloss(1,ipol)=xentrance
			xloss(2,ipol)=0.
			through=.FALSE.
			in=.FALSE.
			GOTO 1000
		END IF
		in=.TRUE.

		rm=REAL(imulti-(nmultiple(isection)+1)/2)+offsetmulti
		xm=-rm*(Wguide(isection,1)+gmultiple(isection))*SINphi
		ym=rm*(Wguide(isection,1)+gmultiple(isection))*COSphi

  100	CONTINUE

		c=cw(isection,2)
		d=dw(isection,2)
		IF (gy.EQ.0.) THEN
			xexit=(by*c+d)/(1.-ay*c)
		ELSE
			F=-gy/(2.*v*v)
			IF (c.EQ.0.) THEN
				xexit=d
				yexit=F*xexit**2+ay*xexit+by
			ELSE
				AA=F*c*c
				BB=2.*F*c*d+ay*c-1.
				CC=F*d*d+ay*d+by
				y1=(-BB-SQRT(BB*BB-4.*AA*CC))/(2.*AA)
				y2=(-BB+SQRT(BB*BB-4.*AA*CC))/(2.*AA)
				yexit=y1	! choose solution with smallest ABS(y)
				IF (ABS(y2).LT.ABS(y1)) yexit=y2
				xexit=c*yexit+d
			END IF
		END IF

		xmin=1.E+10
		imin=0
		DO i=1,4+icavity
			xcross=MAX(xold,xentrance)
			IF (i.LE.2 .OR. i.GE.5) THEN
				ic=icase
				bm=ag(isection,i)*xm+ym
				IF (icase.EQ.0 .OR. icase.EQ.3) THEN	! linear
					a=ag(isection,i)
					b=bg(isection,i)+bm
					c=0.
					d=0.
					ic=0
				ELSE IF (icase.EQ.1) THEN	! curved
					outside=.FALSE.
					IF (R0(i).GT.0. .AND. i.EQ.1) outside=.TRUE.
					IF (R0(i).LT.0. .AND. i.EQ.2) outside=.TRUE.
					a=x0(i)+xm
					b=y0(i)+ym
					c=R0(i)
					d=xexit
				ELSE IF (icase.EQ.2) THEN	! horizontally parabolic
					IF (i.EQ.1) THEN
						gamma=ystart(isection,1)+Wguide(isection,1)*COSphi/2.
     &-xstart(isection,1)*TANphi
					ELSE
						gamma=ystart(isection,2)-Wguide(isection,1)*COSphi/2.
     &-xstart(isection,2)*TANphi
					END IF
					a=alpha
					b=beta
					c=gamma
					d=x00
				ELSE
					PRINT*,'guide_section: Error - icase =',icase
				END IF
				CALL cross(ay,by,gy,v,a,b,c,d,ic,i,outside,xcross)
			ELSE
				if=i
				IF (icase.EQ.3) THEN	! vertically parabolic
					IF (i.EQ.3) THEN
						gamma=zstart(isection,3)+Hguide(isection,1)/2.
					ELSE
						gamma=zstart(isection,4)-Hguide(isection,1)/2.
					END IF
					a=alpha
					b=beta
					c=gamma
					d=x00
					ic=2
					if=i-2
				ELSE
					a=ag(isection,i)
					b=bg(isection,i)
					c=0.
					d=0.
					ic=0
				END IF
				CALL cross(az,bz,gz,v,a,b,c,d,ic,if,outside,xcross)
			END IF
c			IF (show) PRINT*,'i =',i,'   xcross =',REAL(xcross)
			IF (xcross.LT.xmin .AND. xcross.GT.xold) THEN
				IF (show .AND. i.EQ.ilast) PRINT*,'i=',i,' ilast=',ilast,
     &' xold=',xold,' xcross=',xcross
				IF (i.EQ.ilast) THEN
					IF (xcross-xold.LT.0.1) GOTO 150
					slope=-1.
					IF ((i.LE.2 .OR. i.GE.5) 
     &.AND. gy.NE.0.) slope=-gy*xcross/v**2+ay
					IF ((i.EQ.3 .OR. i.EQ.4)
     &.AND. gz.NE.0.) slope=-gz*xcross/v**2+az
					IF (slope.GT.0. .AND. xcross-xold.LT.1.) GOTO 150
				END IF
				xmin=xcross
				imin=i
  150			CONTINUE
			END IF
		END DO

		IF (show) PRINT*,'isection=',isection,' xexit=',REAL(xexit),
     &' xmin=',REAL(xmin)

		IF (xmin.GT.xexit) THEN	! no more reflections
			IF (show) PRINT*,
     &'No further reflections from section',isection
			IF (UCN) THEN
				yexit=ay*xexit+by-0.5*gy*(xexit/v)**2
				zexit=az*xexit+bz-0.5*gz*(xexit/v)**2
				IF (L.LT.0.) THEN
					L=SQRT((xexit-xentrance)**2+(yexit-yentrance)**2
     &+(zexit-zentrance)**2)
				ELSE
					L=L+SQRT((xexit-x)**2+(yexit-y)**2+(zexit-z)**2)
				END IF
			END IF
		ELSE
			x=xmin
			y=ay*x+by-0.5*gy*(x/v)**2
			z=az*x+bz-0.5*gz*(x/v)**2
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
			IF (show) WRITE(*,210) isection, x, y, z
  210		FORMAT(' Neutron hits guide at section',I2,
     &' at (x,y,z)=(',F8.2,',',F7.2,',',F7.2,')')
			IF (show) THEN
				IF (imin.EQ.1) PRINT*,'== right-hand wall'
				IF (imin.EQ.2) PRINT*,'== left-hand wall'
				IF (imin.EQ.3) PRINT*,'== bottom wall'
				IF (imin.EQ.4) PRINT*,'== top wall'
				IF (imin.EQ.5) PRINT*,'== transmission mirror 1'
				IF (imin.EQ.6) PRINT*,'== transmission mirror 2'
			END IF
			IF (pn.AND.through) THEN
				xnew=x
				ynew=y
				znew=z
				CALL pgpanl(1,1)
				CALL pgswin(xmin0,xmax0,ymin0,ymax0)
				IF (gy.EQ.0.) THEN
					CALL pgarro(REAL(xold),REAL(yold),
     &REAL(xnew),REAL(ynew))
					WRITE(*,211) ATAN((ynew-yold)/(xnew-xold))*180./pi
  211				FORMAT(' y angle from plot =',F7.3,' deg')
				ELSE
					dx=(xnew-xold)/100.
					DO j=1,101
						x=xold+REAL(j-1)*dx
						y=ay*x+by-0.5*gy*(x/v)**2
						xp(j)=x
						yp(j)=y
					END DO
					CALL pgline(101,xp,yp)
				END IF
				CALL pgpt1(REAL(xold),REAL(yold),17)
				CALL pgpt1(REAL(xnew),REAL(ynew),17)
				CALL pgpanl(1,2)
				CALL pgswin(xmin0,xmax0,zmin0,zmax0)
				IF (gz.EQ.0.) THEN
					CALL pgarro(REAL(xold),REAL(zold),
     &REAL(xnew),REAL(znew))
					WRITE(*,212) ATAN((znew-zold)/(xnew-xold))*180./pi
  212				FORMAT(' z angle from plot =',F7.3,' deg')
				ELSE
					dx=(xnew-xold)/100.
					DO j=1,101
						x=xold+REAL(j-1)*dx
						z=az*x+bz-0.5*gz*(x/v)**2
						xp(j)=x
						zp(j)=z
					END DO
					CALL pgline(101,xp,zp)
				END IF
				CALL pgpt1(REAL(xold),REAL(zold),17)
				CALL pgpt1(REAL(xnew),REAL(znew),17)
				x=xnew
			END IF
			IF (imin.LE.2 .OR. imin.GE.5) THEN
				slope=-gy*x/v**2+ay
			ELSE
				slope=-gz*x/v**2+az
			END IF
			phi0=ATAN(slope)
			phig0=phig(isection,imin)
			IF (curved .AND. imin.LE.2) THEN
				y=ay*x+by-0.5*gy*(x/v)**2
				phig0=-ATAN((x-x0(imin))/(y-y0(imin)))
			END IF
			IF (hparabolic) THEN
				IF (imin.EQ.1) THEN
					slope=1./(2.*alpha*SQRT(1./(4.*alpha*alpha)
     &+(x00-x)/alpha))+beta
					phig0=ATAN(slope)
				ELSE IF (imin.EQ.2) THEN
					slope=-1./(2.*alpha*SQRT(1./(4.*alpha*alpha)
     &+(x00-x)/alpha))+beta
					phig0=ATAN(slope)
				END IF
			END IF
			IF (vparabolic) THEN
				IF (imin.EQ.3) THEN
					slope=1./(2.*alpha*SQRT(1./(4.*alpha*alpha)
     &+(x00-x)/alpha))+beta
					phig0=ATAN(slope)
				ELSE IF (imin.EQ.4) THEN
					slope=-1./(2.*alpha*SQRT(1./(4.*alpha*alpha)
     &+(x00-x)/alpha))+beta
					phig0=ATAN(slope)
				END IF
			END IF
			IF (Wwav(isection,imin).GT.0.) THEN
				sigma=Wwav(isection,imin)/2.355
  215			dphig=10.*sigma*(random(iseed)-1.)
				IF (EXP(-0.5*(dphig/sigma)**2).LT.random(iseed)) GOTO 215
			ELSE
				dphig=0.
			END IF
			theta=phig0+dphig-phi0
			IF (show) WRITE(*,220) phi0*180./pi, phig0*180./pi, 
     &dphig*180./pi, theta*180./pi
  220		FORMAT(' Incident neutron angle =',F7.3,'  Guide angle =',F7.3,
     &'  Waviness angle =',F7.3,'  => theta=',F7.3,'deg')
			IF (show) WRITE(*,230) 
     &SQRT((x-xold)**2+(y-yold)**2+(z-zold)**2)
  230		FORMAT(' Flight path from last reflection =',F8.2,' mm')

			IF (show) PRINT*,'isection =',isection,'   iside =',imin

			Prefl=refl(isection,imin,ABS(theta),lambda,ipol)
  235		ran0=random(iseed)
c			PRINT*,'ran0=',ran0
c			IF (REAL(ran0).EQ.0.) GOTO 235

c check for chamfer effect
			IF ((imin.LE.2 .OR. imin.GE.5) 
     &.AND. chamfer(isection,imin).GT.0.) THEN
				a=ag(isection,3)
				b=bg(isection,3)
				z3=a*x+b
				a=ag(isection,4)
				b=bg(isection,4)
				z4=a*x+b
				IF (show) WRITE(*,237) z, z3, z4
  237			FORMAT(' Checking chamfer: z=',F8.3,
     &' z(bot)=',F8.3,' z(top)=',F8.3)
				IF (ABS(z-z3).LT.chamfer(isection,imin) .OR. 
     &ABS(z-z4).LT.chamfer(isection,imin)) ran0=2.
			ELSE IF ((imin.EQ.3 .OR. imin.EQ.4) 
     &.AND. chamfer(isection,imin).GT.0.) THEN
				IF (curved) THEN
					x01=x0(1)+xm
					y01=y0(1)+ym
					R01=R0(1)
					y1=y01-SIGN(DBLE(1.),y01)*SQRT(R01**2-(x-x01)**2)
					x02=x0(2)+xm
					y02=y0(2)+ym
					R02=R0(2)
					y2=y02-SIGN(DBLE(1.),y02)*SQRT(R02**2-(x-x02)**2)
				ELSE IF (hparabolic) THEN
					gamma1=ystart(isection,1)+Wguide(isection,1)*COSphi/2.
     &-xstart(isection,1)*TANphi
					gamma2=ystart(isection,2)-Wguide(isection,1)*COSphi/2.
     &-xstart(isection,2)*TANphi
					y1=-SQRT(1./(4.*alpha**2)+(x00-x)/alpha)+beta*x+gamma1
					y2= SQRT(1./(4.*alpha**2)+(x00-x)/alpha)+beta*x+gamma2
				ELSE
					a=ag(isection,1)
					b=bg(isection,1)
					y1=a*x+b
					a=ag(isection,2)
					b=bg(isection,2)
					y2=a*x+b
				END IF
				IF (show) WRITE(*,238) y, y1, y2
  238			FORMAT(' Checking chamfer: y=',F8.3,
     &' y(rhs)=',F8.3,' y(lhs)=',F8.3)
				IF (ABS(y-y1).LT.chamfer(isection,imin) .OR. 
     &ABS(y-y2).LT.chamfer(isection,imin)) ran0=2.
			END IF

			IF (show) WRITE(*,540) theta*180./pi, Prefl, ran0
  540		FORMAT(' theta =',F5.2,' deg,  Refl =',F6.3,
     &',  random(iseed) =',F6.3)

			IF (imin.GE.5) THEN	! transmission mirror
				DO i=2,19
					IF (lambda.LT.Si_lambda(i)) THEN
						Sigma=Si_mu(i-1)+(lambda-Si_lambda(i-1))
     &*(Si_mu(i)-Si_mu(i-1))/(Si_lambda(i)-Si_lambda(i-1))
						GOTO 250
					END IF
				END DO
				PRINT*,'Warning: lambda outside range of tabulated Sigma'
  250			CONTINUE
				Sigma=Sigma/10.	! convert from cm-1 to mm-1
				thick=cthickness(isection)/SIN(ABS(theta))
				T=EXP(-thick*Sigma)	! transmission through the Si substrate
				IF (show) WRITE(*,505) cthickness(isection), thick, Sigma, T
  505			FORMAT(' Si wafer thickness =',F6.3,
     &' -> =',F8.3,' mm  Sigma =',F9.6,' mm-1  T =',F8.5)
				R=Prefl
				IF (nsides(isection).EQ.1) THEN	! single-sided transmission mirror
					Pr=R
					Pt=T*(1-R)
					IF (show) WRITE(*,510) Pr, Pt
  510				FORMAT(' single-sided: Pr=',F6.4,' Pt=',F6.4)
				ELSE							! double-sided transmission mirror
					Pt=T*(1.-R)**2/(1.-(R*T)**2)
					Pr=R*(1.+T*Pt)
					IF (show) WRITE(*,520) Pr, Pt
  520				FORMAT(' double-sided: Pr=',F6.4,' Pt=',F6.4)
				END IF
				IF (ran0.LT.Pr) THEN					! neutron is reflected
					phi0=phi0+2.*theta
					slope=TAN(phi0)
					ay=slope+gy*x/v**2
					by=y+0.5*gy*(x/v)**2-ay*x
					phiy=phi0
					slopez=-gz*x/v**2+az
					phiz=ATAN(slopez)
					nrefl=nrefl+1
					IF (show) WRITE(*,550) phiy*180./pi, phiz*180./pi
  550				FORMAT(' Neutron reflected. New angles: phiy =',F7.3,
     &' deg   phiz =',F7.3,' deg')
				ELSE IF (random(iseed).LT.Pt/(1.-Pr)) THEN	! neutron is transmitted
					IF (show) WRITE(*,560)
  560				FORMAT(' Neutron transmitted')
				ELSE									! neutron is absorbed
					xloss(1,ipol)=x
					xloss(2,ipol)=DBLE(imin)
					through=.FALSE.
					IF (show) PRINT*,' == neutron absorbed =='
					GOTO 1000
				END IF
			ELSE
				IF (ran0.GT.Prefl) THEN
					xloss(1,ipol)=x
					xloss(2,ipol)=DBLE(imin)
					through=.FALSE.
					IF (show) PRINT*,' == neutron absorbed =='
					GOTO 1000
				END IF
				phi0=phi0+2.*theta
				slope=TAN(phi0)
				IF (imin.LE.2) THEN
					ay=slope+gy*x/v**2
					by=y+0.5*gy*(x/v)**2-ay*x
					phiy=phi0
					slopez=-gz*x/v**2+az
					phiz=ATAN(slopez)
				ELSE
					az=slope+gz*x/v**2
					bz=z+0.5*gz*(x/v)**2-az*x
					phiz=phi0
					slopey=-gy*x/v**2+ay
					phiy=ATAN(slopey)
				END IF
				nrefl=nrefl+1
				IF (show) WRITE(*,550) phiy*180./pi, phiz*180./pi
			END IF

			xold=x
			yold=y
			zold=z
			ilast=imin

			IF (pn.AND.show) THEN
				WRITE(*,655)
  655			FORMAT(' Hit return to continue : ',$)
				READ(*,'(A1)') Answer
			END IF

			GOTO 100
		END IF

 1000	RETURN
		END

