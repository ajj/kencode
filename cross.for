		SUBROUTINE cross(an,bn,g,v,a,b,c,d,icase,iface,outside,x)

		INCLUDE 'parameters.inc'
		INCLUDE 'common.inc'

		LOGICAL outside
		CHARACTER*1 Answer

c Find intersection of neutron trajectory with guide face in one dimension:
c
c icase=0 : guide face is a line     y=ag*x+bg
c icase=1 : guide face is a circle   (x-x0)^2+(y-y0)^2=r^2
c icase=2 : guide face is a parabola y=+/-SQRT(1/(4alpha^2)+(x0-x)/alpha)+beta*x+gamma
c
c Neutron trajectory is a parabola: y=-0.5g(x/v)^2+ax+b

		xin=x

		IF (icase.EQ.0) THEN	! straight line
			ag=a
			bg=b
			AA=-0.5*g/v**2
			BB=an-ag
			CC=bn-bg
			Det=BB*BB-4.*AA*CC
			IF (AA.EQ.0.) THEN	! no gravity
				IF (ag.EQ.an) THEN
					x=-1.
					GOTO 1000
				END IF
				x=-CC/BB
				y=an*x+bn
c				IF (show) WRITE(*,110) x, y
  110			FORMAT(' Cross(line): no gravity - (x,y)=(',
     &F9.2,',',F6.2,')')
			ELSE IF (Det.GT.0.) THEN
				x1=(-BB-SQRT(Det))/(2.*AA)
				x2=(-BB+SQRT(Det))/(2.*AA)
				y1=an*x1+bn+AA*x1*x1
				y2=an*x2+bn+AA*x2*x2
c				IF (show) WRITE(*,120) x1, y1, x2, y2
  120			FORMAT(' Cross(line): intersections are (x,y)=(',
     &F9.2,',',F7.2,') and (x,y)=(',F9.2,',',F7.2,')')
				IF (MIN(x1,x2).GT.x) THEN
					x=MIN(x1,x2)
				ELSE
					x=MAX(x1,x2)
				END IF
			ELSE
				x=-1.
c				IF (show) WRITE(*,130) an, bn, g*1000., v, ag, bg
  130			FORMAT(' Cross(line): Determinant is -ve:  an=',F9.5,
     &' bn=',F9.3/,'        g=',F5.2,' v=',F7.1,' ag=',
     &F9.5,' bg=',F9.3)
			END IF
			slopeg=ag
		ELSE IF (icase.EQ.1) THEN	! circle
			x0=a
			y0=b
			r=c
			AA=-0.5*g/v**2
			IF (AA.EQ.0.) THEN	! no gravity
				AA=1.+an*an
				BB=2.*(an*bn-x0-y0*an)
				CC=x0*x0-r*r+(bn-y0)**2
				Det=BB*BB-4.*AA*CC
				IF (Det.GT.0.) THEN
					IF (outside) THEN	! outside mirror
						x=(-BB+SQRT(Det))/(2.*AA)
					ELSE				! inside mirror
						x=(-BB-SQRT(Det))/(2.*AA)
					END IF
					y=an*x+bn
					IF (show) THEN
						PRINT*,'Cross(circle): no gravity - outside=',
     &outside
						x1=(-BB-SQRT(Det))/(2.*AA)
						x2=(-BB+SQRT(Det))/(2.*AA)
						y1=an*x1+bn
						y2=an*x2+bn
c						WRITE(*,310) x1, y1, x2, y2
  310					FORMAT(' Cross(circle): intersections are (x,y)=(',
     &F9.2,',',F7.2,') and (x,y)=(',F9.2,',',F7.2,')')
					END IF
				ELSE
					x=-1.
c					IF (show) WRITE(*,320) an, bn, x0, y0, r, outside
  320				FORMAT(' Cross(circle): Determinant is negative -  an=',
     &F9.5,' bn=',F9.3/'        (x0,y0)=(',F12.2,',',F12.2,
     &') r=',F12.2,' outside=',L1)
				END IF
				y=an*x+bn-0.5*g*(x/v)**2
				slopeg=-(x-x0)/(y-y0)
			ELSE	! with gravity
				IF (show) PRINT*,'Cross(circle): gravity - scan through x'
				xstart=x
				ntries=1
				xexit=d
  340			x=xstart
				yg=y0-SIGN(DBLE(1.),r)*SQRT(r*r-(x-x0)**2)
				yn=an*x+bn+AA*x*x
				test0=SIGN(DBLE(1.),yn-yg)
				IF (show) PRINT*,'x=',x,' yn-yg=',yn-yg
				nx=20
				dx=(xexit-xstart)/REAL(nx)
				GOTO 360
  350			xstart=x-dx
				dx=dx/10.
				IF (dx.LT.0.1) GOTO 370
c				IF (show) PRINT*,'dx=',dx
  360			DO ix=0,nx
					x=xstart+REAL(ix)*dx
					IF (ix.EQ.0.) x=xstart+0.1
					yg=y0-SIGN(DBLE(1.),r)*SQRT(r*r-(x-x0)**2)
					yn=an*x+bn+AA*x*x
c					IF (show) PRINT*,'ix=',ix,' x=',x,' yn-yg=',yn-yg
					test=SIGN(DBLE(1.),yn-yg)
					IF (test.NE.test0) GOTO 350
				END DO
				x=-1.	! no intersection
				GOTO 1000
  370			x2=x
				y2=yn-yg
				x1=x-dx*10.
				yg=y0-SIGN(DBLE(1.),r)*SQRT(r*r-(x1-x0)**2)
				yn=an*x1+bn+AA*x1*x1
				y1=yn-yg
				x=x1-y1*(x2-x1)/(y2-y1)
				IF (show) THEN
c					PRINT*,'x1=',x1,' yn-yg(1)=',y1
c					PRINT*,'x2=',x2,' yn-yg(1)=',y2
					yg=y0-SIGN(DBLE(1.),r)*SQRT(r*r-(x-x0)**2)
					yn=an*x+bn+AA*x*x
c					PRINT*,'x =',x,' yn-yg   =',yn-yg
				END IF
c				x1=xstart
c				IF (show) PRINT*,'Cross(circle):'//
c    &' Newton-Raphson iterations'
c				DO i=1,100
c					dx=dx/10.
c					x2=x1+dx
c					yg1=y0-SIGN(DBLE(1.),r)*SQRT(r*r-(x1-x0)**2)
c					yg2=y0-SIGN(DBLE(1.),r)*SQRT(r*r-(x2-x0)**2)
c					yn1=an*x1+bn+AA*x1*x1
c					yn2=an*x2+bn+AA*x2*x2
c					y1=yn1-yg1
c					y2=yn2-yg2
c					ydash=(y2-y1)/dx
c					IF (ydash.EQ.0.) GOTO 380
c					xnew=x1-y1/ydash
c					IF (show) PRINT*,'i=',i,' xold=',x1,' yn-yg(old)=',yn1-yg1
c					IF (show) PRINT*,'        xnew=',xnew
c					IF (show) READ(*,'(A1)') Answer
c					IF (ABS(xnew-x1).LT.1.E-09) GOTO 380
c					x1=xnew
c				END DO
c				PRINT*,'Cross(circle) - warning: '//
c    &'Newton-Raphson method did not converge'
c 380			x=xnew
c				IF (show) THEN
c					yn=an*x+bn+AA*x*x
c					yg=y0-SIGN(DBLE(1.),r)*SQRT(r*r-(x-x0)**2)
c					PRINT*,'N-R converged: x=',xnew,' yn-yg=',yn-yg
c					IF (show) READ(*,'(A1)') Answer
c				END IF
				y=an*x+bn-0.5*g*(x/v)**2
				slopeg=-(x-x0)/(y-y0)
				slopen=-g*x/v**2+an
				IF (ntries.EQ.2) GOTO 400
				IF (((iface.EQ.1 .OR. iface.EQ.3) .AND. (slopen.GT.slopeg)) 
     &.OR. ((iface.EQ.2 .OR. iface.EQ.4) .AND. (slopen.LT.slopeg))) THEN
					xstart=x+0.01
					ntries=2
					GOTO 340
				END IF
  400			CONTINUE
			END IF
c			IF (show) WRITE(*,390) x, y
  390		FORMAT(' Cross(circle): intersection at (x,y)=(',
     &F9.2,',',F7.2,')')
		ELSE IF (icase.EQ.2) THEN	! parabola
			alpha=a
			beta=b
			gamma=c
			x0=d
			AA=-0.5*g/v**2
			IF (AA.EQ.0.) THEN	! no gravity
				AA=(an-beta)**2
				BB=2.*(an-beta)*(bn-gamma)+1./alpha
				CC=(bn-gamma)**2-1./(4.*alpha*alpha)-x0/alpha
				Det=BB*BB-4.*AA*CC
				IF (show) PRINT*,'Cross(parabola): iface=',iface
				IF (AA.EQ.0.) THEN
					IF (iface.EQ.1) THEN
						IF (b.GT.gamma) THEN
							IF (show) PRINT*,'Cross(parabola): Only one '//
     &'possible soln: No intersection'
							x=-1.
							GOTO 1000
						END IF
					ELSE
						IF (b.LT.gamma) THEN
							IF (show) PRINT*,'Cross(parabola): Only one '//
     &'possible soln: No intersection'
							x=-1.
							GOTO 1000
						END IF
					END IF
					x=x0+1./(4.*alpha)-alpha*(b-gamma)**2
					y=an*x+bn
c					IF (show) WRITE(*,210) x, y
  210				FORMAT(' Cross(parabola): Only one solution (x,y)=(',
     &F9.2,',',F6.2,')')
				ELSE IF (Det.GT.0.) THEN
					x1=(-BB-SQRT(Det))/(2.*AA)
					x2=(-BB+SQRT(Det))/(2.*AA)
					y1=an*x1+bn
					y2=an*x2+bn
c					IF (show) WRITE(*,220) x1, y1, x2, y2
  220				FORMAT(' Cross(parabola): intersections are (x,y)=(',
     &F9.2,',',F7.2,') and (x,y)=(',F9.2,',',F7.2,')')
					argl1=(an-beta)*x1+(bn-gamma)
					argl2=(an-beta)*x2+(bn-gamma)
					argr1=1./(4.*alpha*alpha)+(x0-x1)/alpha
					argr2=1./(4.*alpha*alpha)+(x0-x2)/alpha
c					IF (show) PRINT*,' argl1 =',REAL(argl1),
c    &' argr1 =',REAL(argr1)
c					IF (show) PRINT*,' argl2 =',REAL(argl2),
c    &' argr2 =',REAL(argr2)
					IF (iface.EQ.1)	THEN	! choose soln which gives -ve argl and +ve argr
						IF (argl1.LT.0. .AND. argr1.GT.0.) x=x1
						IF (argl2.LT.0. .AND. argr2.GT.0.) x=x2
						IF (x.LT.0.) GOTO 1000
					ELSE					! choose soln which gives +ve argl and +ve argr
						IF (argl1.GT.0. .AND. argr1.GT.0.) x=x1
						IF (argl2.GT.0. .AND. argr2.GT.0.) x=x2
						IF (x.LT.0.) GOTO 1000
					END IF
					y=an*x+bn
c					IF (show) WRITE(*,230) x, y
  230				FORMAT(' Cross(parabola): chosen intersection at (x,y)=(',
     &F9.2,',',F7.2,')')
				ELSE
c					IF (show) WRITE(*,240) an, bn, x0, iface
  240				FORMAT(' Cross(parabola): Determinant is -ve:  an=',F9.5,
     &' bn=',F9.3/,'        x0=',F12.2,' iface=',I1)
				END IF
				slopeg=-REAL(2*iface-3)/(2.*alpha*SQRT(1./(4.*alpha*alpha)
     &+(x0-x)/alpha))+beta
			ELSE	! with gravity
c				IF (show) PRINT*,'Cross(parabola): gravity - scan through x'
				xstart=x
				ntries=1
				xexit=d
  440			x=xstart
				yg=REAL(2*iface-3)*SQRT(1./(4.*alpha**2)+(x0-x)/alpha)
     &+beta*x+gamma
				yn=an*x+bn+AA*x*x
				test0=SIGN(DBLE(1.),yn-yg)
				nx=20
				dx=(xexit-xstart)/REAL(nx)
				GOTO 460
  450			xstart=x-dx
				dx=dx/10.
				IF (dx.LT.0.1) GOTO 470
  460			DO ix=0,nx
					x=xstart+REAL(ix)*dx
					IF (ix.EQ.0.) x=xstart+0.1
					yg=REAL(2*iface-3)*SQRT(1./(4.*alpha**2)+(x0-x)/alpha)
     &+beta*x+gamma
					yn=an*x+bn+AA*x*x
					test=SIGN(DBLE(1.),yn-yg)
					IF (test.NE.test0) GOTO 450
				END DO
				x=-1.	! no intersection
				GOTO 1000
  470			x2=x
				y2=yn-yg
				x1=x-dx*10.
				yg=REAL(2*iface-3)*SQRT(1./(4.*alpha**2)+(x0-x1)/alpha)
     &+beta*x1+gamma
				yn=an*x1+bn+AA*x1*x1
				y1=yn-yg
				x=x1-y1*(x2-x1)/(y2-y1)
				IF (show) THEN
c					PRINT*,'x1=',x1,' yn-yg(1)=',y1
c					PRINT*,'x2=',x2,' yn-yg(1)=',y2
					yg=REAL(2*iface-3)*SQRT(1./(4.*alpha**2)+(x0-x)/alpha)
     &+beta*x+gamma
					yn=an*x+bn+AA*x*x
c					PRINT*,'x =',x,' yn-yg   =',yn-yg
c					READ(*,'(A1)') Answer
				END IF

c				x1=xstart
c				DO i=1,100
c					dx=dx/10.
c					x2=x1+dx
c					yg1=REAL(2*iface-3)*SQRT(1./(4.*alpha**2)+(x0-x1)/alpha)
c    &+beta*x1+gamma
c					yg2=REAL(2*iface-3)*SQRT(1./(4.*alpha**2)+(x0-x2)/alpha)
c    &+beta*x2+gamma
c					yn1=an*x1+bn+AA*x1*x1
c					yn2=an*x2+bn+AA*x2*x2
c					y1=yn1-yg1
c					y2=yn2-yg2
c					ydash=(y2-y1)/dx
c					IF (ydash.EQ.0.) GOTO 480
c					xnew=x1-y1/ydash
c					IF (ABS(xnew-x1).LT.1.E-09) GOTO 480
c					x1=xnew
c				END DO
c				PRINT*,'Cross(parabola) - warning: '//
c    &'Newton-Raphson method did not converge'
c 480			x=xnew
				slopeg=-REAL(2*iface-3)/(2.*alpha*SQRT(1./(4.*alpha*alpha)
     &+(x0-x)/alpha))+beta
				slopen=-g*x/v**2+an
				IF (ntries.EQ.2) GOTO 500
				IF (((iface.EQ.1 .OR. iface.EQ.3) .AND. (slopen.GT.slopeg)) 
     &.OR. ((iface.EQ.2 .OR. iface.EQ.4) .AND. (slopen.LT.slopeg))) THEN
					xstart=x+0.01
					ntries=2
					GOTO 440
				END IF
  500			CONTINUE
			END IF
		ELSE
			PRINT*,'Cross: Error - icase =',icase
		END IF

		slopen=-g*x/v**2+an
c		IF (show) PRINT*,'neutron slope =',slopen,' guide slope =',slopeg
		IF (iface.EQ.1 .OR. iface.EQ.3) THEN
			IF (slopen.GT.slopeg) x=-1.
		ELSE IF (iface.EQ.2 .OR. iface.EQ.4) THEN
			IF (slopen.LT.slopeg) x=-1.
		END IF

 1000	IF (x.LT.xin) x=-1.

c		IF (show) THEN
c			PRINT*,'Hit return to continue'
c			READ(*,'(A1)') Answer
c		END IF

		RETURN
		END

