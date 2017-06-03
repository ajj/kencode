		SUBROUTINE bin(x,ivar,ipol)
		INCLUDE 'parameters.inc'
		INCLUDE 'common.inc'
		INCLUDE 'c_bin.inc'

		dx=var(ivar,1)-var(ivar,0)
		IF (x.LT.var(ivar,0)) THEN
			DO i=0,-nbin,-1
				IF (x.GT.var(ivar,i)-dx/2.) THEN
					binvar(ivar,ipol,i)=binvar(ivar,ipol,i)+1.
					GOTO 1000
				END IF
			END DO
c			show=.TRUE.
c			PRINT*,'--------------------------------------'
c			WRITE(*,100) ivar, x
  100		FORMAT(' -ve binning variable',I2,
     &' outside binning area: x =',G12.5)
		ELSE
			DO i=0,nbin
				IF (x.LT.var(ivar,i)+dx/2.) THEN
					binvar(ivar,ipol,i)=binvar(ivar,ipol,i)+1.
					GOTO 1000
				END IF
			END DO
c			show=.TRUE.
c			PRINT*,'--------------------------------------'
c			WRITE(*,200) ivar, x
  200		FORMAT(' +ve binning variable',I2,
     &' outside binning area: x =',G12.5)
		END IF
 1000	CONTINUE

		RETURN
		END

