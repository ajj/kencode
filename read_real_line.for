		SUBROUTINE read_real_line(line,array,narray)
		INCLUDE 'parameters.inc'
		INCLUDE 'common.inc'

		DIMENSION array(marray)
		CHARACTER line*(nline)
		LOGICAL space

		IF (show) PRINT*,'read_real_line: line =',line(1:60)
		space=.TRUE.
		narray=0
		DO i=1,nline
			IF (line(i:i).EQ.' ' .OR. line(i:i).EQ.'	') THEN
				space=.TRUE.
			ELSE
				IF (narray.EQ.marray) GOTO 50
				IF (space) narray=narray+1
				space=.FALSE.
			END IF
		END DO

   50	IF (show) PRINT*,'read_real_line: narray =',narray
		DO i=1,marray
			array(i)=0.
		END DO
		IF (narray.NE.0) READ(line,*,ERR=100,END=150) 
     &(array(i), i=1,narray)
		GOTO 200

  100	CONTINUE
c		PRINT*,'read_line: Error in reading from string - narray =',narray
c		PRINT*, line
		DO i=1,nline
			IF (line(i:i).NE.' ') GOTO 110
		END DO
		narray=0
		GOTO 200
  110	line=line(i:nline)
		DO i=1,nline
			IF (line(i:i).EQ.' ') GOTO 120
		END DO
  120	narray=-(i-1)
		GOTO 200

  150	CONTINUE
		narray=100

  200	CONTINUE
		IF (show) PRINT*,'read_real_line: narray =',narray

		RETURN
		END
