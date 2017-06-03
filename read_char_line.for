		SUBROUTINE read_char_line(line,charr,ncharr)
		INCLUDE 'parameters.inc'
		INCLUDE 'common.inc'

		CHARACTER charr(marray)*40, line*(nline), bob*40
		DIMENSION i1(marray), i2(marray)
		LOGICAL space

		IF (show) PRINT*,'read_char_line: line =',line(1:60)
		space=.TRUE.
		ncharr=0
		DO i=1,nline
			IF (line(i:i).EQ.' ' .OR. line(i:i).EQ.'	') THEN
				IF (.NOT.space .AND. ncharr.GE.1) i2(ncharr)=i-1
				space=.TRUE.
			ELSE
				IF (space) THEN
					IF (ncharr.EQ.marray) GOTO 100
					ncharr=ncharr+1
					i1(ncharr)=i
				END IF
				space=.FALSE.
			END IF
		END DO

  100	CONTINUE
		IF (show) PRINT*,'read_char_line: ncharr =',ncharr
		DO i=1,ncharr
c			charr(i)=UPPER(line(i1(i):i2(i)))
			charr(i)=line(i1(i):i2(i))
		END DO
		IF (show) THEN
			DO i=1,ncharr
				PRINT*,'read_char_line: charr = ',charr(i)
			END DO
		END IF

		RETURN
		END
