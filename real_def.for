		FUNCTION real_def(r,ch,OK)
		PARAMETER(length=50)

		REAL	x(length)
		CHARACTER*(length) ch, string
		LOGICAL OK, found

		OK=.TRUE.
		found=.FALSE.
		factor=1.

		IF (ch.EQ.' ') THEN
			real_def=r
			GOTO 1000
		END IF

		j=0
		DO i=1,length-1
			iascii=ICHAR(ch(i:i))
			iascii1=ICHAR(ch(i+1:i+1))
			IF (j.EQ.0 .AND. iascii.EQ.45) THEN
				factor=-1.
			ELSE IF (j.NE.0 .AND. iascii.EQ.32 .AND. iascii1.NE.32) THEN
				OK=.FALSE.
c				PRINT*,'Error -have found non-blank after a blank'
 				GOTO 1000
			ELSE IF (iascii.NE.32) THEN
				j=j+1
				string(j:j)=ch(i:i)
			END IF
		END DO
		n=j

		ndigits=0
		DO i=1,n
			iascii=ICHAR(string(i:i))
			IF (iascii.EQ.46) THEN
				IF (found) THEN
					OK=.FALSE.
c					PRINT*,'Error -have found more than two decimal points'
					GOTO 1000
				END IF
				found=.TRUE.
				rmag=REAL(i-1)
			ELSE IF (iascii.LT.48 .OR. iascii.GT.57) THEN
				OK=.FALSE.
c				PRINT*,'Error -have found non-decimal character'
 				GOTO 1000
			ELSE
				ndigits=ndigits+1
				x(ndigits)=iascii-48
			END IF
		END DO
		IF (.NOT.found) rmag=REAL(n)

		rnum=0.
		DO i=1,ndigits
			rnum=rnum+x(i)*10.**(rmag-REAL(i))
		END DO
		real_def=rnum*factor

 1000	IF (.NOT.OK) WRITE(*,1100)
 1100	FORMAT(' Input must be real number')

		RETURN
		END

