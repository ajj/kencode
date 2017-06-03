		SUBROUTINE makefilename(filename1,lambda,filename2,n)

		INCLUDE 'parameters.inc'

		CHARACTER*50 filename1*50, filename2*50, wstring*5, buf*50

		DO i=1,50
			IF (filename1(i:i).EQ.'.') i2=i-1
			buf(i:i)=' '
		END DO
		DO i=50,1,-1
			IF (filename1(i:i).NE.' ') i1=i
			buf(i:i)=filename1(i:i)
			filename1(i:i)=' '
			filename2(i:i)=' '
		END DO
		n=i2-i1+1
		filename1(1:51-i1)=buf(i1:50)

		IF (lambda.EQ.0.) THEN
			wstring='00000'
		ELSE IF (lambda.LT.0.1) THEN
			WRITE(wstring,21) NINT(100.*lambda)
   21		FORMAT('0000',I1)
		ELSE IF (lambda.LT.1.) THEN
			WRITE(wstring,22) NINT(100.*lambda)
   22		FORMAT('000',I2)
		ELSE IF (lambda.LT.10.) THEN
			WRITE(wstring,23) NINT(100.*lambda)
   23		FORMAT('00',I3)
		ELSE IF (lambda.LT.100.) THEN
			WRITE(wstring,24) NINT(100.*lambda)
   24		FORMAT('0',I4)
		ELSE
			WRITE(wstring,25) NINT(100.*lambda)
   25		FORMAT(I5)
		END IF

		filename2(1:n+6)=filename1(1:n)//'_'//wstring
		n=n+6

		RETURN
		END
