		PROGRAM complicated_guide

c Monte Carlo calculation of guide divergence for arbitrary guide geometry. 
c
c										last modified 26/9/10

		INCLUDE 'parameters.inc'
		INCLUDE 'common.inc'
		INCLUDE 'c_guide.inc'
		INCLUDE 'c_bin.inc'
		INCLUDE 'c_plot.inc'

		REAL*8	Lguide(msections), Lexit, L, Lap, Lend, Lmon, rtime, 
     &nsigmaS, lambdamin, lambdamax, dlambda, lambda_arr(mlambda), 
     &lambda0
		DIMENSION nthrough(msections,2), nmon(2), ncentral(3,2), 
     &dLguide(msections), dvar(mvars), iusevars(mvars), 
     &nUCN(2), PsumUCN(2), xloss(2,2), array(marray), 
     &Rguide0(msections,4), 
     &scan_arr(mscan), iscan_nmultiple(msections), 
     &iscan_gmultiple(msections), iscan_dLguide(msections), 
     &iscan_Wguide(msections,2), iscan_Hguide(msections,2), 
     &iscan_Lguide(msections), iscan_Rguide(msections,4), 
     &iscan_mnumber(msections,5), iscan_refl0(msections,5), 
     &iscan_reflm(msections,5), iscan_Wwav(msections,5), 
     &iscan_tilt(msections,5), iscan_gshift(msections,5), 
     &iscan_chamfer(msections,5) 
		REAL*4	xmin, xmax, ymin, ymax, zmin, zmax, xp(101), yp(101), 
     &zp(101), xplot(2*nbin+1), yplot(2*nbin+1), eplot(2*nbin+1), 
     &yyr, zzr, phiyyr, phizzr, rlambda
		CHARACTER vname(mvars,2)*12, vblurb(mvars)*50, filename*50, 
     &Answer*1, blurb*50, title*80, input_file*50, 
     &string*50, out_file*50, ctime*8, source_file*50, cpol*1, 
     &line*(nline), charr(marray)*40, scanvarname*10, FMT*40, LFMT*120, 
     &line0*(nline)
		LOGICAL	through, binit, printrays, OK, exist, 
     &saveneutrons, batchjob, hitmon, saveresults, 
     &keepgoing, polarised, saveplot, usesavedneutrons, 
     &cylindrical(msections), UCN, intube, looplambda, 
     &saveplots, plot, showplots, circular_source, 
     &cylindrical_source, 
     &collimated_source, in, hparabolic(msections), 
     &vparabolic(msections), curved(msections), refl_file, scan, 
     &wrange, cavity(msections), polarised_in, newneutron, 
     &screeninput, nooutput

		DATA pi/3.141592654/, g/9.82E-03/, const/3956./,	! in units of mm and ms
     &vblurb(2*msections+1)/'Hor pos at exit of last segment'/, 
     &vname(2*msections+1,1)/'_yexit_'/, 
     &dvar(2*msections+1)/1.0/, 
     &vname(2*msections+1,2)/'Neutrons/bin'/, 
     &vblurb(2*msections+2)/'Ver pos at exit of last segment'/,  
     &vname(2*msections+2,1)/'_zexit_'/, 
     &dvar(2*msections+2)/1.0/, 
     &vname(2*msections+2,2)/'Neutrons/bin'/, 
     &vblurb(2*msections+3)/'Horizontal divergence on monitor'/, 
     &vname(2*msections+3,1)/'phiHsam'/, 
     &dvar(2*msections+3)/0.01/, 
     &vname(2*msections+3,2)/'Neutrons/bin'/, 
     &vblurb(2*msections+4)/'Vertical divergence on monitor'/, 
     &vname(2*msections+4,1)/'phiVsam'/, 
     &dvar(2*msections+4)/0.01/, 
     &vname(2*msections+4,2)/'Neutrons/bin'/, 
     &vblurb(2*msections+5)/'Horizontal position at monitor'/, 
     &vname(2*msections+5,1)/'ysample'/, 
     &dvar(2*msections+5)/.50/, 
     &vname(2*msections+5,2)/'Neutrons/bin'/, 
     &vblurb(2*msections+6)/'Vertical position at monitor'/, 
     &vname(2*msections+6,1)/'_zexit_'/, 
     &dvar(2*msections+6)/.50/, 
     &vname(2*msections+6,2)/'Neutrons/bin'/, 
     &vblurb(2*msections+7)/'Total flight path'/, 
     &vname(2*msections+7,1)/'Ltotal_'/, 
     &vname(2*msections+7,2)/'Neutrons/bin'/, 
     &vblurb(2*msections+8)/'Total number of reflections'/, 
     &vname(2*msections+8,1)/'_Nrefl_'/, 
     &dvar(2*msections+8)/1.0/, 
     &vname(2*msections+8,2)/'Neutrons/bin'/, 
     &vblurb(2*msections+9)/'Neutron losses in guide gaps'/, 
     &vname(2*msections+9,1)/'xloss0_'/, 
     &dvar(2*msections+9)/100.0/, 
     &vname(2*msections+9,2)/'Neutrons/bin'/, 
     &vblurb(2*msections+10)/'Neutron losses on face 1'/, 
     &vname(2*msections+10,1)/'xloss1_'/, 
     &dvar(2*msections+10)/100.0/, 
     &vname(2*msections+10,2)/'Neutrons/bin'/, 
     &vblurb(2*msections+11)/'Neutron losses on face 2'/, 
     &vname(2*msections+11,1)/'xloss2_'/, 
     &dvar(2*msections+11)/100.0/, 
     &vname(2*msections+11,2)/'Neutrons/bin'/, 
     &vblurb(2*msections+12)/'Neutron losses on face 3'/, 
     &vname(2*msections+12,1)/'xloss3_'/, 
     &dvar(2*msections+12)/100.0/, 
     &vname(2*msections+12,2)/'Neutrons/bin'/, 
     &vblurb(2*msections+13)/'Neutron losses on face 4'/, 
     &vname(2*msections+13,1)/'xloss4_'/, 
     &dvar(2*msections+13)/100.0/, 
     &vname(2*msections+13,2)/'Neutrons/bin'/, 
     &vblurb(2*msections+14)/'Neutron losses on face 5'/, 
     &vname(2*msections+14,1)/'xloss5_'/, 
     &dvar(2*msections+14)/100.0/, 
     &vname(2*msections+14,2)/'Neutrons/bin'/, 
     &vblurb(2*msections+15)/'Neutron losses on face 6'/, 
     &vname(2*msections+15,1)/'xloss6_'/, 
     &dvar(2*msections+15)/100.0/, 
     &vname(2*msections+15,2)/'Neutrons/bin'/, 
     &Si_lambda/0.0,0.20229,0.40457,0.52230,0.63968,0.90465,
     &1.0813,1.2794,1.4304,1.6517,2.0229,2.8607,3.4192,
     &4.0457,4.5232,5.2230,6.3968,9.0465,100.0/, ! wavelength in Å
     &Si_mu/0.12,0.1125,0.0925,0.0750,0.0625,0.0420,0.0325,
     &0.0265,0.0250,0.0215,0.0200,0.0205,0.0215,0.0250,
     &0.0275,0.0310,0.0365,0.0510,0.5/	! attenuation length in cm^-1

		scanvarname='xxx'
		iscan_Wmod=0
		iscan_Hmod=0
		iscan_Lap=0
		iscan_Wap=0
		iscan_Hap=0
		iscan_Hsourcem=0
		iscan_Vsourcem=0
		DO i=1,msections
			IF (i.LT.10) THEN
				WRITE(vblurb(2*i-1),811) i
  811			FORMAT(' Horizontal divergence at exit of segment ',I1)
				WRITE(vname(2*i-1,1),812) i
  812			FORMAT(' phiH',I1)
				WRITE(vblurb(2*i),813) i
  813			FORMAT(' Vertical divergence at exit of segment ',I1)
				WRITE(vname(2*i,1),814) i
  814			FORMAT(' phiV',I1)
			ELSE
				WRITE(vblurb(2*i-1),821) i
  821			FORMAT(' Horizontal divergence at exit of segment ',I2)
				WRITE(vname(2*i-1,1),822) i
  822			FORMAT(' phiH',I2)
				WRITE(vblurb(2*i),823) i
  823			FORMAT(' Vertical divergence at exit of segment ',I2)
				WRITE(vname(2*i,1),824) i
  824			FORMAT(' phiV',I2)
			END IF
			dvar(2*i-1)=0.02
			vname(2*i-1,2)='Neutrons/bin'
			dvar(2*i)=0.02
			vname(2*i,2)='Neutrons/bin'
			iscan_nmultiple(i)=0
			iscan_gmultiple(i)=0
			iscan_dLguide(i)=0
			iscan_Lguide(i)=0
			DO j=1,2
				iscan_Wguide(i,j)=0
				iscan_Hguide(i,j)=0
				iscan_Rguide(i,2*j-i)=0
			END DO
			DO j=1,5
				iscan_mnumber(i,j)=0
				iscan_refl0(i,j)=0
				iscan_reflm(i,j)=0
				iscan_Wwav(i,j)=0
				iscan_tilt(i,j)=0
				iscan_gshift(i,j)=0
				iscan_chamfer(i,j)=0
			END DO
		END DO


		UCN=.FALSE.
c		UCN=.TRUE.	! calculation for UCN production

c		show=.TRUE.

		polarised_in=.FALSE.
		polarised=.FALSE.
		IF (UCN) THEN
			density=2.185E+22	! He number density (cm^-3)
			sigmaS=3.48E-02		! Total scattering cross-section for 8.9A neutrons (barns)
			sigmaU=5.45E-05		! UCN-creation cross-section for 8.9A neutrons (barns)
			nsigmaS=density*sigmaS*1.E-25	! in mm^-1 (1 barn = 10^-24 cm2)
			dvar(2*msections+7)=5.
		ELSE
			sigmaS=0.
			sigmaU=0.
			nsigmaS=0.
			dvar(2*msections+7)=0.2
		END IF

		screeninput=.FALSE.	! normal running mode
c		screeninput=.TRUE.	! all input/output is from screen

		ierror=0
		IF (screeninput) THEN
			in_unit=5
			input_file='lardon.in'
		ELSE
			in_unit=10
   4		WRITE(*,5)
   5		FORMAT(' Enter input file : ',$)	! filename must end in .in
			READ(*,'(A50)') input_file
			IF (input_file.EQ.' ') input_file='testfile.in'
			OPEN(in_unit,FILE=input_file,STATUS='OLD',ERR=6)
			GOTO 7
    6		PRINT*,'Error - cannot open input file'
			GOTO 4
    7		CONTINUE
		END IF

c****************************************
c Read from input unit

		ierror=ierror+1
   10	READ(in_unit,'(A80)') line
		CALL read_char_line(line,charr,ncharr)
		IF (charr(1)(1:1).EQ.'!') GOTO 10
		IF (charr(1).NE.'SOURCE') GOTO 100

		usesavedneutrons=.FALSE.
		circular_source=.FALSE.
		cylindrical_source=.FALSE.
		collimated_source=.FALSE.
		Hsourcem=0.
		Vsourcem=0.
		Lap=0.
		Wap=0.
		Hap=0.
		halo=0.
		IF (charr(2).EQ.'RECTANGULAR') THEN
			READ(in_unit,'(A80)') line
			CALL read_real_line(line,array,narray)
			IF (narray.EQ.100) GOTO 100
			IF (narray.LT.5 .OR. narray.GT.6) GOTO 100
			Wmod=array(1)
			Hmod=array(2)
			Lap=array(3)
			Wap=array(4)
			Hap=array(5)
			IF (narray.EQ.6) halo=array(6)
		ELSE IF(charr(2).EQ.'CIRCULAR') THEN
			circular_source=.TRUE.
			READ(in_unit,'(A80)') line
			CALL read_real_line(line,array,narray)
			IF (narray.EQ.100) GOTO 100
			IF (narray.LT.4 .OR. narray.GT.5) GOTO 100
			Wmod=array(1)
			Hmod=Wmod
			Lap=array(2)
			Wap=array(3)
			Hap=array(4)
			IF (narray.EQ.5) halo=array(5)
		ELSE IF(charr(2).EQ.'CYLINDRICAL') THEN
			cylindrical_source=.TRUE.
			READ(in_unit,'(A80)') line
			CALL read_real_line(line,array,narray)
			IF (narray.EQ.100) GOTO 100
			IF (narray.LT.5 .OR. narray.GT.6) GOTO 100
			Wmod=array(1)
			Hmod=array(2)
			Lap=array(3)
			Wap=array(4)
			Hap=array(5)
			IF (narray.EQ.6) halo=array(6)
		ELSE IF (charr(2).EQ.'FILE') THEN
			usesavedneutrons=.TRUE.
			READ(in_unit,'(A50)',ERR=100) source_file	! filename must end in .dat
		ELSE IF (charr(2).EQ.'RECTANGULAR_COLLIMATED') THEN
			collimated_source=.TRUE.
			READ(in_unit,'(A80)') line
			CALL read_real_line(line,array,narray)
			IF (narray.EQ.100) GOTO 100
			IF (narray.LT.4 .OR. narray.GT.5) GOTO 100
			Wmod=array(1)
			Hmod=array(2)
			Hsourcem=array(3)	! horizontal and vertical
			Vsourcem=array(4)	! m-number
			IF (narray.EQ.5) halo=array(5)
		ELSE IF (charr(2).EQ.'CIRCULAR_COLLIMATED') THEN
			circular_source=.TRUE.
			collimated_source=.TRUE.
			READ(in_unit,'(A80)') line
			CALL read_real_line(line,array,narray)
			IF (narray.EQ.100) GOTO 100
			IF (narray.LT.3 .OR. narray.GT.4) GOTO 100
			Wmod=array(1)
			Hmod=Wmod
			Hsourcem=array(2)
			Vsourcem=array(3)
			IF (narray.EQ.4) halo=array(4)
		ELSE IF (charr(2).EQ.'CYLINDRICAL_COLLIMATED') THEN
			cylindrical_source=.TRUE.
			collimated_source=.TRUE.
			READ(in_unit,'(A80)') line
			CALL read_real_line(line,array,narray)
			IF (narray.EQ.100) GOTO 100
			IF (narray.LT.4 .OR. narray.GT.5) GOTO 100
			Wmod=array(1)
			Hmod=array(2)
			Hsourcem=array(3)	! horizontal and vertical
			Vsourcem=array(4)	! m-number
			IF (narray.EQ.5) halo=array(5)
		END IF
		IF (Wmod.EQ.-999.) iscan_Wmod=1
		IF (Hmod.EQ.-999.) iscan_Hmod=1
		IF (Lap.EQ.-999.) iscan_Lap=1
		IF (Wap.EQ.-999.) iscan_Wap=1
		IF (Hap.EQ.-999.) iscan_Hap=1
		IF (Hsourcem.EQ.-999.) iscan_Hsourcem=1
		IF (Vsourcem.EQ.-999.) iscan_Vsourcem=1

		nsections=0
		DO i=1,msections

			ierror=ierror+1
   20		READ(in_unit,'(A80)') line
			CALL read_char_line(line,charr,ncharr)
			IF (charr(1)(1:1).EQ.'!') GOTO 20
			IF (charr(1).EQ.'END') GOTO 60	! no more guide sections

			IF (charr(1).EQ.'GUIDE') THEN
				nsections=nsections+1
				curved(i)=.FALSE.
				cylindrical(i)=.FALSE.
				hparabolic(i)=.FALSE.
				vparabolic(i)=.FALSE.
				cavity(i)=.FALSE.
				tilted(i)=.FALSE.
				nmultiple(i)=1
				gmultiple(i)=0.
				ivcavity(i)=0
				READ(in_unit,'(A80)') line
				CALL read_real_line(line,array,narray)
				IF (charr(2).EQ.'STRAIGHT' .OR. 
     &charr(2).EQ.'CONVERGING' .OR. charr(2).EQ.'DIVERGING') THEN
					IF (charr(3).EQ.'MULTIPLE') THEN
						IF (narray.NE.2) GOTO 100
						nmultiple(i)=NINT(array(1))
						gmultiple(i)=array(2)
						READ(in_unit,'(A80)') line
						CALL read_real_line(line,array,narray)
					END IF
					IF (narray.EQ.100 .OR. narray.LT.4) GOTO 100
					dLguide(i)=array(1)
					Wguide(i,1)=array(2)
					Hguide(i,1)=array(3)
					Lguide(i)=array(4)
					Rguide0(i,1)=0.
					IF (narray.EQ.4) THEN
						Wguide(i,2)=Wguide(i,1)
						Hguide(i,2)=Hguide(i,1)
					ELSE IF (narray.EQ.6) THEN
						Wguide(i,2)=array(5)
						Hguide(i,2)=array(6)
					ELSE
						GOTO 100
					END IF
				ELSE IF (charr(2).EQ.'CURVED' .OR.
     &charr(2).EQ.'CURVED_CONVERGING' .OR. 
     &charr(2).EQ.'CONVERGING_CURVED' .OR. 
     &charr(2).EQ.'CURVED_DIVERGING' .OR. 
     &charr(2).EQ.'DIVERGING_CURVED') THEN
					IF (charr(3).EQ.'MULTIPLE') THEN
						IF (narray.NE.2) GOTO 100
						nmultiple(i)=array(1)
						gmultiple(i)=array(2)
						READ(in_unit,'(A80)') line
						CALL read_real_line(line,array,narray)
					END IF
					IF (narray.EQ.100 .OR. narray.LT.5) GOTO 100
					curved(i)=.TRUE.
					dLguide(i)=array(1)
					Wguide(i,1)=array(2)
					Hguide(i,1)=array(3)
					Lguide(i)=array(4)
					Rguide0(i,1)=array(5)
					Rguide0(i,3)=0.
					IF (narray.EQ.5) THEN
						Wguide(i,2)=Wguide(i,1)
						Hguide(i,2)=Hguide(i,1)
					ELSE IF (narray.EQ.7) THEN
						Wguide(i,2)=array(6)
						Hguide(i,2)=array(7)
					ELSE
						GOTO 100
					END IF
				ELSE IF (charr(2).EQ.'CYLINDRICAL') THEN
					IF (charr(3).EQ.'MULTIPLE') THEN
						IF (narray.NE.2) GOTO 100
						nmultiple(i)=array(1)
						gmultiple(i)=array(2)
						READ(in_unit,'(A80)') line
						CALL read_real_line(line,array,narray)
					END IF
					IF (narray.NE.3) GOTO 100
					cylindrical(i)=.TRUE.
					dLguide(i)=array(1)
					Wguide(i,1)=array(2)
					Hguide(i,1)=Wguide(i,1)
					Lguide(i)=array(3)
					Wguide(i,2)=Wguide(i,1)
					Hguide(i,2)=Hguide(i,1)
					Rguide0(i,1)=0.
				ELSE IF (charr(2).EQ.'PARABOLIC' .OR. 
     &charr(2).EQ.'HORIZONTALLY_PARABOLIC' .OR. 
     &charr(2).EQ.'VERTICALLY_PARABOLIC') THEN
					IF (charr(3).EQ.'MULTIPLE') THEN
						IF (narray.NE.2) GOTO 100
						nmultiple(i)=array(1)
						gmultiple(i)=array(2)
						READ(in_unit,'(A80)') line
						CALL read_real_line(line,array,narray)
					END IF
					IF (narray.NE.6) GOTO 100
					dLguide(i)=array(1)
					Wguide(i,1)=array(2)
					Hguide(i,1)=array(3)
					Lguide(i)=array(4)
					Wguide(i,2)=array(5)
					Hguide(i,2)=array(6)
					IF (charr(2).EQ.'VERTICALLY_PARABOLIC') THEN
						vparabolic(i)=.TRUE.
						IF (Hguide(i,1).LE.Hguide(i,2)) THEN
							PRINT*,'Entrance height must be > ',
     &'Exit height for vertically parabolic guide'
							GOTO 100
						END IF
					ELSE
						hparabolic(i)=.TRUE.
						IF (Wguide(i,1).LE.Wguide(i,2)) THEN
							PRINT*,'Entrance width must be > ',
     &'Exit width for horizontally parabolic guide'
							GOTO 100
						END IF
					END IF
				ELSE IF (charr(2).EQ.'TRANSMISSION_MIRROR') THEN
					IF (charr(3).EQ.'MULTIPLE') THEN
						IF (narray.NE.2) GOTO 100
						nmultiple(i)=array(1)
						gmultiple(i)=array(2)
						READ(in_unit,'(A80)') line
						CALL read_real_line(line,array,narray)
					END IF
					IF (narray.EQ.100 .OR. narray.LT.5) GOTO 100
					cavity(i)=.TRUE.
					dLguide(i)=array(1)
					Wguide(i,1)=array(2)
					Hguide(i,1)=array(3)
					Lguide(i)=array(4)
					cthickness(i)=array(5)
					IF (narray.EQ.5) THEN
						Wguide(i,2)=Wguide(i,1)
						Hguide(i,2)=Hguide(i,1)
					ELSE IF (narray.EQ.7) THEN
						Wguide(i,2)=array(6)
						Hguide(i,2)=array(7)
					ELSE
						GOTO 100
					END IF
				ELSE
					GOTO 100
				END IF
			ELSE
				GOTO 100
			END IF

			IF (nmultiple(i).EQ.-999) iscan_nmultiple(i)=1
			IF (gmultiple(i).EQ.-999.) iscan_gmultiple(i)=1
			IF (dLguide(i).EQ.-999.) iscan_dLguide(i)=1
			DO j=1,2
				IF (Wguide(i,j).EQ.-999.) iscan_Wguide(i,j)=1
				IF (Hguide(i,j).EQ.-999.) iscan_Hguide(i,j)=1
				IF (Rguide0(i,2*j-1).EQ.-999.) iscan_Rguide(i,2*j-i)=1
			END DO
			IF (Lguide(i).EQ.-999.) iscan_Lguide(i)=1

			Rguide0(i,1)=Rguide0(i,1)*1.E+03
			Rguide0(i,3)=Rguide0(i,3)*1.E+03

			ierror=ierror+1
   30		READ(in_unit,'(A80)') line
			CALL read_char_line(line,charr,ncharr)
			IF (charr(1)(1:1).EQ.'!') GOTO 30

			DO j=1,4	! default values
				mnumber(i,j)=1.
				refl0(i,j)=1.			! R=1.0 up to m=1, then decreasing
				reflm(i,j)=1.
				Wwav(i,j)=0.
				tilt(i,j)=0.
				gshift(i,j)=0.
				chamfer(i,j)=0.
			END DO

			IF (charr(1).EQ.'FACE') THEN
				READ(charr(2),*,ERR=100) iface
				nfaces=4
				IF (cavity(i)) nfaces=5
			ELSE IF (charr(1).EQ.'ALL' .AND. charr(2).EQ.'FACES') THEN
				nfaces=1
			ELSE
				PRINT*,'Syntax error - faces not specified'
				GOTO 100
			END IF

			DO iface=1,nfaces
				IF (iface.NE.1) THEN
   35				READ(in_unit,'(A80)') line
					CALL read_char_line(line,charr,ncharr)
					IF (charr(1)(1:1).EQ.'!') GOTO 35
					READ(charr(2),*,ERR=100) ifac
					IF (ifac.NE.iface) GOTO 100
				END IF
				refl_file=.FALSE.
				IF (ncharr.GE.3) THEN
					IF (charr(3).EQ.'FILE') refl_file=.TRUE.
				END IF
				IF (iface.EQ.5) THEN	! transmission mirror
					nsides(i)=0
					IF (ncharr.GE.4) THEN
						IF (charr(3).EQ.'SINGLE_SIDED') nsides(i)=1
						IF (charr(3).EQ.'DOUBLE_SIDED') nsides(i)=2
						IF (charr(4).EQ.'V-') ivcavity(i)=-2
						IF (charr(4).EQ.'LEFT') ivcavity(i)=-1
						IF (charr(4).EQ.'RIGHT') ivcavity(i)=1
						IF (charr(4).EQ.'V+') ivcavity(i)=2
					END IF
					IF (nsides(i).EQ.0 .OR. ivcavity(i).EQ.0) THEN
						PRINT*,'Syntax error - face 5 not properly defined'
						GOTO 100
					END IF
					IF (ncharr.GE.5) THEN
						IF (charr(5).EQ.'FILE') refl_file=.TRUE.
					END IF
				END IF
   40			READ(in_unit,'(A80)') line
				CALL read_char_line(line,charr,ncharr)
				IF (charr(1)(1:1).EQ.'!') GOTO 40
				IF (refl_file) THEN
					refl0(i,iface)=-1.
					reflfile(i,iface)=charr(1)
					OPEN(15,FILE=reflfile(i,iface),STATUS='OLD',ERR=41)
					READ(15,'(A80)') line0
					CALL read_real_line(line0,array,narray)
c					Qarr(i,iside,ipol,i)=4.*pi*SIN(array(1)*pi/180.)/7.5	! T3 data
					Qarr(i,iface,1)=array(1)
					Parr(i,iface,1,1)=array(2)
					IF (narray.EQ.2) THEN
						DO iQ=2,mQ
							READ(15,*,END=42) Q, P
							Qarr(i,iface,iQ)=Q
							Parr(i,iface,1,iQ)=P
							Parr(i,iface,2,iQ)=P
						END DO
					ELSE
						polarised=.TRUE.
						mnumber(i,iface)=-1.
						Parr(i,iface,2,1)=array(3)
						DO iQ=2,mQ
							READ(15,*,END=42) Q, P1, P2
							Qarr(i,iface,iQ)=Q
							Parr(i,iface,1,iQ)=P1
							Parr(i,iface,2,iQ)=P2
						END DO
					END IF
					GOTO 42
   41				PRINT*,'Error - cannot open file ',reflfile(i,iface)
					GOTO 1000
   42				nQ(i,iface)=iQ-1
					CLOSE(15)
					DO il=1,nline-1
						IF (line(il:il).NE.' ' .AND.  
     &line(il+1:il+1).EQ.' ') GOTO 45
					END DO
   45				nl=il
					DO il=1,nl
						line(il:il)=' '
					END DO
					CALL read_real_line(line,array,narray)

					IF (narray.GT.4) GOTO 100
					IF (narray.GE.1) Wwav(i,iface)=array(1)
					IF (narray.GE.2) tilt(i,iface)=array(2)
					IF (tilt(i,iface).NE.0.) tilted(i)=.TRUE.
					IF (narray.GE.3) gshift(i,iface)=array(3)
					IF (narray.GE.4) chamfer(i,iface)=array(4)
				ELSE
					CALL read_real_line(line,array,narray)
					IF (narray.GT.7) GOTO 100
					mnumber(i,iface)=array(1)
					IF (narray.GE.2) refl0(i,iface)=array(2)
					IF (narray.GE.3) THEN
						reflm(i,iface)=array(3)
					ELSE
						IF (mnumber(i,iface).GT.1.) 
     &reflm(i,iface)=refl0(i,iface)-0.1*(mnumber(i,iface)-1.)
						IF (mnumber(i,iface).GT.1.119 .AND. 	! make an exception for Ni-58
     &mnumber(i,iface).LT.1.201) reflm(i,iface)=1.0					
					END IF
					IF (narray.GE.4) Wwav(i,iface)=array(4)
					IF (narray.GE.5) tilt(i,iface)=array(5)
					IF (tilt(i,iface).NE.0.) tilted(i)=.TRUE.
					IF (narray.GE.6) gshift(i,iface)=array(6)
					IF (narray.GE.7) chamfer(i,iface)=array(7)
					IF (show) PRINT*,'m =',REAL(mnumber(i,iface)),
     &' R0 =',REAL(refl0(i,iface)),' Rm =',REAL(reflm(i,iface)), 
     &' Wwav =',REAL(Wwav(i,iface)),' tilt=',REAL(tilt(i,iface)), 
     &' shift =',REAL(gshift(i,iface)),
     &' chamfer=',REAL(chamfer(i,iface))
					IF (mnumber(i,iface).LT.0.) polarised=.TRUE.
   50				CONTINUE
				END IF	! reflectivity file or parameter input
				IF (show) PRINT*,'polarised=',polarised
				IF (mnumber(i,iface).EQ.-999.) iscan_mnumber(i,iface)=1
				IF (refl0(i,iface).EQ.-999.) iscan_refl0(i,iface)=1
				IF (reflm(i,iface).EQ.-999.) iscan_reflm(i,iface)=1
				IF (Wwav(i,iface).EQ.-999.) iscan_Wwav(i,iface)=1
				IF (tilt(i,iface).EQ.-999.) iscan_tilt(i,iface)=1
				IF (gshift(i,iface).EQ.-999.) iscan_gshift(i,iface)=1
				IF (chamfer(i,iface).EQ.-999.) iscan_chamfer(i,iface)=1
				Wwav(i,iface)=Wwav(i,iface)*pi/180.
				tilt(i,iface)=tilt(i,iface)*pi/180.

			END DO	! iface=1,nfaces

			IF (nfaces.EQ.1 .AND. .NOT.cylindrical(i)) THEN
				DO iface=2,4
					mnumber(i,iface)=mnumber(i,1)
					IF (iscan_mnumber(i,1).EQ.1) iscan_mnumber(i,iface)=1
					refl0(i,iface)=refl0(i,1)
					IF (iscan_refl0(i,1).EQ.1) iscan_refl0(i,iface)=1
					reflm(i,iface)=reflm(i,1)
					IF (iscan_reflm(i,1).EQ.1) iscan_reflm(i,iface)=1
					Wwav(i,iface)=Wwav(i,1)
					IF (iscan_Wwav(i,1).EQ.1) iscan_Wwav(i,iface)=1
					chamfer(i,iface)=chamfer(i,1)
					IF (iscan_chamfer(i,1).EQ.1) iscan_chamfer(i,iface)=1
					IF (iface.EQ.2) THEN
						tilt(i,iface)=tilt(i,1)
						IF (iscan_tilt(i,1).EQ.1) iscan_tilt(i,iface)=1
						gshift(i,iface)=gshift(i,1)
						IF (iscan_gshift(i,1).EQ.1) iscan_gshift(i,iface)=1
					END IF
					IF (refl0(i,1).LT.0.) THEN	! reflectivity from input file
						refl0(i,iface)=-1.
						nQ(i,iface)=nQ(i,1)
						DO iQ=1,nQ(i,iface)
							Qarr(i,iface,iQ)=Qarr(i,1,iQ)
							Parr(i,iface,1,iQ)=Parr(i,1,1,iQ)
							Parr(i,iface,2,iQ)=Parr(i,1,2,iQ)
						END DO
					END IF
				END DO
			END IF

		END DO	! i=1,msections

   60	CONTINUE
		IF (show) PRINT*,'No more guide sections'
c		PRINT*,'Hit return to continue'
c		READ(*,'(A1)') Answer

		gravity=.FALSE.
		sideways=.FALSE.
		scan=.FALSE.
		wrange=.FALSE.
		lambda=0.
		divlimit=0.50	! +/- divlimit is useful divergence
		ierror=ierror+1
   70	IF (charr(1)(1:3).EQ.'END' .AND. charr(2)(1:5).EQ.'GUIDE') THEN
			READ(in_unit,'(A80)') line
			CALL read_real_line(line,array,narray)
			Lend=array(1)
			Wmon=array(2)
			Hmon=array(3)
			IF (narray.EQ.3) THEN
				Hmonshift=0.
				Vmonshift=0.
			ELSE IF (narray.EQ.5) THEN
				Hmonshift=array(4)
				Vmonshift=array(5)
			ELSE
				PRINT*,'Error - narray should be 3 or 5. It is =',narray
			END IF	
			IF (Lend.EQ.-999.) iscan_Lend=1
			IF (Wmon.EQ.-999.) iscan_Wmon=1
			IF (Hmon.EQ.-999.) iscan_Hmon=1
			IF (Hmonshift.EQ.-999.) iscan_Hmonshift=1
			IF (Vmonshift.EQ.-999.) iscan_Vmonshift=1
		ELSE IF (charr(1).EQ.'GRAVITY' .AND. ncharr.GE.2) THEN
			IF (charr(2).EQ.'SIDEWAYS' .OR. 
     &(ncharr.EQ.3 .AND. charr(3).EQ.'SIDEWAYS')) sideways=.TRUE.
			IF (sideways) gravity=.TRUE.
			IF (charr(2).EQ.'OFF' .OR. 
     &(ncharr.EQ.3 .AND. charr(3).EQ.'OFF')) gravity=.FALSE.
			IF (charr(2).EQ.'ON' .OR. 
     &(ncharr.EQ.3 .AND. charr(3).EQ.'ON')) gravity=.TRUE.
		ELSE IF (charr(1)(1:10).EQ.'WAVELENGTH') THEN
			READ(in_unit,'(A80)') line
			IF (ncharr.GE.2 .AND. charr(2)(1:4).EQ.'FILE') THEN
				INQUIRE(FILE=line,EXIST=exist)
				IF (.NOT.exist) PRINT*,'Error - file does not exist: ',line
				OPEN(11,FILE=line,STATUS='OLD')
				DO ilambda=1,mlambda
					READ(11,*,END=240) lambda_arr(ilambda)
				END DO
  240			CLOSE(11)
				nlambda=ilambda-1
				IF (nlambda.LE.0) PRINT*,
     &'Error - no wavelengths found in file: ',line
				WRITE(*,244) nlambda
  244			FORMAT(I6,' wavelengths :')
				WRITE(*,245) (i, lambda_arr(i), i=1,nlambda)
  245			FORMAT(' Wavelength',I2,' ='F6.2,' A')
				looplambda=.TRUE.
				ilambda=0
			ELSE IF (ncharr.GE.2 .AND. charr(2)(1:4).EQ.'SCAN') THEN
				CALL read_real_line(line,array,narray)
				IF (narray.EQ.3) THEN
					lambdamin=array(1)
					dlambda=array(2)
					lambdamax=array(3)
				ELSE
					PRINT*,'Error: narray should be 3. It is =',narray
				END IF
				nlambda=NINT((lambdamax-lambdamin)/dlambda)+1
				WRITE(*,260) nlambda, lambdamin, dlambda, 
     &lambdamin+REAL(nlambda-1)*dlambda
  260			FORMAT(I6,' wavelengths:',F6.2,' (',F4.2,')',F6.2)
				DO i=1,nlambda
					lambda_arr(i)=lambdamin+REAL(i-1)*dlambda
				END DO
				looplambda=.TRUE.
				ilambda=0
			ELSE IF (ncharr.GE.2 .AND. charr(2)(1:5).EQ.'RANGE') THEN
				CALL read_real_line(line,array,narray)
				IF (narray.EQ.2) THEN
					lambdamin=array(1)
					lambdamax=array(2)
				ELSE
					PRINT*,'Error: narray should be 2. It is =',narray
				END IF
				wrange=.TRUE.
			ELSE
				CALL read_real_line(line,array,narray)
				IF (narray.EQ.1) THEN
					lambda=array(1)
				ELSE
					PRINT*,'Error: narray should be 1. It is =',narray
				END IF
				looplambda=.FALSE.
			END IF
		ELSE IF (charr(1)(1:4).EQ.'SCAN' .AND. ncharr.EQ.2) THEN
			READ(in_unit,'(A80)') line
			IF (ncharr.GE.2 .AND. charr(2)(1:4).EQ.'FILE') THEN
				INQUIRE(FILE=line,EXIST=exist)
				IF (.NOT.exist) PRINT*,'Error - file does not exist: ',line
				OPEN(11,FILE=line,STATUS='OLD')
				DO iscan=1,mscan
					READ(11,*,END=830) scan_arr(iscan)
				END DO
  830			CLOSE(11)
				nscan=iscan-1
				IF (nscan.LE.0) PRINT*,
     &'Error - no scan points found in file: ',line
				IF (.NOT.screeninput) WRITE(*,831) nscan
  831			FORMAT(I6,' scan points :')
				IF (.NOT.screeninput) WRITE(*,832) 
     &(i, scan_arr(i), i=1,nscan)
  832			FORMAT(' Scan point',I3,' ='F9.3)
				scan=.TRUE.
				iscan=0
			ELSE IF (charr(2).EQ.'ON') THEN
				CALL read_real_line(line,array,narray)
				IF (narray.EQ.3) THEN
					scanmin=array(1)
					dscan=array(2)
					scanmax=array(3)
				ELSE
					PRINT*,'Error: narray should be 3. It is =',narray
				END IF
				nscan=NINT((scanmax-scanmin)/dscan)+1
				IF (.NOT.screeninput) WRITE(*,833) nscan, 
     &scanmin, dscan, scanmin+REAL(nscan-1)*dscan
  833			FORMAT(I6,' scan points:',F9.3,' (',F6.3,')',F9.3)
				DO i=1,nscan
					scan_arr(i)=scanmin+REAL(i-1)*dscan
				END DO
				scan=.TRUE.
				iscan=0
			END IF
		ELSE IF (charr(1)(1:8).EQ.'NEUTRONS') THEN
			READ(in_unit,'(A80)') line
			CALL read_real_line(line,array,narray)
			IF (narray.EQ.2) THEN
				neutrons=NINT(array(1))
				iseed=NINT(array(2))
			ELSE
				PRINT*,'Error: narray should be 2. It is =',narray
			END IF
		ELSE IF (charr(1)(1:6).EQ.'OUTPUT') THEN
			READ(in_unit,'(A80)') line
			CALL read_real_line(line,array,narray)
			IF (narray.EQ.2) THEN
				iresults=NINT(array(1))
				idebug=NINT(array(2))
				nooutput=.FALSE.
				IF (idebug.EQ.-2) nooutput=.TRUE.
			ELSE
				PRINT*,'Error: narray should be 2. It is =',narray
			END IF
		ELSE IF (charr(1)(1:3).EQ.'END' .AND. 
     &charr(2)(1:5).EQ.'INPUT') THEN
			GOTO 90
		ELSE IF (charr(1)(1:10).EQ.'ACCEPTANCE') THEN
			READ(in_unit,'(A80)') line
			CALL read_real_line(line,array,narray)
			IF (narray.EQ.1) THEN
				divlimit=array(1)
			ELSE
				PRINT*,'Error: narray should be 1. It is =',narray
			END IF
		END IF

		READ(in_unit,'(A80)') line
		CALL read_char_line(line,charr,ncharr)
		GOTO 70	

   90	IF (.NOT.screeninput) CLOSE(in_unit)

		IF (lambda.EQ.0. .AND. 
     &.NOT.(looplambda.OR.usesavedneutrons.OR.wrange)) THEN
			PRINT*,'Error - wavelength has not been specified'
		END IF

		IF (scan) THEN
			IF (iscan_Wmod.EQ.1) scanvarname='Wmod'
			IF (iscan_Hmod.EQ.1) scanvarname='Hmod'
			IF (iscan_Lap.EQ.1) scanvarname='Lap'
			IF (iscan_Wap.EQ.1) scanvarname='Wap'
			IF (iscan_Hap.EQ.1) scanvarname='Hap'
			IF (iscan_Hsourcem.EQ.1) scanvarname='Hsourcem'
			IF (iscan_Vsourcem.EQ.1) scanvarname='Hsourcem'
			DO i=1,nsections
				IF (iscan_nmultiple(i).EQ.1) scanvarname='nmultiple'
				IF (iscan_gmultiple(i).EQ.1) scanvarname='gmultiple'
				IF (iscan_dLguide(i).EQ.1) scanvarname='dLguide'
				IF (iscan_Lguide(i).EQ.1) scanvarname='Lguide'
				DO j=1,2
					IF (iscan_Rguide(i,2*j-1).EQ.1) scanvarname='Rguide'
					IF (iscan_Wguide(i,j).EQ.1) scanvarname='Wguide'
					IF (iscan_Hguide(i,j).EQ.1) scanvarname='Hguide'
				END DO
				DO iface=1,4
					IF (iscan_mnumber(i,iface).EQ.1) scanvarname='mnumber'
					IF (iscan_refl0(i,iface).EQ.1) scanvarname='refl0'
					IF (iscan_reflm(i,iface).EQ.1) scanvarname='reflm'
					IF (iscan_Wwav(i,iface).EQ.1) scanvarname='Wwav'
					IF (iscan_tilt(i,iface).EQ.1) scanvarname='tilt'
					IF (iscan_gshift(i,iface).EQ.1) scanvarname='gshift'
					IF (iscan_chamfer(i,iface).EQ.1) scanvarname='chamfer'
				END DO
			END DO
			IF (iscan_Lend.EQ.1) scanvarname='Lend'
			IF (iscan_Wmon.EQ.1) scanvarname='Wmon'
			IF (iscan_Hmon.EQ.1) scanvarname='Hmon'
			IF (iscan_Hmonshift.EQ.1) scanvarname='Hmonshift'
			IF (iscan_Vmonshift.EQ.1) scanvarname='Vmonshift'
			IF (scanvarname.EQ.'xxx') THEN
				PRINT*,'Error - scan set up with no scan variable'
				GOTO 1000
			END IF
			PRINT*,'scanvarname=',scanvarname
		END IF

		GOTO 101
  100	PRINT*,'Error in reading from input_file'
		PRINT*,'line =',line
		PRINT*,'ncharr =',ncharr,' narray =',narray,' ierror =',ierror
		CLOSE(10)
		GOTO 1000
  101	CONTINUE

c	step up scan variable if there is one

		IF (scan) THEN
			iscan=iscan+1
			scanvar=scan_arr(iscan)
			IF (.NOT.nooutput) WRITE(*,102) scanvarname, scanvar
  102		FORMAT(' scan variable : ',A10,' =',F9.3)
			IF (iscan_Wmod.EQ.1) Wmod=scanvar
			IF (iscan_Hmod.EQ.1) Hmod=scanvar
			IF (iscan_Lap.EQ.1) Lap=scanvar
			IF (iscan_Wap.EQ.1) Wap=scanvar
			IF (iscan_Hap.EQ.1) Hap=scanvar
			IF (iscan_Hsourcem.EQ.1) Hsourcem=scanvar
			IF (iscan_Vsourcem.EQ.1) Vsourcem=scanvar
			DO i=1,nsections
				IF (iscan_nmultiple(i).EQ.1) nmultiple(i)=NINT(scanvar)
				IF (iscan_gmultiple(i).EQ.1) gmultiple(i)=scanvar
				IF (iscan_dLguide(i).EQ.1) dLguide(i)=scanvar
				IF (iscan_Lguide(i).EQ.1) Lguide(i)=scanvar
				DO j=1,2
					IF (iscan_Rguide(i,2*j-1).EQ.1) 
     &Rguide0(i,2*j-1)=scanvar*1.E+03
					IF (iscan_Wguide(i,j).EQ.1) Wguide(i,j)=scanvar
					IF (iscan_Hguide(i,j).EQ.1) Hguide(i,j)=scanvar
				END DO
				DO iface=1,4
					IF (iscan_mnumber(i,iface).EQ.1) mnumber(i,iface)=scanvar
					IF (iscan_refl0(i,iface).EQ.1) refl0(i,iface)=scanvar
					IF (iscan_reflm(i,iface).EQ.1) reflm(i,iface)=scanvar
					IF (iscan_Wwav(i,iface).EQ.1) Wwav(i,iface)=scanvar*pi/180.
					IF (iscan_tilt(i,iface).EQ.1) tilt(i,iface)=scanvar*pi/180.
					IF (iscan_gshift(i,iface).EQ.1) gshift(i,iface)=scanvar
					IF (iscan_chamfer(i,iface).EQ.1) chamfer(i,iface)=scanvar
				END DO
			END DO
			IF (iscan_Lend.EQ.1) Lend=scanvar
			IF (iscan_Wmon.EQ.1) Wmon=scanvar
			IF (iscan_Hmon.EQ.1) Hmon=scanvar
			IF (iscan_Hmonshift.EQ.1) Hmonshift=scanvar
			IF (iscan_Vmonshift.EQ.1) Vmonshift=scanvar
c			DO iface=1,4	! bodge for D17 S-bender
c				IF (iscan_tilt(3,iface).EQ.1 .AND. 
c     &iscan_tilt(6,iface).EQ.1) tilt(6,iface)=-tilt(3,iface)
c			END DO
		END IF

c****************************************
c Recalculate and show parameters from input file

		IF (.NOT.nooutput) THEN
			IF (usesavedneutrons) THEN
				PRINT*,' Source file : ',source_file
			ELSE IF (collimated_source) THEN
				IF (circular_source) THEN
					WRITE(*,106) NINT(Wmod), NINT(halo), 
     &Hsourcem, Vsourcem
  106				FORMAT(
     &' Diameter of source, halo, Hor and Ver m-number'/,
     &2I5,2F6.2,5X,'- circular collimated source')
				ELSE IF (cylindrical_source) THEN
					WRITE(*,107) NINT(Wmod), NINT(Hmod), NINT(halo), 
     &Hsourcem, Vsourcem
  107				FORMAT(
     &' W, H of source, halo, Hor and Ver m-number'/,
     &3I5,2F6.2,5X,'- cylindrical collimated source')
				ELSE
					WRITE(*,108) NINT(Wmod), NINT(Hmod), NINT(halo), 
     &Hsourcem, Vsourcem
  108				FORMAT(
     &' W, H of source, halo, Hor and Ver m-number'/,
     &3I5,2F6.2,5X,'- rectangular collimated source')
				END IF
			ELSE
				IF (circular_source) THEN
					WRITE(*,109) NINT(Wmod), NINT(halo), 
     &NINT(Lap), NINT(Wap), NINT(Hap)
  109				FORMAT(
     &' Diameter of source, halo, Distance to and W, H',
     &' of virtual aperture'/,2I5,I6,2I5,5X,'- circular source')
				ELSE IF (cylindrical_source) THEN
					WRITE(*,110) NINT(Wmod), NINT(Hmod), NINT(halo), 
     &NINT(Lap), NINT(Wap), NINT(Hap)
  110				FORMAT(
     &' W, H of source, halo, Distance to and W, H',
     &' of virtual aperture'/3I5,I6,2I5,5X,'- cylindrical source')
				ELSE
					WRITE(*,111) NINT(Wmod), NINT(Hmod), NINT(halo), 
     &NINT(Lap), NINT(Wap), NINT(Hap)
  111				FORMAT(
     &' W, H of source, halo, Distance to and W, H',
     &' of virtual aperture'/3I5,I6,2I5,5X,'- rectangular source')
				END IF
			END IF
			WRITE(*,115) nsections
  115		FORMAT(I5,' guide sections')
		END IF

		DO i=1,nsections
			nfaces=4
			string='straight'
			R=0.
			IF (curved(i)) THEN
				string='curved'
				R=Rguide0(i,1)*1.E-03
			ELSE IF (cavity(i)) THEN
				nfaces=5
				string='cavity'
			ELSE IF (cylindrical(i)) THEN
				nfaces=1
				string='cylindrical'
				tilt(i,2)=tilt(i,1)
			ELSE IF (hparabolic(i)) THEN
				string='hparabolic'
			ELSE IF (vparabolic(i)) THEN
				string='vparabolic'
			END IF
			IF (tilted(i)) THEN
				IF (cylindrical(i)) THEN
					Hshift(i)=gshift(i,1)
					Vshift(i)=0.
				ELSE
					IF (tilt(i,1).NE.tilt(i,2)) THEN
						Wguide(i,1)=Wguide(i,1)-gshift(i,1)+gshift(i,2)
						Wguide(i,2)=Wguide(i,1)-Lguide(i)*SIN(tilt(i,1))
     &+Lguide(i)*SIN(tilt(i,2))
					END IF
					IF (tilt(i,3).NE.tilt(i,4)) THEN
						Hguide(i,1)=Hguide(i,1)-gshift(i,3)+gshift(i,4)
						Hguide(i,2)=Hguide(i,1)-Lguide(i)*SIN(tilt(i,3))
     &+Lguide(i)*SIN(tilt(i,4))
					END IF
					Hshift(i)=(gshift(i,1)+gshift(i,2))/2.
					Vshift(i)=(gshift(i,3)+gshift(i,4))/2.
				END IF
			ELSE
				tilt(i,1)=0.
				IF (Lguide(i).GT.0.) 
     &tilt(i,1)=ASIN((Wguide(i,1)-Wguide(i,2))/(2.*Lguide(i)))
				tilt(i,2)=-tilt(i,1)
				tilt(i,3)=0.
				IF (Lguide(i).GT.0.) 
     &tilt(i,3)=ASIN((Hguide(i,1)-Hguide(i,2))/(2.*Lguide(i)))
				tilt(i,4)=-tilt(i,3)
				IF (nfaces.EQ.1) THEN
					Hshift(i)=gshift(i,1)
					Vshift(i)=0.
				ELSE
					Hshift(i)=(gshift(i,1)+gshift(i,2))/2.
					Vshift(i)=(gshift(i,3)+gshift(i,4))/2.
				END IF
			END IF
			IF (.NOT.nooutput) WRITE(*,120) i, NINT(dLguide(i)), 
     &NINT(Wguide(i,1)), NINT(Hguide(i,1)), NINT(Lguide(i)), 
     &NINT(Wguide(i,2)), NINT(Hguide(i,2)), NINT(R), 
     &NINT(Hshift(i)), NINT(Vshift(i)), 
     &string(1:11), nmultiple(i), gmultiple(i)
  120		FORMAT(' Guide',I2,' : Gap  W1   H1  Length  W2   H2', 
     &' Rc[m]  Hshift Vshift   type    nmulti gmulti'/
     &I13,2I5,I7,2I5,I7,I6,I7,3X,A11,I4,F7.3/
     &'    Face Pol Wav[deg] tilt[deg] chamfer m   R1     Rm')
			DO iface=1,nfaces
				string=' no'
				IF (mnumber(i,iface).LT.0.) string='yes'
				IF (mnumber(i,iface).EQ.0.) THEN
					refl0(i,iface)=1.
					reflm(i,iface)=1.
				END IF
				IF (.NOT.nooutput) THEN
					IF (refl0(i,iface).LT.0.) THEN
						WRITE(*,130) iface, string(1:3), 
     &Wwav(i,iface)*180./pi, tilt(i,iface)*180./pi, 
     &chamfer(i,iface), reflfile(i,iface)(1:26)
  130					FORMAT(I7,2X,A3,3F7.3,1X,A26)
					ELSE
						WRITE(*,140) iface, string(1:3), 
     &Wwav(i,iface)*180./pi, tilt(i,iface)*180./pi, 
     &chamfer(i,iface), ABS(mnumber(i,iface)), 
     &refl0(i,iface), reflm(i,iface)
  140					FORMAT(I7,2X,A3,6F7.3)
					END IF
				END IF
			END DO
		END DO

		IF (.NOT.nooutput) THEN
			WRITE(*,150) NINT(Lend), NINT(Wmon), NINT(Hmon), 
     &NINT(Hmonshift), NINT(Vmonshift)
  150		FORMAT(' Distance between end of guide and mon,',
     &' mon width and height, Hor, Ver shift'/,5I6)
			WRITE(*,160) neutrons, iseed, iresults, idebug
  160		FORMAT(I10,' neutrons   iseed =',I10,
     &'     iresults =',I2,'  idebug =',I2)
			IF (gravity) THEN
				IF (sideways) THEN
					PRINT*,'   gravity ON - sideways'
				ELSE
					PRINT*,'   gravity ON'
				END IF
			ELSE
				IF (sideways) THEN
					PRINT*,'   gravity OFF - sideways'
				ELSE
					PRINT*,'   gravity OFF'
				END IF
			END IF
			WRITE(*,170) divlimit
  170		FORMAT(' Useful neutrons fall within +/-',F5.2,' degrees')
		END IF

c****************************************
c set up start parameters

		keepgoing=.FALSE.
		IF (neutrons.LT.0.) keepgoing=.TRUE.
		neutrons=ABS(neutrons)

		IF (iresults.EQ.-1) THEN
			saveneutrons=.TRUE.
			saveresults=.FALSE.
			saveplots=.FALSE.
		ELSE IF (iresults.EQ.0) THEN
			saveneutrons=.FALSE.
			saveresults=.FALSE.
			saveplots=.FALSE.
		ELSE IF (iresults.EQ.1) THEN
			saveneutrons=.FALSE.
			saveresults=.TRUE.
			saveplots=.FALSE.
		ELSE IF (iresults.EQ.2) THEN
			saveneutrons=.FALSE.
			saveresults=.TRUE.
			saveplots=.TRUE.
		ELSE
			PRINT*,'Error - invalid value for iresults =',iresults
		END IF

		indebug=0
c		indebug=223
		IF (idebug.LE.-1) THEN
			show=.FALSE.
			pn=.FALSE.
			printrays=.FALSE.
			batchjob=.TRUE.
			showplots=.FALSE.
		ELSE IF (idebug.EQ.0) THEN
			show=.FALSE.
			pn=.FALSE.
			printrays=.FALSE.
			batchjob=.FALSE.
			showplots=.TRUE.
		ELSE IF (idebug.EQ.1) THEN
			show=.TRUE.
			pn=.TRUE.
			printrays=.FALSE.
			batchjob=.FALSE.
			showplots=.TRUE.
		ELSE IF (idebug.EQ.2) THEN
			show=.FALSE.
			pn=.TRUE.
			printrays=.TRUE.
			batchjob=.FALSE.
			showplots=.FALSE.
		ELSE	! set idebug=0 until neutron number "indebug"
			show=.FALSE.
			pn=.FALSE.
			printrays=.FALSE.
			batchjob=.FALSE.
			showplots=.TRUE.
			indebug=idebug
		END IF

		nvars=2*nsections+nextravars
		DO i=1,2*nsections
			iusevars(i)=i
		END DO
		DO i=2*nsections+1,2*nsections+nextravars
			iusevars(i)=i+2*(msections-nsections)
		END DO

c		IF (iseed.EQ.0) THEN
c			itime=TIME()
c			iseed=4515415+itime*2
c			OPEN(16,FILE='iseed.dat',STATUS='UNKNOWN')
c			WRITE(16,*) iseed
c			CLOSE(16)
c		END IF
		CALL initialise_random(iseed)
c		CALL random_seed()

		IF (printrays) neutrons=20

		Lexit=0.
		DO i=1,nsections
			Lexit=Lexit+dLguide(i)+Lguide(i)
		END DO
		Lmon=Lexit+Lend		! total distance from source to monitor

		npol=1
		IF (polarised) npol=2

c****************************************
c Set up guide parameters

		xstart0=0.
		ystart0=0.
		zstart0=0.
		phiH=0.
		phiV=0.
		DO i=1,nsections
			IF (phiH.NE.0. .AND. phiV.NE.0.) THEN
				PRINT*,'Error - simultaneous hor and ver tilt'
				GOTO 1000
			END IF
			W1=Wguide(i,1)
			W2=Wguide(i,2)
			R=Rguide0(i,1)
			L=Lguide(i)
			dL=dLguide(i)
			H1=Hguide(i,1)
			H2=Hguide(i,2)
			cylindrical(i)=.FALSE.
			IF (H1.LT.0.) THEN	! cylindrical guide
				W2=W1
				H1=W1
				H2=W2
				cylindrical(i)=.TRUE.
				IF (ABS(R).GT.1) THEN
					PRINT*,'Error - cannot produce curved cylindrical guide'
					GOTO 1000
				END IF
			END IF
			xstart0=xstart0+dL*COS(phiH)*COS(phiV)
     &-Hshift(i)*SIN(phiH)-Vshift(i)*SIN(phiV)
			ystart0=ystart0+dL*SIN(phiH)+Hshift(i)*COS(phiH)
			zstart0=zstart0+dL*SIN(phiV)+Vshift(i)*COS(phiV)
			xstart(i,1)=xstart0+W1*SIN(phiH)/2.
			xstart(i,2)=xstart0-W1*SIN(phiH)/2.
			xstart(i,3)=xstart0+H1*SIN(phiV)/2.
			xstart(i,4)=xstart0-H1*SIN(phiV)/2.
			ystart(i,1)=ystart0-W1*COS(phiH)/2.
			ystart(i,2)=ystart0+W1*COS(phiH)/2.
			zstart(i,3)=zstart0-H1*COS(phiV)/2.
			zstart(i,4)=zstart0+H1*COS(phiV)/2.
			phig(i,1)=phiH+tilt(i,1)
			phig(i,2)=phiH+tilt(i,2)
			phig(i,3)=phiV+tilt(i,3)
			phig(i,4)=phiV+tilt(i,4)
			tiltH0=(tilt(i,1)+tilt(i,2))/2.
			tiltV0=(tilt(i,3)+tilt(i,4))/2.
			phiguideH(i)=phiH+tiltH0
			phiguideV(i)=phiV+tiltV0
			IF (.NOT.curved(i)) THEN	! not curved guide
				xend(i,1)=xstart(i,1)+L*COS(tiltH0)*COS(phiV)
				xend(i,2)=xstart(i,2)+L*COS(tiltH0)*COS(phiV)
				xend(i,3)=xstart(i,3)+L*COS(tiltV0)*COS(phiH)
				xend(i,4)=xstart(i,4)+L*COS(tiltV0)*COS(phiH)
				yend(i,1)=ystart(i,1)+L*SIN(phig(i,1))
				yend(i,2)=ystart(i,2)+L*SIN(phig(i,2))
				IF (cylindrical(i)) THEN	! cylindrical guide
					zend(i,3)=zstart(i,3)+L*SIN(phig(i,3))
					zend(i,4)=zstart(i,4)+L*SIN(phig(i,4))
					zend0=(zend(i,3)+zend(i,4))/2.
					ag(i,3)=TAN(phig(i,3))
					bg(i,3)=zstart(i,3)-ag(i,3)*xstart0
					ag(i,4)=TAN(phig(i,4))
					bg(i,4)=zstart(i,4)-ag(i,4)*xstart0
					cw(i,1)=(xstart(i,2)-xstart(i,1))/(ystart(i,2)-ystart(i,1))
					dw(i,1)=xstart(i,1)-cw(i,1)*ystart(i,1)
					cw(i,2)=(xend(i,2)-xend(i,1))/(yend(i,2)-yend(i,1))
					dw(i,2)=xend(i,1)-cw(i,2)*yend(i,1)
					xc1(i)=xstart0
					yc1(i)=ystart0
					xc2(i)=xend0
					yc2(i)=yend0
					bg(i,1)=ystart0-ag(i,1)*xstart0
					ag(i,2)=W1/2.
				ELSE						! straight guide or cavity
					zend(i,3)=zstart(i,3)+L*SIN(phig(i,3))
					zend(i,4)=zstart(i,4)+L*SIN(phig(i,4))
					xend0=(xend(i,1)+xend(i,2)+xend(i,3)+xend(i,4))/4.
					yend0=(yend(i,1)+yend(i,2))/2.
					zend0=(zend(i,3)+zend(i,4))/2.
					ag(i,1)=TAN(phig(i,1))
					bg(i,1)=ystart(i,1)-ag(i,1)*xstart(i,1)
					ag(i,2)=TAN(phig(i,2))
					bg(i,2)=ystart(i,2)-ag(i,2)*xstart(i,2)
					ag(i,3)=TAN(phig(i,3))
					bg(i,3)=zstart(i,3)-ag(i,3)*xstart0
					ag(i,4)=TAN(phig(i,4))
					bg(i,4)=zstart(i,4)-ag(i,4)*xstart0
					cw(i,1)=(xstart(i,2)-xstart(i,1))/(ystart(i,2)-ystart(i,1))
					dw(i,1)=xstart(i,1)-cw(i,1)*ystart(i,1)
					cw(i,2)=(xend(i,2)-xend(i,1))/(yend(i,2)-yend(i,1))
					dw(i,2)=xend(i,1)-cw(i,2)*yend(i,1)
					IF (cavity(i)) THEN
						IF (ivcavity(i).EQ.-2) THEN	! V-shape with apex at exit
							xstart(i,5)=xstart(i,1)
							xstart(i,6)=xstart(i,2)
							ystart(i,5)=ystart(i,1)
							ystart(i,6)=ystart(i,2)
							xend(i,5)=(xend(i,1)+xend(i,2))/2.
							xend(i,6)=xend(i,5)
							yend(i,5)=(yend(i,1)+yend(i,2))/2.
							yend(i,6)=yend(i,5)
							ag(i,5)=(yend(i,5)-ystart(i,5))
     &/(xend(i,5)-xstart(i,5))
							ag(i,6)=(yend(i,6)-ystart(i,6))
     &/(xend(i,6)-xstart(i,6))
							bg(i,5)=ystart(i,5)-ag(i,5)*xstart(i,5)
							bg(i,6)=ystart(i,6)-ag(i,6)*xstart(i,6)
							tilt(i,5)=ATAN(W1/(2.*L))
							tilt(i,6)=-ATAN(W1/(2.*L))
							phig(i,5)=(phig(i,1)+phig(i,2))/2.+tilt(i,5)
							phig(i,6)=(phig(i,1)+phig(i,2))/2.+tilt(i,6)
						ELSE IF (ivcavity(i).EQ.-1) THEN	! single blade, starting from face 2
							xstart(i,5)=xstart(i,2)
							ystart(i,5)=ystart(i,2)
							xend(i,5)=xend(i,1)
							yend(i,5)=yend(i,1)
							ag(i,5)=(yend(i,5)-ystart(i,5))
     &/(xend(i,5)-xstart(i,5))
							bg(i,5)=ystart(i,5)-ag(i,5)*xstart(i,5)
							tilt(i,5)=-ATAN((W1+W2)/(2.*L))
							phig(i,5)=(phig(i,1)+phig(i,2))/2.+tilt(i,5)
						ELSE IF (ivcavity(i).EQ.1) THEN	! single blade, starting from face 1
							xstart(i,5)=xstart(i,1)
							ystart(i,5)=ystart(i,1)
							xend(i,5)=xend(i,2)
							yend(i,5)=yend(i,2)
							ag(i,5)=(yend(i,5)-ystart(i,5))
     &/(xend(i,5)-xstart(i,5))
							bg(i,5)=ystart(i,5)-ag(i,5)*xstart(i,5)
							tilt(i,5)=ATAN((W1+W2)/(2.*L))
							phig(i,5)=(phig(i,1)+phig(i,2))/2.+tilt(i,5)
						ELSE IF (ivcavity(i).EQ.2) THEN	! V-shape with apex at entrance
							xstart(i,5)=(xstart(i,1)+xstart(i,2))/2.
							xstart(i,6)=xstart(i,5)
							ystart(i,5)=(ystart(i,1)+ystart(i,2))/2.
							ystart(i,6)=ystart(i,5)
							xend(i,5)=xend(i,1)
							xend(i,6)=xend(i,2)
							yend(i,5)=yend(i,1)
							yend(i,6)=yend(i,2)
							ag(i,5)=(yend(i,5)-ystart(i,5))
     &/(xend(i,5)-xstart(i,5))
							ag(i,6)=(yend(i,6)-ystart(i,6))
     &/(xend(i,6)-xstart(i,6))
							bg(i,5)=ystart(i,5)-ag(i,5)*xstart(i,5)
							bg(i,6)=ystart(i,6)-ag(i,6)*xstart(i,6)
							tilt(i,5)=-ATAN(W2/(2.*L))
							tilt(i,6)=ATAN(W2/(2.*L))
							phig(i,5)=(phig(i,1)+phig(i,2))/2.+tilt(i,5)
							phig(i,6)=(phig(i,1)+phig(i,2))/2.+tilt(i,6)
						END IF
					END IF
				END IF
				phiH=phiH+tiltH0
				phiV=phiV+tiltV0
				IF (hparabolic(i)) THEN
					a=Lguide(i)/((Wguide(i,1)/2.)**2-(Wguide(i,2)/2.)**2)
					xfocus(i)=a*(Wguide(i,2)/2.)**2+xend0-1./(4.*a)
					apara(i)=a
				END IF
				IF (vparabolic(i)) THEN
					a=Lguide(i)/((Hguide(i,1)/2.)**2-(Hguide(i,2)/2.)**2)
					xfocus(i)=a*(Hguide(i,2)/2.)**2+xend0-1./(4.*a)
					apara(i)=a
				END IF
			ELSE		! curved guide
				phi=phiH
				tilt0=tiltH0
				d1=W1/2.
				d2=W2/2.
				SINphi1=SIN(phi+tilt(i,1))
				COSphi1=COS(phi+tilt(i,1))
				SINphi2=SIN(phi+tilt(i,2))
				COSphi2=COS(phi+tilt(i,2))
				alpha=L/R	! angle subtended by centre-line of curved section
				xend0=xstart0+R*(SIN(phi+tilt0+alpha)-SIN(phi+tilt0))
				dxdy=SIN(phi+tilt0+alpha)/COS(phi+tilt0+alpha)
				yend0=ystart0+R*COS(phi+tilt0)-R/SQRT(dxdy**2+1.)
				xend(i,1)=xend0+d2*SIN(phi+tilt0+alpha)
				xend(i,2)=xend0-d2*SIN(phi+tilt0+alpha)
				yend(i,1)=yend0-d2*COS(phi+tilt0+alpha)
				yend(i,2)=yend0+d2*COS(phi+tilt0+alpha)
				IF (W1.EQ.W2) THEN
					Rguide(i,1)=R
					xc1(i)=xstart0-R*SINphi1+d1*SINphi1
					yc1(i)=ystart0+R*COSphi1-d1*COSphi1
					xc2(i)=xstart0-R*SINphi2-d1*SINphi2
					yc2(i)=ystart0+R*COSphi2+d1*COSphi2
					Rguide(i,2)=R
				ELSE	! converging curved guide
					TANphi=TAN(phi)
					xes1=xend(i,1)-xstart(i,1)
					yes1=yend(i,1)-ystart(i,1)
					xes2=xend(i,2)-xstart(i,2)
					yes2=yend(i,2)-ystart(i,2)
					ycs1=-0.5*(xes1**2+yes1**2)/(xes1*TANphi-yes1)
					yc1(i)=ystart(i,1)+ycs1
					xc1(i)=xstart(i,1)-TANphi*ycs1
					Rguide(i,1)=ycs1*SQRT(1.+TANphi**2)
					ycs2=-0.5*(xes2**2+yes2**2)/(xes2*TANphi-yes2)
					yc2(i)=ystart(i,2)+ycs2
					xc2(i)=xstart(i,2)-TANphi*ycs2
					Rguide(i,2)=ycs2*SQRT(1.+TANphi**2)
				END IF	
				zend(i,3)=zstart(i,3)+L*SIN(phig(i,3))
				zend(i,4)=zstart(i,4)+L*SIN(phig(i,4))
				zend0=(zend(i,3)+zend(i,4))/2.
				ag(i,3)=TAN(phig(i,3))
				bg(i,3)=zstart(i,3)-ag(i,3)*xstart0
				ag(i,4)=TAN(phig(i,4))
				bg(i,4)=zstart(i,4)-ag(i,4)*xstart0
				cw(i,1)=(xstart(i,2)-xstart(i,1))/(ystart(i,2)-ystart(i,1))
				dw(i,1)=xstart(i,1)-cw(i,1)*ystart(i,1)
				cw(i,2)=(xend(i,2)-xend(i,1))/(yend(i,2)-yend(i,1))
				dw(i,2)=xend(i,1)-cw(i,2)*yend(i,1)
				phiH=phi+tilt0+alpha
			END IF
			IF (show) THEN
c				WRITE(*,270) i, xstart0, xend0, ystart0, yend0, phi*180./pi
  270			FORMAT(' Section',I2,':',F7.1,'<x<',F8.1,',',F6.1,'<y',
     &F6.1,' phi(end) =',F6.3,' deg')
c				WRITE(*,275) (j, phig(i,j)*180./pi, ystart(i,j), yend(i,j), j=1,2)
  275			FORMAT(' Face',I2,': angle =',F7.4,'  y =',
     &F7.2,' ->',F7.2,' mm')
c				WRITE(*,280) (j, phig(i,j)*180./pi, zstart(i,j), zend(i,j), j=3,4)
  280			FORMAT(' Face',I2,': angle =',F7.4,'  z =',
     &F7.2,' ->',F7.2,' mm')
			END IF
			xstart0=xend0
			ystart0=yend0
			zstart0=zend0
		END DO
		xmon0=xend0+Lend*COS(phiH)*COS(phiV)
		ymon0=yend0+Lend*SIN(phiH)
		zmon0=zend0+Lend*SIN(phiV)
		xmon1=xmon0+( Wmon/2.-Hmonshift)*SIN(phiH)
     &+( Hmon/2.-Vmonshift)*SIN(phiV)
		xmon2=xmon0+(-Wmon/2.-Hmonshift)*SIN(phiH)
     &+(-Hmon/2.-Hmonshift)*SIN(phiV)
		ymon1=ymon0+(-Wmon/2.+Hmonshift)*COS(phiH)
		ymon2=ymon0+( Wmon/2.+Hmonshift)*COS(phiH)
		ymoncentre=(ymon1+ymon2)/2.
		zmon1=zmon0+(-Hmon/2.+Vmonshift)*COS(phiV)
		zmon2=zmon0+( Hmon/2.+Vmonshift)*COS(phiV)
		zmoncentre=(zmon1+zmon2)/2.
		cw(nsections+1,1)=(xmon2-xmon1)/(ymon2-ymon1)
		dw(nsections+1,1)=xmon1-cw(nsections+1,1)*ymon1

		phiHexit0=phiH	! \ 
		phiVexit0=phiV	!  | parameters at guide exit
		yexit0=yend0	!  |
		zexit0=zend0	! /
		IF (.NOT.batchjob) WRITE(*,290) phiHexit0*180./pi, 
     &phiVexit0*180./pi, yexit0, zexit0, ymon0, zmon0
  290	FORMAT('                           Horizontal  Vertical'/
     &' Direction of guide exit:  ',F7.2,F10.2,' deg'/
     &' centre of guide exit:     ',F7.2,F10.2,' mm'/
     &' centre of monitor:        ',F7.2,F10.2,' mm')

		gy=0.
		gz=0.
		IF (gravity) THEN
			IF (sideways) THEN
				gy=g
			ELSE
				gz=g
			END IF
		END IF

c****************************************
c Plot guide outline

		IF (.NOT.batchjob) THEN
			PRINT*,'Hit return to plot guide outline'
			READ(*,'(A1)') Answer
			CALL pgask(.FALSE.)
			IF (printrays) THEN
				CALL pgbeg(0,'guide_rays.ps/ps',1,2)
			ELSE
				CALL pgbeg(0,'/xw',1,2)
			END IF

			xmin0=-0.1*Lmon
			xmax0=1.1*Lmon
			ymin0=-Wmod/2.
			ymax0= Wmod/2.
			DO i=1,nsections
				ymin0=MIN(ymin0,ystart(i,1))
				ymin0=MIN(ymin0,yend(i,1))
				ymax0=MAX(ymax0,ystart(i,2))
				ymax0=MAX(ymax0,yend(i,2))
			END DO
c			ymin0=ymin0-50.
c			ymax0=ymax0+50.
			ymin0=ymin0-1.
			ymax0=ymax0+1.
			ymin0=MIN(ymin0,-1.1*Hmod/2.)
			ymax0=MAX(ymax0,1.1*Hmod/2.)
			CALL pgpanl(1,1)
			CALL pgvstd
			CALL pgswin(xmin0,xmax0,ymin0,ymax0)
			CALL PGBOX ('BCTN', 0.0, 0, 'BCNST', 0.0, 0)
			CALL pgsci(2)	! plot source in red
			CALL pgsch(0.)	! set arrowhead size to zero
			IF (usesavedneutrons) THEN
				CALL pgsls(2)	! dashed line
				CALL pgarro(0.,ymin0,0.,ymax0)
				CALL pgsls(1)
			ELSE
				CALL pgslw(4)	! set line width to thick
				CALL pgarro(0.,REAL(-Wmod/2.),0.,REAL(Wmod/2.))
				CALL pgslw(1)
			END IF
			CALL pgsci(1)	! set colour back to normal
			CALL pgsch(1.)	! set character size back to normal

c			zmin0=-(ymax0-ymin0)/2.
c			zmax0= (ymax0-ymin0)/2.
			zmin0=-Hmod/2.
			zmax0= Hmod/2.
			DO i=1,nsections
				zmin0=MIN(zmin0,zstart(i,3))
				zmin0=MIN(zmin0,zend(i,3))
				zmax0=MAX(zmax0,zstart(i,4))
				zmax0=MAX(zmax0,zend(i,4))
			END DO
			zmin0=zmin0-50.
			zmax0=zmax0+50.
			zmin0=MIN(zmin0,-1.1*Hmod/2.)
			zmax0=MAX(zmax0,1.1*Hmod/2.)
			CALL pgpanl(1,2)
			CALL pgvstd
			CALL pgswin(xmin0,xmax0,zmin0,zmax0)
			CALL PGBOX ('BCTN', 0.0, 0, 'BCNST', 0.0, 0)
			CALL pgsci(2)	! plot source in red
			CALL pgsch(0.)	! set arrowhead size to zero
			IF (usesavedneutrons) THEN
				CALL pgsls(2)	! dashed line
				CALL pgarro(0.,zmin0,0.,zmax0)
				CALL pgsls(1)
			ELSE
				CALL pgslw(4)	! set line width to thick
				CALL pgarro(0.,REAL(-Hmod/2.),0.,REAL(Hmod/2.))
				CALL pgslw(1)
			END IF
			CALL pgsci(1)	! set colour back to normal
			DO i=1,nsections
				offsetmulti=0.
				IF (2*(nmultiple(i)/2).EQ.nmultiple(i)) offsetmulti=-0.5
				IF (nmultiple(i).EQ.1) CALL pgslw(4)	! set line width to thick
				CALL pgpanl(1,1)	! top view
				CALL pgswin(xmin0,xmax0,ymin0,ymax0)
				IF (curved(i)) THEN	! curved guide
					dx=(xend(i,1)-xstart(i,1))/100.
					SINphi=SIN(phiguideH(i))
					COSphi=COS(phiguideH(i))
					DO imulti=1,nmultiple(i)
						R=Rguide(i,1)
						dir=SIGN(DBLE(1),R)
						rm=REAL(imulti-(nmultiple(i)+1)/2)+offsetmulti
						xm=-rm*(Wguide(i,1)+gmultiple(i))*SINphi
						ym=rm*(Wguide(i,1)+gmultiple(i))*COSphi
						DO j=1,101
							x=xstart(i,1)+REAL(j-1)*dx
							y=yc1(i)-dir*SQRT(R**2-(x-xc1(i))**2)
							xp(j)=x+xm
							yp(j)=y+ym
						END DO
						CALL pgline(101,xp,yp)
						R=Rguide(i,2)
						dir=SIGN(DBLE(1),R)
						dx=(xend(i,2)-xstart(i,2))/100.
						DO j=1,101
							x=xstart(i,2)+REAL(j-1)*dx
							y=yc2(i)-dir*SQRT(R**2-(x-xc2(i))**2)
							xp(j)=x+xm
							yp(j)=y+ym
						END DO
						CALL pgline(101,xp,yp)
					END DO
				ELSE IF (hparabolic(i)) THEN	! parabolic guide
					a=apara(i)
					x0=xfocus(i)
					dx=(xend(i,1)-xstart(i,1))/100.
					SINphi=SIN(phiguideH(i))
					COSphi=COS(phiguideH(i))
					TANphi=TAN(phiguideH(i))
					DO imulti=1,nmultiple(i)
						rm=REAL(imulti-(nmultiple(i)+1)/2)+offsetmulti
						xm=-rm*(Wguide(i,1)+gmultiple(i))*SINphi
						ym=rm*(Wguide(i,1)+gmultiple(i))*COSphi
						DO j=1,101
							x=xstart(i,1)+REAL(j-1)*dx
							y=-SQRT(1./(4.*a*a)+(x0-x)/a)+ystart(i,1)
     &+Wguide(i,1)*COSphi/2.+(x-xstart(i,1))*TANphi
							xp(j)=x+xm
							yp(j)=y+ym
						END DO
						CALL pgline(101,xp,yp)
						dx=(xend(i,2)-xstart(i,2))/100.
						DO j=1,101
							x=xstart(i,2)+REAL(j-1)*dx
							y= SQRT(1./(4.*a*a)+(x0-x)/a)+ystart(i,2)
     &-Wguide(i,1)*COSphi/2.+(x-xstart(i,2))*TANphi
							xp(j)=x+xm
							yp(j)=y+ym
						END DO
						CALL pgline(101,xp,yp)
					END DO
				ELSE					! not curved or parabolic
					SINphi=SIN(phiguideH(i))
					COSphi=COS(phiguideH(i))
					DO imulti=1,nmultiple(i)
						rm=REAL(imulti-(nmultiple(i)+1)/2)+offsetmulti
						xm=-rm*(Wguide(i,1)+gmultiple(i))*SINphi
						ym=rm*(Wguide(i,1)+gmultiple(i))*COSphi
						xmin=xstart(i,1)+xm
						xmax=xend(i,1)+xm
						ymin=ystart(i,1)+ym
						ymax=yend(i,1)+ym
						CALL pgarro(xmin,ymin,xmax,ymax)
						xmin=xstart(i,2)+xm
						xmax=xend(i,2)+xm
						ymin=ystart(i,2)+ym
						ymax=yend(i,2)+ym
						CALL pgarro(xmin,ymin,xmax,ymax)
						IF (cavity(i)) THEN
							CALL pgsci(3)	! set colour to special
							CALL pgslw(1)	! set line width to thin
							DO iface=1,ABS(ivcavity(i))
								xmin=xstart(i,4+iface)+xm
								xmax=xend(i,4+iface)+xm
								ymin=ystart(i,4+iface)+ym
								ymax=yend(i,4+iface)+ym
								CALL pgarro(xmin,ymin,xmax,ymax)
							END DO
							CALL pgsci(1)	! set colour back to normal
						END IF
					END DO
				END IF
				CALL pgpanl(1,2)	! side view
				CALL pgswin(xmin0,xmax0,zmin0,zmax0)
				CALL pgslw(4)	! set line width to thick
				IF (vparabolic(i)) THEN	! parabolic guide
					a=apara(i)
					x0=xfocus(i)
					dx=(xend(i,1)-xstart(i,1))/100.
					SINphi=SIN(phiguideV(i))
					COSphi=COS(phiguideV(i))
					TANphi=TAN(phiguideV(i))
					DO imulti=1,nmultiple(i)
						rm=REAL(imulti-(nmultiple(i)+1)/2)+offsetmulti
						xm=-rm*(Hguide(i,1)+gmultiple(i))*SINphi
						zm=rm*(Hguide(i,1)+gmultiple(i))*COSphi
						DO j=1,101
							x=xstart(i,1)+REAL(j-1)*dx
							z=-SQRT(1./(4.*a*a)+(x0-x)/a)+zstart(i,3)
     &+Hguide(i,1)*COSphi/2.+(x-xstart(i,1))*TANphi
							xp(j)=x+xm
							zp(j)=z+zm
						END DO
						CALL pgline(101,xp,zp)
						dx=(xend(i,2)-xstart(i,2))/100.
						DO j=1,101
							x=xstart(i,2)+REAL(j-1)*dx
							z= SQRT(1./(4.*a*a)+(x0-x)/a)+zstart(i,4)
     &-Hguide(i,1)*COSphi/2.+(x-xstart(i,2))*TANphi
							xp(j)=x+xm
							zp(j)=z+zm
						END DO
						CALL pgline(101,xp,zp)
					END DO
				ELSE
					xmin=(xstart(i,1)+xstart(i,2))/2.
					xmax=(xend(i,1)+xend(i,2))/2.
					zmin=zstart(i,3)
					zmax=zend(i,3)
					CALL pgarro(xmin,zmin,xmax,zmax)
					zmin=zstart(i,4)
					zmax=zend(i,4)
					CALL pgarro(xmin,zmin,xmax,zmax)
				END IF
				CALL pgslw(1)	! set line width to thin
				CALL pgpanl(1,1)
				CALL pgswin(xmin0,xmax0,ymin0,ymax0)
				IF (nmultiple(i).EQ.1) THEN
					xmin=xstart(i,1)
					xmax=xstart(i,2)
					ymin=ystart(i,1)
					ymax=ystart(i,2)
					CALL pgarro(xmin,ymin,xmax,ymax)
					xmin=xend(i,1)
					xmax=xend(i,2)
					ymin=yend(i,1)
					ymax=yend(i,2)
					CALL pgarro(xmin,ymin,xmax,ymax)
				END IF
				CALL pgpanl(1,2)
				CALL pgswin(xmin0,xmax0,zmin0,zmax0)
				xmin=(xstart(i,1)+xstart(i,2))/2.
				xmax=xmin
				ymin=zstart(i,3)
				ymax=zstart(i,4)
				CALL pgarro(xmin,ymin,xmax,ymax)
				xmin=(xend(i,1)+xend(i,2))/2.
				xmax=xmin
				ymin=zend(i,3)
				ymax=zend(i,4)
				CALL pgarro(xmin,ymin,xmax,ymax)
			END DO

			CALL pgsci(2)	! plot monitor in red
			CALL pgslw(4)	! set line width to thick
			CALL pgpanl(1,1)
			CALL pgswin(xmin0,xmax0,ymin0,ymax0)
			CALL pgarro(REAL(xmon1),REAL(ymon1),REAL(xmon2),REAL(ymon2))
			CALL pgpanl(1,2)
			CALL pgswin(xmin0,xmax0,zmin0,zmax0)
			CALL pgarro(REAL(xmon0),REAL(zmon1),REAL(xmon0),REAL(zmon2))
			CALL pgsci(1)
			CALL pgslw(1)
		END IF

		IF (pn) THEN
c			CALL pgsch(1.)	! set arrowhead size
			CALL pgslw(1)	! set line width to thin
c			CALL pgsah(1,180.,0.)
		END IF

c****************************************
c set up and zero arrays

  300	CONTINUE

		IF (looplambda) THEN
			ilambda=ilambda+1
			lambda=lambda_arr(ilambda)
			IF (.NOT.nooutput) WRITE(*,302) ilambda, lambda
  302		FORMAT(' Point',I4,': lambda =',F6.2,' A')
		END IF
		v=const/lambda	! in mm/ms
		iseed=2*iseed*NINT(100.*lambda)+1

		neutron0=0
		DO ipol=1,2
			DO isec=1,nsections
				nthrough(isec,ipol)=0
			END	DO
			nmon(ipol)=0
			ncentral(1,ipol)=0
			ncentral(2,ipol)=0
			ncentral(3,ipol)=0
		END DO
		L=Lmon
		IF (UCN) L=Lguide(nsections)
		DO i=1,nvars
			ivar=iusevars(i)
			step=dvar(ivar)
			DO ibin=-nbin,nbin
				var(ivar,ibin)=step*REAL(ibin)
				IF (ivar.EQ.2*msections+3)	! horizontal divergence at monitor
     &var(ivar,ibin)=var(ivar,ibin)+DBLE(NINT(phiHexit0*180./pi))
				IF (ivar.EQ.2*msections+4)	! vertical divergence at monitor
     &var(ivar,ibin)=var(ivar,ibin)+DBLE(NINT(phiVexit0*180./pi))
				IF (ivar.EQ.2*msections+5)	! y-coordinate at monitor
     &var(ivar,ibin)=var(ivar,ibin)+DBLE(NINT(yexit0))
				IF (ivar.EQ.2*msections+6)	! z-coordinate at monitor
     &var(ivar,ibin)=var(ivar,ibin)+DBLE(NINT(zexit0))
				IF (ivar.EQ.2*msections+7)	! total flight path
     &var(ivar,ibin)=var(ivar,ibin)+DBLE(NINT(L))
				IF (ivar.GE.2*msections+9)	! neutron losses
     &var(ivar,ibin)=var(ivar,ibin)+step*REAL(nbin)
				binvar(ivar,1,ibin)=0.
				binvar(ivar,2,ibin)=0.
			END DO
c			IF (show) WRITE(*,305) vblurb(ivar), var(ivar,-nbin), var(ivar,nbin)
  305		FORMAT(A50,' :',F9.2,' < x <',F9.2)
		END DO

		IF (UCN) THEN
			DO i=1,2
				nUCN(i)=0
				PsumUCN(i)=0.
			END DO
		END IF

		nstartlastsection=0

c*************************************************
c Check and open source input file

		IF (usesavedneutrons) THEN
			IF (looplambda) THEN
				PRINT*,'1:',source_file
				CALL makefilename(source_file,lambda,filename,n)
				PRINT*,'2:',source_file
				filename(1:n+4)=filename(1:n)//'.dat'
				IF (.NOT.nooutput) PRINT*,' From source file ',filename
				OPEN(21,FILE=filename,STATUS='OLD')
				READ(21,'(A50)') string
				READ(21,*) rl, polarised_in
				IF (ABS(rl-lambda).GT.0.01) PRINT*,
     &'Error - wavelength in source file does not match file name'
			ELSE
				filename=source_file
				IF (.NOT.nooutput) PRINT*,' From source file ',filename
				OPEN(21,FILE=filename,STATUS='OLD')
				READ(21,'(A50)') string
				READ(21,*) rl, polarised_in
c---------------------------------------------------------
c				polarised_in=.FALSE.	! pre-2008 format
c				READ(21,*) rl
c---------------------------------------------------------
				lambda=DBLE(rl)
			END IF
			IF (lambda.EQ.0.) wrange=.TRUE.
			neutrons0=1	! set to 0 to check number of neutrons in each file
			IF (neutrons0.EQ.0) THEN
				neutrons1=0
				neutrons2=0
  220			IF (polarised_in) THEN
					IF (wrange) THEN
						READ(21,*,END=222) lambda, y, z, phiy, phiz, ipol
					ELSE
						READ(21,*,END=222) y, z, phiy, phiz, ipol
					END IF
					IF (ipol.EQ.1) neutrons1=neutrons1+1
					IF (ipol.EQ.2) neutrons2=neutrons2+1
				ELSE
					IF (wrange) THEN
						READ(21,*,END=222) lambda, y, z, phiy, phiz
					ELSE
						READ(21,*,END=222) y, z, phiy, phiz
					END IF
					neutrons0=neutrons0+1
				END IF
				GOTO 220
  222			REWIND(21)
				READ(21,*)
				READ(21,*)
				IF (.NOT.nooutput) THEN
					IF (neutrons1.EQ.0 .AND. neutrons2.EQ.0) THEN
						WRITE(*,226) lambda, neutrons0, 
     &MIN(neutrons,neutrons0)
  226					FORMAT('    lambda =',F7.3,' A'/
     &I10,' un-polarised neutrons saved'/
     &I10,' neutrons will be used in the present simulation')
						neutrons=MIN(neutrons,neutrons0)
					ELSE
						WRITE(*,227) lambda, neutrons1, neutrons2, 
     &MIN(neutrons,neutrons1+neutrons2)
  227					FORMAT(12X,'    lambda =',F7.3,' A'/
     &I10,' /',I10,' up/down neutrons saved'/
     &I22,' neutrons will be used in the present simulation')
						neutrons=MIN(neutrons,neutrons1+neutrons2)
					END IF
				END IF
			ELSE
				IF (.NOT.nooutput) WRITE(*,231) filename, lambda, 
     &neutrons
  231			FORMAT('    input file = ',A50/
     &'    lambda =',F7.3,' A'/
     &I10,' neutrons will be used in the present simulation'/
     &10X,' - hopefully there is at least that number in the file')
			END IF
			IF (polarised_in) polarised=.TRUE.
		END IF

c*************************************************
c Set up output file names and trajectory file

		IF (saveneutrons) THEN
			IF (looplambda) THEN
				CALL makefilename(input_file,lambda,filename,n)
				out_file(1:n+4)=filename(1:n)//'.dat'
			ELSE
				DO i=1,50
					IF (input_file(i:i).EQ.'.') n=i-1
				END DO
				out_file(1:n+4)=input_file(1:n)//'.dat'
			END IF
			INQUIRE(FILE=out_file,EXIST=exist)
			IF (.NOT.exist) THEN
				OPEN(12,FILE=out_file,STATUS='NEW')
				WRITE(12,*) input_file
				WRITE(12,*) REAL(lambda), polarised
				neutron00=0
			ELSE
				IF (.NOT.nooutput) PRINT*,'opening ',out_file
				OPEN(12,FILE=out_file,STATUS='OLD')
				READ(12,*)
				READ(12,*)
				neutron00=0
				neutron01=0
				neutron02=0
  340			IF (polarised_in) THEN
					IF (wrange) THEN
						READ(12,*,END=350) rlambda, yyr, zzr, 
     &phiyyr, phizzr, ipol
					ELSE
						READ(12,*,END=350) yyr, zzr, phiyyr, phizzr, ipol
					END IF
					IF (ipol.EQ.1) neutron01=neutron01+1
					IF (ipol.EQ.2) neutron02=neutron02+1
				ELSE
					IF (wrange) THEN
						READ(12,*,END=350) rlambda, yyr, zzr, 
     &phiyyr, phizzr
					ELSE
						READ(12,*,END=350) yyr, zzr, phiyyr, phizzr
					END IF
					neutron00=neutron00+1
				END IF
				GOTO 340
  350			CLOSE(12)
				IF (.NOT.nooutput) THEN
					IF (polarised_in) THEN
						WRITE(*,360) neutron01, neutron02, out_file
  360					FORMAT(I10,' /',I10,
     &' up/down neutrons found in ',A20)
						neutron00=neutron01+neutron02
					ELSE
						WRITE(*,370) neutron00, out_file
  370					FORMAT(I10,
     &' unpolarised neutrons found in ',A20)
					END IF
					IF (neutron00.GE.neutrons) THEN
						PRINT*,'No need to recalculate at this wavelength'
						GOTO 300
					END IF
				END IF
				OPEN(12,FILE=out_file,STATUS='OLD',ACCESS='APPEND')
			END IF
		END IF

c*************************************************
c main loop

		IF (show) THEN
			WRITE(*,390)
  390		FORMAT(' Hit return to start main loop : ',$)
			READ(*,'(A1)') Answer
		END IF

		newneutron=.TRUE.

		DO ineutron=1,neutrons

			IF (ineutron.EQ.indebug) THEN
				show=.TRUE.
				pn=.TRUE.
			END IF

  400		CONTINUE
			IF (show) PRINT*,'********************************************'
			IF (show) PRINT*,'* Neutron',ineutron
			IF (show) PRINT*,'*'

			x=0.
			IF (usesavedneutrons) THEN
  401			IF (newneutron) THEN
  					IF (polarised_in) THEN
						IF (wrange) THEN
							READ(21,*) lambda0, y0, z0, phiy0, phiz0, ipol0
						ELSE
							READ(21,*) y0, z0, phiy0, phiz0, ipol0
						END IF
					ELSE
						ipol0=1
						IF (wrange) THEN
							READ(21,*) lambda0, y0, z0, phiy0, phiz0
						ELSE
							READ(21,*) y0, z0, phiy0, phiz0
						END IF
					END IF
					ipol=ipol0
				END IF
				IF (sideways) THEN
					y=z0
					z=y0
					phiy=phiz0
					phiz=phiy0
				ELSE
					y=y0
					z=z0
					phiy=phiy0
					phiz=phiz0
				END IF
				ay=TAN(phiy)
				az=TAN(phiz)
			ELSE
				IF (newneutron) THEN
					ipol=1
  409				IF (random(iseed).EQ.0.) GOTO 409
					IF (wrange) lambda0=lambdamin
     &+(lambdamax-lambdamin)*random(iseed)
  410				y0=(Wmod+2.*halo)*(random(iseed)-0.5)
					IF (circular_source) THEN
						z0=(Wmod+2.*halo)*(random(iseed)-0.5)
						r=SQRT(y0*y0+z0*z0)
						IF (r.LT.Wmod/2.) THEN
							weight=1.
						ELSE IF (r.LT.Wmod/2.+halo) THEN
							weight=1.-(r-Wmod/2.)/halo
						ELSE
							weight=0.
						END IF
						IF (random(iseed).GT.weight) GOTO 410
					ELSE
						z0=(Hmod+2.*halo)*(random(iseed)-0.5)
						IF (ABS(z0).LT.Hmod/2.) THEN
							weight=1.
						ELSE
							weight=1.-(ABS(z0)-Hmod/2.)/halo
						END IF
						IF (random(iseed).GT.weight) GOTO 410
						IF (cylindrical_source) THEN
							weight=2.*SQRT((Wmod/2.)**2-y0*y0)/Wmod
							IF (random(iseed).GT.weight) GOTO 410
						ELSE
							IF (ABS(y0).LT.Wmod/2.) THEN
								weight=1.
							ELSE
								weight=1.-(ABS(y0)-Wmod/2.)/halo
							END IF
							IF (random(iseed).GT.weight) GOTO 410
						END IF
					END IF
					IF (collimated_source) THEN
						IF (Hsourcem.LT.0.) THEN
  420						phiy0=10.*Hsourcem*(random(iseed)-0.5)
							gauss=EXP(-0.5*(phiy0/(Hsourcem/2.355))**2)
c							PRINT*,'Hsourcem<0: phiy=',phiy,' gauss=',gauss
							IF (gauss.LT.random(iseed)) GOTO 420
							phiy0=phiy0*pi/180.
						ELSE
							phiy0=Hsourcem*0.2*lambda
     &*(random(iseed)-0.5)*pi/180.
						END IF
						IF (Vsourcem.LT.0.) THEN
  430						phiz0=10.*Vsourcem*(random(iseed)-0.5)
							gauss=EXP(-0.5*(phiz0/(Vsourcem/2.355))**2)
							IF (gauss.LT.random(iseed)) GOTO 430
							phiz0=phiz0*pi/180.
						ELSE
							phiz0=Vsourcem*0.2*lambda
     &*(random(iseed)-0.5)*pi/180.
						END IF
						ay0=TAN(phiy0)
						az0=TAN(phiz0)
					ELSE
						x1=Lap
						y1=Wap*(random(iseed)-0.5)
						z1=Hap*(random(iseed)-0.5)
						ay0=(y1+0.5*gy*(x1/v)**2-y0)/x1
						az0=(z1+0.5*gz*(x1/v)**2-z0)/x1
						phiy0=ATAN(ay0)
						phiz0=ATAN(az0)
					END IF
				END IF
				y=y0
				z=z0
				phiy=phiy0
				phiz=phiz0
				ay=ay0
				az=az0
			END IF

			IF (wrange) THEN
				lambda=lambda0
				v=const/lambda	! in mm/ms
			END IF

			IF (pn) CALL pgsci(ipol)

			through=.TRUE.
			by=y
			bz=z
			L=0.
			nrefl=0

			DO isection=1,nsections
				IF (UCN.AND.isection.EQ.nsections) L=-1.
				IF (cylindrical(isection)) THEN
					intube=through
					CALL cylindrical_guide(lambda,x,y,z,ay,by,az,bz,
     &phiy,phiz,in,through,L,nrefl,xloss,
     &isection,UCN,intube,iseed,ipol)
				ELSE
					CALL guide_section(lambda,v,x,y,z,ay,by,gy,az,bz,gz,
     &in,through,L,nrefl,xloss,isection,iseed,
     &curved(isection),hparabolic(isection),
     &vparabolic(isection),ABS(ivcavity(isection)),ipol)
				END IF
				IF (UCN.AND.isection.EQ.nsections.AND.L.GE.0.) THEN
					PUCN=(sigmaU/sigmaS)*(1.-EXP(-nsigmaS*L))
					IF (show) WRITE(*,450) L, (1.-EXP(-nsigmaS*L))*100., 
     &PUCN*100.
  450				FORMAT(' L(in He) =',F7.2,' mm  =>  P(scatt) =',F7.4,
     &' & P(UCN) =',F7.4,' %')
					CALL bin(L,2*msections+7,ipol)		! flight path in tube
					nUCN(ipol)=nUCN(ipol)+1.
					PsumUCN(ipol)=PsumUCN(ipol)+PUCN
				END IF
				IF (xloss(1,ipol).GT.0.) THEN
					IF (show) WRITE(*,460) ipol, NINT(xloss(2,ipol)), 
     &xloss(1,ipol), 2*msections+9+NINT(xloss(2,ipol))
  460				FORMAT('Pol',I1,' Side',I1,': xloss=',F8.2,
     &'mm, ivar=',I3)
					CALL bin(xloss(1,ipol),
     &2*msections+9+NINT(xloss(2,ipol)),ipol) ! neutron loss position
				END IF
				IF (through) THEN
					nthrough(isection,ipol)=nthrough(isection,ipol)+1
					c=cw(isection,2)
					d=dw(isection,2)
					IF (gy.EQ.0.) THEN
						xexit=(by*c+d)/(1.-ay*c)
						yexit=ay*xexit+by
					ELSE
						F=-gy/(2.*v*v)
						IF (c.EQ.0.) THEN
							xexit=d
							yexit=F*xexit**2+ay*xexit+by
						ELSE
							F=-gy/(2.*v*v)
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
					zexit=az*xexit+bz-0.5*gz*(xexit/v)**2
					slopey=-gy*xexit/v**2+ay
					slopez=-gz*xexit/v**2+az
					phiy=ATAN(slopey)
					phiz=ATAN(slopez)
					CALL bin(phiy*180./pi,2*isection-1,ipol) ! hor div at guide exit
					CALL bin(phiz*180./pi,2*isection,ipol)   ! ver div at guide exit
				END IF
				IF (isection.EQ.nsections .AND. in) 
     &nstartlastsection=nstartlastsection+1
c				IF (ABS(phiy*180./pi).GT.10.) PRINT*,'phiy=',phiy,
c    &' ineutron=',ineutron
				IF (.NOT.through) GOTO 500
			END DO

			CALL bin(yexit,2*msections+1,ipol)	! hor pos
			CALL bin(zexit,2*msections+2,ipol)	! ver pos
			IF (nmultiple(nsections).EQ.1) THEN
				IF (ABS(yexit-yexit0).GT.1.1*Wguide(nsections,2)/2. .OR. 
     &ABS(zexit-zexit0).GT.1.1*Hguide(nsections,2)/2.) THEN
					WRITE(*,465) ineutron, yexit, yexit0, Wguide(nsections,2), 
     &zexit, zexit0, Hguide(nsections,2)
  465				FORMAT(' Warning - neutron',I10/
     &' y: neutron exit =',F7.2,' guide centre =',F7.2,' width =',F6.2/
     &' z: neutron exit =',F7.2,' guide centre =',F7.2,
     &' height=',F6.2/
     &' Hit return to continue',$)
c					READ(*,'(A1)') Answer
				END IF
			END IF

			hitmon=.FALSE.
			c=cw(nsections+1,1)
			d=dw(nsections+1,1)
			IF (gy.EQ.0.) THEN
				xmon=(by*c+d)/(1.-ay*c)
				ymon=ay*xmon+by
			ELSE
				F=-gy/(2.*v*v)
				IF (c.EQ.0.) THEN
					xmon=d
					ymon=F*xmon**2+ay*xmon+by
				ELSE
					AA=F*c*c
					BB=2.*F*c*d+ay*c-1.
					CC=F*d*d+ay*d+by
					y1=(-BB-SQRT(BB*BB-4.*AA*CC))/(2.*AA)
					y2=(-BB+SQRT(BB*BB-4.*AA*CC))/(2.*AA)
					ymon=y1	! choose solution with smallest ABS(y)
					IF (ABS(y2).LT.ABS(y1)) ymon=y2
					xmon=c*ymon+d
				END IF
			END IF
			zmon=az*xmon+bz-0.5*gz*(xmon/v)**2

			IF (ymon.GT.ymon1 .AND. ymon.LT.ymon2 .AND. 
     &zmon.GT.zmon1 .AND. zmon.LT.zmon2) hitmon=.TRUE.

			IF (show) WRITE(*,470) ineutron, ymon1, ymon2, zmon1, zmon2, 
     &yend(nsections,1), yend(nsections,2), 
     &zend(nsections,3), zend(nsections,4), 
     &xexit, yexit, zexit, xmon, ymon, zmon
  470		FORMAT(' Neutron',I6,' :'/
     &' monitor y-range :',F8.2,' ->',F8.2,' mm'/
     &' monitor z-range :',F8.2,' ->',F8.2,' mm'/
     &' Guide   y-range :',F8.2,' ->',F8.2,' mm'/
     &' Guide   z-range :',F8.2,' ->',F8.2,' mm'/
     &' Neutron at exit    (x,y,z) = (',F8.2,',',F8.2,',',F8.2,')'/
     &' Neutron at monitor (x,y,z) = (',F8.2,',',F8.2,',',F8.2,')')

			IF (hitmon) THEN
				IF (show) PRINT*,'Neutron hits monitor!'
				IF (.NOT.UCN) L=L+SQRT((xmon-x)**2+(ymon-y)**2+(zmon-z)**2)
				slopey=-gy*xmon/v**2+ay
				slopez=-gz*xmon/v**2+az
				phiy=ATAN(slopey)
				phiz=ATAN(slopez)
				CALL bin(phiy*180./pi,2*msections+3,ipol)	! hor div on monitor
				CALL bin(phiz*180./pi,2*msections+4,ipol) 	! ver div on monitor
				IF (.NOT.UCN) CALL bin(L,2*msections+7,ipol)				! total flight path
				CALL bin(DBLE(nrefl),2*msections+8,ipol)		! total no. of reflections
				nmon(ipol)=nmon(ipol)+1
				IF (ABS(phiy-phiHexit0)*180./pi.LT.divlimit)			! within hor +/- divlimit degrees
     &ncentral(1,ipol)=ncentral(1,ipol)+1
				IF (ABS(phiz-phiVexit0)*180./pi.LT.divlimit)			! within ver +/- divlimit degrees
     &ncentral(2,ipol)=ncentral(2,ipol)+1
				IF (ABS(phiy-phiHexit0)*180./pi.LT.divlimit .AND. 
     &ABS(phiz-phiVexit0)*180./pi.LT.divlimit)			! within hor & ver +/- divlimit degrees
     &ncentral(3,ipol)=ncentral(3,ipol)+1
c------------------------------------------------------------------------
				IF (saveneutrons) THEN
					IF (sideways) THEN
c						yyr=REAL(zmon-zmon0)
c						zzr=REAL(ymon-ymon0)
						yyr=REAL(zmon-zmoncentre)
						zzr=REAL(ymon-ymoncentre)
						phiyyr=REAL(phiz-phiVexit0)
						phizzr=REAL(phiy-phiHexit0)
					ELSE
c						yyr=REAL(ymon-ymon0)
c						zzr=REAL(zmon-zmon0)
						yyr=REAL(ymon-ymoncentre)
						zzr=REAL(zmon-zmoncentre)
						phiyyr=REAL(phiy-phiHexit0)
						phizzr=REAL(phiz-phiVexit0)
					END IF
					IF (polarised) THEN
						IF (wrange) THEN
							WRITE(12,*) REAL(lambda), yyr, zzr, 
     &phiyyr, phizzr, ipol
						ELSE
							WRITE(12,*) yyr, zzr, phiyyr, phizzr, ipol
						END IF
					ELSE
						IF (wrange) THEN
							WRITE(12,*) REAL(lambda), yyr, zzr, 
     &phiyyr, phizzr
						ELSE
							WRITE(12,*) yyr, zzr, phiyyr, phizzr
						END IF
					END IF
c-------------------------------
c dodgy bodge for multi-slit guide for Phil
c					WRITE(12,*) REAL(zmon-zmon0), REAL(phiz-phiVexit0)
c-------------------------------
				END IF
c------------------------------------------------------------------------
			ELSE
				IF (show) PRINT*,'***************************'
				IF (show) PRINT*,'Neutron misses monitor'
				IF (show) PRINT*,'***************************'
			END IF

			IF (through) CALL bin(ymon,2*msections+5,ipol)		! hor pos at monitor
			IF (through) CALL bin(zmon,2*msections+6,ipol)		! ver pos at monitor

			IF (pn) THEN
				xold=x
				yold=y
				zold=z
				CALL pgpanl(1,1)
				CALL pgswin(xmin0,xmax0,ymin0,ymax0)
				IF (gy.EQ.0.) THEN
					CALL pgarro(REAL(xold),REAL(yold),
     &REAL(xmon),REAL(ymon))
					WRITE(*,491) ATAN((ymon-yold)/(xmon-xold))*180./pi
  491				FORMAT(' y angle from plot =',F7.3,' deg')
				ELSE
					dx=(xmon-xold)/100.
					DO j=1,101
						x=xold+REAL(j-1)*dx
						y=ay*x+by-0.5*gy*(x/v)**2
						xp(j)=x
						yp(j)=y
					END DO
					CALL pgline(101,xp,yp)
				END IF
				CALL pgpt1(REAL(xold),REAL(yold),17)
				CALL pgpt1(REAL(xmon),REAL(ymon),17)
				CALL pgpanl(1,2)
				CALL pgswin(xmin0,xmax0,zmin0,zmax0)
				IF (gz.EQ.0.) THEN
					CALL pgarro(REAL(xold),REAL(zold),
     &REAL(xmon),REAL(zmon))
					WRITE(*,492) ATAN((zmon-zold)/(xmon-xold))*180./pi
  492				FORMAT(' z angle from plot =',F7.3,' deg')
				ELSE
					dx=(xmon-xold)/100.
					DO j=1,101
						x=xold+REAL(j-1)*dx
						z=az*x+bz-0.5*gz*(x/v)**2
						xp(j)=x
						zp(j)=z
					END DO
					CALL pgline(101,xp,zp)
				END IF
				CALL pgpt1(REAL(xold),REAL(zold),17)
				CALL pgpt1(REAL(xmon),REAL(zmon),17)
			END IF

  500		CONTINUE

			IF (show) THEN
				WRITE(*,501)
  501			FORMAT(' Hit return for next neutron or C to stop output : ',$)
				READ(*,'(A1)') Answer
				IF (Answer.EQ.'c' .OR. Answer.EQ.'C') THEN
					show=.FALSE.
					pn=.FALSE.
				END IF
			END IF

			IF (.NOT.polarised_in .AND. polarised .AND. ipol.EQ.1) THEN
				ipol=2
				newneutron=.FALSE.
				GOTO 400
			ELSE
				newneutron=.TRUE.
			END IF

			IF (.NOT.nooutput) THEN
				DO i=1,10
					IF (ineutron.EQ.i*neutrons/10) THEN
						PRINT*, ineutron,' neutrons'
					END IF
				END DO
			END IF

			IF (keepgoing) THEN
				neutron0=neutron0+1
				IF (neutrons*(neutron0/neutrons).EQ.neutron0) 
     &PRINT*, neutron0,' neutrons:',nthrough(nsections,1),
     &' made it through'
				IF (nmon(1)+neutron00.LT.neutrons) THEN
					GOTO 400
				ELSE
					GOTO 600
				END IF
			END IF

		END DO	! main Monte Carlo loop

  600	CONTINUE
  
		IF (pn) CALL pgend

		IF (keepgoing) THEN
			n0=neutrons
			neutrons=neutron0
		END IF

		IF (usesavedneutrons) CLOSE(21)

		IF (saveneutrons) CLOSE(12)

		IF (wrange) lambda=0.

		IF (.NOT.polarised) THEN
			Tls=REAL(nthrough(nsections,1))/REAL(nstartlastsection)	! T through last section
		ELSE
			Tls=REAL(nthrough(nsections,1)+nthrough(nsections,2))
     &/REAL(nstartlastsection)	! T through last section
			DO i=1,2
				IF (i.EQ.1) THEN
					P1=nmon(1)
					P2=nmon(2)
				ELSE
					P1=ncentral(3,1)
					P2=ncentral(3,2)
				END IF
				rn=P1-P2
				rd=P1+P2
				IF (rd.EQ.0.) THEN
					P=0.
					dP=0.
				ELSE
					P=rn/rd
					dn=SQRT(ABS(P1-P2))
					dd=SQRT(P1+P2)
					dP=SQRT((dn/rd)**2+(rn*dd/rd**2)**2)
				END IF
				FR=(1.+P)/(1.-P)
				IF (P.LT.0.) FR=-1./FR
				IF (i.EQ.1) THEN
					Pa=P
					dPa=dP
					FRa=FR
				ELSE
					Pb=P
					dPb=dP
					FRb=FR
				END IF
			END DO
		END IF

		IF (.NOT.nooutput) THEN
			IF (scan) WRITE(*,619) scanvarname, scanvar
  619		FORMAT(' Scan variable: ',A10,' =',F9.3)
			IF (.NOT.polarised) THEN
				WRITE(*,630) neutrons, lambda
  630			FORMAT(' Out of',I10,' neutrons at lambda=',F7.3,'A :')
				WRITE(*,640) (nthrough(i,1), i, i=1,nsections)
  640			FORMAT(7X,I10,' made it through section',I2)
				WRITE(*,650) nmon(1), ncentral(1,1), divlimit, 
     &ncentral(2,1), divlimit, ncentral(3,1), divlimit
  650			FORMAT(7X,I10,' hit the monitor'/
     &7X,I10,' hit the monitor within +/-',F5.2,' deg horizontally'/
     &7X,I10,' hit the monitor within +/-',F5.2,' deg vertically'/
     &7X,I10,' hit the monitor within +/-',F5.2,
     &' deg in both directions')
				WRITE(*,651) Tls*100.
  651			FORMAT(' Transmission of last section =',F6.2,'%')
			ELSE
				WRITE(*,660) neutrons, neutrons, lambda
  660			FORMAT(' Out of',I10,'/',I9,' up/down neutrons at lambda=',
     &F5.2,' A :')
				IF (nsections.EQ.1) THEN
					WRITE(*,670) nthrough(1,1), nthrough(1,2), 1
  670				FORMAT(7X,I10,'/',I9,' made it through section',I2)
				ELSE
					WRITE(*,670) (nthrough(i,1), nthrough(i,2),
     &i, i=1,nsections-1)
					IF (UCN) WRITE(*,675) (nUCN(i), i=1,2), nsections
  675				FORMAT(7X,I10,'/',I9,' made it into section',I2)
					WRITE(*,670) (nthrough(nsections,i), i=1,2), nsections
				END IF
				WRITE(*,680) nmon(1), nmon(2), (ncentral(i,1), ncentral(i,2), 
     &divlimit, i=1,3)
  680			FORMAT(7X,I10,'/',I9,' hit the monitor'/
     &7X,I10,'/',I9,' hit the monitor within +/-',F5.2,
     &' deg horizontally'/
     &7X,I10,'/',I9,' hit the monitor within +/-',F5.2,
     &' deg vertically'/
     &7X,I10,'/',I9,' hit the monitor within +/-',F5.2,
     &' deg in both directions')
				WRITE(*,651) Tls*100.
				WRITE(*,690) 100.*Pa, 100.*dPa, FRa
  690			FORMAT(' At monitor,      P =',F6.2,' +/-',F5.2,
     &' %   => FR =',F5.1)
				WRITE(*,695) divlimit, 100.*Pb, 100.*dPb, FRb
  695			FORMAT(' At mon over +/-',F4.2,'deg, P =',F6.2,' +/-',F5.2,
     &' %   => FR =',F5.1)
			END IF
		ELSE
			IF (polarised) THEN
				PRINT*, nmon(1), nmon(2)
			ELSE
				PRINT*, nmon(1)
			END IF
		END IF

		DO iv=3,4	! hor and ver div on monitor
			ivar=2*msections+iv
			sum0=0.
			DO i=-nbin,nbin
				y=binvar(ivar,1,i)
				IF (y.NE.0. .AND. sum0.EQ.0.) i1=i
				IF (y.NE.0.) i2=i
				sum0=sum0+y
			END DO
			IF (sum0.EQ.0.) THEN
				w=0.
			ELSE
c				GOTO 655	! use full 100% width instead of 90% width
				DO istep=1,nbin
					sum=0.
					DO i=i1+istep,i2-istep
						sum=sum+binvar(ivar,1,i)
					END DO
					IF (sum/sum0.LT.0.9) GOTO 652
				END DO
  652			istep=istep-1
				i1=i1+istep
				i2=i2-istep
				DO istep=1,nbin
					sum=0.
					DO i=i1+istep,i2
						sum=sum+binvar(ivar,1,i)
					END DO
					IF (sum/sum0.LT.0.9) GOTO 653
				END DO
  653			istep=istep-1
				i1=i1+istep
				DO istep=1,nbin
					sum=0.
					DO i=i1,i2-istep
						sum=sum+binvar(ivar,1,i)
					END DO
					IF (sum/sum0.LT.0.9) GOTO 654
				END DO
  654			istep=istep-1
				i2=i2-istep
  655			CONTINUE
				xmin=var(ivar,i1)
				xmax=var(ivar,i2)
				w=(xmax-xmin)/2.
			END IF
			IF (.NOT.nooutput) WRITE(*,656) 2.*w, vblurb(ivar)
  656		FORMAT(' > 90% of neutrons within ',F5.2,
     &' deg (FW) of ',A50)
			IF (iv.EQ.3) wH=2.*w
			IF (iv.EQ.4) wV=2.*w
		END DO

		IF (UCN .AND. .NOT.nooutput) WRITE(*,705) 
     &(100.*PsumUCN(i)/REAL(neutrons), i=1,2)
  705	FORMAT(F9.6,' /',F9.6,' % of incident neutrons created a UCN')

		IF (saveresults) THEN
			CALL makefilename(input_file,0.,filename,nin)
			DO ipol=1,npol
				IF (npol.EQ.1) THEN
					n=nin-2
					out_file(1:n)=filename(1:nin-6)//'.out'
				ELSE
					PRINT*,'npol=',npol,' ipol=',ipol,' cpol=',cpol
					WRITE(cpol,'(I1)') ipol-1
					n=nin
					out_file(1:nin)=filename(1:nin-6)//'_'//cpol//'.out'
				END IF
				DO i=n+1,50
					out_file(i:i)=' '
				END DO
				IF (.NOT.nooutput) PRINT*,'Writing results to ',out_file
				INQUIRE(FILE=out_file,EXIST=exist)
				IF (.NOT.exist) THEN
					OPEN(10,FILE=out_file,STATUS='NEW')
					IF (.NOT.scan) scanvarname='Wavelength'
					IF (nsections.LT.10) THEN
						IF (UCN) THEN
							WRITE(LFMT,
     &'("(2X,A10,''  N(start)'',",
     &I1,"('' Section'',I1,1X)",
     &"''   Nmon   N_UCN   P_UCN(%)''",
     &")")') nsections
						ELSE
							WRITE(LFMT,
     &'("(2X,A10,''  N(start)'',",
     &I1,"('' Section'',I1,1X)",
     &"''   Nmon     Ncen(H)   Ncen(V) Ncen(both) HorDiv VerDiv''",
     &")")') nsections
						END IF
					ELSE
						IF (UCN) THEN
							WRITE(LFMT,
     &'("(2X,A10,''  N(start)'',",
     &"9('' Section'',I1,1X)",I2,"('' Section'',I2)",
     &"''   Nmon   N_UCN   P_UCN(%)''",
     &")")') nsections-9
						ELSE
							WRITE(LFMT,
     &'("(2X,A10,''  N(start)'',",
     &"9('' Section'',I1,1X)",I2,"('' Section'',I2)",
     &"''   Nmon     Ncen(H)   Ncen(V) Ncen(both) HorDiv VerDiv''",
     &")")') nsections-9
						END IF
					END IF
					WRITE(10,LFMT) scanvarname, (i, i=1,nsections)
				ELSE
					OPEN(10,FILE=out_file,STATUS='OLD',ACCESS='APPEND')
				END IF
				xvalue=lambda
				IF (scan) xvalue=scanvar
				IF (nsections.LT.10) THEN
					IF (UCN) THEN
						WRITE(FMT,'("(F10.3,I10,",I1,"I10,2I10,F10.6)")') 
     &nsections
					ELSE
						WRITE(FMT,'("(F10.3,I10,",I1,"I10,4I10,2F8.2)")') 
     &nsections
					END IF
				ELSE
					IF (UCN) THEN
						WRITE(FMT,'("(F10.3,I10,",I2,"I10,2I10,F10.6)")') 
     &nsections
					ELSE
						WRITE(FMT,'("(F10.3,I10,",I2,"I10,4I10,2F8.2)")') 
     &nsections
					END IF
				END IF
				IF (UCN) THEN
					WRITE(10,FMT) xvalue, neutrons, 
     &(nthrough(i,ipol), i=1,nsections), nmon(ipol), nUCN(ipol), 
     &100.*PsumUCN(ipol)/REAL(neutrons)
				ELSE
					WRITE(10,FMT) xvalue, neutrons, 
     &(nthrough(i,ipol), i=1,nsections), nmon(ipol), 
     &(ncentral(i,ipol), i=1,3), wH, wV
				END IF
				CLOSE(10)
			END DO
		END IF

		IF (batchjob) THEN
			IF (.NOT.saveplots) GOTO 999
		ELSE
			CALL pgend
		END IF

		IF (showplots) THEN
			WRITE(*,710)
  710		FORMAT(' See plots? (def=Y) ',$)
			READ(*,'(A1)') Answer
			IF (Answer.EQ.'n' .OR. Answer.EQ.'N') THEN
				IF (.NOT.saveplots) GOTO 999
				showplots=.FALSE.
			END IF
		END IF

		npol=1
		IF (polarised) npol=2

		DO iv=1,nvars
			plot=.FALSE.
			saveplot=.FALSE.
			ivar=iusevars(iv)
			IF (showplots) THEN
				WRITE(*,915) vblurb(ivar)
  915			FORMAT(' Plot ',A50,' ? (def=N) : ',$)
				READ(*,'(A1)') Answer
				IF (Answer.EQ.'y' .OR. Answer.EQ.'Y') plot=.TRUE.
			END IF
			IF (saveplots .AND. .NOT.plot) THEN
				saveplot=.TRUE.
				plot=.TRUE.
			END IF
			IF (plot) THEN
  930			IF (saveplot) THEN
					imin=-nbin
					imax=nbin
				ELSE
					IF (ivar.LT.2*msections+9) THEN
						DO i=-nbin,nbin
							IF (binvar(ivar,1,i).GT.0. .OR. 
     &binvar(ivar,2,i).GT.0.) GOTO 940
						END DO
  940					imin=MAX(i-1,-nbin)
						DO i=nbin,-nbin,-1
							IF (binvar(ivar,1,i).GT.0. .OR. 
     &binvar(ivar,2,i).GT.0.) GOTO 945
						END DO
  945					imax=MIN(i+1,nbin)
						IF (imax.LE.imin) THEN
							imin=-1
							imax=1
						END IF
					ELSE
						imin=-nbin
						imax=MIN(NINT(Lmon/dvar(ivar))-nbin+1,nbin)
					END IF
				END IF
				ymin=0.
				ymax=0.
				DO i=imin,imax
					ymax=MAX(ymax,binvar(ivar,1,i))
					ymax=MAX(ymax,binvar(ivar,2,i))
				END DO
				DO ipol=1,npol
					ysum=0.
					DO i=-nbin,nbin
						x=var(ivar,i)
						y=binvar(ivar,ipol,i)
						e=SQRT(y)
						ysum=ysum+y
						xplot(i+nbin+1)=x
						yplot(i+nbin+1)=y
						eplot(i+nbin+1)=e
					END DO
					IF (saveplot) THEN
						CALL makefilename(input_file,scanvar,filename,n)
						IF (scan) THEN
							out_file(1:n+18)=filename(1:n)//scanvarname
     &//'.'//vname(ivar,1)(1:7)
						ELSE
							out_file(1:n+2)=filename(1:n-6)
     &//'.'//vname(ivar,1)(1:7)
						END IF
						PRINT*,'Saving ',out_file
						OPEN(10,FILE=out_file,STATUS='UNKNOWN')
						IF (npol.EQ.1) THEN
							WRITE(10,960) vname(ivar,1), vname(ivar,2), 
     &vname(ivar,2)
  960						FORMAT(2X,A9,8X,A9,8X,'d',A9)
							DO i=imin,imax
								WRITE(10,*) xplot(i+nbin+1), 
     &yplot(i+nbin+1), eplot(i+nbin+1)
							END DO
						ELSE
							WRITE(10,965) vname(ivar,1), vname(ivar,2), 
     &vname(ivar,2), vname(ivar,2), vname(ivar,2)
  965						FORMAT(2X,A9,8X,A9,'_1',8X,'d',A9,'_1',
     &8X,A9,'_2',8X,'d',A9,'_2')
							DO i=imin,imax
								y2=binvar(ivar,2,i)
								WRITE(10,*) xplot(i+nbin+1), 
     &yplot(i+nbin+1), eplot(i+nbin+1), REAL(y2), REAL(SQRT(y2))
							END DO
						END IF
						CLOSE(10)
						GOTO 990
					ELSE
						WRITE(*,970) ipol, NINT(ysum)
  970					FORMAT(' Integral(',I1,') =',I7)
						IF (ipol.EQ.1) THEN
							CALL pgbegin(0,'/xw',1,1)
							CALL pgask(.FALSE.)
							WRITE(title,975) vblurb(ivar), vname(ivar,2), 
     &vname(ivar,1)
  975						FORMAT(A50,': ',A9,' vs ',A9)
							xmin=xplot(imin+nbin+1)
							xmax=xplot(imax+nbin+1)
							CALL pgenv(xmin,xmax,ymin,ymax,0,0)
						END IF
						IF (ipol.EQ.2) CALL pgsci(2)
						im=9+4*(1-ipol)
						CALL pgpt(2*nbin+1,xplot,yplot,im)
						CALL pgerrb(6,2*nbin,xplot,yplot,eplot,0.)
						IF (ipol.EQ.1) CALL pgline(2*nbin+1,xplot,yplot)
						IF (ipol.EQ.1) 
     &CALL pglab(vname(ivar,1),vname(ivar,2),vblurb(ivar))
						IF (ipol.EQ.2) CALL pgsci(1)
					END IF	! if saveplot
				END DO	! do ipol=1,npol
				IF (saveplots) THEN
					saveplot=.TRUE.
				ELSE
					WRITE(*,980)
  980				FORMAT(' Save data in file? (def=N) : ',$)
					READ(*,'(A1)') Answer
					IF (Answer.EQ.'y' .OR. Answer.EQ.'Y') saveplot=.TRUE.
				END IF
				IF (saveplot) GOTO 930
				CALL pgend
			END IF	! if plot
  990		CONTINUE
		END DO	! do ivar=1,nvars

  999	IF (looplambda) THEN
			IF (keepgoing) neutrons=n0
			IF (ilambda.LT.nlambda) GOTO 300
		END IF
		IF (scan) THEN
			IF (keepgoing) neutrons=n0
			IF (iscan.LT.nscan) GOTO 101
		END IF

 1000	STOP
		END

