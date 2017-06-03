		COMMON /c_guide/ Rguide, xc1, yc1, zc1, xc2, yc2, zc2, 
     &ag, bg, cw, dw, phig, xstart, xend, ystart, yend, 
     &zstart, zend, apara, xfocus, mnumber, refl0, reflm, reflfile, 
     &Wwav, tilt, Wguide, Hguide, nmultiple, gmultiple, chamfer, 
     &phiguideH, phiguideV, cthickness, 
     &Qarr, Parr, nQ, Si_lambda, Si_mu

		DIMENSION Rguide(msections,4), refl0(msections,5), 
     &reflm(msections,5), Wwav(msections,5), 
     &xc1(msections), yc1(msections), zc1(msections), 
     &xc2(msections), yc2(msections), zc2(msections), 
     &ag(msections,6), bg(msections,6), 
     &cw(msections+1,2), dw(msections+1,2), 
     &phig(msections,6), gshift(msections,4), 
     &xstart(msections,6), xend(msections,6),
     &ystart(msections,6), yend(msections,6), 
     &zstart(msections,3:4), zend(msections,3:4), 
     &tilt(msections,6), 
     &Hshift(msections), Vshift(msections), 
     &Qarr(msections,5,mQ), nQ(msections,5), 
     &Parr(msections,5,2,mQ), 
     &apara(msections), xfocus(msections), 
     &Wguide(msections,2), Hguide(msections,2), 
     &phiguideH(msections), phiguideV(msections), 
     &nmultiple(msections), gmultiple(msections), 
     &chamfer(msections,5), 
     &cthickness(msections), nsides(msections), 
     &ivcavity(msections), 
     &Si_lambda(19), Si_mu(19)


		REAL*8	mnumber(msections,5)

		CHARACTER*40 reflfile(msections,5)

		LOGICAL tilted(msections)


