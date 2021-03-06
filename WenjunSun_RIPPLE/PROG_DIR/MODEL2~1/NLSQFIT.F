C 8/27/95 Tilted trans-bilayer electron density profile ---- Model 2
C perpendicular to the major side of the ripples.
C
C This program will fit my electron density model of the ripple phase
C to Wack's DMPC x-ray data. It will use the IMSL UNISQ routine to do 
C NLLSQ fitting. 
C
C Needed inputs:
C
C	Ripple lattice parameters: d-spacing, ripple wavelength,oblique angle.
C
C The fitting parameters are: 
C
C	1). Ripple profile: Amplitude A, horizontal projection of the left
C	side x0;
C	2). Trans-membrane: Rhm, the ratio of the headgroup and the methyl
C	term; Zh, the headgroup position from the bilayer center.
C
C The data:
C		Intensities of 20 peaks of DMPC.
C
C---------------------------------------------------------------------------
C
	parameter (n=5,ldfjac=100)
	integer index(ldfjac,2)
	real rinten(ldfjac), bar(ldfjac)
	real x(n), xguess(n), xscale(n), rparam(7)
	real fvec(ldfjac), fjac(ldfjac,n), fscale(ldfjac)
	real sqfjac(n,n),sqinv(n,n)
	real chisq,perr(n)
	integer iparam(6),npoint
	external fcn,unlsf,linrg
C d-spacing,ripple wavelength, oblique angle
	real dsp,rlam,gamma
	common/s0/ dsp,rlam,gamma
	common/s1/ index, rinten, bar
C degree-radian converter
	real pi, radcon

	pi=4.0*atan(1.0)
	radcon=pi/180.0

C Read in ripple lattice parameters and fitting parameters.
	open(unit=19,file='initial.dat')
	read(19,*)dsp,rlam,gamma
	gamma=gamma*radcon
	do j=1, n
		xscale(j)=1.0
		read(19,*) xguess(j)
	enddo
	close(19)

C Read in peak indices and intensities.
	open (unit=21, file='inten.dat')
	npoint=0
	do i=1, 400
		read(21,*,end=321)index(i,1), index(i,2), rinten(i)
		bar(i)=sqrt(rinten(i))
		npoint=npoint+1
	enddo
321	print *,'Data have been successfully read in !'
	print *,npoint,' data points altogether.'
	close(21)
C Data read in.

	do i=1, npoint
		fscale(i)=1.0
	enddo

C Call IMSL NLLSQ subroutine.
	call unlsf(fcn,npoint,n,xguess,xscale,fscale,iparam,
     &           rparam,x,fvec,fjac,ldfjac)
	print *,'Fitting is finished..'

C Fitting is finished. Now save fitting results.
	open(unit=16, file='intenfit.dat')
	do i=1,npoint
		write(16,*)index(i,1),index(i,2),rinten(i),
     &				rinten(i)+fvec(i)*bar(i)
	enddo	
	close(16)

C reduced chi-square.
	chisq=0.0
	do i=1, npoint
		chisq=chisq+fvec(i)*fvec(i)
	enddo
	chisq=chisq/(npoint-n)
	print *,'chisq=',chisq
	print *, 'done....................'

C Calculating errors for fitting parameters.
C First, construct the alpha matrix from Jacobian.
	do i=1,n
		do j=1,n
			sqfjac(i,j)=0.0
			do k=1,npoint
				sqfjac(i,j)=sqfjac(i,j)+fjac(k,i)*fjac(k,j)
			enddo
		enddo
		print *,(sqfjac(i,j),j=1,n)
	enddo
C Second, get the inverse of the alpha matrix.
	call linrg(n,sqfjac,n,sqinv,n)
C Third, fitting errors.
	do i=1,n
		perr(i)=sqrt(sqinv(i,i))
	enddo

C Save fitting parameters to a file.
	open (unit=9, file='intenfit.pa')
	write(9,554) chisq
554	format(1x, 'chi-square=', f10.5)
	write(9,555) (x(j),perr(j), j=1,n)
555	format(1x,f10.5,'	',f10.5)
	close(9)
	stop
	end

C ----------------------------------------------------------------------
C Subroutine to calculate the model intensity.
	subroutine fcn(m,n,x,f)
	real f(m), x(n), rinten(100), bar(100) 
	integer index(100,2)
	real dsp,rlam,gamma
	common/s0/ dsp,rlam,gamma
	common/s1/ index, rinten, bar
	REAL amp,x0,rhm,zh,fact
	REAL alpha
	INTEGER h,k
	REAL qx,qz,w
	REAL qz0,w0,rint0,rint20,rint2
	REAL rint
	real chisq
C degree-radian converter
	real pi

	pi=4.0*atan(1.0)

C Initialize fitting parameters: amp,x0,rhm,zh.
	amp=x(1)
	x0=x(2)
	rhm=x(3)
	zh=x(4)
	fact=x(5)
	print *, (x(i),i=1,n)

C Slope of the major side of the ripples
	alpha=atan(amp/x0)
	print *,alpha*180.0/pi

C Calculate intensity for (10) peak, normalize this to 100.
	qz0=2.0*pi/dsp
	rint20=(rhm*cos(qz0*zh*cos(alpha))-1.0)**2
	w0=0.5*qz0*amp
	rint0=rint20*(sin(w0)/w0)**2

C Calculate intensities for all (h,k) peaks.
	chisq=0.0
	do i=1,m
		h=index(i,1)
		k=index(i,2)
		qx=2.0*pi*k/rlam
		qz=2.0*pi*(h/dsp-k/rlam/tan(gamma))
		rint2=(rhm*cos(qz*zh*cos(alpha)-qx*zh*sin(alpha))-1.0)**2
		w=0.5*(qx*x0+qz*amp)
		rint=sin(w)**2*((k*pi*x0/rlam-w)/(w*(k*pi-w)))**2
		rint=rint*rint2
		f(i)=(fact*rint*100.0/rint0)-rinten(i)
		f(i)=f(i)/bar(i)
		chisq=chisq+f(i)*f(i)
	enddo
	chisq=chisq/(m-n)
	print *,'chisq= ',chisq
	return
	end
