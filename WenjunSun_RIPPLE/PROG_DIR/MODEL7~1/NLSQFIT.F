C 11/27/95 Put in f1, f2 (minor side scaling and delta funtion kinks)
C
C 9/8/95 Tilted trans-bilayer electron density profile ---- Model 7
C 1G hybrid model with bridging function
C fixed Ch, sigh, Cch3, sigch3 to MW DPPC Gel phase 1G model fit results
C T(x,z) tilts freely relative to the major side of the ripples.
C
C This program will fit my electron density model of the ripple phase
C to Wack's DMPC x-ray data. It will use the IMSL UNISQ routine to do 
C NLLSQ fitting. 
C
C Needed inputs:
C
C	Ripple lattice parameters: d-spacing, ripple wavelength,oblique angle.
C
C The fitting parameters (8 altogether) are: 
C
C	1). Ripple profile: Amplitude A, horizontal projection of the left
C	side x0;
C	2). Trans-membrane: Ch, sigh, Xh: headgroup Gaussian parameters;
C	    Cch3, sigch3: methyl trough Gaussian parameters; Overall
C 	    scaling factor
C
C The data:
C		Intensities of 20 peaks of DMPC.
C
C---------------------------------------------------------------------------
C
	parameter (n=7,ldfjac=100)
	integer index(ldfjac,2)
	real rinten(ldfjac), bar(ldfjac)
	real x(n), xguess(n), xscale(n), rparam(7)
	real fvec(ldfjac), fjac(ldfjac,n), fscale(ldfjac)
	real sqfjac(n,n),sqinv(n,n)
	real chisq,perr(n)
	integer iparam(6),npoint
	external fcn,u4lsf,unlsf,linrg
C d-spacing,ripple wavelength, oblique angle
	real dsp,rlam,gamma
	REAL Ch,sigh,Cch3,sigch3
	common/s0/ dsp,rlam,gamma
	common/s1/ index, rinten, bar
	common/s2/ Ch,sigh,Cch3,sigch3
	real pi,radcon

	pi=4.0*atan(1.0)
	radcon=pi/180.0

C Read in ripple lattice parameters and fitting parameters.
	open(unit=19,file='initial.dat')
	read(19,*)dsp,rlam,gamma
	gamma=gamma*radcon
	read(19,*)Ch,sigh
	read(19,*)Cch3,sigch3
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
	call u4lsf(iparam,rparam)
	iparam(3)=10000
	call unlsf(fcn,npoint,n,xguess,xscale,fscale,iparam,
     &           rparam,x,fvec,fjac,ldfjac)
	print *,'Fitting is finished..'

C Fitting is finished. Now save fitting results.
C Save form factor data and fit to formfit.dat
	call formwrt(npoint,n,x)
C Save intensities
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
	REAL Ch,sigh,Cch3,sigch3
	common/s0/ dsp,rlam,gamma
	common/s1/ index, rinten, bar
	common/s2/ Ch,sigh,Cch3,sigch3
	REAL amp,x0
	REAL Xh
	REAL f1,f2
	REAL fact
	REAL psi
	INTEGER h,k
	REAL qx,qz,w
	REAL qzp,yup,yita
	REAL formt,formc,rintt,rintc
	REAL rint0
	REAL rint
	real chisq
C degree-radian converter
	real pi,radcon

	pi=4.0*atan(1.0)
	radcon=pi/180.0

C Initialize fitting parameters: amp,x0,Ch,sigh,Xh,Cch3,sigch3,fact.
	amp=x(1)
	x0=x(2)
	Xh=x(3)
	psi=x(4)
	f1=x(5)
	f2=x(6)
	fact=x(7)
	print *, (x(i),i=1,n)
	psi=psi*radcon
        print *,atan(amp/x0)/radcon

C Calculate intensity for (10) peak, normalize this to 100.
	h=1
	k=0
	qx=2.0*pi*k/rlam
	qz=2.0*pi*(h/dsp-k/rlam/tan(gamma))
	w=0.5*(qx*x0+qz*amp)
	qzp=qz*cos(psi)-qx*sin(psi)
	yup=qzp+pi/2.0/sigh
	yita=qzp-pi/2.0/sigh
C Gaussian transnilayer electron density
	formt=2.0*sqrt(2.0*pi)*Ch*sigh*cos(qzp*Xh)*exp(-0.5*(qzp*sigh)**2)
	formt=formt+sqrt(2.0*pi)*Cch3*sigch3*exp(-0.5*(qzp*sigch3)**2)
	formt=formt-2.0*sin(qzp*(Xh-sigh))/qzp
	formt=formt-2.0*cos(qzp*Xh)*sin(qzp*sigh)/qzp
	formt=formt-sin(qzp*Xh)*(sin(yita*sigh)/yita-sin(yup*sigh)/yup)
	formt=formt*cos(psi)
	rintt=formt*formt
C Contour part, with f1, f2
C Constant major side
	formc=x0*sin(w)/w/rlam
C Constant minor side
	formc=formc+(-1)**k*f1*(rlam-x0)*sin(k*pi-w)/(k*pi-w)/rlam
C Delta funtion kinks
	formc=formc+2*f2*cos(w)/rlam
	rintc=formc*formc
	rint0=rintt*rintc

C Calculate intensities for all (h,k) peaks.
	chisq=0.0
	do i=1,m
	h=index(i,1)
	k=index(i,2)
	qx=2.0*pi*k/rlam
	qz=2.0*pi*(h/dsp-k/rlam/tan(gamma))
	qzp=qz*cos(psi)-qx*sin(psi)
	yup=qzp+pi/2.0/sigh
	yita=qzp-pi/2.0/sigh
C Gaussian transnilayer electron density
	formt=2.0*sqrt(2.0*pi)*Ch*abs(sigh)*cos(qzp*Xh)*exp(-0.5*(qzp*sigh)**2)
	formt=formt+sqrt(2.0*pi)*Cch3*abs(sigch3)*exp(-0.5*(qzp*sigch3)**2)
  	formt=formt-2.0*sin(qzp*(Xh-sigh))/qzp
	formt=formt-2.0*cos(qzp*Xh)*sin(qzp*sigh)/qzp
	formt=formt-sin(qzp*Xh)*(sin(yita*sigh)/yita-sin(yup*sigh)/yup)
	formt=formt*cos(psi)
	rintt=formt*formt
C Contour part, with f1, f2
	w=0.5*(qx*x0+qz*amp)
C Constant major side
        formc=x0*sin(w)/w/rlam
C Constant minor side
        formc=formc+(-1)**k*f1*(rlam-x0)*sin(k*pi-w)/(k*pi-w)/rlam
C Delta funtion kinks
        formc=formc+2*f2*cos(w)/rlam
        rintc=formc*formc
	rint=rintt*rintc
	f(i)=(fact*rint*100.0/rint0)-rinten(i)
	f(i)=f(i)/bar(i)
	chisq=chisq+f(i)*f(i)
	enddo
	chisq=chisq/(m-n)
	print *,'chisq= ',chisq
	return
	end

C ----------------------------------------------------------------------
C Subroutine to calculate the model intensity.
	subroutine formwrt(m,n,x)
	real x(n), rinten(100),bar(100)
	integer index(100,2)
	real dsp,rlam,gamma
	REAL Ch,sigh,Cch3,sigch3
	common/s0/ dsp,rlam,gamma
	common/s1/ index, rinten, bar
	common/s2/ Ch,sigh,Cch3,sigch3
	REAL amp,x0
	REAL Xh
	REAL f1,f2
	REAL fact
	REAL psi
	INTEGER h,k
	REAL qx,qz,w
	REAL qzp,yup,yita
	REAL formt,formc
	REAL fdata,form0,formall
C degree-radian converter
	real pi,radcon

	pi=4.0*atan(1.0)
	radcon=pi/180.0

C Initialize fitting parameters: amp,x0,Ch,sigh,Xh,Cch3,sigch3,fact.
	amp=x(1)
	x0=x(2)
	Xh=x(3)
	psi=x(4)
	f1=x(5)
	f2=x(6)
	fact=sqrt(x(7))
	print *, (x(i),i=1,n)
	psi=psi*radcon
        print *,atan(amp/x0)/radcon

C Calculate intensity for (10) peak, normalize this to 100.
	h=1
	k=0
	qx=2.0*pi*k/rlam
	qz=2.0*pi*(h/dsp-k/rlam/tan(gamma))
	w=0.5*(qx*x0+qz*amp)
	qzp=qz*cos(psi)-qx*sin(psi)
	yup=qzp+pi/2.0/sigh
	yita=qzp-pi/2.0/sigh
C Gaussian transnilayer electron density
	formt=2.0*sqrt(2.0*pi)*Ch*sigh*cos(qzp*Xh)*exp(-0.5*(qzp*sigh)**2)
	formt=formt+sqrt(2.0*pi)*Cch3*sigch3*exp(-0.5*(qzp*sigch3)**2)
	formt=formt-2.0*sin(qzp*(Xh-sigh))/qzp
	formt=formt-2.0*cos(qzp*Xh)*sin(qzp*sigh)/qzp
	formt=formt-sin(qzp*Xh)*(sin(yita*sigh)/yita-sin(yup*sigh)/yup)
	formt=formt*cos(psi)
C Contour part, with f1, f2
C Constant major side
	formc=x0*sin(w)/w/rlam
C Constant minor side
	formc=formc+(-1)**k*f1*(rlam-x0)*sin(k*pi-w)/(k*pi-w)/rlam
C Delta funtion kinks
	formc=formc+2*f2*cos(w)/rlam
	form0=formt*formc

C Calculate intensities for all (h,k) peaks.
	open(unit=17, file='formfit.dat')

	do i=1,m
	h=index(i,1)
	k=index(i,2)
	qx=2.0*pi*k/rlam
	qz=2.0*pi*(h/dsp-k/rlam/tan(gamma))
	w=0.5*(qx*x0+qz*amp)
	qzp=qz*cos(psi)-qx*sin(psi)
	yup=qzp+pi/2.0/sigh
	yita=qzp-pi/2.0/sigh
C Gaussian transnilayer electron density
	formt=2.0*sqrt(2.0*pi)*Ch*abs(sigh)*cos(qzp*Xh)*exp(-0.5*(qzp*sigh)**2)
	formt=formt+sqrt(2.0*pi)*Cch3*abs(sigch3)*exp(-0.5*(qzp*sigch3)**2)
  	formt=formt-2.0*sin(qzp*(Xh-sigh))/qzp
	formt=formt-2.0*cos(qzp*Xh)*sin(qzp*sigh)/qzp
	formt=formt-sin(qzp*Xh)*(sin(yita*sigh)/yita-sin(yup*sigh)/yup)
	formt=formt*cos(psi)
C Contour part, with f1, f2
C Constant major side
        formc=x0*sin(w)/w/rlam
C Constant minor side
        formc=formc+(-1)**k*f1*(rlam-x0)*sin(k*pi-w)/(k*pi-w)/rlam
C Delta funtion kinks
        formc=formc+2*f2*cos(w)/rlam
	formall=formt*formc
	formall=formall*100.0*fact/abs(form0)
	fdata=sqrt(100.0*rinten(i))
	write(17,222)index(i,1),index(i,2),fdata, formall
	enddo

222     format(1x,i3,2x,i3,5x,f6.2,5x,f6.2)
	close(17)
	return
	end
