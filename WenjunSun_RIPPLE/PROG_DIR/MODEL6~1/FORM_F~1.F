C	Program to calculate the first part of the ripple form factor.
C	It is a triangular ripple.
c	Needed input:
c		
c		d,gamma,lambda,amp,xo
c		(h,k)----> the idex of the peak.
c
c	Output:
c		The intensity of the peak.
c
c 2/8/95	The trans-membrane part is added. Delta functions for
c	the headgroups and the methyl groups, also the position of
c	the headgroups.
c
c 9/11/95 Model 7, tilted trans-bilayer electron density profile, 
c free tilt angle for this T(x,z). 1G MW model results, with
c Ch,sigh,Cch3,sigch3 from MW's 1G constrained fit. Xh is a
c free parameter.
c
	REAL dsp,gamma,lambda
	REAL Ch,sigh,Cch3,sigch3
	REAL amp,x0
	REAL alpha,Xh
	INTEGER h,k
	REAL pi, radcon
	REAL qx,qz,w
	REAL qzp,yup,yita
	REAL qz0,w0
	REAL formt,formc
	REAL form0,form
	pi=4.0*atan(1.0)
	radcon=pi/180.0
	print *,'Supply d,lambda,gamma'
	read *,dsp,lambda,gamma
	gamma=gamma*radcon
	print *,'Ch,sigh,Cch3,sigch3'
	read *,Ch,sigh,Cch3,sigch3
	print *,'Supply amp,x0'
	read *,amp,x0
	alpha=atan(amp/x0)
	print *,'Supply Xh'
	read *, Xh
c Calculate the form factor for (10) peak, for normalization
        qz0=2.0*pi/dsp
        qzp=qz0*cos(alpha)
        yup=qzp+pi/2.0/sigh
        yita=qzp-pi/2.0/sigh
	formt=2.0*sqrt(2.0*pi)*Ch*sigh*cos(qzp*Xh)*exp(-0.5*(qzp*sigh)**2)
	formt=formt+sqrt(2.0*pi)*Cch3*sigch3*exp(-0.5*(qzp*sigch3)**2)
	formt=formt-2.0*sin(qzp*(Xh-sigh))/qzp
	formt=formt-2.0*cos(qzp*Xh)*sin(qzp*sigh)/qzp
	formt=formt-sin(qzp*Xh)*(sin(yita*sigh)/yita-sin(yup*sigh)/yup)
	formt=formt*cos(alpha)
	w0=0.5*qz0*amp
	formc=sin(w0)/w0
	form0=formt*formc
	print *,form0

c Calculate form factor for (hk) peak.
10	print *,'(h,k) ?'
	read *,h,k
	if(h .gt. 10)goto 20
	qx=2.0*pi*k/lambda
	qz=2.0*pi*(h/dsp-k/lambda/tan(gamma))
        qzp=qz*cos(alpha)-qx*sin(alpha)
        yup=qzp+pi/2.0/sigh
        yita=qzp-pi/2.0/sigh
	formt=2.0*sqrt(2.0*pi)*Ch*sigh*cos(qzp*Xh)*exp(-0.5*(qzp*sigh)**2)
	formt=formt+sqrt(2.0*pi)*Cch3*sigch3*exp(-0.5*(qzp*sigch3)**2)
	formt=formt-2.0*sin(qzp*(Xh-sigh))/qzp
	formt=formt-2.0*cos(qzp*Xh)*sin(qzp*sigh)/qzp
	formt=formt-sin(qzp*Xh)*(sin(yita*sigh)/yita-sin(yup*sigh)/yup)
	formt=formt*cos(alpha)
	w=0.5*(qx*x0+qz*amp)
	formc=sin(w)*((k*pi*x0/lambda-w)/(w*(k*pi-w)))
	form=formt*formc
	print *,h,k,form*100.0/abs(form0)
	goto 10
20	stop
	END

