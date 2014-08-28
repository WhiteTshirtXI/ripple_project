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
c 8/27/95 Tilted trans-bilayer electron density profile, with the tilt
c angle as a free parameter.
c
	REAL d,gamma,lambda,amp,xo
	REAL rhm, z0
	REAL alpha
	INTEGER h,k
	REAL pi, radcon
	REAL qx,qz,w
	REAL qz0,w0,rint0,rint20,rint2
	REAL rint
	pi=4.0*atan(1.0)
	radcon=pi/180.0
	print *,'Supply d,lambda,gamma'
	read *,d,lambda,gamma
	gamma=gamma*radcon
	print *,'Supply amp,xo'
	read *,amp,xo
	print *,'Supply alpha:'
	read *,alpha
	alpha=alpha*radcon
	print *,'rhm, z0'
	read *,rhm, z0
c I(10), for normalization
	qz0=2.0*pi/d
	rint20=(rhm*cos(qz0*z0*cos(alpha))-1.0)**2
	w0=0.5*qz0*amp
	rint0=rint20*(sin(w0)/w0)**2
10	print *,'(h,k) ?'
	read *,h,k
	if(h .gt. 10)goto 20
	qx=2.0*pi*k/lambda
	qz=2.0*pi*(h/d-k/lambda/tan(gamma))
	rint2=(rhm*cos(qz*z0*cos(alpha)-qx*z0*sin(alpha))-1.0)**2
	w=0.5*(qx*xo+qz*amp)
	rint=sin(w)**2*((k*pi*xo/lambda-w)/(w*(k*pi-w)))**2
	rint=rint*rint2
	print *,h,k,rint*100.0/rint0
	goto 10
20	stop
	END

