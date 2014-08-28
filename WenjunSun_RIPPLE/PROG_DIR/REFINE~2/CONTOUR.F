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
	REAL d,gamma,lambda,amp,xo
	INTEGER h,k
	REAL pi, radcon
	REAL qx,qz,w
	REAL qz0,w0,rint0
	REAL rint
	pi=4.0*atan(1.0)
	radcon=pi/180.0
	print *,'Supply d,lambda,gamma'
	read *,d,lambda,gamma
	gamma=gamma*radcon
	print *,'Supply amp,xo'
	read *,amp,xo
	qz0=2.0*pi/d
	w0=0.5*qz0*amp
	rint0=(sin(w0)/w0)**2
10	print *,'(h,k) ?'
	read *,h,k
	if(h .gt. 10)goto 20
	qx=2.0*pi*k/lambda
	qz=2.0*pi*(h/d-k/lambda/tan(gamma))
	w=0.5*(qx*xo+qz*amp)
	rint=sin(w)**2*((k*pi*xo/lambda-w)/(w*(k*pi-w)))**2
	print *,h,k,rint*100.0/rint0
	goto 10
20	stop
	END

