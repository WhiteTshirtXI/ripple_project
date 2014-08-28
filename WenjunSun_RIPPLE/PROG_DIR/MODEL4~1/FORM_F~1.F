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
c 8/27/95 Model 2, tilted trans-bilayer electron density profile, 
c perpendicular to the major side of the ripples.
c
	REAL d,gamma,lambda,amp,xo
	REAL rhm, zh
	REAL alpha
	INTEGER h,k
	REAL pi, radcon
	REAL qx,qz,w
	REAL qz0,w0,form0,form20,form2
	REAL form
	pi=4.0*atan(1.0)
	radcon=pi/180.0
	print *,'Supply d,lambda,gamma'
	read *,d,lambda,gamma
	gamma=gamma*radcon
	print *,'Supply amp,xo'
	read *,amp,xo
	alpha=atan(amp/xo)
	print *,'rhm, zh'
	read *,rhm, zh
c Calculate the form factor for (10) peak, for normalization
	qz0=2.0*pi/d
	form20=(rhm*cos(qz0*zh*cos(alpha))-1.0)
	w0=0.5*qz0*amp
	form0=form20*(sin(w0)/w0)
c Calculate form factor for (hk) peak.
10	print *,'(h,k) ?'
	read *,h,k
	if(h .gt. 10)goto 20
	qx=2.0*pi*k/lambda
	qz=2.0*pi*(h/d-k/lambda/tan(gamma))
	form2=(rhm*cos(qz*zh*cos(alpha)-qx*zh*sin(alpha))-1.0)
	w=0.5*(qx*xo+qz*amp)
	form=sin(w)*((k*pi*xo/lambda-w)/(w*(k*pi-w)))
	form=form*form2
	print *,h,k,form*100.0/abs(form0)
	goto 10
20	stop
	END

