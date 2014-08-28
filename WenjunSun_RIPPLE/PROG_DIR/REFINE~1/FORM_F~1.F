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
c 8/27/95 Model 2, tilted trans-bilayer electron density profile, with
c tilt angle as a free parameter.
c
c 11/13/95 Refined model. f1: minor/major electron density ratio;
c			  f2: kink delta function height.
c
	REAL d,gamma,lambda
	REAL amp,xo
	REAL psi
	REAL rhm, zh
	REAL factor
	REAL f1,f2
	INTEGER h,k
	REAL pi, radcon
	REAL qx,qz,q,w
	REAL form0,form1,form2
	REAL form
	pi=4.0*atan(1.0)
	radcon=pi/180.0
	print *,'Supply d,lambda,gamma'
	read *,d,lambda,gamma
	gamma=gamma*radcon
	print *,'Supply amp,xo'
	read *,amp,xo
	print *,'Supply psi:'
	read *,psi
	psi=psi*radcon
	print *,'rhm, zh'
	read *,rhm, zh
	print *,'factor:'
	read *,factor
	print *,'f1, f2:'
	read *,f1, f2
c Calculate the form factor for (10) peak, for normalization
	h=1
	k=0
	qx=2.0*pi*k/lambda
	qz=2.0*pi*(h/d-k/lambda/tan(gamma))
	form2=rhm*cos(qz*zh*cos(psi)-qx*zh*sin(psi))-1.0
	w=0.5*(qx*xo+qz*amp)
	form1=xo*sin(w)/w/lambda
	form1=form1+(-1)**k*f1*(lambda-xo)*sin(k*pi-w)/(k*pi-w)/lambda
	form1=form1+2*f2*cos(w)/lambda
	form0=form1*form2
	print *,'F(1,0)=',form0
c Calculate form factor for (hk) peak.
10	print *,'(h,k) ?'
	read *,h,k
	if(h .gt. 10)goto 20
	qx=2.0*pi*k/lambda
	qz=2.0*pi*(h/d-k/lambda/tan(gamma))
	q=sqrt(qx*qx+qz*qz)
	form2=rhm*cos(qz*zh*cos(psi)-qx*zh*sin(psi))-1.0
	w=0.5*(qx*xo+qz*amp)
	form1=xo*sin(w)/w/lambda
	form1=form1+(-1)**k*f1*(lambda-xo)*sin(k*pi-w)/(k*pi-w)/lambda
	form1=form1+2*f2*cos(w)/lambda
	form=form1*form2
	print *,h,k,q*100.0,factor*form*100.0/abs(form0)
	goto 10
20	stop
	END

