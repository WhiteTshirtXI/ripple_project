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
	REAL d,gamma,lambda,amp,xo
	REAL rhm, z0
	INTEGER h,k
	REAL pi, radcon
	REAL qx,qz,w
	REAL qz0,w0,rint0,rint20,rint2
	REAL rint
	pi=4.0*atan(1.0)
	radcon=pi/180.0
        open(unit=16, file='model1_form.dat')
	print *,'Supply d,lambda,gamma'
	read *,d,lambda,gamma
	gamma=gamma*radcon
	print *,'Supply amp,xo'
	read *,amp,xo
	print *,'rhm, z0'
	read *,rhm, z0
	qz0=2.0*pi/d
	rint20=(rhm*cos(qz0*z0)-1.0)
	w0=0.5*qz0*amp
	rint0=rint20*(sin(w0)/w0)
	write(16,*)d,lambda,gamma/radcon
	write(16,*)(21*21-1)
	do h=0,20
		do k=-10,10
			if( h .eq. 0 .and. k .eq. 0) goto 10
			qx=2.0*pi*k/lambda
			qz=2.0*pi*(h/d-k/lambda/tan(gamma))
			rint2=(rhm*cos(qz*z0)-1.0)
			w=0.5*(qx*xo+qz*amp)
			rint=sin(w)*((k*pi*xo/lambda-w)/(w*(k*pi-w)))
			rint=rint*rint2
			write(16, *)h,k,rint*100.0/abs(rint0)
10		enddo
	enddo
	close(16)
	stop
	END

