C	Program to calculate the first part of the ripple form factor.
C	It is a triangular ripple.
c	Needed input:
c		
c		d,gamma,rlam,amp,xo
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
	REAL dsp,gamma,rlam
	REAL Ch,sigh,Cch3,sigch3
	REAL amp,x0
	REAL psi,Xh
	REAL f1, f2
	REAL fact
	INTEGER h,k
	REAL pi, radcon
	REAL qx,qz,q,w
	REAL qzp,yup,yita
	REAL formt,formc
	REAL form0,formall

	pi=4.0*atan(1.0)
	radcon=pi/180.0
	print *,'Supply d,rlam,gamma'
	read *,dsp,rlam,gamma
	gamma=gamma*radcon
	print *,'Ch,sigh,Cch3,sigch3'
	read *,Ch,sigh,Cch3,sigch3
	print *,'Supply amp,x0'
	read *,amp,x0
	print *,'Supply psi'
	read *,psi
	psi=psi*radcon
	print *,'Supply Xh'
	read *, Xh
	print *,'f1, f2?'
	read *, f1, f2
	print *,'factor?'
	read *,fact
c Calculate the form factor for (10) peak, for normalization
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

c Calculate form factor for (hk) peak.
10	print *,'(h,k) ?'
	read *,h,k
	if(h .gt. 10)goto 20
	qx=2.0*pi*k/rlam
	qz=2.0*pi*(h/dsp-k/rlam/tan(gamma))
	q=sqrt(qx*qx+qz*qz)
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
        formall=formt*formc
        formall=formall*100.0*fact/abs(form0)
	print *,h,k,q*100.0,formall
	goto 10
20	stop
	END

