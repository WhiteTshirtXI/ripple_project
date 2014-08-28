	real a,b,tilt,alpha
	integer m,n
	real pi, radcon
	real qx,qy,qz
	real sigma,qr,q
	pi=4.0*atan(1.0)
	radcon=pi/180.0
	print *,'Supply a,tilt:'
	read *,a,tilt
	tilt=tilt*radcon
	b=a/sqrt(3.0)/cos(tilt)
10	print *,'Supply alpha(-1000 to stop):'
	read *,alpha
	if(alpha .eq. -1000)go to 30
	alpha=alpha*radcon
20	print *,'Supply M,N (M=-10 to go back):'
	read *,m,n
	if(m .eq. -10)go to 10
	qx=2.0*pi*m*cos(alpha)/a-2.0*pi*n*tan(tilt)*sin(alpha)
	qy=2.0*pi*n/b
	qz=2.0*pi*m*sin(alpha)/a+2.0*pi*n*tan(tilt)*cos(alpha)/b
	qr=sqrt(qx*qx+qy*qy)
	q=sqrt(qr*qr+qz*qz)
	sigma=acos(qr/q)/radcon
	print *,'(',m,',',n,')-------->',sigma
	go to 20
30	stop
	END

