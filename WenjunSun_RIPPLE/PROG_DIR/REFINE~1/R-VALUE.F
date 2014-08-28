	integer num
	real fo, fc
	real sumfo
	real sumdiff
	real rvalue

	open(unit=19,file='mdf.dat')
	read(19,*)num
	sumfo=0.0
	sumdiff=0.0
	do i=1,num
		read(19,*)fo,fc
		sumdiff=sumdiff+abs(abs(fo)-abs(fc))
		sumfo=sumfo+abs(fo)
	enddo
	rvalue=sumdiff/sumfo
	print *,'R-value=',rvalue
	stop
	end
