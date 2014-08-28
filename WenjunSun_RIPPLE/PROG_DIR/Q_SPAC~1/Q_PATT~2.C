/*    
**	For Model 2 (9/6/95)
**
**	This program is to calculate the q-space pattern of the ripple phase.
**	It will read from a file that contains the 2-D oblique lattice
**	parameters and the intensities of the detected peaks.
**	Then it will generate a PS q-space pattern using grey scales to
**	depict each peak.
**
**	Usage:
**
**		%pattern name.in name.ps
*/
#include <stdio.h>
#include <math.h>

#define INCH 1.0
#define PROFILE_GRAY 0.2
#define LINEWIDTH 2
#define RADIUS 7

#define MIN(a,b) ((a > b)? b:a)
#define MAX(a,b) ((a > b)? a:b)

struct unitcell {
	float llx,lly;
	float ulx,uly;
	float urx,ury;
	float lrx,lry;
};

void main(argc, argv)
int argc; 
char *argv[];
{
	char author[30];	/* Author and experimental conditions */
	char sample[30],water[30],temp[30];
	float d, rlam, gamma;	/* D-spacing,ripple wavelength,oblique angle */
	float Amp,x0,Rhm,Zh;	/* Model parameters */
	float Amp_con,x0_con,Zh_con;
	float alpha, beta;	/* Ripple slopes */
	float psi;		/* T(x,z) psi */
	int h, k;		/* Miller indices */
	int hmax,kmax;		/* Max h, k, for plotting use */
	float thre,rint;	/* Grayscale threshold and intensity */
	float reci_con, qx, qz;	/* Q-space coordinates */
	float real_con;		/* User-space to PS page converting factor */
	struct unitcell real_sp;
	struct unitcell reci_sp;
	int c, i;
	float gs;		/* Grayscale of scattering peak */
	float radcon;		/* Radial to degree converter */
	float x,y,slope,length;
	FILE *ofpi, *ofpo, *fopen();

	radcon = M_PI/180.0;

	if(argc < 3){
		printf("Usage: %s inputfile outputfile\n",argv[0]);
		exit(-1);
	}
	ofpi = fopen(argv[1], "r");
	ofpo = fopen(argv[2], "w");

	printf("Start...........\n");
	printf("Input the grayscale threshold: ");
	scanf("%f", &thre);

	fscanf(ofpi,"%s %s %s %s",author,sample,water,temp);
	fscanf(ofpi,"%f %f %f", &d, &rlam, &gamma);
	fscanf(ofpi,"%f %f %f %f %f",&Amp,&x0,&Rhm,&Zh,&psi);
	gamma = gamma*radcon;
	psi=psi*radcon;
	real_con = 0.6*INCH / d;	/* d -----> 0.6 inch */
	Amp_con=Amp*real_con;
	x0_con=x0*real_con;
	Zh_con=Zh*real_con;
	reci_con = 0.6*INCH*d/2.0/M_PI; /* 2*pi/d -----> 0.6 inch */

	/* Real space unit cell box */
	real_sp.llx = -0.5*rlam*real_con;
	real_sp.lly = -2.25*INCH;
	real_sp.ulx = real_sp.llx + d/tan(gamma)*real_con;
	real_sp.uly = real_sp.lly + d*real_con;
	real_sp.urx = real_sp.ulx + rlam*real_con;
	real_sp.ury = real_sp.uly;
	real_sp.lrx = real_sp.llx + rlam*real_con;
	real_sp.lry = real_sp.lly;

	/* Reciprocal space unit cell box */
	reci_sp.llx = 0.0;	/* (0,0) */
	reci_sp.lly = 0.0;
	reci_sp.ulx = 0.0;	/* (1,0) */
	reci_sp.uly = reci_con*(2.0*M_PI/d);
	reci_sp.urx = reci_con*(2.0*M_PI/rlam);	/* (1,1) */
	reci_sp.ury = reci_con*(2.0*M_PI/d-2.0*M_PI/rlam/tan(gamma));
	reci_sp.lrx = reci_sp.urx;	/* (0,1) */
	reci_sp.lry = reci_con*(-2.0*M_PI/rlam/tan(gamma));

	/* Write into a PostScript file */
	ps_header(ofpo);

	hmax = 0;
	kmax = 0;
	while ((c = getc(ofpi)) != EOF){
		ungetc(c,ofpi);
		fscanf(ofpi,"%d %d %f", &h, &k, &rint);
		hmax=MAX(hmax,h);
		kmax=MAX(kmax,_ABS(k));
		qx = 2.0*M_PI*k/rlam;
		qx = qx*reci_con;
		qz = 2.0*M_PI*h/d - 2.0*M_PI*k/rlam/tan(gamma);
		qz = qz*reci_con;
		if(rint > thre) rint = thre;
		gs = 1.0-rint/thre; /* gs=0 -------> darkest */
		fprintf(ofpo,"%.3f i %.3f i %.3f peak\n",qx,qz,gs);
	}

	fprintf(ofpo,"grestore\n\n");

	/* Print indices */
	for(h=0;h<=hmax;++h){	/* lamellar indices */	
		qx=reci_con*(2.0*M_PI*(kmax+1)/rlam);
		qz=reci_con*(2.0*M_PI*h/d - 2.0*M_PI*(kmax+1)/rlam/tan(gamma));
		fprintf(ofpo,"%.3f i %.3f i %d sub moveto (%d) show\n",
			qx,qz,(int)RADIUS/2,h);
	}
	qx=reci_con*(2.0*M_PI*(kmax+2)/rlam);
	qz=reci_con*(2.0*M_PI*(hmax/2.0)/d - 2.0*M_PI*(kmax+2)/rlam/tan(gamma));
	fprintf(ofpo,"%.3f i %.3f i %d sub moveto (h) show\n\n",
		qx,qz,(int)RADIUS/2);

	for(k=-kmax;k<=kmax;++k){	/* ripple indices */	
		qx=reci_con*(2.0*M_PI*k/rlam);
		qz=reci_con*(2.0*M_PI*(hmax+0.5)/d-2.0*M_PI*k/rlam/tan(gamma));
		fprintf(ofpo,"%.3f i %d sub %.3f i moveto (%d) show\n",
			qx,RADIUS,qz,k);
	}
	qx=reci_con*(2.0*M_PI*(-kmax/2.0)/rlam);
	qz=reci_con*(2.0*M_PI*(hmax+1)/d-2.0*M_PI*(-kmax/2.0)/rlam/tan(gamma));
	fprintf(ofpo,"%.3f i %.3f i moveto (k) show\n\n",qx,qz);

	/* Draw reciprocal space unit cell */
	fprintf(ofpo,"newpath\n");
	fprintf(ofpo,"%.3f i %.3f i moveto\n",reci_sp.llx,reci_sp.lly);
	fprintf(ofpo,"%.3f i %.3f i lineto\n",reci_sp.ulx,reci_sp.uly);
	fprintf(ofpo,"%.3f i %.3f i lineto\n",reci_sp.urx,reci_sp.ury);
	fprintf(ofpo,"%.3f i %.3f i lineto\n",reci_sp.lrx,reci_sp.lry);
	fprintf(ofpo,"closepath %d setlinewidth stroke\n\n",LINEWIDTH);

	if(Zh > 1.0) { /* Draw bilayer */
		alpha=atan(Amp/x0);
		beta=atan(Amp/(rlam-x0));

		x=0.5*(real_sp.llx+real_sp.ulx);
		y=0.5*(real_sp.lly+real_sp.uly);
		bilayer(ofpo,x,y,rlam,x0,real_con,x0_con,Amp_con,Zh_con,psi);

		x=rlam*real_con+0.5*(real_sp.llx+real_sp.ulx);
		y=0.5*(real_sp.lly+real_sp.uly);
		bilayer(ofpo,x,y,rlam,x0,real_con,x0_con,Amp_con,Zh_con,psi);

		x=-rlam*real_con+0.5*(real_sp.llx+real_sp.ulx);
		y=0.5*(real_sp.lly+real_sp.uly);
		bilayer(ofpo,x,y,rlam,x0,real_con,x0_con,Amp_con,Zh_con,psi);

		x=1.5*real_sp.llx-0.5*real_sp.ulx;
		y=1.5*real_sp.lly-0.5*real_sp.uly;
		bilayer(ofpo,x,y,rlam,x0,real_con,x0_con,Amp_con,Zh_con,psi);
	
		x=rlam*real_con+1.5*real_sp.llx-0.5*real_sp.ulx;
		y=1.5*real_sp.lly-0.5*real_sp.uly;
		bilayer(ofpo,x,y,rlam,x0,real_con,x0_con,Amp_con,Zh_con,psi);

		x=-rlam*real_con+1.5*real_sp.llx-0.5*real_sp.ulx;
		y=1.5*real_sp.lly-0.5*real_sp.uly;
		bilayer(ofpo,x,y,rlam,x0,real_con,x0_con,Amp_con,Zh_con,psi);
	}

	/* Draw the ripple profile within the real space unit cell */
	if(Amp > 1.0){
		x=0.5*(real_sp.llx+real_sp.ulx);
		y=0.5*(real_sp.lly+real_sp.uly);
		profile(ofpo,x,y,rlam,x0,real_con,x0_con,Amp_con,Zh_con,psi);
	
		x=rlam*real_con+0.5*(real_sp.llx+real_sp.ulx);
		y=0.5*(real_sp.lly+real_sp.uly);
		profile(ofpo,x,y,rlam,x0,real_con,x0_con,Amp_con,Zh_con,psi);

		x=-rlam*real_con+0.5*(real_sp.llx+real_sp.ulx);
		y=0.5*(real_sp.lly+real_sp.uly);
		profile(ofpo,x,y,rlam,x0,real_con,x0_con,Amp_con,Zh_con,psi);

		x=1.5*real_sp.llx-0.5*real_sp.ulx;
		y=1.5*real_sp.lly-0.5*real_sp.uly;
		profile(ofpo,x,y,rlam,x0,real_con,x0_con,Amp_con,Zh_con,psi);
	
		x=rlam*real_con+1.5*real_sp.llx-0.5*real_sp.ulx;
		y=1.5*real_sp.lly-0.5*real_sp.uly;
		profile(ofpo,x,y,rlam,x0,real_con,x0_con,Amp_con,Zh_con,psi);

		x=-rlam*real_con+1.5*real_sp.llx-0.5*real_sp.ulx;
		y=1.5*real_sp.lly-0.5*real_sp.uly;
		profile(ofpo,x,y,rlam,x0,real_con,x0_con,Amp_con,Zh_con,psi);

	/* Draw normals of ripple sides */
		alpha=atan(Amp/x0);
		beta=atan(Amp/(rlam-x0));
	/* left side */
		x=0.5*(real_sp.llx+real_sp.urx);
		y=0.5*(real_sp.lly+real_sp.ury);
		slope=alpha/radcon;
		length=1.8;
		fprintf(ofpo,"%.3g i %.3g i %.3g %.3g i normal\n",
			x,y,slope,length);
	/* right side */
		x=0.5*(real_sp.lrx+real_sp.urx);
		y=0.5*(real_sp.lry+real_sp.ury);
		slope=-beta/radcon;
		length=2.0;
		fprintf(ofpo,"%.3g i %.3g i %.3g %.3g i normal\n",
			x,y,slope,length);
	/* left in reciprocal space */
		x=0.0;
		y=0.0;
		slope=alpha/radcon;
		qz=reci_con*(2.0*M_PI*hmax/d);
		length=qz;
		fprintf(ofpo,"%.3g i %.3g i %.3g %.3g i normal\n",
			x,y,slope,length);
	/* right in reciprocal space */
		x=0.0;
		y=0.0;
		slope=-beta/radcon;
		qx=reci_con*(2.0*M_PI*kmax/rlam);
		length=qx/tan(beta);
		fprintf(ofpo,"%.3g i %.3g i %.3g %.3g i normal\n\n",
			x,y,slope,length);
	}

	/* Draw real space unit cell */
	fprintf(ofpo,"newpath\n");
	fprintf(ofpo,"%.3f i %.3f i moveto\n",
		real_sp.llx,real_sp.lly);
	fprintf(ofpo,"%.3f i %.3f i lineto\n",
		real_sp.ulx,real_sp.uly);
	fprintf(ofpo,"%.3f i %.3f i lineto\n",
		real_sp.urx,real_sp.ury);
	fprintf(ofpo,"%.3f i %.3f i lineto\n",
		real_sp.lrx,real_sp.lry);
	fprintf(ofpo,"closepath %d setlinewidth stroke\n\n",LINEWIDTH);

	/* Print lattice parameters and legends on output page */
	fprintf(ofpo,"LM -3.5 i moveto\n");
	fprintf(ofpo,"getCou (d = %gA, ) show\n",d);
	fprintf(ofpo,"getSym (l) show getCou ( = %gA, ) show\n",rlam);
	fprintf(ofpo,"getSym (g) show getCou ( = %g deg) show\n",gamma/radcon);
	fprintf(ofpo,"newline\n");
	if (Amp > 1.0) {
		fprintf(ofpo,"getCou (A = %gA, x0 = %gA) show\n",Amp,x0);
		fprintf(ofpo,"newline\n");
	}
	if (Zh > 1.0) {
		fprintf(ofpo,"getCou (Rhm = %g, Zh = %gA, Beta = %gdeg) show\n"
			,Rhm,Zh,psi/radcon);
		fprintf(ofpo,"newline\n");
	}
	fprintf(ofpo,"((%s %s %s %s)) show\n",author,sample,water,temp);

	ps_end(ofpo);

	printf("Done...........\n");
	fclose(ofpi);
	fclose(ofpo);
	exit(0);
}

ps_header (ofp)
FILE *ofp;
{
	fprintf(ofp,"%%!PS-Adobe-2.0\n");
	fprintf(ofp,"%%%%Creator: Wenjun Sun\n");
	fprintf(ofp,"%%%%Title: Ripple pattern\n");
	fprintf(ofp,"%%%%Pages: 1\n");
	fprintf(ofp,"%%%%BoundingBox: 0 0 612 792\n");
	fprintf(ofp,"%%%%EndComments\n");
	fprintf(ofp,"%%-------------Variables and procedures------------\n");
	fprintf(ofp,"/i {72 mul} def\n");
	fprintf(ofp,"/width 8.5 i def\n");
	fprintf(ofp,"/height 11 i def\n");
	fprintf(ofp,"/LM -1.5 i def\n");
	fprintf(ofp,"/center {width 2 div height 2 div translate} def\n");
	fprintf(ofp,"/getHel {/Helvetica-Bold findfont 15 scalefont ");
	fprintf(ofp,"setfont} def\n");
	fprintf(ofp,"/getCou {/Courier findfont 12 scalefont setfont} def\n");
	fprintf(ofp,"/getSym {/Symbol findfont 12 scalefont setfont} def\n");
	fprintf(ofp,"/newline{currentpoint 14 sub exch pop LM ");
	fprintf(ofp,"exch moveto} def\n");
	fprintf(ofp,"/label {/subscript exch def\n");
	fprintf(ofp," currentpoint moveto ( q) show\n");
	fprintf(ofp," currentpoint 4 sub moveto subscript show} def\n");
	fprintf(ofp,"/peak \n");
	fprintf(ofp," {/gs exch def /y exch def /x exch def\n");
	fprintf(ofp," x y %d 360 0 arcn\n",RADIUS);
	fprintf(ofp," closepath gs setgray fill} def\n");
	fprintf(ofp,"/normal{/length exch def /slope exch def\n");
	fprintf(ofp,"/y exch def /x exch def\n");
	fprintf(ofp,"gsave newpath x y translate slope rotate 0 0 moveto\n");
	fprintf(ofp,"0 length rlineto [6 4 2 4] 0 setdash 1 setlinewidth\n");
	fprintf(ofp,"0.5 setgray stroke grestore} def\n");
	fprintf(ofp,"%%------------End variables and procedures---------\n");
	fprintf(ofp,"%%%%EndProlog\n");
	fprintf(ofp,"%%%%Page: 1 1\n");
	fprintf(ofp,"%%------------Begin program------------------------\n");
	fprintf(ofp,"center\n\n");
	fprintf(ofp,"getHel\n");
	fprintf(ofp,"-2 i 0 i moveto\n");
	fprintf(ofp,"2 i 0 i lineto\n");
	fprintf(ofp,"(x) label\n");
	fprintf(ofp,"0 i -1 i moveto\n");
	fprintf(ofp,"0 i 4.5 i lineto\n");
	fprintf(ofp,"(z) label\n");
	fprintf(ofp,"1 setlinewidth\n");
	fprintf(ofp,"stroke\n\n");
	fprintf(ofp,"gsave\n");
	return(1);
}

ps_end (ofp)
FILE *ofp;
{
	fprintf(ofp,"\n");
	fprintf(ofp,"showpage\n");
	fprintf(ofp,"%%%%Trailer\n");
	return(1);
}

bilayer (ofp,x,y,rlam,x0,real_con,x0_con,Amp_con,Zh_con,psi)
FILE *ofp;
float x,y;
float rlam,x0,real_con;
float x0_con,Amp_con,Zh_con;
float psi;
{
	fprintf(ofp,"gsave\n");
	fprintf(ofp,"newpath\n");
	fprintf(ofp,"%.3g i %.3g i moveto\n",x,y);
	fprintf(ofp,"%.3g i %.3g i rlineto\n",-Zh_con*sin(psi),
		Zh_con*cos(psi));
	fprintf(ofp,"%.3g i %.3g i rlineto\n",0.5*(rlam-x0)*real_con,
		-0.5*Amp_con);
	fprintf(ofp,"%.3g i %.3g i rlineto\n",x0_con,Amp_con);
	fprintf(ofp,"%.3g i %.3g i rlineto\n",0.5*(rlam-x0)*real_con,
		-0.5*Amp_con);
	fprintf(ofp,"%.3g i %.3g i rlineto\n",2.0*Zh_con*sin(psi),
		-2.0*Zh_con*cos(psi));
	fprintf(ofp,"%.3g i %.3g i rlineto\n",-0.5*(rlam-x0)*real_con,
		0.5*Amp_con);
	fprintf(ofp,"%.3g i %.3g i rlineto\n",-x0_con,-Amp_con);
	fprintf(ofp,"%.3g i %.3g i rlineto\n",-0.5*(rlam-x0)*real_con,
		0.5*Amp_con);
	fprintf(ofp,"closepath 0.75 setgray fill\n");
	fprintf(ofp,"grestore\n\n");
	return(1);
}

profile (ofp,x,y,rlam,x0,real_con,x0_con,Amp_con,Zh_con,psi)
FILE *ofp;
float x,y;
float rlam,x0,real_con;
float x0_con,Amp_con,Zh_con;
float psi;
{
	fprintf(ofp,"gsave\n");
	fprintf(ofp,"newpath\n");
	fprintf(ofp,"%.3g i %.3g i moveto\n",x,y);
	fprintf(ofp,"%.3f i %.3f i rlineto\n",0.5*(rlam-x0)*real_con,
		-0.5*Amp_con);
	fprintf(ofp,"%.3f i %.3f i rlineto\n",x0_con,Amp_con);
	fprintf(ofp,"%.3f i %.3f i rlineto\n",0.5*(rlam-x0)*real_con,
		-0.5*Amp_con);
	fprintf(ofp,"%g setgray %d setlinewidth stroke\n",
		PROFILE_GRAY,2*LINEWIDTH);
	fprintf(ofp,"grestore\n\n");
	return(1);
}
