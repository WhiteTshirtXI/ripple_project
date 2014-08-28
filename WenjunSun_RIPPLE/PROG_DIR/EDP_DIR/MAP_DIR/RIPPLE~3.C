/* 
**	Program: ripple_edp.c
**
**		By Wenjun Sun, 8/14/95, Carnegie Mellon University
**
**	Purpose: To reconstruct the 2-D electron density
**	map using the experimental amplitudes of the
**	form factors, and the phases obtained from the
**	model proposed by Wenjun Sun.
**
**	Needed inputs:
**			Ripple wavelength, lamellar spacing
**			Oblique angle of the ripple unit cell
**			indices and form factors of ripple
**			scattering peaks.
**
**	Output:
**		A ready-for-print PS file showing the 2D electron
**		density distribution in gray scales from 0 to
**		256.	
**	Usage:
**		%ripple_edp form_factor_data electron_density_map
**
*/

#include <stdio.h>
#include <math.h>

void main(argc, argv)
int argc; 
char *argv[];
{
	float d, rlam, gamma;
	float radcon;
	int Nmax;
	int h, k;
	float Fhk;
	float x, z;
	float qx, qz;
	float **p;  /* pointer to a 2-D array to store the pixel matrix */

	int i, j, m;
	float multiple;
	float xinc;
	float gsmin, gsmax, size;
	float bar, barlen;
	float amp,x0;
	int pres, gs;
	char *infile, *outfile;
	FILE *ifp, *ofp, *fopen();

	radcon = M_PI/180.0;

	printf("Start...........\n");

	printf("Size of picture (inch) and Resolution (in pts):\n");
	scanf("%f %d", &size, &pres);
	printf("Number of ripple wavelength in the picture: ");
	scanf("%f",&multiple);
	printf("Dimension bar length (A)? ");
	scanf("%f",&bar);
	printf("amp,x0? ");
	scanf("%f %f",&amp,&x0);

	p = (float **) malloc(pres*sizeof(float *));
	for (i=0; i < pres; i++) {
		p[i] = (float *) malloc(pres*sizeof(float));
		for (j=0; j < pres; j++) {
			p[i][j]=0.0;
		}
	}

	infile = argv[1];
	ifp = fopen(infile,"r");

	fscanf(ifp,"%f %f %f",&d,&rlam,&gamma);
	xinc=(multiple*rlam)/pres;
	gamma *= radcon;
	barlen=size*bar/(multiple*rlam);

	fscanf(ifp,"%d",&Nmax);
	
	for(m=0;m<Nmax;++m){
		fscanf(ifp,"%d %d %f",&h,&k,&Fhk);
		qx=2.0*M_PI*k/rlam;
		qz=2.0*M_PI*(h/d-k/rlam/tan(gamma));
		for(i=0;i<pres;++i){
			x = i*xinc;
			for(j=0;j<pres;++j){
				z = j*xinc;
				p[i][j] += Fhk*cos(qx*(x-0.5*x0)+qz*(z-0.5*amp));
			}
		}
	}
	fclose(ifp);

	outfile = argv[2];
	ofp = fopen(outfile, "w");

	gsmax=0.0;
	gsmin=10000.0;
	for (i=0; i < pres; i++) {
		for (j=0; j < pres; j++) {
			if ( gsmax <= p[i][j] ) gsmax=p[i][j];
			if ( gsmin >= p[i][j] ) gsmin=p[i][j];
		}
	}
	printf("%g %g\n",gsmax,gsmin);

	ps_header(ofp, pres, size, bar);

	for (i=0; i < pres; i++) {
		for (m=0; m < (pres-24); m=m+25) { 
			for (j=m; j < (m+25); j++) { 
				gs = 255.0*(p[j][i]-gsmin)/(gsmax-gsmin);
				fprintf(ofp, "%02x",gs);
			}
			fprintf(ofp,"\n"); 
		} 
	}

	ps_end(ofp,size,barlen); 

	printf("Number of files: %d\n\n",argc);
	printf("Done...........\n");
	fclose(ofp);
	exit (0);
}

ps_header (ofp, num, size, bar)
FILE *ofp;
int num;
float size, bar;
{
	fprintf(ofp,"%%!PS-Adobe-2.0\n");
	fprintf(ofp,"%%%%Creator: Wenjun Sun\n");
	fprintf(ofp,"%%%%Title: Ripple EDP\n");
	fprintf(ofp,"%%%%Pages: 1\n");
	fprintf(ofp,"%%%%BoundingBox: 0 0 612 792\n");
	fprintf(ofp,"%%%%EndComments\n");
	fprintf(ofp,"%%-------------Variables and procedures------------\n");
	fprintf(ofp,"/i {72 mul} def\n");
	fprintf(ofp,"/width 8.5 i def\n");
	fprintf(ofp,"/height 11 i def\n");
	fprintf(ofp,"/center {width 2 div height 2 div translate} def\n");
fprintf(ofp,"/barbox{/barlen exch def /size exch def\n");
fprintf(ofp,"/half size 2 div i def /quarter size 4 div i def\n");
fprintf(ofp,"quarter 10 sub half neg 10 add translate\n");
fprintf(ofp,"newpath 0 0 moveto quarter 0 rlineto 0 quarter 2 div rlineto\n");
fprintf(ofp,"quarter neg 0 rlineto closepath 1 setgray fill\n");
fprintf(ofp,"5 quarter 4 div moveto barlen i 0 rlineto\n");
fprintf(ofp,"5 setlinewidth 0 setgray stroke\n");
fprintf(ofp,"/Helvetica-Bold findfont 18 scalefont setfont\n");
fprintf(ofp,"quarter 2 div 8 add quarter 4 div -6 add moveto\n");
fprintf(ofp,"(%.0fA) show currentpoint 14 add exch -11 add exch moveto\n",bar);
fprintf(ofp,"/Helvetica-Bold findfont 14 scalefont setfont\n");
fprintf(ofp,"(o) show} def\n");
	fprintf(ofp,"/pix %d string def\n",num);
	fprintf(ofp,"%%------------End variables and procedures---------\n");
	fprintf(ofp,"%%%%EndProlog\n");
	fprintf(ofp,"%%%%Page: 1 1\n");
	fprintf(ofp,"%%------------Begin program------------------------\n");
	fprintf(ofp,"center\n");
	fprintf(ofp,"gsave\n");
	fprintf(ofp,"-%.1f i dup neg translate\n", size/2.0);
	fprintf(ofp,"%.1f i dup scale\n", size);
	fprintf(ofp,"%d %d 8 [%d 0 0 %d 0 %d]\n", num, num, num, num, num);
	fprintf(ofp,"{currentfile pix readhexstring pop} image\n\n");
	return(1);
}

ps_end (ofp,size,barlen)
FILE *ofp;
float size,barlen;
{
	fprintf(ofp,"\n");
	fprintf(ofp,"grestore\n");
	fprintf(ofp,"%.1f %.3f barbox\n",size,barlen);
	fprintf(ofp,"showpage\n");
	fprintf(ofp,"%%%%Trailer\n");
	return(1);
}
