/* 
**	Program: ripple_patterson.c
**
**		By Wenjun Sun, 8/14/95, Carnegie Mellon University
**
**	Purpose: To reconstruct the 2-D electron density
**	patterson function using the experimental amplitudes of the
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
**		%ripple_patterson form_factor_data electron_density_map
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
	int pres, gs;
	char *infile, *outfile;
	FILE *ifp, *ofp, *fopen();

	radcon = M_PI/180.0;

	printf("Start...........\n");

	printf("Size of picture (inch) and Resolution (in pts):\n");
	scanf("%f %d", &size, &pres);
	printf("Number of ripple wavelength in the picture: ");
	scanf("%f",&multiple);

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

	fscanf(ifp,"%d",&Nmax);
	
	for(m=0;m<Nmax;++m){
		fscanf(ifp,"%d %d %f",&h,&k,&Fhk);
		qx=2.0*M_PI*k/rlam;
		qz=2.0*M_PI*(h/d-k/rlam/tan(gamma));
		for(i=0;i<pres;++i){
			x = i*xinc;
			for(j=0;j<pres;++j){
				z = j*xinc;
				p[i][j] += Fhk*Fhk*cos(qx*x+qz*z);
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

	ps_header(ofp, pres, size);

	for (i=0; i < pres; i++) {
		for (m=0; m < (pres-24); m=m+25) { 
			for (j=m; j < (m+25); j++) { 
				gs = 255.0*(p[j][i]-gsmin)/(gsmax-gsmin);
				fprintf(ofp, "%02x",gs);
			}
			fprintf(ofp,"\n"); 
		} 
	}

	ps_end(ofp); 

	printf("Number of files: %d\n\n",argc);
	printf("Done...........\n");
	fclose(ofp);
	exit (0);
}

ps_header (ofp, num, size)
FILE *ofp;
int num;
float size;
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
	fprintf(ofp,"/pix %d string def\n",num);
	fprintf(ofp,"%%------------End variables and procedures---------\n");
	fprintf(ofp,"%%%%EndProlog\n");
	fprintf(ofp,"%%%%Page: 1 1\n");
	fprintf(ofp,"%%------------Begin program------------------------\n");
	fprintf(ofp,"gsave\n");
	fprintf(ofp,"center\n");
	fprintf(ofp,"-%.1f i dup neg translate\n", size/2.0);
	fprintf(ofp,"%.1f i dup scale\n", size);
	fprintf(ofp,"%d %d 8 [%d 0 0 %d 0 %d]\n", num, num, num, num, num);
	fprintf(ofp,"{currentfile pix readhexstring pop} image\n\n");
	return(1);
}

ps_end (ofp)
FILE *ofp;
{
	fprintf(ofp,"\n");
	fprintf(ofp,"grestore\n");
	fprintf(ofp,"showpage\n");
	fprintf(ofp,"%%%%Trailer\n");
	return(1);
}
