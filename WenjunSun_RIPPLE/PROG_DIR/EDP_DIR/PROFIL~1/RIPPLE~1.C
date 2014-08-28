/* 
**	Program: ripple_profile.c
**
**		By Wenjun Sun, 8/14/95, Carnegie Mellon University
**
**		Modified on 9/19/95
**		Modefied on 11/20/95
**
**	Purpose: To obtain the trans-bilayer electron density
**	profile at any given direction from the reconstructed 
**	2D electron density map using the experimental amplitudes 
**	of the form factors, and the phases obtained from the
**	model proposed by Wenjun Sun.
**
**	Needed inputs:
**			Ripple wavelength, lamellar spacing
**			Oblique angle of the ripple unit cell
**			indices and form factors of ripple
**			scattering peaks.
**
**	Output:
**		A ready-for-print PS file showing the electron
**		density distribution in gray scales from 0 to
**		256.	
**	Usage:
**		%ripple_edp form_factor_data electron_density_profile
**
*/
 

#include <stdio.h>
#include <math.h>

void main(argc, argv)
int argc; 
char *argv[];
{
	float d, rlam, gamma;
	float x0, alpha, beta;
	float radcon;
	int Nmax;
	int *h, *k;
	float *Fhk;
	int ans;
	float x_o,z_o;
	float x, z, z0;
	float qx, qz;
	float p;  

	int i, j, m;
	float multiple;
	float xinc, zinc;
	float zpmax,znmax,pmax;
	int pres;
	char *infile, *outfile;
	FILE *ifp, *ofp, *fopen();

	radcon = M_PI/180.0;

	printf("Start...........\n");

	printf("Resolution (in pts)? ");
	scanf("%d", &pres);
	printf("Number of ripple wavelength? ");
	scanf("%f",&multiple);
	printf("Origin (1,2,3,4)? ");
	scanf("%d",&ans);
	printf("x0, alpha and beta? ");
	scanf("%f %f %f",&x0,&alpha,&beta);
	alpha=alpha*radcon;
	beta=beta*radcon;

	infile = argv[1];
	ifp = fopen(infile,"r");

	fscanf(ifp,"%f %f %f",&d,&rlam,&gamma);
	xinc=multiple*rlam/pres;
	zinc=0.5*d/pres;
	gamma *= radcon;
	switch(ans){
		case 1:
			x_o=0.0;
			z_o=0.0;
			break;
		case 2:
			x_o=0.5*rlam;
			z_o=0.0;
			break;
		case 3:
			x_o=0.5*d/tan(gamma);
			z_o=0.5*d;
			break;
		case 4:
			x_o=0.5*rlam+0.5*d/tan(gamma);
			z_o=0.5*d;
			break;
		default:
			printf("Origin (x_o,z_o)? ");
			scanf("%f %f",&x_o,x_o,z_o);
			break;
	}

	fscanf(ifp,"%d",&Nmax);
	h = (int *) malloc(Nmax*sizeof(int));
	k = (int *) malloc(Nmax*sizeof(int));
	Fhk = (float *) malloc(Nmax*sizeof(float));
	
	for(m=0;m<Nmax;++m){
		fscanf(ifp,"%d %d %f",&h[m],&k[m],&Fhk[m]);
	}
	fclose(ifp);

	outfile = argv[2];
	ofp = fopen(outfile, "w");

	for(i=0;i<pres;++i){
		x = i*xinc;
		if( x>=0.0 && x<=(0.5*x0)){
			z0=x*tan(alpha);
			z0 += d;
			pmax=0.0;
			for(j=0;j<pres;++j){
				z = z0+j*zinc;
				p = 0.0;
				for(m=0;m<Nmax;++m){
					qx=2.0*M_PI*k[m]/rlam;
					qz=2.0*M_PI*(h[m]/d-k[m]/rlam
						/tan(gamma));
					p += Fhk[m]*cos(qx*(x-x_o)+qz*(z-z_o));
				}
				if(p>=pmax) {pmax=p; zpmax=z;};
			}
			pmax=0.0;
			for(j=0;j<pres;++j){
				z = z0-j*zinc;
				p = 0.0;
				for(m=0;m<Nmax;++m){
					qx=2.0*M_PI*k[m]/rlam;
					qz=2.0*M_PI*(h[m]/d-k[m]/rlam
						/tan(gamma));
					p += Fhk[m]*cos(qx*(x-x_o)+qz*(z-z_o));
				}
				if(p>=pmax) {pmax=p; znmax=z;};
			}
		fprintf(ofp,"%.2f\t%.3f\t%.3f\n",x,zpmax,znmax);
		}
		if( x>(0.5*x0) && x<=(rlam-0.5*x0)){
			z0=0.5*x0*tan(alpha)-(x-0.5*x0)*tan(beta);
			z0 += d;
			pmax=0.0;
			for(j=0;j<pres;++j){
				z = z0+j*zinc;
				p = 0.0;
				for(m=0;m<Nmax;++m){
					qx=2.0*M_PI*k[m]/rlam;
					qz=2.0*M_PI*(h[m]/d-k[m]/rlam
						/tan(gamma));
					p += Fhk[m]*cos(qx*(x-x_o)+qz*(z-z_o));
				}
				if(p>=pmax) {pmax=p; zpmax=z;};
			}
			pmax=0.0;
			for(j=0;j<pres;++j){
				z = z0-j*zinc;
				p = 0.0;
				for(m=0;m<Nmax;++m){
					qx=2.0*M_PI*k[m]/rlam;
					qz=2.0*M_PI*(h[m]/d-k[m]/rlam
						/tan(gamma));
					p += Fhk[m]*cos(qx*(x-x_o)+qz*(z-z_o));
				}
				if(p>=pmax) {pmax=p; znmax=z;};
			}
		fprintf(ofp,"%.2f\t%.3f\t%.3f\n",x,zpmax,znmax);
		}
		if( x>(rlam-0.5*x0) && x<=(rlam+0.5*x0)){
			z0=(x-rlam)*tan(alpha);
			z0 += d;
			pmax=0.0;
			for(j=0;j<pres;++j){
				z = z0+j*zinc;
				p = 0.0;
				for(m=0;m<Nmax;++m){
					qx=2.0*M_PI*k[m]/rlam;
					qz=2.0*M_PI*(h[m]/d-k[m]/rlam
						/tan(gamma));
					p += Fhk[m]*cos(qx*(x-x_o)+qz*(z-z_o));
				}
				if(p>=pmax) {pmax=p; zpmax=z;};
			}
			pmax=0.0;
			for(j=0;j<pres;++j){
				z = z0-j*zinc;
				p = 0.0;
				for(m=0;m<Nmax;++m){
					qx=2.0*M_PI*k[m]/rlam;
					qz=2.0*M_PI*(h[m]/d-k[m]/rlam
						/tan(gamma));
					p += Fhk[m]*cos(qx*(x-x_o)+qz*(z-z_o));
				}
				if(p>=pmax) {pmax=p; znmax=z;};
			}
		fprintf(ofp,"%.2f\t%.3f\t%.3f\n",x,zpmax,znmax);
		}

		if( x>(rlam+0.5*x0) && x<=(2.0*rlam-0.5*x0)){
			z0=0.5*x0*tan(alpha)-(x-(rlam+0.5*x0))*tan(beta);
			z0 += d;
			pmax=0.0;
			for(j=0;j<pres;++j){
				z = z0+j*zinc;
				p = 0.0;
				for(m=0;m<Nmax;++m){
					qx=2.0*M_PI*k[m]/rlam;
					qz=2.0*M_PI*(h[m]/d-k[m]/rlam
						/tan(gamma));
					p += Fhk[m]*cos(qx*(x-x_o)+qz*(z-z_o));
				}
				if(p>=pmax) {pmax=p; zpmax=z;};
			}
			pmax=0.0;
			for(j=0;j<pres;++j){
				z = z0-j*zinc;
				p = 0.0;
				for(m=0;m<Nmax;++m){
					qx=2.0*M_PI*k[m]/rlam;
					qz=2.0*M_PI*(h[m]/d-k[m]/rlam
						/tan(gamma));
					p += Fhk[m]*cos(qx*(x-x_o)+qz*(z-z_o));
				}
				if(p>=pmax) {pmax=p; znmax=z;};
			}
		fprintf(ofp,"%.2f\t%.3f\t%.3f\n",x,zpmax,znmax);
		}

		if( x>(2.0*rlam-0.5*x0) && x<=(2.0*rlam)){
			z0=(x-2.0*rlam)*tan(alpha);
			z0 += d;
			pmax=0.0;
			for(j=0;j<pres;++j){
				z = z0+j*zinc;
				p = 0.0;
				for(m=0;m<Nmax;++m){
					qx=2.0*M_PI*k[m]/rlam;
					qz=2.0*M_PI*(h[m]/d-k[m]/rlam
						/tan(gamma));
					p += Fhk[m]*cos(qx*(x-x_o)+qz*(z-z_o));
				}
				if(p>=pmax) {pmax=p; zpmax=z;};
			}
			pmax=0.0;
			for(j=0;j<pres;++j){
				z = z0-j*zinc;
				p = 0.0;
				for(m=0;m<Nmax;++m){
					qx=2.0*M_PI*k[m]/rlam;
					qz=2.0*M_PI*(h[m]/d-k[m]/rlam
						/tan(gamma));
					p += Fhk[m]*cos(qx*(x-x_o)+qz*(z-z_o));
				}
				if(p>=pmax) {pmax=p; znmax=z;};
			}
		fprintf(ofp,"%.2f\t%.3f\t%.3f\n",x,zpmax,znmax);
		}
	}

	printf("Number of files: %d\n\n",argc);
	printf("Done...........\n");
	fclose(ofp);
	exit (0);
}
