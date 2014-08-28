/* 
**	Program: map_profile.c
**
**		By Wenjun Sun, 11/29/95, Carnegie Mellon University
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
**		A data file and a PS file for the ripple profile.
**	Usage:
**		%map_profile form_factor_data profile.dat profile.ps
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
	float n_bilayer;
	float size;
	float factor;
	float xinc, zinc;
	float *zpmax,*znmax,pmax;
	int pres;
	char *infile, *outfile1, *outfile2;
	FILE *ifp, *ofp1, *ofp2, *fopen();

	radcon = M_PI/180.0;

	printf("Start...........\n");

	printf("Resolution (in pts)? ");
	scanf("%d", &pres);
	printf("Size of picture (inch)? ");
	scanf("%f",&size);
	printf("Number of ripple wavelength? ");
	scanf("%f",&multiple);
	printf("Which bilayer? ");
	scanf("%f",&n_bilayer);
	printf("Origin (1,2,3,4)? ");
	scanf("%d",&ans);
	printf("x0, alpha and beta? ");
	scanf("%f %f %f",&x0,&alpha,&beta);
	alpha=alpha*radcon;
	beta=beta*radcon;

	infile = argv[1];
	ifp = fopen(infile,"r");

	zpmax = (float *)malloc(pres*sizeof(float));
	znmax = (float *)malloc(pres*sizeof(float));

	fscanf(ifp,"%f %f %f",&d,&rlam,&gamma);
	factor=72.0*size/(2.0*rlam);
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

	outfile1 = argv[2];
	ofp1 = fopen(outfile1, "w");

	for(i=0;i<pres;++i){
		x = i*xinc;
		if( x>=0.0 && x<=(0.5*x0)){
			z0=x*tan(alpha);
			z0 += n_bilayer*d;
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
				if(p>=pmax) {pmax=p; zpmax[i]=z;};
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
				if(p>=pmax) {pmax=p; znmax[i]=z;};
			}
		fprintf(ofp1,"%.2f\t%.3f\t%.3f\n",x,zpmax[i],znmax[i]);
		}
		if( x>(0.5*x0) && x<=(rlam-0.5*x0)){
			z0=0.5*x0*tan(alpha)-(x-0.5*x0)*tan(beta);
			z0 += n_bilayer*d;
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
				if(p>=pmax) {pmax=p; zpmax[i]=z;};
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
				if(p>=pmax) {pmax=p; znmax[i]=z;};
			}
		fprintf(ofp1,"%.2f\t%.3f\t%.3f\n",x,zpmax[i],znmax[i]);
		}
		if( x>(rlam-0.5*x0) && x<=(rlam+0.5*x0)){
			z0=(x-rlam)*tan(alpha);
			z0 += n_bilayer*d;
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
				if(p>=pmax) {pmax=p; zpmax[i]=z;};
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
				if(p>=pmax) {pmax=p; znmax[i]=z;};
			}
		fprintf(ofp1,"%.2f\t%.3f\t%.3f\n",x,zpmax[i],znmax[i]);
		}

		if( x>(rlam+0.5*x0) && x<=(2.0*rlam-0.5*x0)){
			z0=0.5*x0*tan(alpha)-(x-(rlam+0.5*x0))*tan(beta);
			z0 += n_bilayer*d;
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
				if(p>=pmax) {pmax=p; zpmax[i]=z;};
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
				if(p>=pmax) {pmax=p; znmax[i]=z;};
			}
		fprintf(ofp1,"%.2f\t%.3f\t%.3f\n",x,zpmax[i],znmax[i]);
		}

		if( x>(2.0*rlam-0.5*x0) && x<=(2.0*rlam)){
			z0=(x-2.0*rlam)*tan(alpha);
			z0 += n_bilayer*d;
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
				if(p>=pmax) {pmax=p; zpmax[i]=z;};
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
				if(p>=pmax) {pmax=p; znmax[i]=z;};
			}
		fprintf(ofp1,"%.2f\t%.3f\t%.3f\n",x,zpmax[i],znmax[i]);
		}
	}
	fclose(ofp1);

	outfile2 = argv[3];
	ofp2 = fopen(outfile2, "w");

	
	fprintf(ofp2,"%.1f dup translate\n",-0.5*size*72.0);
	fprintf(ofp2,"newpath\n");
	if(n_bilayer == 0.0){	
		x=0.0;
		zpmax[0] *= factor;
		fprintf(ofp2,"%.1f\t%.1f moveto\n",x,zpmax[0]);
		for(i=1;i<pres;++i){
			x = i*xinc*factor;
			zpmax[i] *= factor;
			fprintf(ofp2,"%.1f\t%.1f lineto\n",x,zpmax[i]);
		}
		fprintf(ofp2,"0 setgray 2 setlinewidth stroke\n");
	}
	if(n_bilayer >= 1.0){	
                x=0.0;
                zpmax[0] *= factor;
                fprintf(ofp2,"%.1f\t%.1f moveto\n",x,zpmax[0]);
                for(i=1;i<pres;++i){
                        x = i*xinc*factor;
                        zpmax[i] *= factor;
                        fprintf(ofp2,"%.1f\t%.1f lineto\n",x,zpmax[i]);
                }
                fprintf(ofp2,"0 setgray 2 setlinewidth stroke\n");
                x=0.0;
                znmax[1] *= factor;
                fprintf(ofp2,"%.1f\t%.1f moveto\n",x,znmax[0]);
                for(i=1;i<pres;++i){
                        x = i*xinc*factor;
                        znmax[i] *= factor;
                        fprintf(ofp2,"%.1f\t%.1f lineto\n",x,znmax[i]);
                }
                fprintf(ofp2,"0 setgray 2 setlinewidth stroke\n");
	}
	fclose(ofp2);

	printf("Done...........\n");
	exit (0);
}
