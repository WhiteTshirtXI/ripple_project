/* 
**	Program: ripple_projection.c
**
**		By Wenjun Sun, 8/14/95, Carnegie Mellon University
**
**		Modified on 9/19/95
**
**		Modefied on 11/17/95
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
	float ds;
	float x0,A;
	float slp_mj,slp_mn;
	float Lmj,Lmn,Lp;
	float alpha;
	float radcon;
	int Nmax;
	int *h, *k;
	float *Fhk;
	int ans;
	float xp;
	float qx, qz;
	float qxp, qzp;
	float p;  

	int i, m;
	float xinc;
	int pres;
	char *infile, *outfile;
	FILE *ifp, *ofp, *fopen();

	radcon = M_PI/180.0;

	printf("Start...........\n");

	printf("Resolution (in pts)? ");
	scanf("%d", &pres);
	printf("x0 and A? ");
	scanf("%f %f",&x0,&A);
	printf("Direction angle (in deg)? ");
	scanf("%f",&alpha);
	alpha=alpha*radcon;

	infile = argv[1];
	ifp = fopen(infile,"r");

	fscanf(ifp,"%f %f %f",&d,&rlam,&gamma);
	gamma *= radcon;
	ds=d/sin(gamma);
	slp_mj=atan(A/x0);
	slp_mn=atan(A/(rlam-x0));
	printf("%.2f %.2f\n",slp_mj/radcon,slp_mn/radcon);
	Lmj=x0/cos(slp_mj);
	Lmn=(rlam-x0)/cos(slp_mn);
	printf("%.2f %.2f\n",Lmj,Lmn);
	Lp=Lmj+Lmn*cos(slp_mj+slp_mn);
	printf("%.2f\n",Lp);
	xinc=0.5*Lp/pres;

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

	for(i=(pres-1);i>0;--i){
		xp = -i*xinc;
		p = 0.0;
		for(m=0;m<Nmax;++m){
			qx=2.0*M_PI*k[m]/rlam;
			qz=2.0*M_PI*(h[m]/d-k[m]/rlam/tan(gamma));
			qxp=qx*cos(alpha)+qz*sin(alpha);
			qzp=-qx*sin(alpha)+qz*cos(alpha);
			p += 2.0*Fhk[m]*cos(qxp*xp)*sin(qzp*ds/2.0)/qzp/ds;
		}
		fprintf(ofp,"%.2f\t%.3f\n",xp,p);
		fprintf(ofp,"%.2f\t%.3f\n",xp+Lp,p);
	}

	for(i=0;i<pres;++i){
		xp = i*xinc;
		p = 0.0;
		for(m=0;m<Nmax;++m){
			qx=2.0*M_PI*k[m]/rlam;
			qz=2.0*M_PI*(h[m]/d-k[m]/rlam/tan(gamma));
			qxp=qx*cos(alpha)+qz*sin(alpha);
			qzp=-qx*sin(alpha)+qz*cos(alpha);
			p += 2.0*Fhk[m]*cos(qxp*xp)*sin(qzp*ds/2.0)/qzp/ds;
		}
		fprintf(ofp,"%.2f\t%.3f\n",xp,p);
		fprintf(ofp,"%.2f\t%.3f\n",xp+Lp,p);
	}

	printf("Number of files: %d\n\n",argc);
	printf("Done...........\n");
	fclose(ofp);
	exit (0);
}
