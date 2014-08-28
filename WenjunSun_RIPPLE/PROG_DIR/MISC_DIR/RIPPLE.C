/* Given a (primary ripple repeat), b(secondary ripple repeat) and c (
** bilayer stacking repeat) and gamma (monoclinic angle), calculate
** d(hkl) and compare with Hatta et. al.'s results. 
**
**	Usage: ripple output_filename
*/

#include <stdio.h>
#include <math.h>

#define PI 3.1415926

void main(argc, argv)
int argc; char *argv[];
{
  float a, b, c, r, d;
  int h, k, l;
  float radcon; 
  FILE *ofpo, *fopen();

  radcon=PI/180.0;

  ofpo = fopen(argv[1], "w");
  printf("Start...........\n");
  while ( getchar() != EOF){
  	printf("Read in a, b, c (in A), gamma (in deg): \n");
  	scanf("%f %f %f %f", &a, &b, &c, &r);
  	r=r*radcon;
	
	fprintf(ofpo,"\na= %.1f, b= %.1f, c= %.1f, gamma= %.1f\n\n",a,b,c,r/radcon);
	fprintf(ofpo,"h\tk\tl\td(hkl)\n\n");

	for(h=0; h<5; h++){
		for(k=0; k<5; k++){
			for(l=0; l<5; l++){
				if( (h+k+l) != NULL){
					d=h*h/a/a/sin(r)/sin(r);
					d=d+k*k/b/b;
					d=d+l*l/c/c/sin(r)/sin(r);
					d=d-2.*h*l*cos(r)/a/c/sin(r)/sin(r);
					d=1.0/sqrt(d);
					printf("%d\t%d\t%d\t%.1f\n",h,k,l,d);
					fprintf(ofpo,"%d\t%d\t%d\t%.1f\n",h,k,l,d);
				}
			}
		}
	}
  }

  printf("Done...........\n");
  fclose(ofpo);
  exit (0);
}
