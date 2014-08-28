/*
**	Program: average.c
**		By Wenjun Sun (8/16/95)
**
**	Purpose: Get the average of a few indentical scans
**
**	Usage: % average raw_data_files smooth_data_file
**
*/

#include <stdio.h>

void main(argc, argv)
int argc; 
char *argv[];
{
	char c;
	int i;
	FILE **ifp, *ofp, *fopen();

	int bin,count,monitor;
	float tth,temp,pressure;
	int average_cnt;

	ifp = (FILE **) malloc((argc-2)*sizeof(FILE *));
	for(i=0;i<(argc-2);i++){
		ifp[i] = fopen(argv[i+1],"r");
	}
	ofp = fopen(argv[argc-1],"w");

	while((c=getc(ifp[0])) != EOF){
		ungetc(c,ifp[0]);
		average_cnt = 0;
		for(i=0;i<(argc-2);i++){
			fscanf(ifp[i],"%d %d %d %f %f %f",&bin,&count,	
				&monitor,&tth,&temp,&pressure);
			average_cnt += count;
		}
		fprintf(ofp,"%.5f %d\n",tth,average_cnt);
	}

        for(i=0;i<(argc-2);i++){
                fclose(ifp[i]);
        }
	fclose(ofp);
	exit(0);
}
