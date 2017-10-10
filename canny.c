#include <stdio.h>                          /* Sobel.c */
#include <math.h>
#include <stdlib.h>

#define  PICSIZE 256
#define  MAXMASK 100
#define  SIGONE 1
#define  SIGTWO 2

void imageToFile(char*, int[256][256]);
void candidates(double image[][PICSIZE]);


int    pic[PICSIZE][PICSIZE];
int outpicx[PICSIZE][PICSIZE];
int outpicy[PICSIZE][PICSIZE];
int    edgeflag[PICSIZE][PICSIZE];
double xmask[MAXMASK][MAXMASK], ymask[MAXMASK][MAXMASK];
double conv[PICSIZE][PICSIZE];
double ival[256][256];
int lowOutPic[256][256];
int highOutPic[256][256];
double hThresh = 110.0, lThresh = 40.0;

int main(int argc, char** argv){
	int     i,j,p,q,s,t,mr,centx,centy;
	double  xmaskval,ymaskval,sum1,sum2,sig,maxival,minival,maxval,ZEROTOL;
	FILE    *fo1, *fo2,*fp1, *fopen();
	char    *foobar, *pgmType;
	int picLength, picWidth;

	argc--; argv++;
	foobar = *argv;
	fp1=fopen(foobar,"rb");


	fo1=fopen("cannyMag.pgm","wb");


	fo2=fopen("canny2output.pgm","wb");

	argc--; argv++;
	foobar = *argv;
	sig = atof(foobar);

	mr = (int)(sig * 3);
	centx = (MAXMASK / 2);
	centy = (MAXMASK / 2);

	fscanf(fp1, "%s", pgmType);
	fscanf(fp1, "%d", &picLength);
	fscanf(fp1, "%d", &picWidth);
	for (i=0;i<256;i++)
	{ for (j=0;j<256;j++)
					{
						pic[i][j]  =  getc (fp1);
					}
	}
// convolution to create the x and y mask
	for (p=-mr;p<=mr;p++)
	{  for (q=-mr;q<=mr;q++)
		 {
				xmaskval = (p * (exp(-1*(((p*p)+(q*q))/(2*(sig*sig))))));
				ymaskval = (q * (exp(-1*(((p*p)+(q*q))/(2*(sig*sig))))));
				(xmask[p+centy][q+centx]) = xmaskval;
				(ymask[p+centy][q+centx]) = ymaskval;
		 }
	}
	for (i=mr;i<=255-mr;i++)
	{ for (j=mr;j<=255-mr;j++)
		{
			 sum1 = 0;
			 sum2 = 0;
			 for (p=-mr;p<=mr;p++)
			 {
					for (q=-mr;q<=mr;q++)
					{
						 sum1 += pic[i+p][j+q] * ymask[p+centy][q+centx];
						 sum2 += pic[i+p][j+q] * xmask[p+centy][q+centx];
					}
			 }
			 outpicx[i][j] = (int)sum1;
			 outpicy[i][j] = (int)sum2;
		}
	}
// sobel code
maxival = 0.0; // value for scaling
// getting magnitude for each pixel
for (i=mr;i<256-mr;i++)
{ for (j=mr;j<256-mr;j++)
	{
		 ival[i][j]=sqrt((double)((outpicx[i][j]*outpicx[i][j]) +
															(outpicy[i][j]*outpicy[i][j])));
		 if (ival[i][j] > maxival)
				maxival = ival[i][j];

	 }
}

candidates(ival);
// printing magnitude imageToFile

fprintf(fo1,"P5\n256 256\n255\n");  /* Print out the .pgm header */
for (i=0;i<256;i++)
{ for (j=0;j<256;j++)
		{
		 ival[i][j] = (ival[i][j] / maxival) * 255;
		 // applying high threshold
		 if(ival[i][j] > hThresh) highOutPic[i][j] = 255;
		 else highOutPic[i][j] = 0;
		 // applying low threshold
		 if(ival[i][j] > lThresh) lowOutPic[i][j] = 255;
		 else lowOutPic[i][j] = 0;
		 fprintf(fo1,"%c",(char)((int)(ival[i][j])));

		}
}


imageToFile("cannyXGrad.pgm", outpicx);
imageToFile("cannyYGrad.pgm", outpicy);

 return 0;
} // end main

void candidates(double image[][PICSIZE])
{
	int i, j;
	for(i = 0; i < PICSIZE; i++)
	{	for(j = 0; j < PICSIZE; j++)
		{
			if(image[i][j] <= lThresh)
			{
				edgeflag[i][j] = 0;
			}
			else if(image[i][j] >= hThresh)
			{
				edgeflag[i][j] = 1;
			}
			else{ //
				edgeflag[i][j] = 1;
			}
		}
	}
} // end candidates

void imageToFile(char *file, int output[256][256]){
	int i, j;
	FILE *fOut = fopen(file, "wb");
	fprintf(fOut,"P5\n256 256\n255\n");  /* Print out the .pgm header */
	for (i=0;i<256;i++)
	{ for (j=0;j<256;j++)
			{
			 fprintf(fOut,"%c", ((char)output[i][j]));
			}
	}
} // end imageToFile
