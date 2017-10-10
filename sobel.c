#include <stdio.h>                          /* Sobel.c */
#include <math.h>
#include <stdlib.h>

void imageToFile(char*, int[256][256]);

int pic[256][256];
int outpicx[256][256];
int outpicy[256][256];
int lowOutPic[256][256];
int highOutPic[256][256];
int maskx[3][3] = {{-1,0,1},{-2,0,2},{-1,0,1}};
int masky[3][3] = {{1,2,1},{0,0,0},{-1,-2,-1}};
double ival[256][256],maxival;


int main(argc,argv)
int argc;
char **argv;
{
        int i,j,p,q,mr,sum1,sum2;
        double hThresh = 110.0, lThresh = 40.0;
        FILE *fo1, *fo2, *fp1, *fopen();
        char *foobar, *pgmType;
				int picLength, picWidth;
				// input image file to read as fp1
        argc--; argv++;
        foobar = *argv;
        fp1=fopen(foobar,"rb");
				// original sobel output as fo1
				fo1=fopen("sobelMag.pgm","wb");

				fscanf(fp1, "%s", pgmType);
				fscanf(fp1, "%d", &picLength);
				fscanf(fp1, "%d", &picWidth);
				// processing image bytes
        for (i=0;i<256;i++)
        { for (j=0;j<256;j++)
                {
                  pic[i][j]  =  getc(fp1);
                  pic[i][j]  &= 0377;
                }
        }
				fclose(fp1);
        mr = 1;
        for (i=mr;i<256-mr;i++)
        { for (j=mr;j<256-mr;j++)
          {
             sum1 = 0;
             sum2 = 0;

             for (p=-mr;p<=mr;p++)
             {
                for (q=-mr;q<=mr;q++)
                {	// applying sobel filter
                   sum1 += pic[i+p][j+q] * maskx[p+mr][q+mr];
                   sum2 += pic[i+p][j+q] * masky[p+mr][q+mr];
                }
             }
						 // applying x/y gradient to pixel maps
             outpicx[i][j] = sum1;
             outpicy[i][j] = sum2;
          }
        }
				// Calculating the gradient magnitude for the image
        maxival = 0;
        for (i=mr;i<256-mr;i++)
        { for (j=mr;j<256-mr;j++)
          {
             ival[i][j]=sqrt((double)((outpicx[i][j]*outpicx[i][j]) +
                                      (outpicy[i][j]*outpicy[i][j])));
             if (ival[i][j] > maxival) // getting max value for
                maxival = ival[i][j];

           }
        }

				// printing magnitude of gradient to image
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
				fclose(fo1);
				imageToFile("sobelXGrad.pgm", outpicx);
				imageToFile("sobelYGrad.pgm", outpicy);
				imageToFile("sobelHiThresh.pgm", highOutPic);
				imageToFile("sobelLoThresh.pgm", lowOutPic);

				return 0;
} // end main

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
	fclose(fOut);
} // end imageToFile
