#include <stdio.h>
#include <math.h>
#define  PICSIZE 256
#define  MAXMASK 100
#define  ON      255
#define  OFF     0


           int    pic[PICSIZE][PICSIZE];

main(int argc,char **argv)
{

        int i,j,k,p,q,mr,centx,centy,hi,lo,moved,histocount,pic2mod10;
        int maxvotes,index;

        int maghisto[256],labelindex,votewinindex;
        double maskvalx,maskvaly,sig,maxival,sum1,sum2,tandir,threshper,label,angle;

         int marker1, max_marker1, marker2, max_marker2;
         int delta_angle;
         double i_angle, shift_start;

        char throwaway[80];

        FILE    *fo1, *fo2,*fo3,*fo4,*fo5,*fo6,*fo7,*fo8,*fp1, *fopen();
      char    *foobar;
        argc--; argv++;
        foobar = *argv;
        fp1=fopen(foobar,"rb");

        argc--; argv++;
        foobar = *argv;
        fo1=fopen(foobar,"wb");


        fgets(throwaway, 80, fp1);
        fgets(throwaway, 80, fp1);
        fgets(throwaway, 80, fp1);
if ( !( (throwaway[0]=='2') && (throwaway[1]=='5') && (throwaway[2]=='5')))
        fgets(throwaway, 80, fp1);


          for (i=0;i<256;i++)
        {
                for (j=0;j<256;j++)
                {
                  pic[i][j]  =  getc (fp1);
                  pic[i][j]  &= 0377;
                  printf("%c", (char)((int)(pic[i][j])));
                }
        }

}
