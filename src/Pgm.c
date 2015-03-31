#include <pgm.h>
#include <stdio.h>

//gcc -I /usr/local/netpbm/include/ -L /usr/local/netpbm/lib -o Pgm Pgm.c -lnetpbm -lm

int main(int argc, char *argv[])
{
    //2d array that contains image's data
    gray **image;

    //Maximum value of our input image, probably
    gray max;

    //Num cols,row,matrix's indices
    int cols, rows,y,x;

    //Initialize libpgm
    pgm_init(&argc, argv);
   
    //Read the image
    FILE *f=fopen("img/imm.pgm","r");
    image = pgm_readpgm(f, &cols, &rows, &max);

    for (y=0; y<rows; y++)
      {
        for (x=0; x<cols; x++)
	  {
            image[y][x] = image[y][x]/2;
	    printf("%d ",image[y][x]);
	  }
	printf("\n");
      }

    //Write the modified image to another file */
    FILE *fout=fopen("img/imm2.pgm","w");
    pgm_writepgm(fout, image, cols, rows, max, 1);


    /* cleanup */
    pgm_freearray(image, rows);
    return 0;
}

