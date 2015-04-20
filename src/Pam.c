#include <pam.h>
#include <stdio.h>

//gcc -I /usr/local/netpbm/include/ -L /usr/local/netpbm/lib -o Pam Pam.c -lnetpbm -lm
//gcc -o Pam Pam.c -I /usr/local/include/netpbm/ -lnetpbm -lm

int main(int argc, char *argv[]) {

	// The two pam, one for input and one for output,structures
	// containing infos about the image
	struct pam inpam, outpam;

	// The 3D array containing the values
	tuple ** imageArray;

	// Indexes to browse the array
	unsigned int row, column, plane;

	// Initialize libpgm
	pgm_init(&argc, argv);

	// Read the input image and pass
	inpam.file = fopen("../img/imm.pgm", "r");
	imageArray = pnm_readpam(inpam.file, &inpam, PAM_STRUCT_SIZE(tuple_type));

	// Prepare the struct for the output image
	outpam = inpam;
	outpam.file = fopen("../img/imm2.pgm", "w");
	outpam.plainformat = 1; // Force plain, so that we, useless meatbags, can read what is being written

	// Loops!
	for (row = 0; row < inpam.height; row++) {

		for (column = 0; column < inpam.width; ++column) {
			for (plane = 0; plane < inpam.depth; ++plane) {
				// We don't actually use the 3rd dimension... I miss the 90s!
				imageArray[row][column][plane] = imageArray[row][column][plane]/2;
			}
		}
	}

	// Write all the resulting values into the file
	pnm_writepam(&outpam, imageArray);

	// Free space
	pnm_freepamarray(imageArray, &inpam);

	return 0;
}

