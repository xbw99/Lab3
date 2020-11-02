#include <stdio.h>
#include <string.h>
import "image_queue";


behavior Monitor(i_image_receiver ch) {
	img IMAGE;
	int cols, rows;
	char *outfilename;
	char *comment;
	unsigned char* image;
	int maxval;
	FILE *fp;
	void main(void) {
		ch.receive(&IMAGE);
		
		outfilename = (char*)"ba.pgm";
		comment = (char*)"";
		maxval = 255;
		image = IMAGE;
		
		if (outfilename == NULL) fp = stdout;
		else {
			if ((fp = fopen(outfilename, "w")) == NULL) {
				fprintf(stderr, "Error writing the file %s in write_pgm_image().\n",
					outfilename);
				return;
			}
		}
		
		
	
		fprintf(fp, "P5\n%d %d\n", cols, rows);
		if (comment != NULL)
			if (strlen(comment) <= 70) fprintf(fp, "# %s\n", comment);
		fprintf(fp, "%d\n", maxval);

	
		if ((unsigned)rows != fwrite(image, cols, rows, fp)) {
			fprintf(stderr, "Error writing the image data in write_pgm_image().\n");
			if (fp != stdout) fclose(fp);
			return;
		}

		if (fp != stdout) fclose(fp);
		return;
		}
	
	
};
