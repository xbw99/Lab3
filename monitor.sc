#include <stdio.h>
#include <string.h>
import "image_queue";
#define _CRT_SECURE_NO_WARNINGS
#define ROWS 240
#define COLS 320

behavior Monitor(i_image_receiver ch) {
	img IMAGE;
	int cols, rows;
	char *outfilename;
	char outnamebuf[128];
	unsigned char* image;
	char *comment;
	int maxval;
	FILE *fp;
	void main(void) {
		int i = 0, stop = 1;
		while(i < stop) 
		{
			ch.receive(&IMAGE);			
			sprintf(outnamebuf, "ba%d.pgm", i);
		
			outfilename = outnamebuf;
			image = IMAGE;
			rows = ROWS;
			cols = COLS;
			comment = (char*)"";
			maxval = 255;
		
			if (outfilename == NULL) fp = stdout;
			else {
				if ((fp = fopen(outfilename, "w")) == NULL) {
					fprintf(stderr, "Error writing the file %s in write_pgm_image().\n",
						outfilename);
					exit(1);
				}
			}
			
			//fprintf(fp, "P5\n");
			
			fprintf(fp, "P5\n%d %d\n", cols, rows);
			if (comment != NULL)
				if (strlen(comment) <= 70) fprintf(fp, "# %s\n", comment);
			fprintf(fp, "%d\n", maxval);
			
			if ((unsigned)rows != fwrite(image, cols, rows, fp)) {
				fprintf(stderr, "Error writing the image data in write_pgm_image().\n");
				if (fp != stdout) fclose(fp);
				exit(1);
			}

			if (fp != stdout) fclose(fp);
			i++;
		}
			exit(0);
	}	
};
