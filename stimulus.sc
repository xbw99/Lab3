#include <stdio.h>
#include <string.h>
import "image_queue";


behavior Stimulus(i_image_sender ch, in char* infilename) {
	img IMAGE;	
	unsigned char *image = 0;
	int rows, cols;
	
	void main(void) {
		
		FILE *fp;
		char buf[71];

	
		if (infilename == NULL) fp = stdin;
		else {
			if ((fp = fopen(infilename, "r")) == NULL) {
				fprintf(stderr, "Error reading the file %s in read_pgm_image().\n",
					infilename);
				return;
			}
		}
		
		fgets(buf, 70, fp);
		if (strncmp(buf, "P5", 2) != 0) {
			fprintf(stderr, "The file %s is not in PGM format in ", infilename);
			fprintf(stderr, "read_pgm_image().\n");
			if (fp != stdin) fclose(fp);
			return;
		}
		
		do { fgets(buf, 70, fp); } while (buf[0] == '#');  
		sscanf(buf, "%d %d", &cols, &rows);
		do { fgets(buf, 70, fp); } while (buf[0] == '#');  
		
		
		image = IMAGE;
	
		if ((unsigned)rows != fread((*(&image)), cols, rows, fp)) {
			fprintf(stderr, "Error reading the image data in read_pgm_image().\n");
			if (fp != stdin) fclose(fp);
			return;
		}

		if (fp != stdin) fclose(fp);
		
		ch.send(IMAGE);
		return;
	}
	
};
