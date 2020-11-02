#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define _USE_MATH_DEFINES
#define VERBOSE 0
#define NOEDGE 255
#define POSSIBLE_EDGE 128
#define EDGE 0
#define BOOSTBLURFACTOR 90.0

import "image_queue";

behavior DataIn(i_image_receiver ch1, i_image_sender ch2) {
	img IMAGE;	
	
	void main(void) {
		while(true) {
			ch1.receive(&IMAGE);
			ch2.send(IMAGE);
		}
	}
};

behavior DataOut(i_image_receiver ch1, i_image_sender ch2) {
	int counter = 0;
	img IMAGE;
	
	void main(void) {
		int stop = 10;
		while(true) {
			ch1.receive(&IMAGE);
			ch2.send(IMAGE);
			if (counter++ == stop) {
				exit(1);
			}
		}
	}
};

behavior DUT(i_image_receiver ch1, i_image_sender ch2) {
	img IMAGE;
	void main(void) {
		ch1.receive(&IMAGE);
		ch2.send(IMAGE);
	}
};

behavior Platform(i_image_receiver toPlatform, i_image_sender outPlatform) {
	c_image_queue inDut(2ul);
	c_image_queue outDut(2ul);
	DataIn din(toPlatform, inDut);
	DUT canny(inDut, outDut);
	DataOut dout(outDut, outPlatform);
	
	void main(void) {
		par{ din; canny; dout; };
	}	
};
