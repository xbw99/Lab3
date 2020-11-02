#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


#define VERBOSE 0

#define NOEDGE 255
#define POSSIBLE_EDGE 128
#define EDGE 0

#define BOOSTBLURFACTOR 90.0

import "image_queue"

behavior DataIn(c_image_receiver ch1, c_image_sender ch2) {
	img IMAGE;	
	
	void main(void) {
		while(true) {
			ch1.receive(&IMAGE);
			ch2.send(IMAGE);
		}
	}
};

behavior DataOut(c_image_receiver ch1, c_image_sender ch2) {
	int counter = 0;
	img IMAGE;
	
	void main(void) {
		int stop = 5;
		while(true) {
			ch1.receive(&IMAGE);
			ch2.send(IMAGE);
			if (counter++ == stop) {
				exit(1);
			}
		}
};

behavior DUT(c_image_receiver ch1, c_image_sender ch2) {
	img IMAGE;
	void main(void) {
		ch1.receive(&IMAGE);
		ch2.send(IMAGE);
	}
};

behavior Platform(c_image_receiver toPlatform, c_image_sender outPlatform) {
	c_image_queue inDut(2ul);
	c_image_queue outDut(2ul);
	DataIn din(toPlatform, inDut);
	DUT canny(inDut, outDut);
	DataOut dout(outDut, outPlatform);
	
	void main(void) {
		par( din; canny; dout; };
	}	
};
