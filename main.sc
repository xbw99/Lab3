import "image_queue";
import "stimulus";
import "platform";
import "monitor";

behavior Main() {
	c_image_queue toPlatform(2ul);
	c_image_queue outPlatform(2ul);
	char *infilename;
	Stimulus stim(toPlatform, infilename);
	Platform plat(toPlatform, outPlatform);
	Monitor mon(outPlatform);

	int main(int argc, char *argv[]) {
		infilename = argv[1];
		
		par { stim; plat; mon; };
		return 0;
	}
};
