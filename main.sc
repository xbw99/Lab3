import "image_queue"
import "stimulus"
import "platform"

behavior Main() {
	c_image_queue toPlatform(2ul)
	char *infilename;
	Stimulus stim(toPlatform, infilename);
	Platform plat(toPlatform);

	int main(int argc, char *argv[]) {
		infilename = argv[1];
		
		par { stim; plat; };
		return 0;
	}
};
