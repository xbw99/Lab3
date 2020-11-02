import "image_queue";

behavior Platform(i_image_receiver ch) {
	img IMAGE;
	void main(void) {
		ch.receive(&IMAGE);
	}
};
