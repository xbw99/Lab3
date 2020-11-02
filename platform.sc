import "image_queue"

behavior Platform(i_image_receiver ch) {
	img IMAGE;
	void main() {
		ch.receive(&IMAGE);
	}
};