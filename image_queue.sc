#include <stdio.h>
#include <stdlib.h>
#include <c_typed_queue.sh>

#define ROWS 240
#define COLS 320
#define SIZE ROWS * COLS

typedef unsigned char img[SIZE];
DEFINE_IC_TYPED_QUEUE(image, img)