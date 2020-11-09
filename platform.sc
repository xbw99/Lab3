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

#define ROWS 240
#define COLS 320
#define SIGMA 0.6
#define TLOW 0.3
#define THIGH 0.8
#define SIZE ROWS*COLS
#define WINSIZE 21


import "image_queue";

behavior DataIn(i_image_receiver ch1, i_image_sender ch2) {
	img IMAGE;	
	
	void main(void) {
		ch1.receive(&IMAGE);	
		while(true) {			
			ch2.send(IMAGE);
		}
	}
};

behavior DataOut(i_image_receiver ch1, i_image_sender ch2) {
	img IMAGE;
	
	void main(void) {
		while(true) {
			ch1.receive(&IMAGE);
			ch2.send(IMAGE);			
		}
	}
};

behavior DUT(i_image_receiver ch1, i_image_sender ch2) {
	img IMAGE;
	img EDGES;
	unsigned char NMS[SIZE];
short int SMOOTHEDIM[SIZE];
short int DELTA_X[SIZE];
short int DELTA_Y[SIZE];
short int MAGNITUDE[SIZE];

void canny(unsigned char *image, unsigned char **edge);
void gaussian_smooth(unsigned char *image, short int **smoothedim);
void make_gaussian_kernel(float **kernel);
void derrivative_x_y(short int *smoothedim, short int **delta_x, short int **delta_y);
void magnitude_x_y(short int *delta_x, short int *delta_y, short int **magnitude);
void apply_hysteresis(short int *mag, unsigned char *nms, unsigned char *edge);

int non_max_supp(short int *magnitude, short int *delta_x, short int *delta_y, unsigned char* nms);

void canny(unsigned char *image, unsigned char **edge)
{
	unsigned char *nms = 0;        /* Points that are local maximal magnitude. */
	short int *smoothedim = 0,     /* The image after gaussian smoothing.      */
		*delta_x = 0,        /* The first devivative image, x-direction. */
		*delta_y = 0,        /* The first derivative image, y-direction. */
		*magnitude = 0;      /* The magnitude of the gadient image.      */

	
	nms = NMS;
	smoothedim = SMOOTHEDIM;
	delta_x = DELTA_X;
	delta_y = DELTA_Y;
	magnitude = MAGNITUDE;

								 /****************************************************************************
								 * Perform gaussian smoothing on the image using the input standard
								 * deviation.
								 ****************************************************************************/

	if (VERBOSE) printf("Smoothing the image using a gaussian kernel.\n");
	gaussian_smooth(image, &smoothedim);

	/****************************************************************************
	* Compute the first derivative in the x and y directions.
	****************************************************************************/
	if (VERBOSE) printf("Computing the X and Y first derivatives.\n");
	derrivative_x_y(smoothedim, &delta_x, &delta_y);

	/****************************************************************************
	* Compute the magnitude of the gradient.
	****************************************************************************/
	if (VERBOSE) printf("Computing the magnitude of the gradient.\n");
	magnitude_x_y(delta_x, delta_y, &magnitude);

	/****************************************************************************
	* Perform non-maximal suppression.
	****************************************************************************/
	if (VERBOSE) printf("Doing the non-maximal suppression.\n");
	

	non_max_supp(magnitude, delta_x, delta_y, nms);

	/****************************************************************************
	* Use hysteresis to mark the edge pixels.
	****************************************************************************/
	if (VERBOSE) printf("Doing hysteresis thresholding.\n");
	

	apply_hysteresis(magnitude, nms, *edge);

	/****************************************************************************
	* Free all of the memory that we allocated except for the edge image that
	* is still being used to store out result.
	****************************************************************************/
	
	
	
}


/*******************************************************************************
* PROCEDURE: magnitude_x_y
* PURPOSE: Compute the magnitude of the gradient. This is the square root of
* the sum of the squared derivative values.
* NAME: Mike Heath
* DATE: 2/15/96
*******************************************************************************/
void magnitude_x_y(short int *delta_x, short int *delta_y, short int **magnitude)
{
	int r, c, pos, sq1, sq2, rows, cols;
	rows = ROWS, cols = COLS;

	/****************************************************************************
	* Allocate an image to store the magnitude of the gradient.
	****************************************************************************/
	

	for (r = 0, pos = 0; r<rows; r++) {
		for (c = 0; c<cols; c++, pos++) {
			sq1 = (int)delta_x[pos] * (int)delta_x[pos];
			sq2 = (int)delta_y[pos] * (int)delta_y[pos];
			(*magnitude)[pos] = (short)(0.5 + sqrt((float)sq1 + (float)sq2));
		}
	}

}

/*******************************************************************************
* PROCEDURE: derrivative_x_y
* PURPOSE: Compute the first derivative of the image in both the x any y
* directions. The differential filters that are used are:
*
*                                          -1
*         dx =  -1 0 +1     and       dy =  0
*                                          +1
*
* NAME: Mike Heath
* DATE: 2/15/96
*******************************************************************************/
void derrivative_x_y(short int *smoothedim, short int **delta_x, short int **delta_y)
{
	int r, c, pos, rows, cols;
	rows = ROWS, cols = COLS;

	/****************************************************************************
	* Allocate images to store the derivatives.
	****************************************************************************/
	

	/****************************************************************************
	* Compute the x-derivative. Adjust the derivative at the borders to avoid
	* losing pixels.
	****************************************************************************/
	if (VERBOSE) printf("   Computing the X-direction derivative.\n");
	for (r = 0; r<rows; r++) {
		pos = r * cols;
		(*delta_x)[pos] = smoothedim[pos + 1] - smoothedim[pos];
		pos++;
		for (c = 1; c<(cols - 1); c++, pos++) {
			(*delta_x)[pos] = smoothedim[pos + 1] - smoothedim[pos - 1];
		}
		(*delta_x)[pos] = smoothedim[pos] - smoothedim[pos - 1];
	}

	/****************************************************************************
	* Compute the y-derivative. Adjust the derivative at the borders to avoid
	* losing pixels.
	****************************************************************************/
	if (VERBOSE) printf("   Computing the Y-direction derivative.\n");
	for (c = 0; c<cols; c++) {
		pos = c;
		(*delta_y)[pos] = smoothedim[pos + cols] - smoothedim[pos];
		pos += cols;
		for (r = 1; r<(rows - 1); r++, pos += cols) {
			(*delta_y)[pos] = smoothedim[pos + cols] - smoothedim[pos - cols];
		}
		(*delta_y)[pos] = smoothedim[pos] - smoothedim[pos - cols];
	}
}

/*******************************************************************************
* PROCEDURE: gaussian_smooth
* PURPOSE: Blur an image with a gaussian filter.
* NAME: Mike Heath
* DATE: 2/15/96
*******************************************************************************/
void gaussian_smooth(unsigned char *image, short int **smoothedim)
{
	int r, c, rr, cc,rows, cols,     /* Counter variables. */
		windowsize,        /* Dimension of the gaussian kernel. */
		center;            /* Half of the windowsize. */
	float TEMPIM[SIZE], KERNEL[WINSIZE];
	float *tempim = 0,        /* Buffer for separable filter gaussian smoothing. */
		*kernel = 0,        /* A one dimensional gaussian kernel. */
		dot,            /* Dot product summing variable. */
		sum;            /* Sum of the kernel weights variable. */
	
	tempim = TEMPIM;
	kernel = KERNEL;
	rows = ROWS, cols = COLS;
	windowsize = WINSIZE;

						/****************************************************************************
						* Create a 1-dimensional gaussian smoothing kernel.
						****************************************************************************/

	if (VERBOSE) printf("   Computing the gaussian smoothing kernel.\n");
	make_gaussian_kernel(&kernel);
	center = windowsize / 2;

	/****************************************************************************
	* Allocate a temporary buffer image and the smoothed image.
	****************************************************************************/
	

	/****************************************************************************
	* Blur in the x - direction.
	****************************************************************************/
	if (VERBOSE) printf("   Bluring the image in the X-direction.\n");
	for (r = 0; r<rows; r++) {
		for (c = 0; c<cols; c++) {
			dot = 0.0;
			sum = 0.0;
			for (cc = (-center); cc <= center; cc++) {
				if (((c + cc) >= 0) && ((c + cc) < cols)) {
					dot += (float)image[r*cols + (c + cc)] * kernel[center + cc];
					sum += kernel[center + cc];
				}
			}
			tempim[r*cols + c] = dot / sum;
		}
	}

	/****************************************************************************
	* Blur in the y - direction.
	****************************************************************************/
	if (VERBOSE) printf("   Bluring the image in the Y-direction.\n");
	for (c = 0; c<cols; c++) {
		for (r = 0; r<rows; r++) {
			sum = 0.0;
			dot = 0.0;
			for (rr = (-center); rr <= center; rr++) {
				if (((r + rr) >= 0) && ((r + rr) < rows)) {
					dot += tempim[(r + rr)*cols + c] * kernel[center + rr];
					sum += kernel[center + rr];
				}
			}
			(*smoothedim)[r*cols + c] = (short int)(dot*BOOSTBLURFACTOR / sum + 0.5);
		}
	}

	
}

/*******************************************************************************
* PROCEDURE: make_gaussian_kernel
* PURPOSE: Create a one dimensional gaussian kernel.
* NAME: Mike Heath
* DATE: 2/15/96
*******************************************************************************/
void make_gaussian_kernel(float **kernel)
{
	int i, center, windowsize;
	float x, fx, sum = 0.0, sigma = SIGMA;

	windowsize = WINSIZE;
	center = windowsize / 2;

	if (VERBOSE) printf("      The kernel has %d elements.\n", windowsize);
	

	for (i = 0; i<windowsize; i++) {
		x = (float)(i - center);
		fx = pow(2.71828, -0.5*x*x / (sigma*sigma)) / (sigma * sqrt(6.2831853));
		(*kernel)[i] = fx;
		sum += fx;
	}

	for (i = 0; i<windowsize; i++) (*kernel)[i] /= sum;

	if (VERBOSE) {
		printf("The filter coefficients are:\n");
		for (i = 0; i<windowsize; i++)
			printf("kernel[%d] = %f\n", i, (*kernel)[i]);
	}
}

int follow_edges(unsigned char *edgemapptr, short *edgemagptr, short lowval,
	int cols)
{
	short *tempmagptr;
	unsigned char *tempmapptr;
	int i;
	int x[8] = { 1,1,0,-1,-1,-1,0,1 },
		y[8] = { 0,1,1,1,0,-1,-1,-1 };

	for (i = 0; i<8; i++) {
		tempmapptr = edgemapptr - y[i] * cols + x[i];
		tempmagptr = edgemagptr - y[i] * cols + x[i];

		if ((*tempmapptr == POSSIBLE_EDGE) && (*tempmagptr > lowval)) {
			*tempmapptr = (unsigned char)EDGE;
			follow_edges(tempmapptr, tempmagptr, lowval, cols);
		}
	}
	return 0;
}

/*******************************************************************************
* PROCEDURE: apply_hysteresis
* PURPOSE: This routine finds edges that are above some high threshhold or
* are connected to a high pixel by a path of pixels greater than a low
* threshold.
* NAME: Mike Heath
* DATE: 2/15/96
*******************************************************************************/
void apply_hysteresis(short int *mag, unsigned char *nms, unsigned char *edge)
{
	int r, c, pos, numedges, highcount, lowthreshold, highthreshold,
		hist[32768], rows, cols;
	short int maximum_mag;
	float tlow, thigh;
	rows = ROWS, cols = COLS, tlow = TLOW, thigh = THIGH;

	/****************************************************************************
	* Initialize the edge map to possible edges everywhere the non-maximal
	* suppression suggested there could be an edge except for the border. At
	* the border we say there can not be an edge because it makes the
	* follow_edges algorithm more efficient to not worry about tracking an
	* edge off the side of the image.
	****************************************************************************/
	for (r = 0, pos = 0; r<rows; r++) {
		for (c = 0; c<cols; c++, pos++) {
			if (nms[pos] == POSSIBLE_EDGE) edge[pos] = POSSIBLE_EDGE;
			else edge[pos] = NOEDGE;
		}
	}

	for (r = 0, pos = 0; r<rows; r++, pos += cols) {
		edge[pos] = NOEDGE;
		edge[pos + cols - 1] = NOEDGE;
	}
	pos = (rows - 1) * cols;
	for (c = 0; c<cols; c++, pos++) {
		edge[c] = NOEDGE;
		edge[pos] = NOEDGE;
	}

	/****************************************************************************
	* Compute the histogram of the magnitude image. Then use the histogram to
	* compute hysteresis thresholds.
	****************************************************************************/
	for (r = 0; r<32768; r++) hist[r] = 0;
	for (r = 0, pos = 0; r<rows; r++) {
		for (c = 0; c<cols; c++, pos++) {
			if (edge[pos] == POSSIBLE_EDGE) hist[mag[pos]]++;
		}
	}

	/****************************************************************************
	* Compute the number of pixels that passed the nonmaximal suppression.
	****************************************************************************/
	for (r = 1, numedges = 0; r<32768; r++) {
		if (hist[r] != 0) maximum_mag = r;
		numedges += hist[r];
	}

	highcount = (int)(numedges * thigh + 0.5);

	/****************************************************************************
	* Compute the high threshold value as the (100 * thigh) percentage point
	* in the magnitude of the gradient histogram of all the pixels that passes
	* non-maximal suppression. Then calculate the low threshold as a fraction
	* of the computed high threshold value. John Canny said in his paper
	* "A Computational Approach to Edge Detection" that "The ratio of the
	* high to low threshold in the implementation is in the range two or three
	* to one." That means that in terms of this implementation, we should
	* choose tlow ~= 0.5 or 0.33333.
	****************************************************************************/
	r = 1;
	numedges = hist[1];
	while ((r<(maximum_mag - 1)) && (numedges < highcount)) {
		r++;
		numedges += hist[r];
	}
	highthreshold = r;
	lowthreshold = (int)(highthreshold * tlow + 0.5);

	if (VERBOSE) {
		printf("The input low and high fractions of %f and %f computed to\n",
			tlow, thigh);
		printf("magnitude of the gradient threshold values of: %d %d\n",
			lowthreshold, highthreshold);
	}

	/****************************************************************************
	* This loop looks for pixels above the highthreshold to locate edges and
	* then calls follow_edges to continue the edge.
	****************************************************************************/
	for (r = 0, pos = 0; r<rows; r++) {
		for (c = 0; c<cols; c++, pos++) {
			if ((edge[pos] == POSSIBLE_EDGE) && (mag[pos] >= highthreshold)) {
				edge[pos] = EDGE;
				follow_edges((edge + pos), (mag + pos), lowthreshold, cols);
			}
		}
	}

	/****************************************************************************
	* Set all the remaining possible edges to non-edges.
	****************************************************************************/
	for (r = 0, pos = 0; r<rows; r++) {
		for (c = 0; c<cols; c++, pos++) if (edge[pos] != EDGE) edge[pos] = NOEDGE;
	}
}

/*******************************************************************************
* PROCEDURE: non_max_supp
* PURPOSE: This routine applies non-maximal suppression to the magnitude of
* the gradient image.
* NAME: Mike Heath
* DATE: 2/15/96
*******************************************************************************/

int non_max_supp(short int *mag, short int *gradx, short int *grady, unsigned char *result)
{
	int rowcount, colcount, count, nrows, ncols;
	short *magrowptr, *magptr;
	short *gxrowptr, *gxptr;
	short *gyrowptr, *gyptr, z1, z2;
	short m00, gx, gy;
	float mag1, mag2, xperp, yperp;
	unsigned char *resultrowptr, *resultptr;

	nrows = ROWS, ncols = COLS;
	/****************************************************************************
	* Zero the edges of the result image.
	****************************************************************************/
	for (count = 0, resultrowptr = result, resultptr = result + ncols*(nrows - 1);
		count<ncols; resultptr++, resultrowptr++, count++) {
		*resultrowptr = *resultptr = (unsigned char)0;
	}

	for (count = 0, resultptr = result, resultrowptr = result + ncols - 1;
		count<nrows; count++, resultptr += ncols, resultrowptr += ncols) {
		*resultptr = *resultrowptr = (unsigned char)0;
	}

	/****************************************************************************
	* Suppress non-maximum points.
	****************************************************************************/
	for (rowcount = 1, magrowptr = mag + ncols + 1, gxrowptr = gradx + ncols + 1,
		gyrowptr = grady + ncols + 1, resultrowptr = result + ncols + 1;
		rowcount<=nrows - 2;
		rowcount++, magrowptr += ncols, gyrowptr += ncols, gxrowptr += ncols,
		resultrowptr += ncols) {
		for (colcount = 1, magptr = magrowptr, gxptr = gxrowptr, gyptr = gyrowptr,
			resultptr = resultrowptr; colcount<=ncols - 2;
			colcount++, magptr++, gxptr++, gyptr++, resultptr++) {
			m00 = *magptr;
			if (m00 == 0) {
				*resultptr = (unsigned char)NOEDGE;
			}
			else {
				xperp = -(gx = *gxptr) / ((float)m00);
				yperp = (gy = *gyptr) / ((float)m00);
			}

			if (gx >= 0) {
				if (gy >= 0) {
					if (gx >= gy)
					{
						/* 111 */
						/* Left point */
						z1 = *(magptr - 1);
						z2 = *(magptr - ncols - 1);

						mag1 = (m00 - z1)*xperp + (z2 - z1)*yperp;

						/* Right point */
						z1 = *(magptr + 1);
						z2 = *(magptr + ncols + 1);

						mag2 = (m00 - z1)*xperp + (z2 - z1)*yperp;
					}
					else
					{
						/* 110 */
						/* Left point */
						z1 = *(magptr - ncols);
						z2 = *(magptr - ncols - 1);

						mag1 = (z1 - z2)*xperp + (z1 - m00)*yperp;

						/* Right point */
						z1 = *(magptr + ncols);
						z2 = *(magptr + ncols + 1);

						mag2 = (z1 - z2)*xperp + (z1 - m00)*yperp;
					}
				}
				else
				{
					if (gx >= -gy)
					{
						/* 101 */
						/* Left point */
						z1 = *(magptr - 1);
						z2 = *(magptr + ncols - 1);

						mag1 = (m00 - z1)*xperp + (z1 - z2)*yperp;

						/* Right point */
						z1 = *(magptr + 1);
						z2 = *(magptr - ncols + 1);

						mag2 = (m00 - z1)*xperp + (z1 - z2)*yperp;
					}
					else
					{
						/* 100 */
						/* Left point */
						z1 = *(magptr + ncols);
						z2 = *(magptr + ncols - 1);

						mag1 = (z1 - z2)*xperp + (m00 - z1)*yperp;

						/* Right point */
						z1 = *(magptr - ncols);
						z2 = *(magptr - ncols + 1);

						mag2 = (z1 - z2)*xperp + (m00 - z1)*yperp;
					}
				}
			}
			else
			{
				if ((gy = *gyptr) >= 0)
				{
					if (-gx >= gy)
					{
						/* 011 */
						/* Left point */
						z1 = *(magptr + 1);
						z2 = *(magptr - ncols + 1);

						mag1 = (z1 - m00)*xperp + (z2 - z1)*yperp;

						/* Right point */
						z1 = *(magptr - 1);
						z2 = *(magptr + ncols - 1);

						mag2 = (z1 - m00)*xperp + (z2 - z1)*yperp;
					}
					else
					{
						/* 010 */
						/* Left point */
						z1 = *(magptr - ncols);
						z2 = *(magptr - ncols + 1);

						mag1 = (z2 - z1)*xperp + (z1 - m00)*yperp;

						/* Right point */
						z1 = *(magptr + ncols);
						z2 = *(magptr + ncols - 1);

						mag2 = (z2 - z1)*xperp + (z1 - m00)*yperp;
					}
				}
				else
				{
					if (-gx > -gy)
					{
						/* 001 */
						/* Left point */
						z1 = *(magptr + 1);
						z2 = *(magptr + ncols + 1);

						mag1 = (z1 - m00)*xperp + (z1 - z2)*yperp;

						/* Right point */
						z1 = *(magptr - 1);
						z2 = *(magptr - ncols - 1);

						mag2 = (z1 - m00)*xperp + (z1 - z2)*yperp;
					}
					else
					{
						/* 000 */
						/* Left point */
						z1 = *(magptr + ncols);
						z2 = *(magptr + ncols + 1);

						mag1 = (z2 - z1)*xperp + (m00 - z1)*yperp;

						/* Right point */
						z1 = *(magptr - ncols);
						z2 = *(magptr - ncols - 1);

						mag2 = (z2 - z1)*xperp + (m00 - z1)*yperp;
					}
				}
			}

			/* Now determine if the current point is a maximum point */

			if ((mag1 > 0.0) || (mag2 > 0.0))
			{
				*resultptr = (unsigned char)NOEDGE;
			}
			else
			{
				if (mag2 == 0.0)
					*resultptr = (unsigned char)NOEDGE;
				else
					*resultptr = (unsigned char)POSSIBLE_EDGE;
			}
		}
	}
	return 0;
}
	
	unsigned char *image = 0;  
	unsigned char *edge = 0;
	void main(void) {
		while(true) {
			ch1.receive(&IMAGE);
			image = IMAGE;
			edge = EDGES;
			canny(image, &edge);
			ch2.send(EDGES);			
		}
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
