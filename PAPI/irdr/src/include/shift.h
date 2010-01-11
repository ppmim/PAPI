/* shift.h -- header file for shift.c */

extern float *shift_image_int(float *img, int nx, int ny, int border,
                              float xshift, float yshift);

extern float *shift_image(float *img, float *wimg, int nx, int ny, int border,
                          float xshift, float yshift, float **wimgout);

extern float *new_shift_image(float *img, float *wimg, int nx, int ny, 
		int xbelow, int xabove, int ybelow, int yabove,
                          float xshift, float yshift, float **wimgout);

extern float *unshift_image(float *img, float *wimg, int nx, int ny, int border,
                            float xshift, float yshift, float **wunshift);

extern int get_border(float *xshift, float *yshift, int n);

extern int new_get_border(float *xshift, float *yshift, int n,
		int *xbelow, int *xabove, int *ybelow, int *yabove);
