/* This is a collection of all function prototypes for OMEGA_pipeline.c */

 "flat_correction.c"		    		
/* Module 1 */
int Flat(float *p_inframe, float *p_flat, float *p_outframe, int total_pix, 
	 int image_no_out, char *name_inframe, char *flat_frame, int dark_flag, float *p_dark, 
	 char *dark_name, double dark_itime, int inframe_imno) 


"subframe_statistics.c"			
/* Module 2 */
int Statistics(float *p_inframe, int image_no)


 "normalize_frame.c"				
/* Module 3 */
int Normalize(float *p_inframe, int image_no, int tot_pix)


 "extract_name.c"				
/* Module 4 */
int Name(char *p_in_name, char *p_out_name)


 "load_stack.c"				
/* Module 5 */
int Load(float *p_input, int image_no, char frame_name[], int sky_mode)


 "get_sky.c"					
/* Module 6 */
int Sky()


 "save_sky.c"					
/* Module 7 */
int Savesky(int sky_imno, float *p_skyframe)


 "subtract_sky.c"				
/* Module 8 */
int Subtract(int sky_mode, int n_smallest, float kappa)


 "copy_frame.c"				
/* Module 9 */
int Copy(float *p_inframe, float *p_outframe, int total_pix, int image_no_out, char *name_inframe)


"original_level.c"
/* Module 10 */
int Original(float *p_inframe, int image_no, int tot_pix)


"minimum_sky.c"
/* Module 11 */
int Minimum(int n_smallest)


"sky_clipping.c"
/* Module 12 */
int Clipping(float kappa, int end_flag)


"frame_info.c"
/* Module 13 */
int Info(char *p_name, int image_no)


"sum_preparation.c"
/* Module 14 */
int Preparation(float *p_frame, int n_average)


"make_cleansum.c"
/* Module 15 */
int Cleansum(int n_average)


"dirty_sum.h"
/* Module 16 */
int Dirty(float *p_inframe)


"put_tosum.c"
/* Module 17 */
int Tosum(int sum_save_flag, int n_average, int n_real)


"save_dirty.c"
/* Module 18 */
int Savedirty()


"find_object.c"
/* Module 19 */
int Object(int n_average)


"cut_frame.c"
/* Module 20 */
int Cut(int flag, float *p_outframe, int *p_xstart, int *p_xend, int *p_ystart, int *p_yend)


"badpixel_correction.c"
/* Moduel 21 */
int Badpixel()


"elaborate_sum.c"
/* Module 22 */
int Elaborate(int n_average, float kappa_sum)
