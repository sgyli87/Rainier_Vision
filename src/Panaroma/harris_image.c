#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#include "matrix.h"
#include <time.h>
#define TWOPI 6.2831853

// Frees an array of descriptors.
// descriptor *d: the array.
// int n: number of elements in array.
void free_descriptors(descriptor *d, int n)
{
    int i;
    for(i = 0; i < n; ++i){
        free(d[i].data);
    }
    free(d);
}

// Create a feature descriptor for an index in an image.
// image im: source image.
// int i: index in image for the pixel we want to describe.
// returns: descriptor for that index.
descriptor describe_index(image im, int i)
{
    int w = 5;
    descriptor d;
    d.p.x = i%im.w;
    d.p.y = i/im.w;
    d.data = calloc(w*w*im.c, sizeof(float));
    d.n = w*w*im.c;
    int c, dx, dy;
    int count = 0;
    // If you want you can experiment with other descriptors
    // This subtracts the central value from neighbors
    // to compensate some for exposure/lighting changes.
    for(c = 0; c < im.c; ++c){
        float cval = im.data[c*im.w*im.h + i];
        for(dx = -w/2; dx < (w+1)/2; ++dx){
            for(dy = -w/2; dy < (w+1)/2; ++dy){
                float val = get_pixel(im, i%im.w+dx, i/im.w+dy, c);
                d.data[count++] = cval - val;
            }
        }
    }
    return d;
}

// Marks the spot of a point in an image.
// image im: image to mark.
// ponit p: spot to mark in the image.
void mark_spot(image im, point p)
{
    int x = p.x;
    int y = p.y;
    int i;
    for(i = -9; i < 10; ++i){
        set_pixel(im, x+i, y, 0, 1);
        set_pixel(im, x, y+i, 0, 1);
        set_pixel(im, x+i, y, 1, 0);
        set_pixel(im, x, y+i, 1, 0);
        set_pixel(im, x+i, y, 2, 1);
        set_pixel(im, x, y+i, 2, 1);
    }
}

// Marks corners denoted by an array of descriptors.
// image im: image to mark.
// descriptor *d: corners in the image.
// int n: number of descriptors to mark.
void mark_corners(image im, descriptor *d, int n)
{
    int i;
    for(i = 0; i < n; ++i){
        mark_spot(im, d[i].p);
    }
}

// Creates a 1d Gaussian filter.
// float sigma: standard deviation of Gaussian.
// returns: single row image of the filter.
image make_1d_gaussian(float sigma)
{
    int w = (int)roundf(sigma*6) + 1;
    if(w%2==0) w++;
    float sigma_sqr = pow(sigma, 2);
    int x,y;
    image filter = make_image(w, 1, 1);
    for(int i = 0; i < w; i++){
      for(int j = 0; j < 1; j++){
        x = pow(i-w/2, 2);
        y = pow(j-w/2, 2);
        float v = 1.0/(TWOPI*sigma_sqr)*exp(-(x+y)/(2*sigma_sqr));
        set_pixel(filter, i, j, 0, v);
      }
    }
    l1_normalize(filter);
    return filter;
}

// Smooths an image using separable Gaussian filter.
// image im: image to smooth.
// float sigma: std dev. for Gaussian.
// returns: smoothed image.
image smooth_image(image im, float sigma)
{
    image filter = make_1d_gaussian(sigma);
    image im_1 = convolve_image(im, filter, 1);
    
    image transpose_filter = make_image(1, filter.w, 1);
    for(int i = 0; i < filter.w; i++){
        set_pixel(transpose_filter, 0, i, 0, get_pixel(filter, i, 0, 0));
    }
    image im_2 = convolve_image(im_1, transpose_filter, 1);
    return im_2;
}

// Calculate the structure matrix of an image.
// image im: the input image.
// float sigma: std dev. to use for weighted sum.
// returns: structure matrix. 1st channel is Ix^2, 2nd channel is Iy^2,
//          third channel is IxIy.
image structure_matrix(image im, float sigma)
{
    assert(im.c == 1 || im.c == 3);
    image im_to_process;
   
    im_to_process = copy_image(im);
    //else im_to_process = rgb_to_grayscale(im);

    image S = make_image(im_to_process.w, im_to_process.h, 3);

    image Ix_im = convolve_image(im_to_process, make_gx_filter(),0);
    image Iy_im = convolve_image(im_to_process, make_gy_filter(),0);

    matrix Ix = make_matrix(im_to_process.w, im_to_process.h);
    matrix Iy = make_matrix(im_to_process.w, im_to_process.h);

    for(int i = 0; i < im_to_process.w; i++){
        for(int j = 0; j < im_to_process.h; j++){
            Ix.data[i][j] = get_pixel(Ix_im, i, j, 0);
            Iy.data[i][j] = get_pixel(Iy_im, i, j, 0);
        }
    }
    
    matrix IxIx = matrix_elmult_matrix(Ix,Ix);
    matrix IyIy = matrix_elmult_matrix(Iy,Iy);
    matrix IxIy = matrix_elmult_matrix(Ix,Iy);

    for(int i = 0; i < im_to_process.w; i++){
        for(int j = 0; j < im_to_process.h; j++){
            set_pixel(S, i, j, 0, IxIx.data[i][j]);
            set_pixel(S, i, j, 1, IyIy.data[i][j]);
            set_pixel(S, i, j, 2, IxIy.data[i][j]);
        }
    }   
    return smooth_image(S,sigma);
}

// Estimate the cornerness of each pixel given a structure matrix S.
// image S: structure matrix for an image.
// returns: a response map of cornerness calculations.
image cornerness_response(image S)
{
    image R = make_image(S.w, S.h, 1);
    for(int w = 0; w < S.w; w++){
        for(int h = 0; h < S.h; h++){
            float IxIx = get_pixel(S,w,h,0);
            float IyIy = get_pixel(S,w,h,1);
            float IxIy = get_pixel(S,w,h,2);

            set_pixel(R, w, h, 0, (IxIx*IyIy) - powf(IxIy,2) - 0.06 * powf(IxIx+IyIy,2));
        }
    }   
    return R;
}

// Perform non-max supression on an image of feature responses.
// image im: 1-channel image of feature responses.
// int w: distance to look for larger responses.
// returns: image with only local-maxima responses within w pixels.
image nms_image(image im, int w)
{
    image r = copy_image(im);
    for(int x = 0; x < im.w; x++){
        for(int y = 0; y < im.h; y++){
            for(int x_adj = x - w; x_adj <= x + w; x_adj++){
                for(int y_adj = y - w; y_adj <= y + w; y_adj++){
                    if(get_pixel(im,x,y,0)<get_pixel(im,x_adj,y_adj,0)){
                        set_pixel(r, x, y, 0, -INFINITY);
                    }
                }
            }
        }
    }
    return r;
}

// Perform harris corner detection and extract features from the corners.
// image im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
// int *n: pointer to number of corners detected, should fill in.
// returns: array of descriptors of the corners in the image.
descriptor *harris_corner_detector(image im, float sigma, float thresh, int nms, int *n)
{
    // Calculate structure matrix
    image S = structure_matrix(im, sigma);
    // Estimate cornerness
    image R = cornerness_response(S);
    // Run NMS on the responses
    image Rnms = nms_image(R, nms);

    //TODO: count number of responses over threshold
    int count = 0; // change this
    
    for(int x = 0; x < im.w; x++){
        for(int y = 0; y < im.h; y++){
            if(get_pixel(Rnms,x,y,0)>thresh){
                count++;
            }
        }
    }

    *n = count; // <- set *n equal to number of corners in image.

    descriptor *d = calloc(count, sizeof(descriptor));
    //TODO: fill in array *d with descriptors of corners, use describe_index.
    int i = 0;
    for(int x = 0; x < im.w; x++){
        for(int y = 0; y < im.h; y++){
            if(get_pixel(Rnms,x,y,0)>thresh){
                d[i] = describe_index(im, im.w*y + x );
                i++;
            }
        }
    }
    free_image(S);
    free_image(R);
    free_image(Rnms);
    return d;
}

// Find and draw corners on an image.
// image im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
void detect_and_draw_corners(image im, float sigma, float thresh, int nms)
{
    int n = 0;
    descriptor *d = harris_corner_detector(im, sigma, thresh, nms, &n);
    mark_corners(im, d, n);
}
