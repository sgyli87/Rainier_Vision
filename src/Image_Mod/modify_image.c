#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#define TWOPI 6.2831853

/******************************** Resizing *****************************
  To resize we'll need some interpolation methods and a function to create
  a new image and fill it in with our interpolation methods.
************************************************************************/

float nn_interpolate(image im, float x, float y, int c)
{
    // TODO
    /***********************************************************************
      This function performs nearest-neighbor interpolation on image "im"
      given a floating column value "x", row value "y" and integer channel "c",
      and returns the interpolated value.
    ************************************************************************/
    assert(c < im.c && c>=0);    
    x = roundf(x);
    y = roundf(y);

    return get_pixel(im, (int)x, (int)y, c);
}

image nn_resize(image im, int w, int h)
{
    // TODO Fill in (also fix the return line)
    /***********************************************************************
      This function uses nearest-neighbor interpolation on image "im" to a new
      image of size "w x h"
    ************************************************************************/
    image re = make_image(w, h, im.c);

    float aw = (float)im.w / (float)w;
    float bw = 0.5*(aw-1);
    float ah = (float)im.h / (float)h;
    float bh = 0.5*(ah-1);
    float x,y;

    for(int i = 0; i < w; i++){
      for(int j = 0; j < h; j++){
          x = i*aw + bw;
          y = j*ah + bh;
          for(int c = 0; c < im.c; c++){
            set_pixel(re, i, j, c, nn_interpolate(im,x,y,c));
          }
      }
    }
    
    return re;
}

float bilinear_interpolate(image im, float x, float y, int c)
{
    // TODO
    /***********************************************************************
      This function performs bilinear interpolation on image "im" given
      a floating column value "x", row value "y" and integer channel "c".
      It interpolates and returns the interpolated value.
    ************************************************************************/
    assert(c>=0 && c<im.c);
    int x1 = floor(x);
    int x2 = x1 + 1;
    int y1 = floor(y);
    int y2 = y1 + 1;
    float V1 = get_pixel(im,x1, y2, c);
    float V2 = get_pixel(im,x2, y2, c);
    float V3 = get_pixel(im,x1, y1, c);
    float V4 = get_pixel(im,x2, y1, c);
    float A1 = (x2-x) * (y-y1);
    float A2 = (x-x1) * (y-y1);
    float A3 = (x2-x) * (y2-y);
    float A4 = (x-x1) * (y2-y);
    return V1*A1 + V2*A2 + V3*A3 + V4*A4;
}

image bilinear_resize(image im, int w, int h)
{
    // TODO
    /***********************************************************************
      This function uses bilinear interpolation on image "im" to a new image
      of size "w x h". Algorithm is same as nearest-neighbor interpolation.
    ************************************************************************/
    image re = make_image(w,h,im.c);

    float aw = (float)im.w / (float)w;
    float bw = 0.5*(aw-1);
    float ah = (float)im.h / (float)h;
    float bh = 0.5*(ah-1);
    float x,y;

    for(int i = 0; i < w; i++){
      for(int j = 0; j < h; j++){
          x = i*aw + bw;
          y = j*ah + bh;
          for(int c = 0; c < im.c; c++){
            set_pixel(re, i, j, c, bilinear_interpolate(im,x,y,c));
          }
      }
    }
    
    return re;
}


/********************** Filtering: Box filter ***************************
  We want to create a box filter. We will only use square box filters.
************************************************************************/

void l1_normalize(image im)
{
    // TODO
    /***********************************************************************
      This function divides each value in image "im" by the sum of all the
      values in the image and modifies the image in place.
    ************************************************************************/
   float s = 0;
   for(int c = 0; c < im.c; c++){
     s = 0;
     for(int i = c*im.w*im.h; i < (c+1)*im.w*im.h; i++){
      s += im.data[i];
     }
     for(int i = c*im.w*im.h; i < (c+1)*im.w*im.h; i++){
      if(s>0) im.data[i]/=s;
      else im.data[i] = 1.0/im.w/im.h;
     }
   }
}

image make_box_filter(int w)
{
    // TODO
    /***********************************************************************
      This function makes a square filter of size "w x w". Make an image of
      width = height = w and number of channels = 1, with all entries equal
      to 1. Then use "l1_normalize" to normalize your filter.
    ************************************************************************/
    assert(w%2);
    image filter = make_image(w,w,1);
    for(int i = 0; i < w*w; i++){
      filter.data[i] = 1;
    }
    l1_normalize(filter);
    return filter;
}

image convolve_image(image im, image filter, int preserve)
{
    // TODO
    /***********************************************************************
      This function convolves the image "im" with the "filter". The value
      of preserve is 1 if the number of input image channels need to be 
      preserved. Check the detailed algorithm given in the README.  
    ************************************************************************/
    assert(filter.c==1);
    image re;

    float val;
    if(preserve){
      re = make_image(im.w, im.h, im.c);
      for(int x = 0; x < im.w; x++){
        for(int y = 0; y < im.h; y++){
          for(int c = 0; c < im.c; c++){
            val = 0;
            for(int w = 0; w < filter.w; w++){
              for(int h = 0; h < filter.h; h++){
                val += get_pixel(filter, w, h, 0) * get_pixel(im,x-filter.w/2+w, y-filter.h/2+h, c);
              }
            }
            set_pixel(re, x, y, c, val);
          }
        }
      }
    }
    else{
      re = make_image(im.w, im.h, 1);
      for(int x = 0; x < im.w; x++){
        for(int y = 0; y < im.h; y++){
          val = 0;
          for(int c = 0; c < im.c; c++){
            for(int w = 0; w < filter.w; w++){
              for(int h = 0; h < filter.h; h++){
                val += get_pixel(filter, w, h, 0) * get_pixel(im,x-filter.w/2+w, y-filter.h/2+h, c);
              }
            }
          }
          set_pixel(re, x, y, 0, val);
        }
      }
    }
    
    return re;
}

image make_highpass_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 filter with highpass filter values using image.data[]
    ************************************************************************/
    image filter = make_image(3,3,1);
    set_pixel(filter,0,0,0,0);
    set_pixel(filter,0,1,0,-1);
    set_pixel(filter,0,2,0,0);
    set_pixel(filter,1,0,0,-1);
    set_pixel(filter,1,1,0,4);
    set_pixel(filter,1,2,0,-1);
    set_pixel(filter,2,0,0,0);
    set_pixel(filter,2,1,0,-1);
    set_pixel(filter,2,2,0,0);

    return filter;
}

image make_sharpen_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 filter with sharpen filter values using image.data[]
    ************************************************************************/
    image filter = make_highpass_filter();
    set_pixel(filter, 1, 1, 0, 5);
    return filter;
}

image make_emboss_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 filter with emboss filter values using image.data[]
    ************************************************************************/
    image filter = make_image(3,3,1);
    set_pixel(filter,0,0,0,-2);
    set_pixel(filter,0,1,0,-1);
    set_pixel(filter,0,2,0,0);
    set_pixel(filter,1,0,0,-1);
    set_pixel(filter,1,1,0,1);
    set_pixel(filter,1,2,0,1);
    set_pixel(filter,2,0,0,0);
    set_pixel(filter,2,1,0,1);
    set_pixel(filter,2,2,0,2);

    return filter;
}

// Question 2.3.1: Which of these filters should we use preserve when we run our
// convolution and which ones should we not? Why?
// Answer: we should use preserve when we are using sharpen and emboss filters because they are applied to
// all three bands, while for highpass filter we should not because it applies on graytone images.

// Question 2.3.2: Do we have to do any post-processing for the above filters? 
//Which ones and why?
// Answer: For all filters above we have to do some post-processing such as smoothing because 
// images are noisy and they amplifies noise.

image make_gaussian_filter(float sigma)
{
    // TODO
    /***********************************************************************
      sigma: a float number for the Gaussian.
      Create a Gaussian filter with the given sigma. Note that the kernel size 
      is the next highest odd integer from 6 x sigma. Return the Gaussian filter.
    ************************************************************************/
    int w = (int)roundf(sigma*6) + 1;
    if(w%2==0) w++;
    float sigma_sqr = pow(sigma, 2);
    int x,y;
    image filter = make_image(w, w, 1);
    for(int i = 0; i < w; i++){
      for(int j = 0; j < w; j++){
        x = pow(i-w/2, 2);
        y = pow(j-w/2, 2);
        float v = 1.0/(TWOPI*sigma_sqr)*exp(-(x+y)/(2*sigma_sqr));
        set_pixel(filter, i, j, 0, v);
      }
    }
    l1_normalize(filter);
    return filter;
}

image add_image(image a, image b)
{
    // TODO
    /***********************************************************************
      The input images a and image b have the same height, width, and channels.
      Sum the given two images and return the result, which should also have
      the same height, width, and channels as the inputs. Do necessary checks.
    ************************************************************************/
    assert(a.w==b.w && a.h==b.h && a.c==b.c);

    image re = make_image(a.w, a.h, a.c);

    for(int i = 0; i < a.w*a.h*a.c; i++){
      re.data[i] = a.data[i] + b.data[i];
    }
    
    return re;
}

image sub_image(image a, image b)
{
    // TODO
    /***********************************************************************
      The input image a and image b have the same height, width, and channels.
      Subtract the given two images and return the result, which should have
      the same height, width, and channels as the inputs. Do necessary checks.
    ************************************************************************/
    assert(a.w==b.w && a.h==b.h && a.c==b.c);

    image re = make_image(a.w, a.h, a.c);

    for(int i = 0; i < a.w*a.h*a.c; i++){
      re.data[i] = a.data[i] - b.data[i];
    }
    
    return re;
}

image make_gx_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 Sobel Gx filter and return it
    ************************************************************************/
    image filter = make_image(3,3,1);
    set_pixel(filter,0,0,0,-1);
    set_pixel(filter,0,1,0,-2);
    set_pixel(filter,0,2,0,-1);
    set_pixel(filter,1,0,0,0);
    set_pixel(filter,1,1,0,0);
    set_pixel(filter,1,2,0,0);
    set_pixel(filter,2,0,0,1);
    set_pixel(filter,2,1,0,2);
    set_pixel(filter,2,2,0,1);
    
    return filter;
}

image make_gy_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 Sobel Gy filter and return it
    ************************************************************************/
    image filter = make_image(3,3,1);
    set_pixel(filter,0,0,0,-1);
    set_pixel(filter,0,1,0,0);
    set_pixel(filter,0,2,0,1);
    set_pixel(filter,1,0,0,-2);
    set_pixel(filter,1,1,0,0);
    set_pixel(filter,1,2,0,2);
    set_pixel(filter,2,0,0,-1);
    set_pixel(filter,2,1,0,0);
    set_pixel(filter,2,2,0,1);
    
    return filter;
}

void feature_normalize(image im)
{
    // TODO
    /***********************************************************************
      Calculate minimum and maximum pixel values. Normalize the image by
      subtracting the minimum and dividing by the max-min difference.
    ************************************************************************/
  assert(im.w>0&&im.h>0);
  float val_max = im.data[0], val_min=im.data[0];
  for(int i = 0; i < im.w*im.h*im.c; i++){
    if(im.data[i] > val_max) val_max = im.data[i];
    if(im.data[i] < val_min) val_min = im.data[i];
  }
  float diff = val_max - val_min;
  int zero = diff == 0;
  for(int i = 0; i < im.w*im.h*im.c; i++){
    if(zero){
      im.data[i] = 0;
    }
    else{
      im.data[i] = (im.data[i]-val_min)/diff;
    }
  }
}

image *sobel_image(image im)
{
    // TODO
    /***********************************************************************
      Apply Sobel filter to the input image "im", get the magnitude as sobelimg[0]
      and gradient as sobelimg[1], and return the result.
    ************************************************************************/
    image *sobelimg = calloc(2, sizeof(image));
    
    image gx = make_gx_filter();
    image gy = make_gy_filter();

    image re_x = convolve_image(im, gx, 0);
    image re_y = convolve_image(im, gy, 0);

    image mag = make_image(im.w, im.h, 1);
    image theta = make_image(im.w, im.h, 1);
    float x, y;

    for(int i = 0; i < mag.h*mag.w*mag.c; i++){
      x = re_x.data[i];
      y = re_y.data[i];
      mag.data[i] = sqrt(pow(x,2)+pow(y,2));
      theta.data[i] = atan2f(y,x);
    }
    
    sobelimg[0] = mag;
    sobelimg[1] = theta;

    return sobelimg;
}

image colorize_sobel(image im)
{
  // TODO
  /***********************************************************************
    Create a colorized version of the edges in image "im" using the 
    algorithm described in the README.
  ************************************************************************/
  image filter = make_gaussian_filter(4);
  image im_c = convolve_image(im, filter, 1);

  image *sobel = sobel_image(im_c);
  image mag = sobel[0];
  image theta = sobel[1];

  image re = make_image(im.w, im.h, 3);
  feature_normalize(mag);
  
  for(int i = 0; i < mag.w; i++){
    for(int j = 0; j < mag.h; j++){
      set_pixel(re, i, j, 0, get_pixel(theta, i, j, 0)/TWOPI + 0.5);
      set_pixel(re, i, j, 1, get_pixel(mag, i, j, 0));
      set_pixel(re, i, j, 2, get_pixel(mag, i, j, 0));
    }
  }
  hsv_to_rgb(re);
  return re;
}

// EXTRA CREDIT: Median filter

/*
image apply_median_filter(image im, int kernel_size)
{
  return make_image(1,1,1);
}
*/

// SUPER EXTRA CREDIT: Bilateral filter

/*
image apply_bilateral_filter(image im, float sigma1, float sigma2)
{
  return make_image(1,1,1);
}
*/