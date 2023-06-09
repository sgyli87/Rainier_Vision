#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "image.h"

float get_pixel(image im, int x, int y, int c)
{
  assert( c<im.c && c>=0);
  x = x >= 0 ? x : 0;
  y = y >= 0 ? y : 0;
  x = x < im.w ? x : im.w-1;
  y = y < im.h ? y : im.h-1;
  int idx = c*im.w*im.h + im.w*y + x;
  return im.data[idx];
}

void set_pixel(image im, int x, int y, int c, float v)
{
  if( x<0 || x >= im.w || y < 0 || y >= im.h || c < 0 || c >= im.c){
    return;
  }
  int idx = c*im.w*im.h + im.w*y + x;
  im.data[idx] = v;
}

image copy_image(image im)
{
    image copy = make_image(im.w, im.h, im.c);
    memcpy(copy.data, im.data, im.w*im.h*im.c*sizeof(float));
    return copy;
}

image rgb_to_grayscale(image im)
{
    assert(im.c == 3);
    image gray = make_image(im.w, im.h, 1);

    for(int x = 0; x < im.w; x++){
      for(int y = 0; y < im.h; y++){
	float gray_v = 0.299*get_pixel(im,x,y,0)+0.587*get_pixel(im,x,y,1)+0.114*get_pixel(im,x,y,2);
	set_pixel(gray,x,y,0,gray_v);
      }
    }
    return gray;
}

void shift_image(image im, int c, float v)
{
    assert(c>=0 && c<im.c);
	
	for(int x = 0; x < im.w; x++){
		for(int y = 0; y < im.h; y++){
		
			float curr_v = get_pixel(im,x,y,c) + v;
			curr_v = curr_v > 0 ? curr_v:0;
			curr_v = curr_v <= 1.0 ? curr_v:1.0;
			set_pixel(im,x,y,c,curr_v);
		}
	}
}

//extra credit 
void scale_image(image im, int c, float v){
	assert(c>=0&&c<im.c);
	
	for(int x = 0; x < im.w; x++){
		for(int y = 0; y < im.h; y++){
			float curr_v = get_pixel(im,x,y,c)*(1+v);
			curr_v = curr_v > 0 ? curr_v:0;
			curr_v = curr_v <= 1.0 ? curr_v:1.0;
			set_pixel(im, x, y, c, curr_v);
		}
	}
}

void clamp_image(image im)
{
   for(int i = 0; i < im.w*im.h*im.c; i++){
		im.data[i] = im.data[i] >= 0 ? im.data[i]:0;
		im.data[i] = im.data[i] <= 1 ? im.data[i]:1;
   }
}


// These might be handy
float three_way_max(float a, float b, float c)
{
    return (a > b) ? ( (a > c) ? a : c) : ( (b > c) ? b : c) ;
}

float three_way_min(float a, float b, float c)
{
    return (a < b) ? ( (a < c) ? a : c) : ( (b < c) ? b : c) ;
}

void rgb_to_hsv(image im)
{
    assert(im.c == 3);
	
	float R, G, B, S, V, m, C, H;
	
	for(int x = 0; x < im.w; x++){
		for(int y = 0; y < im.h; y++){
			
			R = get_pixel(im,x,y,0);
			G = get_pixel(im,x,y,1);
			B = get_pixel(im,x,y,2);
			
			V = three_way_max(R,G,B);
			m = three_way_min(R,G,B);
			
			C = V - m;
			
			if(V==0) S = 0;
			else S = C/V;
			
			if(C==0) H = 0;
			else{
				if(V==R){
					H = (G-B)/C;
				}
				else if(V==G){
					H = (B-R)/C + 2;
				}
				else{
					H = (R-G)/C + 4;
				}
			}
			
			if(H<0){
				H = H/6+1;
			}
			else{
				H = H/6;
			}
			
			while(H<0){
				H += 1;
			}
			
			set_pixel(im,x,y,0,H);
			set_pixel(im,x,y,1,S);
			set_pixel(im,x,y,2,V);	
		}
	}
}

void hsv_to_rgb(image im)
{
  assert(im.c==3);
	float Hi, F, H, S, V, R, G, B;
	float P, Q, T;
	
	for(int x = 0; x < im.w; x++){
		for(int y = 0; y < im.h; y++){
		
		H = get_pixel(im, x, y, 0);
		S = get_pixel(im, x, y, 1);
		V = get_pixel(im, x, y, 2);
		
		H = H*6;
		Hi = floor(H);
		
		F = H - Hi;

		P = V * (1-S);
		Q = V * (1-F*S);
		T = V * (1-(1-F)*S);

		if(Hi == 0){
			R = V;
			G = T;
			B = P;
		}
		else if(Hi == 1){
			R = Q;
			G = V;
			B = P;
		}
		else if(Hi == 2){
			R = P;
			G = V;
			B = T;
		}
		else if(Hi == 3){
			R = P;
			G = Q;
			B = V;
		}
		else if(Hi == 4){
			R = T;
			G = P;
			B = V;
		}
		else if(Hi == 5){
			R = V;
			G = P;
			B = Q;
		}
		
		set_pixel(im,x,y,0,R);
		set_pixel(im,x,y,1,G);
		set_pixel(im,x,y,2,B);

		}
	}
}