#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#include "matrix.h"

// Comparator for matches
// const void *a, *b: pointers to the matches to compare.
// returns: result of comparison, 0 if same, 1 if a > b, -1 if a < b.
int match_compare(const void *a, const void *b)
{
    match *ra = (match *)a;
    match *rb = (match *)b;
    if (ra->distance < rb->distance) return -1;
    else if (ra->distance > rb->distance) return  1;
    else return 0;
}

// Helper function to create 2d points.
// float x, y: coordinates of point.
// returns: the point.
point make_point(float x, float y)
{
    point p;
    p.x = x; p.y = y;
    return p;
}

// Place two images side by side on canvas, for drawing matching pixels.
// image a, b: images to place.
// returns: image with both a and b side-by-side.
image both_images(image a, image b)
{
    image both = make_image(a.w + b.w, a.h > b.h ? a.h : b.h, a.c > b.c ? a.c : b.c);
    int i,j,k;
    for(k = 0; k < a.c; ++k){
        for(j = 0; j < a.h; ++j){
            for(i = 0; i < a.w; ++i){
                set_pixel(both, i, j, k, get_pixel(a, i, j, k));
            }
        }
    }
    for(k = 0; k < b.c; ++k){
        for(j = 0; j < b.h; ++j){
            for(i = 0; i < b.w; ++i){
                set_pixel(both, i+a.w, j, k, get_pixel(b, i, j, k));
            }
        }
    }
    return both;
}

// Draws lines between matching pixels in two images.
// image a, b: two images that have matches.
// match *matches: array of matches between a and b.
// int n: number of matches.
// int inliers: number of inliers at beginning of matches, drawn in green.
// returns: image with matches drawn between a and b on same canvas.
image draw_matches(image a, image b, match *matches, int n, int inliers)
{
    image both = both_images(a, b);
    int i,j;
    for(i = 0; i < n; ++i){
        int bx = matches[i].p.x; 
        int ex = matches[i].q.x; 
        int by = matches[i].p.y;
        int ey = matches[i].q.y;
        for(j = bx; j < ex + a.w; ++j){
            int r = (float)(j-bx)/(ex+a.w - bx)*(ey - by) + by;
            set_pixel(both, j, r, 0, i<inliers?0:1);
            set_pixel(both, j, r, 1, i<inliers?1:0);
            set_pixel(both, j, r, 2, 0);
        }
    }
    return both;
}

// Draw the matches with inliers in green between two images.
// image a, b: two images to match.
// matches *
image draw_inliers(image a, image b, matrix H, match *m, int n, float thresh)
{
    int inliers = model_inliers(H, m, n, thresh);
    image lines = draw_matches(a, b, m, n, inliers);
    return lines;
}

// Find corners, match them, and draw them between two images.
// image a, b: images to match.
// float sigma: gaussian for harris corner detector. Typical: 2
// float thresh: threshold for corner/no corner. Typical: 1-5
// int nms: window to perform nms on. Typical: 3
image find_and_draw_matches(image a, image b, float sigma, float thresh, int nms)
{
    int an = 0;
    int bn = 0;
    int mn = 0;
    descriptor *ad = harris_corner_detector(a, sigma, thresh, nms, &an);
    descriptor *bd = harris_corner_detector(b, sigma, thresh, nms, &bn);
    match *m = match_descriptors(ad, an, bd, bn, &mn);

    mark_corners(a, ad, an);
    mark_corners(b, bd, bn);
    image lines = draw_matches(a, b, m, mn, 0);

    free_descriptors(ad, an);
    free_descriptors(bd, bn);
    free(m);
    return lines;
}

// Calculates L1 distance between to floating point arrays.
// float *a, *b: arrays to compare.
// int n: number of values in each array.
// returns: l1 distance between arrays (sum of absolute differences).
float l1_distance(float *a, float *b, int n)
{
    float diff = 0.0;
    for(int i = 0; i < n; i++){
        diff += fabs(a[i]-b[i]);
    }
    return diff;
}

// Finds best matches between descriptors of two images.
// descriptor *a, *b: array of descriptors for pixels in two images.
// int an, bn: number of descriptors in arrays a and b.
// int *mn: pointer to number of matches found, to be filled in by function.
// returns: best matches found. each descriptor in a should match with at most
//          one other descriptor in b.
match *match_descriptors(descriptor *a, int an, descriptor *b, int bn, int *mn)
{
    int i,j,k;
    // We will have at most an matches.
    *mn = an;
    match *m = calloc(an, sizeof(match));
    
    for(j = 0; j < an; ++j){
        // TODO: for every descriptor in a, find best match in b.
        // record ai as the index in *a and bi as the index in *b.
        int bind = 0; // <- find the best match
        float best_distance = 1e10f;
        
        for(i = 0; i < bn; i++){
            float dist = l1_distance(a[j].data, b[i].data, b[i].n);
            //printf("%f",dist);
            //printf("\n");
            if(dist < best_distance){
                best_distance = dist;
                bind = i;
            }
        }
        m[j].ai = j;
        m[j].bi = bind; // <- should be index in b.
        m[j].p = a[j].p;
        m[j].q = b[bind].p;
        m[j].distance = best_distance; // <- should be the smallest L1 distance!
    }

    int count = 0;
    int *seen = calloc(bn, sizeof(int));
    // TODO: we want matches to be injective (one-to-one).
    // Sort matches based on distance using match_compare and qsort.
    // Then throw out matches to the same element in b. Use seen to keep track.
    // Each point should only be a part of one match.
    // Some points will not be in a match.
    // In practice just bring good matches to front of list, set *mn.
    qsort(m, an, sizeof(match), match_compare);
    
    for(i=0;i<an;i++){
        int b_no = m[i].bi;
        if(seen[b_no]!=1){
            seen[b_no]=1;
            count++;
        }
        else{
            for(j = i; j < an-1; j++){
                m[j] = m[j+1];
            }
            an--;
            i--;
        }
    }

    *mn = count;
    free(seen);
    return m;
}

// Apply a projective transformation to a point.
// matrix H: homography to project point.
// point p: point to project.
// returns: point projected using the homography.
point project_point(matrix H, point p)
{
    matrix c = make_matrix(3, 1);
    // TODO: project point p with homography H.
    // Remember that homogeneous coordinates are equivalent up to scalar.
    // Have to divide by.... something...
    c.data[0][0] = p.x;
    c.data[1][0] = p.y;
    c.data[2][0] = 1;

    matrix c_prime = make_matrix(3,1);
    c_prime = matrix_mult_matrix(H, c);

    float xprime = c_prime.data[0][0]/c_prime.data[2][0];
    float yprime = c_prime.data[1][0]/c_prime.data[2][0];
    
    point q = make_point(xprime, yprime);
    return q;
}

// Calculate L2 distance between two points.
// point p, q: points.
// returns: L2 distance between them.
float point_distance(point p, point q)
{
    float x2 = pow(p.x-q.x,2);
    float y2 = pow(p.y-q.y,2);
    return sqrt(x2+y2);
}

// Count number of inliers in a set of matches. Should also bring inliers
// to the front of the array.
// matrix H: homography between coordinate systems.
// match *m: matches to compute inlier/outlier.
// int n: number of matches in m.
// float thresh: threshold to be an inlier.
// returns: number of inliers whose projected point falls within thresh of
//          their match in the other image. Should also rearrange matches
//          so that the inliers are first in the array. For drawing.
int model_inliers(matrix H, match *m, int n, float thresh)
{
    int i;
    int count = 0;
    // TODO: count number of matches that are inliers
    // i.e. distance(H*p, q) < thresh
    // Also, sort the matches m so the inliers are the first 'count' elements.
    int j = 0;
    for(i = 0; i < n; i++){
        if(point_distance(project_point(H,m[i].p),m[i].q)<thresh){
            count++; 
            match temp = m[j];
            m[j] = m[i];
            m[i] = temp;
            j++;
        }
    }
    return count;
}

// Randomly shuffle matches for RANSAC.
// match *m: matches to shuffle in place.
// int n: number of elements in matches.
void randomize_matches(match *m, int n)
{
    // TODO: implement Fisher-Yates to shuffle the array.
    for(int i = 0; i < n-1; i++){
        int j = (rand()%(n-i))+i;
        match temp = m[i];
        m[i] = m[j];
        m[j] = temp;
    }
}

// Computes homography between two images given matching pixels.
// match *matches: matching points between images.
// int n: number of matches to use in calculating homography.
// returns: matrix representing homography H that maps image a to image b.
matrix compute_homography(match *matches, int n)
{
    matrix M = make_matrix(n*2, 8);
    matrix b = make_matrix(n*2, 1);

    int i;
    for(i = 0; i < n; ++i){
        double x  = matches[i].p.x; //mx
        double xp = matches[i].q.x; //nx
        double y  = matches[i].p.y; //my
        double yp = matches[i].q.y; //ny
        // TODO: fill in the matrices M and b.
        int row1 = 2*i;
        int row2 = row1 + 1;

        M.data[row1][0] = x;
        M.data[row1][1] = y;
        M.data[row1][2] = 1;
        M.data[row1][3] = 0;
        M.data[row1][4] = 0;
        M.data[row1][5] = 0;
        M.data[row1][6] = -xp*x;
        M.data[row1][7] = -xp*y;
        b.data[row1][0] = xp;

        M.data[row2][0] = 0;
        M.data[row2][1] = 0;
        M.data[row2][2] = 0;
        M.data[row2][3] = x;
        M.data[row2][4] = y;
        M.data[row2][5] = 1;
        M.data[row2][6] = -yp*x;
        M.data[row2][7] = -yp*y;
        b.data[row2][0] = yp;

    }
    matrix a = solve_system(M, b);
    free_matrix(M); free_matrix(b); 

    // If a solution can't be found, return empty matrix;
    matrix none = {0};
    if(!a.data) return none;

    matrix H = make_matrix(3, 3);
    // TODO: fill in the homography H based on the result in a.
    H.data[0][0] = a.data[0][0];
    H.data[0][1] = a.data[1][0];
    H.data[0][2] = a.data[2][0];
    H.data[1][0] = a.data[3][0];
    H.data[1][1] = a.data[4][0];
    H.data[1][2] = a.data[5][0];
    H.data[2][0] = a.data[6][0];
    H.data[2][1] = a.data[7][0];
    H.data[2][2] = 1;

    free_matrix(a);
    return H;
}

// Perform RANdom SAmple Consensus to calculate homography for noisy matches.
// match *m: set of matches.
// int n: number of matches.
// float thresh: inlier/outlier distance threshold.
// int k: number of iterations to run.
// int cutoff: inlier cutoff to exit early.
// returns: matrix representing most common homography between matches.
matrix RANSAC(match *m, int n, float thresh, int k, int cutoff)
{
    int best = 0;
    matrix Hb = make_translation_homography(256, 0);
    // TODO: fill in RANSAC algorithm.
    // for k iterations:
    //     shuffle the matches
    //     compute a homography with a few matches (how many??)
    //     if new homography is better than old (how can you tell?):
    //         compute updated homography using all inliers
    //         remember it and how good it is
    //         if it's better than the cutoff:
    //             return it immediately
    // if we get to the end return the best homography
    for(int i = 0; i < k; i++){
        randomize_matches(m, n);
        
        match *few_matches = calloc(4, sizeof(match));//at least 4 points for RANSAC
        for(int j = 0; j < 4; j++){
            few_matches[j] = m[j];
        }
        matrix homography = compute_homography(few_matches,4);

        int inliers_size = model_inliers(homography, m, n,thresh);
        if(inliers_size > best){
            Hb = compute_homography(m,inliers_size);
            best = inliers_size;
            if(best>cutoff){
                return Hb;
            }
        }
        
    }  
    return Hb;
}

// Stitches two images together using a projective transformation.
// image a, b: images to stitch.
// matrix H: homography from image a coordinates to image b coordinates.
// returns: combined image stitched together.
image combine_images(image a, image b, matrix H)
{
    matrix Hinv = matrix_invert(H);

    // Project the corners of image b into image a coordinates.
    point c1 = project_point(Hinv, make_point(0,0));
    point c2 = project_point(Hinv, make_point(b.w-1, 0));
    point c3 = project_point(Hinv, make_point(0, b.h-1));
    point c4 = project_point(Hinv, make_point(b.w-1, b.h-1));

    // Find top left and bottom right corners of image b warped into image a.
    point topleft, botright;
    botright.x = MAX(c1.x, MAX(c2.x, MAX(c3.x, c4.x)));
    botright.y = MAX(c1.y, MAX(c2.y, MAX(c3.y, c4.y)));
    topleft.x = MIN(c1.x, MIN(c2.x, MIN(c3.x, c4.x)));
    topleft.y = MIN(c1.y, MIN(c2.y, MIN(c3.y, c4.y)));

    // Find how big our new image should be and the offsets from image a.
    int dx = MIN(0, topleft.x);
    int dy = MIN(0, topleft.y);
    int w = MAX(a.w, botright.x) - dx;
    int h = MAX(a.h, botright.y) - dy;

    // Can disable this if you are making very big panoramas.
    // Usually this means there was an error in calculating H.
    if(w > 7000 || h > 7000){
         fprintf(stderr, "output too big, stopping\n");
         return copy_image(a);
    }

    int i,j,k;
    image c = make_image(w, h, a.c);
    
    // Paste image a into the new image offset by dx and dy.
    for(k = 0; k < a.c; ++k){
        for(j = 0; j < a.h; ++j){
            for(i = 0; i < a.w; ++i){
                set_pixel(c, i - dx, j - dy, k, get_pixel(a,i,j,k));
            }
        }
    }
    // TODO: Paste in image b as well.
    // You should loop over some points in the new image (which? all?)
    // and see if their projection from a coordinates to b coordinates falls
    // inside of the bounds of image b. If so, use bilinear interpolation to
    // estimate the value of b at that projection, then fill in image c.
    for(int x = topleft.x; x <= botright.x; x++){
        for(int y = topleft.y; y <= botright.y; y++){
            point p = make_point(x,y);
            point project_p = project_point(H, p);
            
            if((project_p.x>-0.5f && project_p.x<b.w-0.5f 
            && project_p.y>-0.5f && project_p.y<b.h-0.5f)&& x-dx >= 0 
            && x-dx < c.w && y-dy >= 0 
            && y-dy < c.h ){

            int b_nonempty = 1;

            int c1=0;
            for(int ch=0;ch<b.c;ch++)if(get_pixel(b,project_p.x,project_p.y,ch))c1++;
            if(c1==0)b_nonempty=0;
   
                for(int ch = 0; ch < c.c; ch++){
                    if(b_nonempty){
                        set_pixel(c, x-dx, y-dy, ch, bilinear_interpolate(b,project_p.x,project_p.y,ch));
                    }
                }
            }
        }
    }
    return c;
}

// Create a panoramam between two images.
// image a, b: images to stitch together.
// float sigma: gaussian for harris corner detector. Typical: 2
// float thresh: threshold for corner/no corner. Typical: 1-5
// int nms: window to perform nms on. Typical: 3
// float inlier_thresh: threshold for RANSAC inliers. Typical: 2-5
// int iters: number of RANSAC iterations. Typical: 1,000-50,000
// int cutoff: RANSAC inlier cutoff. Typical: 10-100
// int draw: flag to draw inliers.
image panorama_image(image a, image b, float sigma, float thresh, int nms, float inlier_thresh, int iters, int cutoff, int draw)
{
    srand(10);
    int an = 0;
    int bn = 0;
    int mn = 0;
    
    // Calculate corners and descriptors
    descriptor *ad = harris_corner_detector(a, sigma, thresh, nms, &an);
    descriptor *bd = harris_corner_detector(b, sigma, thresh, nms, &bn);

    // Find matches
    match *m = match_descriptors(ad, an, bd, bn, &mn);

    // Run RANSAC to find the homography
    matrix H = RANSAC(m, mn, inlier_thresh, iters, cutoff);

    if(draw){
        // Mark corners and matches between images
        mark_corners(a, ad, an);
        mark_corners(b, bd, bn);
        image inlier_matches = draw_inliers(a, b, H, m, mn, inlier_thresh);
        save_image(inlier_matches, "output/inliers");
    }

    free_descriptors(ad, an);
    free_descriptors(bd, bn);
    free(m);

    // Stitch the images together with the homography
    image comb = combine_images(a, b, H);
    return comb;
}

// Project an image onto a cylinder.
// image im: image to project.
// float f: focal length used to take image (in pixels).
// returns: image projected onto cylinder, then flattened.
image cylindrical_project(image im, float f)
{
    //TODO: project image onto a cylinder
    image c = make_image(im.w, im.h,im.c);
    
    float x_c = im.w/2.0;
    float y_c = im.h/2.0;

    for(int x = 0; x < c.w; x++){
        for(int y = 0; y < c.h; y++){
            float theta = (x-x_c)/f;
            float h = (y-y_c)/f;

            float x_p = sinf(theta);
            float y_p = h;
            float z_p = cosf(theta);

            float nx = (f*x_p/z_p) + x_c;
            float ny = (f*y_p/z_p) + y_c;

            if(nx>-0.5f && nx<im.w-0.5f && ny>-0.5f && ny<im.h-0.5f){
                for(int ch = 0; ch < c.c; ch++){
                    set_pixel(c,x,y,ch,bilinear_interpolate(im, nx, ny,ch));
                }
            }
        }
    }
    return c;
}
