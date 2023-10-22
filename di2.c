#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#define debug(...)

//#define R(t,x,y) _t=x;x-=t*y;y+=t*_t;_t=(3-x*x-y*y)/2;x*=_t;y*=_t;
#define R(s,x,y) x-=(y>>s);y+=(x>>s)
const int dz = 5, r1 = 1, r2 = 2;

// CORDIC algorithm to find magnitude of |x,y| in 8.8 fixed point
int length_cordic(int16_t x, int16_t y) {
  if (x<0) x=-x;
  //printf("%d %d (%f) ->", x, y, sqrt((float)x*x + (float)y*y));
  for (int i = 0; i < 8; i++) {
    int t=x;
    if (y<0) {
      x-=y>>i;
      y+=t>>i;
    } else {
      x+=y>>i;
      y-=t>>i;
    }
  }
  //printf(" %d %d -> %d\n", x, y, x * 155 >> 8);
   // * 155 >> 8; // correct for the 0.607 factor
  return (x >> 1) + (x >> 3); // divide by 0.625 as a cheap approximation to the 0.607 factor
}

void main() {
  /*
  float sA = 0, cA = 1,  // sines and cosines of
        sB = 0, cB = 1,  // donut rotation angles A and B
        */
  int16_t sBf = 0, cBf = 16384;
  int16_t sAf = 11583, cAf = 11583;
  for (;;) {
    int sA = sAf>>6, cA = cAf>>6;
    int sB = sBf>>6, cB = cBf>>6;
    int x0_16 = cB*sA;
    int x1_16 = cA*cB;

    int p0x = dz*sB;
    int p0y = dz*x0_16 >> 8;
    int p0z = -dz*x1_16 >> 8;

    int lxi = sB;
    int lyi = -cA + (cB*sA>>8);
    int lzi = (-cA*cB>>8) - sA;  // original lighting
    //float lx = cB, ly = -sA*sB, lz = cA*sB;  // lit from right side

    debug("sA=%d cA=%d sB=%d cB=%d x0=%d x1=%d p0x=%d p0y=%d p0z=%d lxi=%d lyi=%d lzi=%d\n", sA, cA, sB, cB, x0_16, x1_16, p0x, p0y, p0z, lxi, lyi, lzi);

    const int r1i = r1*256;
    const int r2i = r2*256;

    int niters = 0;
    int nnormals = 0;
    for (int j = 0; j < 23; j++) {
      float y = (j-11)*.05;
      for (int i = 0; i < 79; i++) {
        float x = (i-39.5)*.025;
        int x2 = sB*x;
        int vxi = cB*x - sB;
        int vyi = cA*y - (sA*x2 + x0_16 >> 8);
        int vzi = (x1_16 + cA*x2 >> 8) + sA*y;
        //int t = (int) (256 * dz) - r2i - r1i;
        int t = 512;
        debug("[%2d,%2d] vxyz (%+4d,%+4d,%+4d) p0xyz (%+4d,%+4d,%+4d) t=%4d ", j, i, vxi, vyi, vzi, p0x, p0y, p0z, t);

        int px = p0x + (vxi << 1); // assuming t = 512, t*vxi>>8 == vxi<<1
        int py = p0y + (vyi << 1);
        int pz = p0z + (vzi << 1);
        debug("pxyz (%+4d,%+4d,%+4d)\n", px, py, pz);
        for (;;) {
          int t0, t1, t2, d;
          t0 = length_cordic(px, py);
          t1 = t0 - r2i;
          t2 = length_cordic(pz, t1);
          d = t2 - r1i;
          debug("[%2d,%2d] (pz, t1) = (%d, %d) -> %d, d=%d, t0=%d\n", j, i, pz, t1, t2, d, t);
          t += d;
          debug("p(%d %d %d) -%d (%d, %d, %d)-> ", px, py, pz, d, vxi, vyi, vzi);
          px += d*vxi >> 8;
          py += d*vyi >> 8;
          pz += d*vzi >> 8;
          debug("p(%d %d %d); t=%d\n", px, py, pz, t);
          if (t > 8*256) {
            putchar(' ');
            break;
          } else if (d < 2) {
            int x4 = t0*t2;
            debug("lxi*px=%d x4=%d lyi*py=%d lzi*pz=%d t2=%d\n", lxi*px, x4, lyi*py, lzi*pz, t2);
            debug("N = %d + %d + %d = %d, ", lxi*px/x4, lyi*py/x4, lzi*pz/t2, lxi*px/x4 + lyi*py/x4 + lzi*pz/t2);
            int N = (lxi*px*t1/x4 + lyi*py*t1/x4 + lzi*pz/t2) >> 5;
            debug(">>5 = %d\n", N);
            putchar(".,-~:;!*=$@#"[N > 0 ? N < 12 ? N : 11 : 0]);
            nnormals++;
            break;
          }
          niters++;
        }
      }
      puts("");
    }
    printf("%d iterations %d lit pixels", niters, nnormals);
    fflush(stdout);
    R(5, cAf, sAf);
    R(6, cBf, sBf);
    usleep(15000);
    printf("\r\x1b[23A");
  }
}
