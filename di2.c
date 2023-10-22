#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#define R(t,x,y) _t=x;x-=t*y;y+=t*_t;_t=(3-x*x-y*y)/2;x*=_t;y*=_t;
const float dz = 5, r1 = 1, r2 = 2;

// CORDIC algorithm to find magnitude of |x,y| in 8.8 fixed point, effectively multiplied by 0.607
int li(int x, int y) {
  if (x<0) x=-x;
  //printf("%d %d (%f) ->", x, y, sqrt((float)x*x + (float)y*y));
  for (int i = 0; i < 10; i++) {
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
  return x; // * 155 >> 8; // correct for the 0.607 factor
}

void main() {
  float sA = 0, cA = 1,  // sines and cosines of
        sB = 0, cB = 1,  // donut rotation angles A and B
        _t;
  for (;;) {
    float x0 = cB*sA,
          x1 = cA*cB,
          p0x = dz*sB,
          p0y = dz*x0,
          p0z = -dz*x1;
    float lx = sB, ly = -cA + cB*sA, lz = -cA*cB - sA;  // original lighting
    //float lx = cB, ly = -sA*sB, lz = cA*sB;  // lit from right side

    int lxi = 256*lx;
    int lyi = 256*ly;
    int lzi = 256*lz;

    int niters = 0;
    int nnormals = 0;
    for (int j = 0; j < 23; j++) {
      float y = (j-11)*.05;
      for (int i = 0; i < 79; i++) {
        float x = (i-39.5)*.025,
              x2 = sB*x,
              vx = cB*x - sB,
              vy = cA*y - sA*x2 - x0,
              vz = cA*x2 + sA*y + x1;
        int t = (int) (256 * (dz - r2 - r1));

        int px = (int)((p0x*256 + t * vx));
        int py = (int)((p0y*256 + t * vy));
        int pz = (int)((p0z*256 + t * vz));
        int vxi = (int)(vx*256);
        int vyi = (int)(vy*256);
        int vzi = (int)(vz*256);
        for (;;) {
          int t0, t1, t2, d;
          t0 = li(px, py) * 155 >> 8;
          t1 = t0 - 256*r2;
          t2 = li(pz, t1) * 155 >> 8;
          d = t2 - r1*256;
          //printf("[%d,%d] (pz, t1) = (%d, %d) -> %d, d=%d, t0=%d\n", i, j, pz, t1, t2, d, t);
          t += d;
          //printf("%d %d %d -%d (%0.3f, %0.3f, %0.3f)-> ", px, py, pz, d, vx, vy, vz);
          px += d*vxi >> 8;
          py += d*vyi >> 8;
          pz += d*vzi >> 8;
          //printf("[%d,%d] %d %d %d; t=%d\n", i, j, px, py, pz, t);
          if (t > 8*256) {
            putchar(' ');
            break;
          } else if (d < 2) {
            int x4 = t0*t2;
            //printf("lxi*px=%d x4=%d lyi*py=%d lzi*pz=%d t2=%d\n", lxi*px, x4, lyi*py, lzi*pz, t2);
            //printf("N = %d + %d + %d = %d, ", lxi*px/x4, lyi*py/x4, lzi*pz/t2, lxi*px/x4 + lyi*py/x4 + lzi*pz/t2);
            int N = (lxi*px*t1/x4 + lyi*py*t1/x4 + lzi*pz/t2) >> 5;
            //printf(">>5 = %d\n", N);
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
    R(.03, cA, sA);
    R(.011, cB, sB);
    usleep(15000);
    printf("\r\x1b[23A");
  }
}
