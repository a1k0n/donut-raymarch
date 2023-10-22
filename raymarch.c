#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#define R(t,x,y) _t=x;x-=t*y;y+=t*_t;_t=(3-x*x-y*y)/2;x*=_t;y*=_t;
const float dz = 5, r1 = 1, r2 = 2;


float l(float x, float y) {
  // instead of using sqrt(x*x+y*y) and linking a math library:
  float a = x*x + y*y;
  // three Newton steps of f(x) = x^2 - a
  // initial guess is 1; first two steps are unrolled into one expression
  x = a/(a+1)+a/4+.25;
  // and then we return the third step
  return a/(2*x) + x*.5;

  // this algorithm will be very easy to convert to fixed point, also.
}

void main() {
  float sA = 0, cA = 1,  // sines and cosines of
        sB = 0, cB = 1,  // donut rotation angles A and B
        sC = 0, cC = 1,
        _t;
  for (;;) {
    float x0 = cB*sA,
          x1 = cA*cB,
          p0x = dz*sB,
          p0y = dz*x0,
          p0z = -dz*x1;
    float lx = sB, ly = -cA + cB*sA, lz = -cA*cB - sA;  // original lighting
    // float lx = cB, ly = -sA*sB, lz = cA*sB;  // lit from right side

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
        float t = dz - (r2+0.5*sC) - r1;
        float d = t;

        for (;;) {
          float px = p0x+t*vx,
                py = p0y+t*vy,
                pz = p0z+t*vz;
          float t0 = l(px, py),
                t1 = -(r2+0.5*sC) + t0,
                t2 = l(pz, t1),
                d1 = t2 - r1;
          //printf("t=%f d=%f t0=%f t1=%f t2=%f d1=%f\n", t, d, t0, t1, t2, d1);
          if (t > 8) {
            putchar(' ');
            break;
          } else if (d1 < .01) {
            float x3 = 1.0/t2;
            float x4 = t1*x3/t0;
            int N = 8*(lx*px*x4 + ly*py*x4 + lz*pz*x3);
            putchar(".,-~:;!*=$@#"[N > 0 ? N < 12 ? N : 11 : 0]);
            nnormals++;
            break;
          }
          d = d1;
          t += d;
          niters++;
        }
      }
      puts("");
    }
    printf("%d iterations %d lit pixels", niters, nnormals);
    fflush(stdout);
    R(.011, cA, sA);
    R(.03, cB, sB);
    R(.003, cC, sC);
    usleep(15000);
    printf("\r\x1b[23A");
  }
}
