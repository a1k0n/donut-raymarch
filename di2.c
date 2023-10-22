#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#define debug(...)
//#define debug printf

#define R(s,x,y) x-=(y>>s);y+=(x>>s)
const int dz = 5, r1 = 1, r2 = 2;

// CORDIC algorithm to find magnitude of |x,y|
// also bring vector (x2,y2) along for the ride, and write back to x2
int length_cordic(int16_t x, int16_t y, int16_t *x2_, int16_t y2) {
  int x2 = *x2_;
  if (x < 0) { // start in right half-plane
    x = -x;
    x2 = -x2;
  }
  for (int i = 0; i < 8; i++) {
    int t = x;
    int t2 = x2;
    if (y < 0) {
      x -= y >> i;
      y += t >> i;
      x2 -= y2 >> i;
      y2 += t2 >> i;
    } else {
      x += y >> i;
      y -= t >> i;
      x2 += y2 >> i;
      y2 -= t2 >> i;
    }
  }
  // divide by 0.625 as a cheap approximation to the 0.607 scaling factor factor
  // introduced by this algorithm (see https://en.wikipedia.org/wiki/CORDIC)
  *x2_ = (x2 >> 1) + (x2 >> 3);
  return (x >> 1) + (x >> 3);
}

void main() {
  // high-precision rotation directions
  int16_t sB = 0, cB = 16384;
  int16_t sA = 11583, cA = 11583;
  int16_t sAsB = 0, cAsB = 0;
  int16_t sAcB = 11583, cAcB = 11583;
  for (;;) {
    int x1_16 = cAcB << 2;

    int p0x = dz * sB >> 6;
    int p0y = dz * sAcB >> 6;
    int p0z = -dz * cAcB >> 6;

    const int r1i = r1*256;
    const int r2i = r2*256;

    int niters = 0;
    int nnormals = 0;
    int yincC = (cA >> 4) + (cA >> 3);      // 12*cA >> 6;
    int yincS = (sA >> 4) + (sA >> 3);      // 12*sA >> 6;
    int xincX = (cB >> 5) + (cB >> 4);      // 6*cB >> 6;
    int xincY = (sAsB >> 5) + (sAsB >> 4);  // 6*sAsB >> 6;
    int xincZ = (cAsB >> 5) + (cAsB >> 4);  // 6*cAsB >> 6;
    int ycA = -((cA << 1) + (cA >> 2));     // -12 * yinc1 = -9*cA >> 2;
    int ysA = -((sA << 1) + (sA >> 2));     // -12 * yinc2 = -9*sA >> 2;
    for (int j = 0; j < 23; j++, ycA += yincC, ysA += yincS) {
      //int xsAsB = -40*xincY;
      int xsAsB = (sAsB >> 2) - (sAsB << 2);  // -40*xincY
      int xcAsB = (cAsB >> 2) - (cAsB << 2);  // -40*xincZ;

      int vxi16 = (cB >> 2) - (cB << 2) - (sB << 2); // -40*xincX - (sB << 2);
      int vyi16 = ycA - (sAcB<<2) - xsAsB;
      int vzi16 = ysA + (cAcB<<2) + xcAsB;

      for (int i = 0; i < 79; i++, vxi16 += xincX, vyi16 -= xincY, vzi16 += xincZ) {
        int t = 512; // (256 * dz) - r2i - r1i;

        int16_t px = p0x + (vxi16 >> 7); // assuming t = 512, t*vxi>>8 == vxi<<1
        int16_t py = p0y + (vyi16 >> 7);
        int16_t pz = p0z + (vzi16 >> 7);
        debug("pxyz (%+4d,%+4d,%+4d)\n", px, py, pz);
        for (;;) {
          int t0, t1, t2, d;
          int16_t lx = sB >> 2;
          int16_t ly = sAcB - cA >> 2;
          int16_t lz = -cAcB - sA >> 2;
          debug("[%2d,%2d] (px, py) = (%d, %d), (lx, ly) = (%d, %d) -> ", j, i, px, py, lx, ly);
          t0 = length_cordic(px, py, &lx, ly);
          debug("t0=%d (lx', ly') = (%d, %d)\n", t0, lx, ly);
          t1 = t0 - r2i;
          t2 = length_cordic(pz, t1, &lz, lx);
          d = t2 - r1i;
          t += d;
          px += d*vxi16 >> 16;
          py += d*vyi16 >> 16;
          pz += d*vzi16 >> 16;
          if (t > 8*256) {
            putchar(' ');
            break;
          } else if (d < 2) {
            int N = lz >> 9;
            putchar(".,-~:;!*=#$@"[N > 0 ? N < 12 ? N : 11 : 0]);
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
    R(5, cA, sA);
    R(5, cAsB, sAsB);
    R(5, cAcB, sAcB);
    R(6, cB, sB);
    R(6, cAcB, cAsB);
    R(6, sAcB, sAsB);
    usleep(15000);
    printf("\r\x1b[23A");
  }
}
