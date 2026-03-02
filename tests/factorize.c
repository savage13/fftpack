#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define Treal double
#define NSPECIAL 4
#define MAXFAC 13
void factorize(int n, int ifac[MAXFAC+2], const int ntryh[NSPECIAL])
  /* Factorize n in factors in ntryh and rest. On exit,
ifac[0] contains n and ifac[1] contains number of factors,
the factors start from ifac[2]. */
  {
    int ntry=3, i, j=0, ib, nf=0, nl=n, nq, nr;
startloop:
    if (j < NSPECIAL)
      ntry = ntryh[j];
    else
      ntry+= 2;
    j++;
    //printf("%d %d\n", j, ntry);
    do {
      nq = nl / ntry;
      nr = nl - ntry*nq;
      //printf("   %d %d\n", nq, nr);
      if (nr != 0) goto startloop;
      nf++;
      ifac[nf + 1] = ntry;
      nl = nq;
      if (ntry == 2 && nf != 1) {
        for (i=2; i<=nf; i++) {
          ib = nf - i + 2;
          ifac[ib + 1] = ifac[ib];
        }
        ifac[2] = 2;
      }
    } while (nl != 1);
    ifac[0] = n;
    ifac[1] = nf;
  }


static void cffti1(int n, Treal wa[], int ifac[MAXFAC+2])
  {
    static const Treal twopi = 2.0 * M_PI; //= 6.28318530717959;
    Treal arg, argh, argld, fi;
    int idot, i, j;
    int i1, k1, l1, l2;
    int ld, ii, nf, ip;
    int ido, ipm;

    static const int ntryh[NSPECIAL] = {
      3,4,2,5    }; /* Do not change the order of these. */

    factorize(n,ifac,ntryh);
    nf = ifac[1];
    //printf("---------------------------------------\n");
    //printf("NF: %d N: %d \n", nf, n);
    argh = twopi/(Treal)n;
    i = 1;
    l1 = 1;
    //printf("%d\n", nf);
    for (k1=1; k1<=nf; k1++) {
      ip = ifac[k1+1];
      ld = 0;
      //printf("ip: %d ld: %d\n", ip, ld);
      l2 = l1*ip;
      ido = n / l2;
      idot = ido + ido + 2;
      ipm = ip - 1;
      //printf("%d %d\n", k1, ipm);
      for (j=1; j<=ipm; j++) {
        i1 = i;
        wa[i-1] = 1;
        wa[i] = 0;
        //printf("wa[%d,%d]_: %f %f (ld: %d fi: %f IPM: %d/%d)\n", i-1,i, wa[i-1],wa[i], ld, fi, j,ipm);
        ld += l1;
        fi = 0;
        argld = ld*argh;
        //printf("  %d\n",idot);
        for (ii=4; ii<=idot; ii+=2) {
          i+= 2;
          fi+= 1;
          arg = fi*argld;
          wa[i-1] = cos(arg);
          wa[i] = sin(arg);
          //printf("wa[%d,%d] : %f %f (ld: %d fi: %f l1,l2: %d %d)\n", i-1,i, wa[i-1],wa[i], ld, fi, l1,l2);
        }
        if (ip > 5) {
          wa[i1-1] = wa[i-1];
          wa[i1] = wa[i];
          //printf("wa[%d,%d]*: %f %f\n", i1-1,i1, wa[i-1],wa[i]);
        }
      }
      l1 = l2;
    }
    //for(ii = 0; ii <= i; ii++) {
    //    printf("%d/%d: %f\n", ii,i, wa[ii]);
    //}
  } /* cffti1 */


void cffti(int n, Treal wsave[]){
    int iw1, iw2;
    if (n == 1) return;
    iw1 = 2*n;
    iw2 = iw1 + 2*n;
    cffti1(n, wsave+iw1, (int*)(wsave+iw2));
  } /* cffti */

int main(int argc, char *argv[]) {
    int n = 6;
    int n1 = 6;

    int i = 0;
    Treal wsave[100000];
    Treal *wa = &wsave[2*n];
    int *ifac = (int *)(&wsave[2*n+2*n]);
    //int ifac[MAXFAC+2] = {0,0,0,0,0,0,0,0,0,0,0,0,0};
    int ntryh[NSPECIAL] = {3,4,2,5};
    factorize(36, ifac, ntryh);
    while(ifac[i] != 0) {
        //printf("%d: %d\n", i, ifac[i]);
        i++;
    }
    if(argc < 2) {
        return 0;
    }
    n1 = atoi(argv[1]);
    printf("import { cffti1 } from './fftpack.ts'\n");
    printf("import { deepEqualTolerance } from '../tests/deepEqualTolerance.js'\n");
    printf("function eq(a, b, tol = 2e-15) { return deepEqualTolerance(a, b, tol) }\n");
    for(n = 2; n <= n1; n++) {
        wa = &wsave[2*n];
        ifac = (int *)(&wsave[2*n+2*n]);

        cffti(n, wsave);
        printf("eq(cffti1(%d), { wa: [\n", n);
        for(i = 0; i < 2 * n; i+=2) {
            printf("  %.16f, %.16f,\n", wa[i], wa[i+1]);
        }
        printf("],\n");
        printf("ifac: [%d, %d, ", ifac[0], ifac[1]);
        for(i = 0; i < ifac[1]; i++) {
            printf("%d, ", ifac[2+i]);
        }
        printf("]\n");
        printf("})\n");
    }
    return 0;
}
