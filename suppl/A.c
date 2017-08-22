#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>

/*  a tool to prepare the transformation matrices
 *  $A_{mm'}^l(\vec{z})$
 *  of real spherical harmonics upon spatial rotation
 */

/*  notation:
 *   Y_{l,m}(r) -- complex spherical harmonic
 *   y_{l,m}(r) -- real    spherical harmonic
 *
 *  options:
 *   ./A c  $m1 $m2        prints <Y_{l,m1}|y_{l,m2}>, but wrong for <Y_{l,0}|y_{l,0}>
 *   ./A d  $l $m1 $m2     prints d^{l}_{m1,m2}, an element of Wigner's small d-matrix;
 *                         l \in [0,2], cb = cos(β), sb = sin(β)
 *   ./A r  $l $m1 $m2     prints <y_{l,m}(r')|y_{l,m}_(r)> as a function of the Euler angles α, β, γ
 *   ./A t                 prints a list of available <y(r')|y(r)>
 *
 *   ($x denotes an integer argument)
 */

complex double Yy(int m1, int m2){
  /* Condon-Shortley phase convention */
  if(abs(m1)!=abs(m2)){
    abort();
  }
  if(m2<0){
    if(m1>0) return m2%2 ? I*M_SQRT1_2 : -I*M_SQRT1_2;
    if(m1<0) return I*M_SQRT1_2;
  }
  else if(m2==0){
    return 0.5;
  }
  else{
    if(m1>0) return m2%2 ? -M_SQRT1_2 : M_SQRT1_2;
    if(m1<0) return M_SQRT1_2;
  }
}

void smalld(int l, int m1, int m2){
  /* Varshalovich, Moskalev, Khersonsky
   * "Quantum Theory Of Angular Momemtum"
   * World Scientific, 1988
   * page 119
   */
  if(l == 0) { printf("  1.0\n"); return; }
  if(l == 1){
    if     ( (m1 ==  1) && (m2 ==  1) ) { printf("  0.5*(1.0+cb)\n");  return; }
    else if( (m1 ==  1) && (m2 ==  0) ) { printf("  -sb*M_SQRT1_2\n"); return; }
    else if( (m1 ==  1) && (m2 == -1) ) { printf("  0.5*(1.0-cb)\n");  return; }
    else if( (m1 ==  0) && (m2 ==  1) ) { printf("  sb*M_SQRT1_2\n");  return; }
    else if( (m1 ==  0) && (m2 ==  0) ) { printf("  cb\n");            return; }
    else if( (m1 ==  0) && (m2 == -1) ) { printf("  -sb*M_SQRT1_2\n"); return; }
    else if( (m1 == -1) && (m2 ==  1) ) { printf("  0.5*(1.0-cb)\n");  return; }
    else if( (m1 == -1) && (m2 ==  0) ) { printf("  sb*M_SQRT1_2\n");  return; }
    else if( (m1 == -1) && (m2 == -1) ) { printf("  0.5*(1.0+cb)\n");  return; }
  }
  else if(l == 2){
    /* sqrt3_2 = sqrt(3.0/2.0) */
    /* sqrt3_8 = sqrt(3.0/8.0) */
    if     ( (m1 ==  2) && (m2 ==  2) ) { printf("  0.25*(1.0+cb)*(1.0+cb)\n"); return ;}
    else if( (m1 ==  2) && (m2 ==  1) ) { printf("  -0.5*sb*(1.0+cb)\n");       return ;}
    else if( (m1 ==  2) && (m2 ==  0) ) { printf("  sqrt3_8*sb*sb\n");          return ;}
    else if( (m1 ==  2) && (m2 == -1) ) { printf("  -0.5*sb*(1.0-cb)\n");       return ;}
    else if( (m1 ==  2) && (m2 == -2) ) { printf("  0.25*(1.0-cb)*(1.0-cb)\n"); return ;}
    else if( (m1 ==  1) && (m2 ==  2) ) { printf("  0.5*sb*(1.0+cb)\n");        return ;}
    else if( (m1 ==  1) && (m2 ==  1) ) { printf("  ( cb*cb+0.5*cb-0.5)\n");    return ;}
    else if( (m1 ==  1) && (m2 ==  0) ) { printf("  -sqrt3_2*sb*cb\n");         return ;}
    else if( (m1 ==  1) && (m2 == -1) ) { printf("  (-cb*cb+0.5*cb+0.5)\n");    return ;}
    else if( (m1 ==  1) && (m2 == -2) ) { printf("  -0.5*sb*(1.0-cb)\n");       return ;}
    else if( (m1 ==  0) && (m2 ==  2) ) { printf("  sqrt3_8*sb*sb\n");          return ;}
    else if( (m1 ==  0) && (m2 ==  1) ) { printf("  sqrt3_2*sb*cb\n");          return ;}
    else if( (m1 ==  0) && (m2 ==  0) ) { printf("  (1.5*cb*cb-0.5)\n");        return ;}
    else if( (m1 ==  0) && (m2 == -1) ) { printf("  -sqrt3_2*sb*cb\n");         return ;}
    else if( (m1 ==  0) && (m2 == -2) ) { printf("  sqrt3_8*sb*sb\n");          return ;}
    else if( (m1 == -1) && (m2 ==  2) ) { printf("  0.5*sb*(1.0-cb)\n");        return ;}
    else if( (m1 == -1) && (m2 ==  1) ) { printf("  (-cb*cb+0.5*cb+0.5)\n");    return ;}
    else if( (m1 == -1) && (m2 ==  0) ) { printf("  sqrt3_2*sb*cb\n");          return ;}
    else if( (m1 == -1) && (m2 == -1) ) { printf("  ( cb*cb+0.5*cb-0.5)\n");    return ;}
    else if( (m1 == -1) && (m2 == -2) ) { printf("  -0.5*sb*(1.0+cb)\n");       return ;}
    else if( (m1 == -2) && (m2 ==  2) ) { printf("  0.25*(1.0-cb)*(1.0-cb)\n"); return ;}
    else if( (m1 == -2) && (m2 ==  1) ) { printf("  0.5*sb*(1.0-cb)\n");        return ;}
    else if( (m1 == -2) && (m2 ==  0) ) { printf("  sqrt3_8*sb*sb\n");          return ;}
    else if( (m1 == -2) && (m2 == -1) ) { printf("  0.5*sb*(1.0+cb)\n");        return ;}
    else if( (m1 == -2) && (m2 == -2) ) { printf("  0.25*(1.0+cb)*(1.0+cb)\n"); return ;}
  }
  abort();
}

void yy(int l, int m, int m1){

  double complex c1 = conj(Yy( m,m)) * Yy( m1, m1);
  double complex c2 = conj(Yy(-m,m)) * Yy(-m1, m1);
  double complex c3 = conj(Yy( m,m)) * Yy(-m1, m1);
  double complex c4 = conj(Yy(-m,m)) * Yy( m1, m1);
  double complex c5 = conj(Yy( m,m)) * Yy( m1, m1) * ( I);
  double complex c6 = conj(Yy(-m,m)) * Yy(-m1, m1) * (-I);
  double complex c7 = conj(Yy( m,m)) * Yy(-m1, m1) * (-I);
  double complex c8 = conj(Yy(-m,m)) * Yy( m1, m1) * ( I);

  printf(" + (% lf + % lfi) * ", creal(c1), cimag(c1)); printf("cos(% d*alpha+(% d)*gamma) * ", m1, m); smalld(l, m1, m);
  printf(" + (% lf + % lfi) * ", creal(c2), cimag(c2)); printf("cos(% d*alpha+(% d)*gamma) * ", m1, m); smalld(l,-m1,-m);
  printf(" + (% lf + % lfi) * ", creal(c3), cimag(c3)); printf("cos(% d*alpha-(% d)*gamma) * ", m1, m); smalld(l,-m1, m);
  printf(" + (% lf + % lfi) * ", creal(c4), cimag(c4)); printf("cos(% d*alpha-(% d)*gamma) * ", m1, m); smalld(l, m1,-m);
  printf(" + (% lf + % lfi) * ", creal(c5), cimag(c5)); printf("sin(% d*alpha+(% d)*gamma) * ", m1, m); smalld(l, m1, m);
  printf(" + (% lf + % lfi) * ", creal(c6), cimag(c6)); printf("sin(% d*alpha+(% d)*gamma) * ", m1, m); smalld(l,-m1,-m);
  printf(" + (% lf + % lfi) * ", creal(c7), cimag(c7)); printf("sin(% d*alpha-(% d)*gamma) * ", m1, m); smalld(l,-m1, m);
  printf(" + (% lf + % lfi) * ", creal(c8), cimag(c8)); printf("sin(% d*alpha-(% d)*gamma) * ", m1, m); smalld(l, m1,-m);

}

int main(int argc, char * argv[]){
  if(argv[1][0] == 'c'){
    int m1 = atoi(argv[2]);
    int m2 = atoi(argv[3]);
    complex double c = Yy(m1,m2);
    printf("<Y(x,%d)|y(x,%d)> = % lf + % lfi\n", m1,m2, creal(c), cimag(c));
  }
  else if(argv[1][0] == 'r'){
    int l  = atoi(argv[2]);
    int m  = atoi(argv[3]);
    int m1 = atoi(argv[4]);
    yy(l, m, m1);
  }
  else if(argv[1][0] == 'd'){
    int l  = atoi(argv[2]);
    int m1 = atoi(argv[3]);
    int m2 = atoi(argv[4]);
    smalld(l, m1, m2);
  }
  else if(argv[1][0] == 't'){
    for(int l=0; l<3; l++){
      for(int m1=-l; m1<=l; m1++){
        for(int m=-l; m<=l; m++){
          printf("<y'(%d,% d)|y(%d,% d)> = \n", l,m,l,m1);
          yy(l, m, m1);
          printf("\n");
        }
      }
    }
  }
  return 0;
}

