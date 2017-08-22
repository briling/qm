#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>

/*  a tool to prepare the triple products
 *  $B_{mm'm''}^{l\;\;l'\;\;l''}$
 *  of real spherical harmonics
 */

/*  notation:
 *   Y_{l,m}(r) -- complex spherical harmonic
 *   y_{l,m}(r) -- real    spherical harmonic
 *
 *  options:
 *   ./B c $m                          prints y(x,m) as a sum of Y(x,m) and Y(x,-m)
 *   ./B s $j1 $j2 $j3 $m1 $m2 $m3     prints Wigner's 3-jm symbol
 *   ./B r $j1 $j2 $j3 $m1 $m2 $m3     prints its analogue for real spherical harmonics
 *   ./B i $j1 $j2 $j3 $m1 $m2 $m3     prints \iint y_{l1,m1} y_{l2,m2} y_{l3,m3} \sin\theta d\phi d\theta
 *   ./B t                             prints non-zero integrals
 *   ./B f                             prints a C code
 *
 *   ($x denotes an integer argument)
 */


#define SQRT1_PI 0.56418958354775628695

const double f[] = {
  [ 0] =         1.0,  [ 1] =          1.0,   [ 2] =           2.0,  [ 3] =        6.0,
  [ 4] =        24.0,  [ 5] =        120.0,   [ 6] =         720.0,  [ 7] =     5040.0,
  [ 8] =     40320.0,  [ 9] =     362880.0,   [10] =     3628800.0,  [11] = 39916800.0,
  [12] = 479001600.0,  [13] = 6227020800.0,   [14] = 87178291200.0};
#define MAXF 14

const double f_sqrt[] = {
  [ 0] =     1.0, [ 1] =     1.0, [ 2] =     M_SQRT2,
  [ 3] =     2.44948974278317809819, [ 4] =     4.89897948556635619638,
  [ 5] =    10.95445115010332226913, [ 6] =    26.83281572999747635691,
  [ 7] =    70.99295739719539251080, [ 8] =   200.79840636817813151476,
  [ 9] =   602.39521910453439454428, [10] =  1904.94094396650522516116};
#define MAXFSQRT 10

double symbol(int j1, int j2, int j3, int m1, int m2, int m3){
/* Numerical Tables for Angular Correlation Computations:
   3j-, 6j-, 9j-Symbols, F- and Gamma-Coefficients
   pages 8, 25
 */
  if(m1+m2+m3){
    return 0.0;
  }
  if( (!m1) && (!m2) && (!m3) ){
    if((j1+j2+j3)%2){
      return 0.0;
    }
  }
  int i3 =  j1+j2-j3;
  int i2 =  j1-j2+j3;
  int i1 = -j1+j2+j3;
  if( (i1<0) || (i2<0) || (i3<0) ){
    return 0.0;
  }
  if( (i1>MAXFSQRT) || (i2>MAXFSQRT) || (i3>MAXFSQRT) ){
    abort();
  }
  if(j1+j2+j3+1 > MAXFSQRT){
    abort();
  }
  double t1  = f_sqrt[i3]*f_sqrt[i2]*f_sqrt[i1]/f_sqrt[j1+j2+j3+1];
  double t2  = f_sqrt[j1+m1]*f_sqrt[j1-m1]*f_sqrt[j2+m2]*f_sqrt[j2-m2]*f_sqrt[j3+m3]*f_sqrt[j3-m3];
  double t12 = t1*t2;
  if((j1-j2-m3)%2){
    t12 = -t12;
  }
  double s = 0.0;

  int k1 = 2*j1+1;
  int k2 = 2*j2+1;
  int k3 = j1+j2-j3+1;
  int K;
  K = k1>k2?k1:k2;
  if(k3>K){
    K = k3;
  }
  if(K>MAXF){
    abort();
  }

  for(int k=0; k<K; k++){
    int i1 = j1+j2-j3-k;
    int i2 = j1-m1-k;
    int i3 = j2+m2-k;
    int i4 = j3-j2+m1+k;
    int i5 = j3-j1-m2+k;
    if( (i1<0) || (i2<0) || (i3<0) || (i4<0) || (i5<0) ){
      continue;
    }
    double t = 1.0/(f[k]*f[i1]*f[i2]*f[i3]*f[i4]*f[i5]);
    s += ((k%2)?(-t):(t));
  }
  return t12*s;
}

void coeff(int m, complex double c[2]){
  /* Condon-Shortley phase convention */
  double complex c1;  // +|m|
  double complex c2;  // -|m|
  if(m<0){
    c1 = m%2 ? I*M_SQRT1_2 : -I*M_SQRT1_2;
    c2 = I*M_SQRT1_2;
  }
  else if(m==0){
    c1 = 0.5;
    c2 = 0.5;
  }
  else{
    c1 = m%2 ? -M_SQRT1_2 : M_SQRT1_2;
    c2 = M_SQRT1_2;
  }
  c[0] = c1;
  c[1] = c2;
  return;
}

double complex real_naive(int j1, int j2, int j3, int m1, int m2, int m3){

  double complex c1[2];
  double complex c2[2];
  double complex c3[2];

  coeff(m1, c1);
  coeff(m2, c2);
  coeff(m3, c3);

  double complex r =
    (c1[0]*c2[0]*c3[0])  *  symbol(j1,j2,j3,  abs(m1),  abs(m2),  abs(m3))  +
    (c1[1]*c2[0]*c3[0])  *  symbol(j1,j2,j3, -abs(m1),  abs(m2),  abs(m3))  +
    (c1[0]*c2[1]*c3[0])  *  symbol(j1,j2,j3,  abs(m1), -abs(m2),  abs(m3))  +
    (c1[1]*c2[1]*c3[0])  *  symbol(j1,j2,j3, -abs(m1), -abs(m2),  abs(m3))  +
    (c1[0]*c2[0]*c3[1])  *  symbol(j1,j2,j3,  abs(m1),  abs(m2), -abs(m3))  +
    (c1[1]*c2[0]*c3[1])  *  symbol(j1,j2,j3, -abs(m1),  abs(m2), -abs(m3))  +
    (c1[0]*c2[1]*c3[1])  *  symbol(j1,j2,j3,  abs(m1), -abs(m2), -abs(m3))  +
    (c1[1]*c2[1]*c3[1])  *  symbol(j1,j2,j3, -abs(m1), -abs(m2), -abs(m3))  ;

  return r;
}

double complex real(int j1, int j2, int j3, int m1, int m2, int m3){

  double complex r;
  if( (!m1)&&(!m2)&&(!m3) ){
    r = symbol(j1,j2,j3, m1,m2,m3);
  }

  else if( (m1)&&(!m2)&&(!m3) || (!m1)&&(m2)&&(!m3) || (!m1)&&(!m2)&&(m3) ){
    r = 0.0;
  }

  else if( (!m1)&&(m2)&&(m3) ){
    if(abs(m2)==abs(m3)){
      double complex c2[2];
      double complex c3[2];
      coeff(m2, c2);
      coeff(m3, c3);
      r =
        (c2[1]*c3[0])  *  symbol(j1,j2,j3, m1, -abs(m2),  abs(m3))  +
        (c2[0]*c3[1])  *  symbol(j1,j2,j3, m1,  abs(m2), -abs(m3))  ;
    }
    else{
      r = 0.0;
    }
  }
  else if( (m1)&&(!m2)&&(m3) ){
    if(abs(m1)==abs(m3)){
      double complex c1[2];
      double complex c3[2];
      coeff(m1, c1);
      coeff(m3, c3);
      r =
        (c1[1]*c3[0])  *  symbol(j1,j2,j3, -abs(m1), m2,  abs(m3))  +
        (c1[0]*c3[1])  *  symbol(j1,j2,j3,  abs(m1), m2, -abs(m3))  ;
    }
    else{
      r = 0.0;
    }
  }
  else if( (m1)&&(m2)&&(!m3) ){
    if(abs(m1)==abs(m2)){
      double complex c1[2];
      double complex c2[2];
      coeff(m1, c1);
      coeff(m2, c2);
      r =
        (c1[1]*c2[0])  *  symbol(j1,j2,j3, -abs(m1),  abs(m2), m3)  +
        (c1[0]*c2[1])  *  symbol(j1,j2,j3,  abs(m1), -abs(m2), m3)  ;
    }
    else{
      r = 0.0;
    }
  }

  else{
    double complex c1[2];
    double complex c2[2];
    double complex c3[2];
    coeff(m1, c1);
    coeff(m2, c2);
    coeff(m3, c3);
    int am1 = abs(m1);
    int am2 = abs(m2);
    int am3 = abs(m3);

    if(am1+am2==am3){
      r =
       (c1[0]*c2[0]*c3[1])  *  symbol(j1,j2,j3,  am1,  am2, -am3)  +
       (c1[1]*c2[1]*c3[0])  *  symbol(j1,j2,j3, -am1, -am2,  am3)  ;
   }
    else if(am1+am3==am2){
      r =
        (c1[0]*c2[1]*c3[0])  *  symbol(j1,j2,j3,  am1, -am2,  am3)  +
        (c1[1]*c2[0]*c3[1])  *  symbol(j1,j2,j3, -am1,  am2, -am3)  ;
    }
    else if(am2+am3==am1){
      r =
        (c1[1]*c2[0]*c3[0])  *  symbol(j1,j2,j3, -am1,  am2,  am3)  +
        (c1[0]*c2[1]*c3[1])  *  symbol(j1,j2,j3,  am1, -am2, -am3)  ;
    }
    else{
      r = 0.0;
    }
  }
  return r;
}

double integ_real(int j1, int j2, int j3, int m1, int m2, int m3){
  /* Varshalovich, Moskalev, Khersonsky
   * "Quantum Theory Of Angular Momemtum"
   * World Scientific, 1988
   * page 148
   */
  if((j1+j2+j3)%2){
    return 0.0;
  }

  if( (j1+j2<j3) || (j1+j3<j2) || (j2+j3<j1) ){
    return 0.0;
  }

  int am1 = abs(m1);
  int am2 = abs(m2);
  int am3 = abs(m3);
  if( (( m1)&&(!m2)&&(!m3)) ||
      ((!m1)&&( m2)&&(!m3)) ||
      ((!m1)&&(!m2)&&( m3)) ){
    return 0.0;
  }

  else if( (!m1)&&( m2)&&( m3) ){
    if( (am2!=am3) || (m2==-m3) ){
      return 0.0;
    }
  }
  else if( ( m1)&&(!m2)&&( m3) ){
    if( (am1!=am3) || (m1==-m3) ){
      return 0.0;
    }
  }
  else if( ( m1)&&( m2)&&(!m3) ){
    if( (am1!=am2) || (m1==-m2) ){
      return 0.0;
    }
  }

  else{
    if( (am1+am2!=am3) && (am1+am3!=am2) && (am2+am3!=am1) ){
    return 0.0;
    }
  }

  double t1 = sqrt((2.0*j1+1.0)*(2.0*j2+1.0)*(2.0*j3+1.0))*0.5*SQRT1_PI;
  double t2 = symbol(j1,j2,j3, 0,0,0);
  complex double t3 = real(j1,j2,j3, m1,m2,m3);
  complex double ret = t1*t2*t3;
  if(fabs(cimag(ret))>1e-15){
    printf("% d % d % d    % d % d % d     % e +% ei\n", j1,j2,j3, m1,m2,m3, creal(ret), cimag(ret));
    fflush(stdout);
    abort();
  }
  return creal(ret);
}

double integ_real_norm(int j1, int j2, int j3, int m1, int m2, int m3){
  /* computes integ_real(j1,j2,j3,m1,m2,m3)/integ_real(0,0,0,0,0,0) = integ_real(j1,j2,j3,m1,m2,m3)*(2.0*SQRTPI) */
  if((j1+j2+j3)%2){
    return 0.0;
  }

  if( (j1+j2<j3) || (j1+j3<j2) || (j2+j3<j1) ){
    return 0.0;
  }

  int am1 = abs(m1);
  int am2 = abs(m2);
  int am3 = abs(m3);
  if( (( m1)&&(!m2)&&(!m3)) ||
      ((!m1)&&( m2)&&(!m3)) ||
      ((!m1)&&(!m2)&&( m3)) ){
    return 0.0;
  }

  else if( (!m1)&&( m2)&&( m3) ){
    if( (am2!=am3) || (m2==-m3) ){
      return 0.0;
    }
  }
  else if( ( m1)&&(!m2)&&( m3) ){
    if( (am1!=am3) || (m1==-m3) ){
      return 0.0;
    }
  }
  else if( ( m1)&&( m2)&&(!m3) ){
    if( (am1!=am2) || (m1==-m2) ){
      return 0.0;
    }
  }

  else{
    if( (am1+am2!=am3) && (am1+am3!=am2) && (am2+am3!=am1) ){
    return 0.0;
    }
  }

  double t1 = sqrt((2.0*j1+1.0)*(2.0*j2+1.0)*(2.0*j3+1.0));
  double t2 = symbol(j1,j2,j3, 0,0,0);
  complex double t3 = real(j1,j2,j3, m1,m2,m3);
  complex double ret = t1*t2*t3;
  if(fabs(cimag(ret))>1e-15){
    printf("% d % d % d    % d % d % d     % e +% ei\n", j1,j2,j3, m1,m2,m3, creal(ret), cimag(ret));
    fflush(stdout);
    abort();
  }
  return creal(ret);
}

int main(int argc, char * argv[]){
  if(argv[1][0] == 's'){
    int j1 = atoi(argv[2]);
    int j2 = atoi(argv[3]);
    int j3 = atoi(argv[4]);
    int m1 = atoi(argv[5]);
    int m2 = atoi(argv[6]);
    int m3 = atoi(argv[7]);
    printf("% d % d % d    % d % d % d   % 20.15lf\n", j1,j2,j3, m1,m2,m3, symbol(j1,j2,j3,m1,m2,m3));
  }
  else if(argv[1][0] == 'c'){
    int m = atoi(argv[2]);
    complex double c[2];
    coeff(m, c);
    printf("y(x,%d)  =    (%lf+%lfi) Y(x,%d) + (%lf+%lfi) Y(x,%d) \n", m, creal(c[0]), cimag(c[0]), abs(m), creal(c[1]), cimag(c[1]), -abs(m)  );
  }
  else if(argv[1][0] == 'r'){
    int j1 = atoi(argv[2]);
    int j2 = atoi(argv[3]);
    int j3 = atoi(argv[4]);
    int m1 = atoi(argv[5]);
    int m2 = atoi(argv[6]);
    int m3 = atoi(argv[7]);
    complex double r = real(j1,j2,j3,m1,m2,m3);
    printf("% d % d % d    % d % d % d   % 20.15lf +% 20.15lfi\n", j1,j2,j3, m1,m2,m3, creal(r), cimag(r));
  }
  else if(argv[1][0] == 'i'){
    int j1 = atoi(argv[2]);
    int j2 = atoi(argv[3]);
    int j3 = atoi(argv[4]);
    int m1 = atoi(argv[5]);
    int m2 = atoi(argv[6]);
    int m3 = atoi(argv[7]);
    double r = integ_real(j1,j2,j3,m1,m2,m3);
    printf("% d % d % d    % d % d % d   % 20.15lf\n", j1,j2,j3, m1,m2,m3, r);
  }
  else if(argv[1][0] == 't'){
    for(int j1=0; j1<4; j1++){
      for(int j2=0; j2<=j1; j2++){
        for(int j3=0; j3<=j2; j3++){
          for(int m1=-j1; m1<=j1; m1++){
            for(int m2=-j2; m2<=j2; m2++){
              for(int m3=-j3; m3<=j3; m3++){
                double r = integ_real(j1,j2,j3,m1,m2,m3);
                if(fabs(r) < 1e-15){
                  continue;
                }
                printf("% d % d % d    % d % d % d   % 20.15lf\n", j1,j2,j3, m1,m2,m3, r);
              }
            }
          }
        }
      }
    }
  }
  else if(argv[1][0] == 'f'){
    for(int j1=0; j1<4; j1++){
      printf("if(l1 == %d){\n", j1);
      for(int j2=0; j2<=j1; j2++){
        printf("  if(l2 == %d){\n", j2);
        for(int j3=0; j3<=j2; j3++){
          printf("    if(l3 == %d){\n", j3);
          if( !(((j1+j2+j3)%2) || (j1+j2<j3) || (j1+j3<j2) || (j2+j3<j1)) ){
            for(int m1=-j1; m1<=j1; m1++){
              printf("      if(m1 == %d){\n", m1);
              for(int m2=-j2; m2<=j2; m2++){
                printf("        if(m2 == %d){\n", m2);
                for(int m3=-j3; m3<=j3; m3++){

                  int am1 = abs(m1);
                  int am2 = abs(m2);
                  int am3 = abs(m3);
                  if( (( m1)&&(!m2)&&(!m3)) ||
                      ((!m1)&&( m2)&&(!m3)) ||
                      ((!m1)&&(!m2)&&( m3)) ){
                    continue;
                  }
                  else if( (!m1)&&( m2)&&( m3) ){
                    if( (am2!=am3) || (m2==-m3) ){
                      continue;
                    }
                  }
                  else if( ( m1)&&(!m2)&&( m3) ){
                    if( (am1!=am3) || (m1==-m3) ){
                      continue;
                    }
                  }
                  else if( ( m1)&&( m2)&&(!m3) ){
                    if( (am1!=am2) || (m1==-m2) ){
                      continue;
                    }
                  }

                  else{
                    if( (am1+am2!=am3) && (am1+am3!=am2) && (am2+am3!=am1) ){
                      continue;
                    }
                  }
#if 0
                  double r = integ_real(j1,j2,j3,m1,m2,m3);
#else
                  double r = integ_real_norm(j1,j2,j3,m1,m2,m3);
#endif
                  printf("          if(m3 == %d){\n", m3);
                  printf("            return % 20.15lf;\n", r);
                  printf("          }\n");
                }
                printf("        }\n");
              }
              printf("      }\n");
            }
          }
          printf("    }\n");
        }
        printf("  }\n");
      }
      printf("}\n");
    }
  }

  return 0;
}

