#include "eq.h"
#include "vec3.h"
#include "matrix.h"

static int find_min_proj(double z[3]){
  int i_min = 0;
  for(int i=1; i<3; i++){
    if( fabs(z[i]) < fabs(z[i_min]) ){
      i_min = i;
    }
  }
  return i_min;
}

static void new_xy_axis(int i, double * xx1, double * xx2, double x[3], double y[3], double z[3]){
  r3cpsc( x, z, -z[i] );
  x[i] += 1.0;
  double s2 = 1.0 / (1.0-z[i]*z[i]);
  double s1 = sqrt(s2);
  r3scal(x, s1);
  r3x(y, z, x);
  *xx1 = s1;
  *xx2 = s2;
  return;
}

double A_full(int l, int m1, int m2, double z[3]){
  if(l == 0){
    return 1.0;
  }
  int i_min = find_min_proj(z);
  double xx1, xx2;
  axis xyz;
  r3cp(xyz.z, z);
  new_xy_axis(i_min, &xx1, &xx2, xyz.x, xyz.y, z);
  return A(l, m1, m2, &xyz);
}

double A(int l, int m1, int m2, axis * xyz){

  double * z = xyz->z;
  double * y = xyz->y;
  double * x = xyz->x;

  if(l == 0){
    return 1.0;
  }

  if(l == 1){
    if((m2==-1)&&(m1==-1)) return  y[1]  ;
    if((m2== 0)&&(m1==-1)) return  y[2]  ;
    if((m2== 1)&&(m1==-1)) return  y[0]  ;
    if((m2==-1)&&(m1== 0)) return  z[1]  ;
    if((m2== 0)&&(m1== 0)) return  z[2]  ;
    if((m2== 1)&&(m1== 0)) return  z[0]  ;
    if((m2==-1)&&(m1== 1)) return  x[1]  ;
    if((m2== 0)&&(m1== 1)) return  x[2]  ;
    if((m2== 1)&&(m1== 1)) return  x[0]  ;
  }

  if(l == 2){
    if(m1 == -2){
      if( m2==-2) return  x[0]*y[1]+x[1]*y[0]  ;
      if( m2==-1) return  x[1]*y[2]+x[2]*y[1]  ;
      if( m2== 0) return  x[2]*y[2] * SQRT3    ;
      if( m2== 1) return  x[0]*y[2]+x[2]*y[0]  ;
      if( m2== 2) return  x[0]*y[0]-x[1]*y[1]  ;
    }

    else if(m1 == -1){
      if( m2==-2) return  y[0]*z[1]+y[1]*z[0]  ;
      if( m2==-1) return  y[1]*z[2]+y[2]*z[1]  ;
      if( m2== 0) return  y[2]*z[2] * SQRT3    ;
      if( m2== 1) return  y[0]*z[2]+y[2]*z[0]  ;
      if( m2== 2) return  y[0]*z[0]-y[1]*z[1]  ;
    }

    else if(m1 == 0){
      if( m2==-2) return  z[0]*z[1] * SQRT3                    ;
      if( m2==-1) return  z[1]*z[2] * SQRT3                    ;
      if( m2== 0) return  1.5*z[2]*z[2] - 0.5                  ;
      if( m2== 1) return  z[0]*z[2] * SQRT3                    ;
      if( m2== 2) return  (z[0]*z[0]-z[1]*z[1]) * 0.5 * SQRT3  ;
    }

    else if(m1 == 1){
      if( m2==-2) return  x[0]*z[1]+x[1]*z[0]  ;
      if( m2==-1) return  x[1]*z[2]+x[2]*z[1]  ;
      if( m2== 0) return  x[2]*z[2] * SQRT3    ;
      if( m2== 1) return  x[0]*z[2]+x[2]*z[0]  ;
      if( m2== 2) return  x[0]*z[0]-x[1]*z[1]  ;
    }

    else if(m1 == 2){
      if( m2==-2) return  x[0]*x[1]-y[0]*y[1]                              ;
      if( m2==-1) return  x[1]*x[2]-y[1]*y[2]                              ;
      if( m2== 0) return  (x[2]*x[2]-y[2]*y[2]) * 0.5 * SQRT3              ;
      if( m2== 1) return  x[0]*x[2]-y[0]*y[2]                              ;
      if( m2== 2) return  (x[0]*x[0]-x[1]*x[1]+y[1]*y[1]-y[0]*y[0]) * 0.5  ;
    }
  }
  GOTOHELL;
}

void distang(double * rij, axis * xyz, mol * m) {
  for(int i=0; i<m->n; i++){
    for(int j=i+1; j<m->n; j++){
      double z[3];
      r3diff(z, m->r+i*3, m->r+j*3);
      double r = sqrt(r3dot(z,z));
      r3scal(z, 1.0/r);
      int ij = i*m->n+j;
      int ji = j*m->n+i;
      xyz[ij].r = r;
      xyz[ji].r = r;
      r3cp  (xyz[ij].z, z);
      r3cpsc(xyz[ji].z, z, -1);

      int i_min = find_min_proj(z);
      double xx1,xx2;
      new_xy_axis(i_min, &xx1, &xx2, xyz[ij].x, xyz[ij].y, xyz[ij].z);
      new_xy_axis(i_min, &xx1, &xx2, xyz[ji].x, xyz[ji].y, xyz[ji].z);

      rij[mpos(i,j)] = r;
    }
  }
  return;
}

