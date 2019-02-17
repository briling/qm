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

static void dx_dz_fill(int i, double xx1, double xx2, double z[3], double dx_dz[9]){
  /* dx_dz(3*i+j) = dx[j] / dz[i] */
  vecset(9, dx_dz, 0.0);
  dx_dz[3*i+i] = -z[i]*xx1;
  for(int j=0; j<3; j++){
    if(j==i) continue;
    dx_dz[3*i+j] = -z[j]*xx1*xx2;
    dx_dz[3*j+j] = -z[i]*xx1;
  }
  return;
}

static void dy_dz_fill( double x[3], double z[3], double dx_dz[3], double dy_dz[9]){
  double dy_dx[9] = {
    + 0.0,  z[2], -z[1],
    -z[2],   0.0,  z[0],
    +z[1], -z[0],   0.0
  };
  double dy_dz_expl[9] = {
    + 0.0,  -x[2],  x[1],
    +x[2],    0.0, -x[0],
    -x[1],   x[0],   0.0
  };
  mx_multmx(3,3,3, dy_dz, dx_dz, dy_dx);
  vecadd(9, dy_dz, dy_dz_expl);
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

double A_grad_z(int l, int m1, int m2, double g[3], double z[3]){

  r3set(g, 0.0);

  if(l == 0){
    return 1.0;
  }

  int i_min = find_min_proj(z);
  double xx1, xx2, x[3], y[3];
  new_xy_axis(i_min, &xx1, &xx2, x, y, z);

  double A=0, dx_dz[9], dy_dz[9];
  dx_dz_fill(i_min, xx1, xx2, z, dx_dz);
  dy_dz_fill(x, z, dx_dz, dy_dz);

  if(l == 1){

    if(m1==-1){
      double dAdy[3]={};
      if(m2==-1)  dAdy[0] = 0.0, dAdy[1] = 1.0, dAdy[2] = 0.0, A = y[1];
      if(m2== 0)  dAdy[0] = 0.0, dAdy[1] = 0.0, dAdy[2] = 1.0, A = y[2];
      if(m2== 1)  dAdy[0] = 1.0, dAdy[1] = 0.0, dAdy[2] = 0.0, A = y[0];
      r3mx(g, dAdy, dy_dz);
    }

    else if(m1==0){
      if(m2==-1)  g[0] =  0.0, g[1] =  1.0, g[2] =  0.0, A = z[1];
      if(m2== 0)  g[0] =  0.0, g[1] =  0.0, g[2] =  1.0, A = z[2];
      if(m2== 1)  g[0] =  1.0, g[1] =  0.0, g[2] =  0.0, A = z[0];
    }

    else if(m1==1){
      double dAdx[3]={};
      if(m2==-1)  dAdx[0] = 0.0, dAdx[1] = 1.0, dAdx[2] = 0.0, A = x[1];
      if(m2== 0)  dAdx[0] = 0.0, dAdx[1] = 0.0, dAdx[2] = 1.0, A = x[2];
      if(m2== 1)  dAdx[0] = 1.0, dAdx[1] = 0.0, dAdx[2] = 0.0, A = x[0];
      r3mx(g, dAdx, dx_dz);
    }

  }

  else if(l == 2){

    if(m1==-2){
      double dAdx[3]={}, dAdy[3]={}, tmp[3];
      if(m2==-2)  dAdx[0] = y[1],  dAdx[1] =  y[0],  dAdx[2] =  0.0,        dAdy[0] = x[1],  dAdy[1] =  x[0],   dAdy[2] =  0.0       ,  A = x[0]*y[1]+x[1]*y[0];
      if(m2==-1)  dAdx[0] =  0.0,  dAdx[1] =  y[2],  dAdx[2] = y[1],        dAdy[0] =  0.0,  dAdy[1] =  x[2],   dAdy[2] = x[1]       ,  A = x[1]*y[2]+x[2]*y[1];
      if(m2== 0)  dAdx[0] =  0.0,  dAdx[1] =   0.0,  dAdx[2] = y[2]*SQRT3,  dAdy[0] =  0.0,  dAdy[1] =   0.0,   dAdy[2] = x[2]*SQRT3 ,  A = x[2]*y[2] * SQRT3  ;
      if(m2== 1)  dAdx[0] = y[2],  dAdx[1] =   0.0,  dAdx[2] = y[0],        dAdy[0] = x[2],  dAdy[1] =   0.0,   dAdy[2] = x[0]       ,  A = x[0]*y[2]+x[2]*y[0];
      if(m2== 2)  dAdx[0] = y[0],  dAdx[1] = -y[1],  dAdx[2] =  0.0,        dAdy[0] = x[0],  dAdy[1] = -x[1],   dAdy[2] =  0.0       ,  A = x[0]*y[0]-x[1]*y[1];
      r3mx(g, dAdy, dy_dz);
      r3mx(tmp, dAdx, dx_dz);
      r3add(g, tmp);
    }

    else if(m1==-1){
      double dAdy[3]={}, tmp[3];
      if(m2==-2) dAdy[0] = z[1], dAdy[1] = z[0],  dAdy[2] =  0.0      ,  g[0] = y[1], g[1] = y[0] , g[2] =  0.0,       A = y[0]*z[1]+y[1]*z[0] ;
      if(m2==-1) dAdy[0] =  0.0, dAdy[1] = z[2],  dAdy[2] = z[1]      ,  g[0] =  0.0, g[1] = y[2] , g[2] = y[1],       A = y[1]*z[2]+y[2]*z[1] ;
      if(m2== 0) dAdy[0] =  0.0, dAdy[1] =  0.0,  dAdy[2] = z[2]*SQRT3,  g[0] =  0.0, g[1] =  0.0 , g[2] = y[2]*SQRT3, A = y[2]*z[2] * SQRT3   ;
      if(m2== 1) dAdy[0] = z[2], dAdy[1] =  0.0,  dAdy[2] = z[0]      ,  g[0] = y[2], g[1] =  0.0 , g[2] = y[0],       A = y[0]*z[2]+y[2]*z[0] ;
      if(m2== 2) dAdy[0] = z[0], dAdy[1] =-z[1],  dAdy[2] =  0.0      ,  g[0] = y[0], g[1] =-y[1] , g[2] =  0.0,       A = y[0]*z[0]-y[1]*z[1] ;
      r3mx(tmp, dAdy, dy_dz);
      r3add(g, tmp);
    }

    else if(m1==0){
      if(m2==-2)  g[0] =  SQRT3 * z[1] ,  g[1] =  SQRT3 * z[0]     ,  g[2] =  0.0          , A = SQRT3 * z[0] * z[1];
      if(m2==-1)  g[0] =  0.0          ,  g[1] =  SQRT3 * z[2]     ,  g[2] =  SQRT3 * z[1] , A = SQRT3 * z[1] * z[2];
      if(m2== 0)  g[0] =  0.0          ,  g[1] =  0.0              ,  g[2] =  3.0*z[2]     , A = 1.5*z[2]*z[2]-0.5;
      if(m2== 1)  g[0] =  SQRT3 * z[2] ,  g[1] =  0.0              ,  g[2] =  SQRT3 * z[0] , A = SQRT3 * z[0]*z[2];
      if(m2== 2)  g[0] =  0.0          ,  g[1] = -2.0*SQRT3 * z[1] ,  g[2] = -SQRT3 * z[2] , A = SQRT3 * ((1.0-z[2]*z[2])*0.5-z[1]*z[1]);
    }

    else if(m1==1){
      double dAdx[3]={}, tmp[3];
      if(m2==-2)  dAdx[0] = z[1], dAdx[1] = z[0], dAdx[2] =  0.0,       g[0] = x[1], g[1] = x[0], g[2] =  0.0,        A = x[0]*z[1]+x[1]*z[0] ;
      if(m2==-1)  dAdx[0] =  0.0, dAdx[1] = z[2], dAdx[2] = z[1],       g[0] =  0.0, g[1] = x[2], g[2] = x[1],        A = x[1]*z[2]+x[2]*z[1] ;
      if(m2== 0)  dAdx[0] =  0.0, dAdx[1] =  0.0, dAdx[2] = z[2]*SQRT3, g[0] =  0.0, g[1] =  0.0, g[2] = x[2]*SQRT3,  A = x[2]*z[2] * SQRT3   ;
      if(m2== 1)  dAdx[0] = z[2], dAdx[1] =  0.0, dAdx[2] = z[0],       g[0] = x[2], g[1] =  0.0, g[2] = x[0],        A = x[0]*z[2]+x[2]*z[0] ;
      if(m2== 2)  dAdx[0] = z[0], dAdx[1] =-z[1], dAdx[2] =  0.0,       g[0] = x[0], g[1] =-x[1], g[2] =  0.0,        A = x[0]*z[0]-x[1]*z[1] ;
      r3mx(tmp, dAdx, dx_dz);
      r3add(g, tmp);
    }

    else if(m1==2){
      double dAdx[3]={}, dAdy[3]={}, tmp[3];
      if( m2==-2) dAdx[0] =  x[1], dAdx[1] =  x[0], dAdx[2] =   0.0,       dAdy[0] = -y[1], dAdy[1] = -y[0], dAdy[2] =   0.0,       A = x[0]*x[1]-y[0]*y[1];
      if( m2==-1) dAdx[0] =   0.0, dAdx[1] =  x[2], dAdx[2] =  x[1],       dAdy[0] =   0.0, dAdy[1] = -y[2], dAdy[2] = -y[1],       A = x[1]*x[2]-y[1]*y[2];
      if( m2==-0) dAdx[0] =   0.0, dAdx[1] =   0.0, dAdx[2] =  x[2]*SQRT3, dAdy[0] =   0.0, dAdy[1] =   0.0, dAdy[2] = -y[2]*SQRT3, A = (x[2]*x[2]-y[2]*y[2]) * 0.5 * SQRT3;
      if( m2== 1) dAdx[0] =  x[2], dAdx[1] =   0.0, dAdx[2] =  x[0],       dAdy[0] = -y[2], dAdy[1] =   0.0, dAdy[2] = -y[0],       A = x[0]*x[2]-y[0]*y[2];
      if( m2== 2) dAdx[0] =  x[0], dAdx[1] = -x[1], dAdx[2] =   0.0,       dAdy[0] = -y[0], dAdy[1] =  y[1], dAdy[2] =   0.0,       A = (x[0]*x[0]-x[1]*x[1]+y[1]*y[1]-y[0]*y[0]) * 0.5;
      r3mx(g, dAdy, dy_dz);
      r3mx(tmp, dAdx, dx_dz);
      r3add(g, tmp);
    }
  }
  return A;
}

void A_grad_z2r(double g[3], double z[3], double r1){
  double s = r3dot(z, g);
  r3adds(g, z, -s);
  r3scal(g, r1);
  return;
}

void A_grad_test(){
  double r[3] = {2,3,4};
  double d = 1e-4;
  for(int l=2; l<=2; l++){
    for(int m1=-l; m1<=l; m1++){
      for(int m2=-l; m2<=l; m2++){
        double z[3], ga[3], gn[3], dg[3];
        double r1 = 1.0/sqrt(r3dot(r,r));
        r3cpsc(z, r, r1);
        A_grad_z(l, m1, m2, ga, z);
        A_grad_z2r(ga, z, r1);
        for(int i=0; i<3; i++){
          double t = r[i];
          r[i] = t+d;
          r1 = 1.0/sqrt(r3dot(r,r));
          r3cpsc(z, r, r1);
          double A1 = A_full(l, m1, m2, z);
          r[i] = t-d;
          r1 = 1.0/sqrt(r3dot(r,r));
          r3cpsc(z, r, r1);
          double A2 = A_full(l, m1, m2, z);
          r[i] = t;
          gn[i] = (A1-A2)*0.5/d;
        }
        r3diff(dg, gn, ga);
        printf("l = %d  m1 =% d  m2 =% d%s\n", l, m1, m2, r3dot(dg,dg)>1e-15?" \e[1;31mfail\e[0m":"");
        r3print(ga, stdout);
        r3print(gn, stdout);
        r3print(dg, stdout);
      }
    }
  }
  return;
}

void distang(double * rij, axis * xyz, mol * m) { //TODO
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

