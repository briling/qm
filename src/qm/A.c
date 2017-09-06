#include "eq.h"
#include "vec3.h"
#include "matrix.h"

double A(int l, int m1, int m2, double z[3]){

  if(l == 0){
    return 1.0;
  }

  double cb,sb,cg,sg; // ca = 1.0; sa = 0.0

  cb = z[2];
  double acb = fabs(cb);
  if(acb < 1e-15){
    sb =  1.0;
    cg = -z[0];
    sg =  z[1];
  }
  else if(fabs(acb-1.0)>1e-15){
    sb = sqrt(1.0 - cb*cb);
    double sb1 = 1.0/sb;
    cg = -z[0]*sb1;
    sg =  z[1]*sb1;
  }
  else{
    sb = 0.0;
    cg = 1.0;
    sg = 0.0;
  }

  if(l == 1){
    if((m2==-1)&&(m1==-1)) return    cg    ;
    if((m2== 0)&&(m1==-1)) return    0.0   ;
    if((m2== 1)&&(m1==-1)) return    sg    ;
    if((m2==-1)&&(m1== 0)) return    sg*sb ;
    if((m2== 0)&&(m1== 0)) return    cb    ;
    if((m2== 1)&&(m1== 0)) return   -sb*cg ;
    if((m2==-1)&&(m1== 1)) return   -sg*cb ;
    if((m2== 0)&&(m1== 1)) return    sb    ;
    if((m2== 1)&&(m1== 1)) return    cg*cb ;
  }

  if(l == 2){

    if((m2==-2)&&(m1==-2)) return    (1.0-2.0*sg*sg)*cb;
    if((m2==-1)&&(m1==-2)) return    cg*sb;
    if((m2== 0)&&(m1==-2)) return    0.0;
    if((m2== 1)&&(m1==-2)) return    sg*sb;
    if((m2== 2)&&(m1==-2)) return    2.0*sg*cg*cb;

    if((m2==-2)&&(m1==-1)) return    sb*(2.0*sg*sg-1.0);
    if((m2==-1)&&(m1==-1)) return    cg*cb;
    if((m2== 0)&&(m1==-1)) return    0.0;
    if((m2== 1)&&(m1==-1)) return    sg*cb;
    if((m2== 2)&&(m1==-1)) return   -2.0*sb*sg*cg;

    if((m2==-2)&&(m1== 0)) return   -SQRT3 * sg*cg*sb*sb;
    if((m2==-1)&&(m1== 0)) return    SQRT3 * sg*sb*cb;
    if((m2== 0)&&(m1== 0)) return    1.5*cb*cb-0.5;
    if((m2== 1)&&(m1== 0)) return   -SQRT3 * cg*sb*cb;
    if((m2== 2)&&(m1== 0)) return    SQRT3 * sb*sb*(0.5-sg*sg);

    if((m2==-2)&&(m1== 1)) return    2.0*sb*sg*cg*cb;
    if((m2==-1)&&(m1== 1)) return    sg*(1.0-2.0*cb*cb);
    if((m2== 0)&&(m1== 1)) return    SQRT3 * sb*cb;
    if((m2== 1)&&(m1== 1)) return    cg*(2.0*cb*cb-1.0);
    if((m2== 2)&&(m1== 1)) return    sb*cb*(2.0*sg*sg-1.0);

    if((m2==-2)&&(m1== 2)) return -sg*cg*(1.0+cb*cb);
    if((m2==-1)&&(m1== 2)) return -cb*sb*sg;
    if((m2== 0)&&(m1== 2)) return SQRT3*0.5*sb*sb;
    if((m2== 1)&&(m1== 2)) return sb*cg*cb;
    if((m2== 2)&&(m1== 2)) return (0.5-sg*sg)*(1.0+cb*cb);
  }
  abort();
}

static euler z2eu(double z[3]){
  double cb,sb,cg,sg; // ca = 1.0; sa = 0.0
  cb = z[2];
  double acb = fabs(cb);
  if(acb < 1e-15){
    sb =  1.0;
    cg = -z[0];
    sg =  z[1];
  }
  else if(fabs(acb-1.0)>1e-15){
    sb = sqrt(1.0 - cb*cb);
    double sb1 = 1.0/sb;
    cg = -z[0]*sb1;
    sg =  z[1]*sb1;
#if 1
    if(sg > 1.0){
      sg = 1.0;
      cg = 0.0;
    }
    else if(sg < -1.0){
      sg = -1.0;
      cg = 0.0;
    }
#endif
  }
  else{
    sb = 0.0;
    cg = 1.0;
    sg = 0.0;
  }
  euler eu;
  eu.cos_b = cb;
  eu.sin_b = sb;
  eu.cos_g = cg;
  eu.sin_g = sg;
  return eu;
}

void distang(double * rij, euler * eu, mol * m) {
  for(int i=0; i<m->n; i++){
    for(int j=i+1; j<m->n; j++){
      double dij[3];
      double dji[3];
      r3diff(dij, m->r+i*3, m->r+j*3);
      double r = sqrt(r3dot(dij,dij));
      r3scal(dij, 1.0/r);
      r3cpsc(dji, dij, -1);
      int ij = (i*m->n+j);
      int ji = (j*m->n+i);
      eu[ij] = z2eu(dij);
      eu[ji] = z2eu(dji);
      eu[ij].r = r;
      eu[ji].r = r;
      rij[mpos(i,j)] = r;
    }
  }
  return;
}

double A_new(int l, int m1, int m2, euler * eu){

  if(l == 0){
    return 1.0;
  }

  double cb = eu->cos_b;
  double sb = eu->sin_b;
  double cg = eu->cos_g;
  double sg = eu->sin_g;

  if(l == 1){
    if((m2==-1)&&(m1==-1)) return    cg    ;
    if((m2== 0)&&(m1==-1)) return    0.0   ;
    if((m2== 1)&&(m1==-1)) return    sg    ;
    if((m2==-1)&&(m1== 0)) return    sg*sb ;
    if((m2== 0)&&(m1== 0)) return    cb    ;
    if((m2== 1)&&(m1== 0)) return   -sb*cg ;
    if((m2==-1)&&(m1== 1)) return   -sg*cb ;
    if((m2== 0)&&(m1== 1)) return    sb    ;
    if((m2== 1)&&(m1== 1)) return    cg*cb ;
  }

  if(l == 2){

    if((m2==-2)&&(m1==-2)) return    (1.0-2.0*sg*sg)*cb;
    if((m2==-1)&&(m1==-2)) return    cg*sb;
    if((m2== 0)&&(m1==-2)) return    0.0;
    if((m2== 1)&&(m1==-2)) return    sg*sb;
    if((m2== 2)&&(m1==-2)) return    2.0*sg*cg*cb;

    if((m2==-2)&&(m1==-1)) return    sb*(2.0*sg*sg-1.0);
    if((m2==-1)&&(m1==-1)) return    cg*cb;
    if((m2== 0)&&(m1==-1)) return    0.0;
    if((m2== 1)&&(m1==-1)) return    sg*cb;
    if((m2== 2)&&(m1==-1)) return   -2.0*sb*sg*cg;

    if((m2==-2)&&(m1== 0)) return   -SQRT3 * sg*cg*sb*sb;
    if((m2==-1)&&(m1== 0)) return    SQRT3 * sg*sb*cb;
    if((m2== 0)&&(m1== 0)) return    1.5*cb*cb-0.5;
    if((m2== 1)&&(m1== 0)) return   -SQRT3 * cg*sb*cb;
    if((m2== 2)&&(m1== 0)) return    SQRT3 * sb*sb*(0.5-sg*sg);

    if((m2==-2)&&(m1== 1)) return    2.0*sb*sg*cg*cb;
    if((m2==-1)&&(m1== 1)) return    sg*(1.0-2.0*cb*cb);
    if((m2== 0)&&(m1== 1)) return    SQRT3 * sb*cb;
    if((m2== 1)&&(m1== 1)) return    cg*(2.0*cb*cb-1.0);
    if((m2== 2)&&(m1== 1)) return    sb*cb*(2.0*sg*sg-1.0);

    if((m2==-2)&&(m1== 2)) return -sg*cg*(1.0+cb*cb);
    if((m2==-1)&&(m1== 2)) return -cb*sb*sg;
    if((m2== 0)&&(m1== 2)) return SQRT3*0.5*sb*sb;
    if((m2== 1)&&(m1== 2)) return sb*cg*cb;
    if((m2== 2)&&(m1== 2)) return (0.5-sg*sg)*(1.0+cb*cb);
  }
  abort();
}

