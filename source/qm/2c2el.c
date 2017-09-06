#include "eq.h"
#include "common.h"
#include "vec3.h"

#define EPS 1e-15

static int qlll_mm(int qu, int lu, int lv, int l, double * q, qmdata * qmd){
  if(l==0){
    *q = 1.0;
    return 1;
  }
  if(lu>lv){
    int t;
    SWAP(lu,lv,t);
  }
  int bra = qmd->q_list[qu-1].qa;
  int ket = qmd->q_list[qu  ].qa;
  for(int i=bra; i<ket; i++){
    if ( (lu == qmd->qa[i].lu) &&
         (lv == qmd->qa[i].lv) &&
         (l  == qmd->qa[i].l)){
      *q = qmd->qa[i].qa;
      return 1;
    }
  }
  return 0;
}

static int qlll_pm(int qu, int la, int lv, int l, double * q, qmdata * qmd){
  int bra = qmd->q_list[qu-1].q1a;
  int ket = qmd->q_list[qu  ].q1a;
  for(int i=bra; i<ket; i++){
    if ( (la == qmd->q1a[i].lu) &&
         (lv == qmd->q1a[i].lv) &&
         (l  == qmd->q1a[i].l)){
      *q = qmd->q1a[i].qa;
      return 1;
    }
  }
  return 0;
}

static double fundconst(int l, int l1, int m){
  const double q0l[] = {
    [0] = 1.0,
    [1] = SQRT3,
    [2] = 6.70820393249936908920
  };
  const double q11[] = {
    [0] = 6.0,
    [1] = 3.0
  };
  const double q12[] = {
    [0] = 34.85685011586675196661,
    [1] = 20.12461179749810726768
  };
  const double q22[] = {
    [0] = 270.0,
    [1] = 180.0,
    [2] =  45.0
  };
  if(l>l1){
    int t;
    SWAP(l,l1,t);
  }
  if(l==0){
    return q0l[l1];
  }
  unsigned int am = abs(m);
  if( l==1 && l1==1 ){
    return q11[am];
  }
  else if( l==2 && l1==2 ){
    return q22[am];
  }
  else if( l==1 && l1==2){
    return q12[am];
  }
  GOTOHELL;
}

double G_eq52_mmmm(int m, int l, int l1, int lu, int lv, int lu1, int lv1, int qu, int qu1, double r, qmdata * qmd){
  double q1,q2,q3;
  if( ! qlll_mm(qu,  lu, lv, l,  &q1, qmd)) GOTOHELL;
  if( ! qlll_mm(qu1, lu1,lv1,l1, &q2, qmd)) GOTOHELL;
  q3 = fundconst(l, l1, m);
  double au  = qmd->agb[qu*qmd->nLo+lu];
  double av  = qmd->agb[qu*qmd->nLo+lv];
  double au1 = qmd->agb[qu1*qmd->nLo+lu1];
  double av1 = qmd->agb[qu1*qmd->nLo+lv1];
  double q   = q1*q2*q3;
  //хорошо бы вычислять их один раз за цикл, но только если нужны TODO
  double aa  = au+av;
  double aa1 = au1+av1;
  double a   = aa*aa1/(aa+aa1);
  if( (l==0) && (l1==0)){
    return g_eq43_l0(r, q, a, -0.5);
  }
  else{
    int aph = (l+m)%2?-1:1;
    return aph*g_eq43_c0(l+l1, r, q, a);
  }
}

double G_eq52_mmmp(int m, int l, int l1, int lu, int lv, int lu1, int lv1, int qu, int qu1, double q1, double q2, double r, qmdata * qmd){
  double q3 = fundconst(l, l1, m);
  double au  = qmd->ag1b[qu *qmd->nLv+lu ];
  double av  = qmd->ag0b[qu *qmd->nLo+lv ];
  double au1 = qmd->ag0b[qu1*qmd->nLo+lu1];
  double av1 = qmd->ag0b[qu1*qmd->nLo+lv1];
  double q   = q1*q2*q3;
  double aa  = au+av;
  double aa1 = au1+av1;
  double a   = aa*aa1/(aa+aa1);
  int    aph = (l+m)%2?-1:1;
  return aph*g_eq43_c0(l+l1, r, q, a);
}

double R0_eq39_mmmm_old(int mu, int mv, int mu1, int mv1, int lu, int lv, int lu1, int lv1, int qu, int qu1, double ruu1[3], qmdata * qmd){
  double z[3];
  double r = sqrt(r3dot(ruu1, ruu1));  //eq29
  r3cpsc(z, ruu1, 1.0/r);              //eq29
  double ret = 0.0;
  for(int m_=-lu;m_<=lu;m_++){
    double A1 = A(lu,m_,mu,z);
    if(fabs(A1)<EPS){
      continue;
    }
    for(int m__=-lv;m__<=lv;m__++){
      double AA = A1 * A(lv,m__,mv,z);
      if(fabs(AA)<EPS){
        continue;
      }
      for(int m1_=-lu1;m1_<=lu1;m1_++){
        double AAA = AA * A(lu1,m1_,mu1,z);
        if(fabs(AAA)<EPS){
          continue;
        }
        for(int m1__=-lv1;m1__<=lv1;m1__++){
          double AAAA = AAA * A(lv1,m1__,mv1,z);
          if(fabs(AAAA)<EPS){
            continue;
          }
          double BBG  = 0.0;
          for(int l=abs(lu-lv); l<=lu+lv; l++){
            for(int l1=abs(lu1-lv1); l1<=lu1+lv1; l1++){
              int lm = MIN(l,l1);
              for(int m=-lm; m<=lm; m++){
                double B1 = B(l,  lu,  lv,  m, m_,  m__);
                if(fabs(B1)<EPS){
                  continue;
                }
                double B2 = B(l1, lu1, lv1, m, m1_, m1__);
                if(fabs(B2)<EPS){
                  continue;
                }
                double G  = G_eq52_mmmm(m,l,l1,lu,lv,lu1,lv1,qu,qu1,r,qmd);
                BBG += B1*B2*G;
              }
            }
          }
          ret += AAAA*BBG;
        }
      }
    }
  }
  return ret;
}

double R0_eq39_mmmp_old(int mu, int mv, int mu1, int mv1, int lu, int lv, int lu1, int lv1, int qu, int qu1, double ruu1[3], qmdata * qmd){
  double z[3];
  double r = sqrt(r3dot(ruu1, ruu1));  //eq29
  r3cpsc(z, ruu1, 1.0/r);              //eq29
  double ret = 0.0;
  for(int m_=-lu;m_<=lu;m_++){
    double A1 = A(lu,m_,mu,z);
    if(fabs(A1)<EPS){
      continue;
    }
    for(int m__=-lv;m__<=lv;m__++){
      double AA = A1 * A(lv,m__,mv,z);
      if(fabs(AA)<EPS){
        continue;
      }
      for(int m1_=-lu1;m1_<=lu1;m1_++){
        double AAA = AA * A(lu1,m1_,mu1,z);
        if(fabs(AAA)<EPS){
          continue;
        }
        for(int m1__=-lv1;m1__<=lv1;m1__++){
          double AAAA = AAA * A(lv1,m1__,mv1,z);
          if(fabs(AAAA)<EPS){
            continue;
          }
          double BBG  = 0.0;
          for(int l=abs(lu-lv); l<=lu+lv; l++){
            double q1;
            if( ! qlll_pm(qu,  lu,  lv,  l,  &q1, qmd)){
              continue;
            }
            for(int l1=abs(lu1-lv1); l1<=lu1+lv1; l1++){
              double q2;
              if( ! qlll_mm(qu1, lu1, lv1, l1, &q2, qmd)){
                continue;
              }
              int lm = MIN(l,l1);
              for(int m=-lm; m<=lm; m++){
                double B1 = B(l,  lu,  lv,  m, m_,  m__);
                if(fabs(B1)<EPS){
                  continue;
                }
                double BB = B1 * B(l1, lu1, lv1, m, m1_, m1__);
                if(fabs(BB)<EPS){
                  continue;
                }
                double G  = G_eq52_mmmp(m,l,l1,lu,lv,lu1,lv1,qu,qu1,q1,q2,r,qmd);
                BBG += BB*G;
              }
            }
          }
          ret += AAAA*BBG;
        }
      }
    }
  }
  return ret;
}

double R0_eq39_mmmm(int u, int v, int u1, int v1, int qu, int qu1, euler * z, basis * bo, qmdata * qmd){
  int mu  = bo->m[u ];
  int lu  = bo->l[u ];
  int mv  = bo->m[v ];
  int lv  = bo->l[v ];
  int mu1 = bo->m[u1];
  int lu1 = bo->l[u1];
  int mv1 = bo->m[v1];
  int lv1 = bo->l[v1];
  double r = z->r;

  double ret = 0.0;
  for(int m_=-lu;m_<=lu;m_++){
    double A1 = A_new(lu,m_,mu,z);
    if(fabs(A1)<EPS){
      continue;
    }
    for(int m__=-lv;m__<=lv;m__++){
      double AA = A1 * A_new(lv,m__,mv,z);
      if(fabs(AA)<EPS){
        continue;
      }
      for(int m1_=-lu1;m1_<=lu1;m1_++){
        double AAA = AA * A_new(lu1,m1_,mu1,z);
        if(fabs(AAA)<EPS){
          continue;
        }
        for(int m1__=-lv1;m1__<=lv1;m1__++){
          double AAAA = AAA * A_new(lv1,m1__,mv1,z);
          if(fabs(AAAA)<EPS){
            continue;
          }
          double BBG  = 0.0;
          for(int l=abs(lu-lv); l<=lu+lv; l++){
            for(int l1=abs(lu1-lv1); l1<=lu1+lv1; l1++){
              int lm = MIN(l,l1);
              for(int m=-lm; m<=lm; m++){
                double B1 = B(l,  lu,  lv,  m, m_,  m__);
                if(fabs(B1)<EPS){
                  continue;
                }
                double B2 = B(l1, lu1, lv1, m, m1_, m1__);
                if(fabs(B2)<EPS){
                  continue;
                }
                double G  = G_eq52_mmmm(m,l,l1,lu,lv,lu1,lv1,qu,qu1,r,qmd);
                BBG += B1*B2*G;
              }
            }
          }
          ret += AAAA*BBG;
        }
      }
    }
  }
  return ret;
}

double R0_eq39_mmmp(int u, int v, int u1, int v1, int qu, int qu1, euler * z, basis * bo, basis * bv, qmdata * qmd){

  int mu  = bv->m[u ];
  int lu  = bv->l[u ];
  int mv  = bo->m[v ];
  int lv  = bo->l[v ];
  int mu1 = bo->m[u1];
  int lu1 = bo->l[u1];
  int mv1 = bo->m[v1];
  int lv1 = bo->l[v1];
  double r = z->r;

  double ret = 0.0;
  for(int m_=-lu;m_<=lu;m_++){
    double A1 = A_new(lu,m_,mu,z);
    if(fabs(A1)<EPS){
      continue;
    }
    for(int m__=-lv;m__<=lv;m__++){
      double AA = A1 * A_new(lv,m__,mv,z);
      if(fabs(AA)<EPS){
        continue;
      }
      for(int m1_=-lu1;m1_<=lu1;m1_++){
        double AAA = AA * A_new(lu1,m1_,mu1,z);
        if(fabs(AAA)<EPS){
          continue;
        }
        for(int m1__=-lv1;m1__<=lv1;m1__++){
          double AAAA = AAA * A_new(lv1,m1__,mv1,z);
          if(fabs(AAAA)<EPS){
            continue;
          }
          double BBG  = 0.0;
          for(int l=abs(lu-lv); l<=lu+lv; l++){
            double q1;
            if( ! qlll_pm(qu,  lu,  lv,  l,  &q1, qmd)){
              continue;
            }
            for(int l1=abs(lu1-lv1); l1<=lu1+lv1; l1++){
              double q2;
              if( ! qlll_mm(qu1, lu1, lv1, l1, &q2, qmd)){
                continue;
              }
              int lm = MIN(l,l1);
              for(int m=-lm; m<=lm; m++){
                double B1 = B(l,  lu,  lv,  m, m_,  m__);
                if(fabs(B1)<EPS){
                  continue;
                }
                double BB = B1 * B(l1, lu1, lv1, m, m1_, m1__);
                if(fabs(BB)<EPS){
                  continue;
                }
                double G  = G_eq52_mmmp(m,l,l1,lu,lv,lu1,lv1,qu,qu1,q1,q2,r,qmd);
                BBG += BB*G;
              }
            }
          }
          ret += AAAA*BBG;
        }
      }
    }
  }
  return ret;
}

