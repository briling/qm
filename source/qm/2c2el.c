#include "eq.h"
#include "common.h"
#include "vec3.h"

#define EPS 1e-15

double R0_eq39_mmmm_old(int mu, int mv, int mu1, int mv1, int lu, int lv, int lu1, int lv1, int qu, int qu1, double ruu1[3], qmdata * qmd){
  double z[3];
  double r = sqrt(r3dot(ruu1, ruu1));  //eq29
  r3cpsc(z, ruu1, 1.0/r);              //eq29
  double ret = 0.0;
  for(int m_=-lu;m_<=lu;m_++){
    double A1 = A_full(lu,m_,mu,z);
    if(fabs(A1)<EPS){
      continue;
    }
    for(int m__=-lv;m__<=lv;m__++){
      double AA = A1 * A_full(lv,m__,mv,z);
      if(fabs(AA)<EPS){
        continue;
      }
      for(int m1_=-lu1;m1_<=lu1;m1_++){
        double AAA = AA * A_full(lu1,m1_,mu1,z);
        if(fabs(AAA)<EPS){
          continue;
        }
        for(int m1__=-lv1;m1__<=lv1;m1__++){
          double AAAA = AAA * A_full(lv1,m1__,mv1,z);
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
    double A1 = A_full(lu,m_,mu,z);
    if(fabs(A1)<EPS){
      continue;
    }
    for(int m__=-lv;m__<=lv;m__++){
      double AA = A1 * A_full(lv,m__,mv,z);
      if(fabs(AA)<EPS){
        continue;
      }
      for(int m1_=-lu1;m1_<=lu1;m1_++){
        double AAA = AA * A_full(lu1,m1_,mu1,z);
        if(fabs(AAA)<EPS){
          continue;
        }
        for(int m1__=-lv1;m1__<=lv1;m1__++){
          double AAAA = AAA * A_full(lv1,m1__,mv1,z);
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

double R0_eq39_mmmm(int u, int v, int u1, int v1, axis * z, basis * bo, qmdata * qmd){
  int qu  = bo->Q[u ];
  int qu1 = bo->Q[u1];
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

double R0_eq39_mmmp(int a, int v, int u1, int v1, axis * z, basis * bo, basis * bv, qmdata * qmd){

  int qa  = bv->Q[a ];
  int qu1 = bo->Q[u1];
  int ma  = bv->m[a ];
  int la  = bv->l[a ];
  int mv  = bo->m[v ];
  int lv  = bo->l[v ];
  int mu1 = bo->m[u1];
  int lu1 = bo->l[u1];
  int mv1 = bo->m[v1];
  int lv1 = bo->l[v1];
  double r = z->r;

  double ret = 0.0;
  for(int m_=-la;m_<=la;m_++){
    double A1 = A(la,m_,ma,z);
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
          for(int l=abs(la-lv); l<=la+lv; l++){
            double q1;
            if( ! qlll_pm(qa,  la,  lv,  l,  &q1, qmd)){
              continue;
            }
            for(int l1=abs(lu1-lv1); l1<=lu1+lv1; l1++){
              double q2;
              if( ! qlll_mm(qu1, lu1, lv1, l1, &q2, qmd)){
                continue;
              }
              int lm = MIN(l,l1);
              for(int m=-lm; m<=lm; m++){
                double B1 = B(l,  la,  lv,  m, m_,  m__);
                if(fabs(B1)<EPS){
                  continue;
                }
                double BB = B1 * B(l1, lu1, lv1, m, m1_, m1__);
                if(fabs(BB)<EPS){
                  continue;
                }
                double G  = G_eq52_mmmp(m,l,l1,la,lv,lu1,lv1,qa,qu1,q1,q2,r,qmd);
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

