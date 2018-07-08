#include "qm.h"
#include "eq.h"
#include "matrix.h"

static double Ddiff(int M, double * Da, double * Db, double * oldD){
  double dD = 0.0;
  for(int u=0; u<M; u++){
    for(int v=u; v<M; v++){
      int uv = mpos(u,v);
      double ab = Da[uv]+Db[uv];
      double d  = oldD[uv] - ab;
      dD += d*d*(u==v?1:2);
      oldD[uv] = ab;
    }
  }
  return sqrt(dD)/M;
}

static void diagF(int M, int use_old_c, double * F, double * C, double * V){
  if(use_old_c){
    /* C1 -- coefficients from the previous step
       transform       F' = C1^T * F * C1
       solve           F'B  = BV
       transform back  C^T = B^T C1^T (simultaneously with solving)
     */
    mx_BHBt_sym(M, F, C);
  }
  else{
    mx_id(M, C);
  }
  jacobi(F, C, V, M, 1e-15, 20, NULL);
  eigensort(M, V, C);
  return;
}

static inline void errormx(int M, double * R, double * F, double * D){
  /* R := FD-DF */
  mx_symmultsymmx(M, R, F, D);
  mx_antisym(M, R);
  return;
}

static void iter_print(int k, double E, double dE, double dD, FILE * fo){
  if(!k){
    fprintf(fo, " it %3d     E = % 17.10lf\n", k+1, E);
  }
  else{
    fprintf(fo, " it %3d     E = % 17.10lf    dE = % 17.10lf    dD = % 5.2e\n", k+1, E, dE, dD);
  }
  return;
}

static void iter_diis_print(int k, double E, double dE, double dD, double Rmax, FILE * fo){
  if(!k){
    fprintf(fo, " it %3d     E = % 17.10lf                                                Rmax = % 5.2e\n", k+1, E, Rmax);
  }
  else{
    fprintf(fo, " it %3d     E = % 17.10lf    dE = % 17.10lf    dD = % 5.2e    Rmax = % 5.2e\n", k+1, E, dE, dD, Rmax);
  }
  return;
}

static void E012_print(double E, double E0, double E1, double E2, FILE * fo){
  fprintf(fo, "\n");
  fprintf(fo, " (E0   = %20.10lf)\n", E0);
  fprintf(fo, " (E0+1 = %20.10lf)\n", E0+E1);
  fprintf(fo, " (E2   = %20.10lf)\n", E2);
  fprintf(fo, "  E    = %20.10lf\n",  E);
  return;
}

double scf(int Na, int Nb, double E0,
    double * Ca, double * Cb, double * Va, double * Vb,
    double * Da, double * Db, double * Dmp,
    int maxit, double dDmax, int * alo, int * alv,
    double * H, double * Hmp, double * mmmm, double * pmmm,
    basis * bo, basis * bv, mol * m, qmdata * qmd, FILE * fo){

  int Mo = bo->M;
  int Mv = bv->M;

  double * Fa    = malloc(sizeof(double)*symsize(Mo));
  double * Fb    = malloc(sizeof(double)*symsize(Mo));
  double * oldD  = malloc(sizeof(double)*symsize(Mo));
  double * Fw    = malloc(sizeof(double)*symsize(Mo));
  double * Fmpa  = malloc(sizeof(double)*Mo*Mv);
  double * Fmpb  = malloc(sizeof(double)*Mo*Mv);
  double * Xa    = malloc(sizeof(double)*Mo*Mv);
  double * Xb    = malloc(sizeof(double)*Mo*Mv);
  double * FaXa  = malloc(sizeof(double)*Mo*Mo);
  double * FbXb  = malloc(sizeof(double)*Mo*Mo);
  double * sa    = malloc(sizeof(double)*Mo);
  double * sb    = malloc(sizeof(double)*Mo);
  double * F2a   = malloc(sizeof(double)*symsize(Mo));
  double * F2b   = malloc(sizeof(double)*symsize(Mo));
  double * FA    = malloc(sizeof(double)*symsize(Mo));
  double * FB    = malloc(sizeof(double)*symsize(Mo));
  double * dEdFa = malloc(sizeof(double)*Mo*Mv);
  double * dEdFb = malloc(sizeof(double)*Mo*Mv);

  double E1 = 0.0;
  double E2 = 0.0;
  double E  = E0+E1+E2;
  vecset(symsize(Mo), oldD, 0.0);
  int k = 0;

  while(k < maxit){
    double oldE = E;
    D_eq9 (Na, Mo, Ca, Da);
    D_eq9 (Nb, Mo, Cb, Db);
    F_eq4 (Da, Db, H, Fa, Fb, alo, mmmm, bo, m, qmd);
    F2_8_7_14_15_6(Da, Db, Hmp, Fmpa, Fmpb, Xa, Xb, FaXa, FbXb, sa, sb, F2a, F2b, alo, alv, pmmm, bo, bv, m, qmd);
    dEdF(Da, Db, Xa, Xb, FaXa, FbXb, sa, sb, Fmpa, Fmpb, dEdFa, dEdFb, alo, bo, bv, qmd);
    Heff(dEdFa, dEdFb, Fa, Fb, F2a, F2b, FA, FB, alo, alv, pmmm, bo, bv, m, qmd);

    E1 = E1_eq3(Mo, H, Da, Db, Fa, Fb);
    E2 = E2_eq5(Mo, Da, Db, F2a, F2b);
    E  = E0+E1+E2;

    double dD = Ddiff(Mo, Da, Db, oldD);
    double dE = E-oldE;
    if(fo) iter_print(k, E, dE, dD, fo);
    if(dD < dDmax){
      if(fo) fprintf(fo, "converged\n");
      break;
    }

    veccp(symsize(Mo), Fw, FA);
    diagF(Mo, 1, Fw, Ca, Va);
    veccp(symsize(Mo), Fw, FB);
    diagF(Mo, 1, Fw, Cb, Vb);

    k++;
  }

  if(fo) E012_print(E, E0, E1, E2, fo);
  vecsum(Mo*Mv, Dmp, dEdFa, dEdFb);

  free(Fa);
  free(Fb);
  free(oldD);
  free(Fw);
  free(Fmpa);
  free(Fmpb);
  free(Xa);
  free(Xb);
  free(FaXa);
  free(FbXb);
  free(sa);
  free(sb);
  free(F2a);
  free(F2b);
  free(FA);
  free(FB);
  free(dEdFa);
  free(dEdFb);

  return E;
}

double scf_r(int N, double E0,
    double * C, double * V, double * D, double * Dmp,
    int maxit, double dDmax, int * alo, int * alv,
    double * H, double * Hmp, double * mmmm, double * pmmm,
    basis * bo, basis * bv, mol * m, qmdata * qmd, FILE * fo){

  int Mo = bo->M;
  int Mv = bv->M;

  double * F    = malloc(sizeof(double)*symsize(Mo));
  double * oldD = malloc(sizeof(double)*symsize(Mo));
  double * Fw   = malloc(sizeof(double)*symsize(Mo));
  double * Fmp  = malloc(sizeof(double)*Mo*Mv);
  double * X    = malloc(sizeof(double)*Mo*Mv);
  double * FX   = malloc(sizeof(double)*Mo*Mo);
  double * s    = malloc(sizeof(double)*Mo);
  double * F2   = malloc(sizeof(double)*symsize(Mo));
  double * Feff = malloc(sizeof(double)*symsize(Mo));

  double E1 = 0.0;
  double E2 = 0.0;
  double E  = E0+E1+E2;
  vecset(symsize(Mo), oldD, 0.0);
  int k = 0;

  while(k < maxit){
    double oldE = E;
    D_eq9(N/2, Mo, C, D);
    F_eq4_r(D, H, F, alo, mmmm, bo, m, qmd);
    F2_8_7_14_15_6_r(D, Hmp, Fmp, X, FX, s, F2, alo, alv, pmmm, bo, bv, m, qmd);
    dEdF_r(D, X, FX, s, Fmp, Dmp, alo, bo, bv, qmd);
    Heff_r(Dmp, F, F2, Feff, alo, alv, pmmm, bo, bv, m, qmd);

    E1 = E1_eq3_r(Mo, H, D, F);
    E2 = E2_eq5_r(Mo, D, F2);
    E  = E0+E1+E2;

    double dD = Ddiff(Mo, D, D, oldD);
    double dE = E-oldE;
    if(fo) iter_print(k, E, dE, dD, fo);
    if(dD < dDmax){
      if(fo) fprintf(fo, "converged\n");
      break;
    }

    veccp(symsize(Mo), Fw, Feff);
    diagF(Mo, 1, Fw, C, V);
    k++;
  }

  if(fo) E012_print(E, E0, E1, E2, fo);

  free(oldD);
  free(s);
  free(F);
  free(Fmp);
  free(F2);
  free(Feff);
  free(Fw);
  free(FX);
  free(X);
  return E;
}

/*---------------------------------------------------------------------------*/

static inline void addB(int k, double * B, int M, double ** allR){
  for(int i=0; i<=k; i++){
    B[mpos(i,k)] += mx_multtrmx_tr(M, allR[i], allR[k]);
  }
  return;
}

static void makeb(int k, int l, double * b, double * B){
  for(int i=0; i<k-l; i++){
    b[i*(k+1-l)+i] = B[mpos(i+l,i+l)];
    for(int j=i+1; j<k-l; j++){
      b[i*(k+1-l)+j] = b[j*(k+1-l)+i] = B[mpos(i+l,j+l)];
    }
  }
  return;
}

static int bols(int k, int l, double * cf, double * b, double * B){
  if(k-l<1){
    GOTOHELL;
  }
  vecset((k+1)*(k+1), b, -1.0);
  makeb(k, l, b, B);
  b[(k-l)*(k+1-l)+k-l] = 0.0;
  vecset(k+1, cf, 0.0);
  cf[k-l] = -1.0;
  return mx_inv(k+1-l, 1, cf, b, 1e-15);
}

static int fcoef(int k, int l, double * cf, double * B){
  double * b = malloc(sizeof(double)*(k+1)*(k+1));
  int p;
  do{
    p = bols(k, l++, cf, b, B);
  } while(p);
  free(b);
  return l-1;
}

static inline void lincomb(int n, int m, double * p, double ** v, double * c){
  vecset(n, p, 0.0);
  for(int i=0; i<m; i++){
    vecadds(n, p, v[i], c[i]);
  }
  return;
}

double scf_diis(int Na, int Nb, double E0,
    double * Ca, double * Cb, double * Va, double * Vb,
    double * Da, double * Db, double * Dmp,
    int maxit, int memit, double dDmax, int * alo, int * alv,
    double * H, double * Hmp, double * mmmm, double * pmmm,
    basis * bo, basis * bv, mol * m, qmdata * qmd, FILE * fo){

  int Mo = bo->M;
  int Mv = bv->M;

  double * Fa    = malloc(sizeof(double)*symsize(Mo));
  double * Fb    = malloc(sizeof(double)*symsize(Mo));
  double * oldD  = malloc(sizeof(double)*symsize(Mo));
  double * Fmpa  = malloc(sizeof(double)*Mo*Mv);
  double * Fmpb  = malloc(sizeof(double)*Mo*Mv);
  double * Xa    = malloc(sizeof(double)*Mo*Mv);
  double * Xb    = malloc(sizeof(double)*Mo*Mv);
  double * FaXa  = malloc(sizeof(double)*Mo*Mo);
  double * FbXb  = malloc(sizeof(double)*Mo*Mo);
  double * sa    = malloc(sizeof(double)*Mo);
  double * sb    = malloc(sizeof(double)*Mo);
  double * F2a   = malloc(sizeof(double)*symsize(Mo));
  double * F2b   = malloc(sizeof(double)*symsize(Mo));
  double * FA    = malloc(sizeof(double)*symsize(Mo));
  double * FB    = malloc(sizeof(double)*symsize(Mo));
  double * dEdFa = malloc(sizeof(double)*Mo*Mv);
  double * dEdFb = malloc(sizeof(double)*Mo*Mv);

  double * B  = calloc(sizeof(double)*(symsize(maxit+1)+(maxit+1)),1);
  double * cf = B + symsize(maxit+1);
  double ** allFa = malloc(4*maxit*sizeof(double *));
  double ** allFb = allFa + maxit;
  double ** allRa = allFb + maxit;
  double ** allRb = allRa + maxit;
  int K1 = MIN(maxit,memit);
  double * Fs = malloc(sizeof(double)*K1*2*symsize(Mo));
  double * Rs = malloc(sizeof(double)*K1*2*Mo*Mo);
  for(int i=0; i<K1; i++){
    allFa[i] = Fs + (2*i  ) * symsize(Mo);
    allFb[i] = Fs + (2*i+1) * symsize(Mo);
    allRa[i] = Rs + (2*i  ) * Mo*Mo;
    allRb[i] = Rs + (2*i+1) * Mo*Mo;
  }

  double E1 = 0.0;
  double E2 = 0.0;
  double E  = E0+E1+E2;
  vecset(symsize(Mo), oldD, 0.0);
  int k = 0, l = 0;

  while(k < maxit){

    double oldE = E;
    D_eq9 (Na, Mo, Ca, Da);
    D_eq9 (Nb, Mo, Cb, Db);

    if(k>=K1){ // out of memory
      allFa[k] = allFa[k-K1];
      allFb[k] = allFb[k-K1];
      allRa[k] = allRa[k-K1];
      allRb[k] = allRb[k-K1];
      if(l==k-K1){
        l++;
      }
    }

    F_eq4 (Da, Db, H, Fa, Fb, alo, mmmm, bo, m, qmd);
    F2_8_7_14_15_6(Da, Db, Hmp, Fmpa, Fmpb, Xa, Xb, FaXa, FbXb, sa, sb, F2a, F2b, alo, alv, pmmm, bo, bv, m, qmd);
    dEdF(Da, Db, Xa, Xb, FaXa, FbXb, sa, sb, Fmpa, Fmpb, dEdFa, dEdFb, alo, bo, bv, qmd);
    Heff(dEdFa, dEdFb, Fa, Fb, F2a, F2b, allFa[k], allFb[k], alo, alv, pmmm, bo, bv, m, qmd);

    errormx(Mo, allRa[k], allFa[k], Da);
    errormx(Mo, allRb[k], allFb[k], Db);
    addB(k, B, Mo, allRa);
    addB(k, B, Mo, allRb);
    l = fcoef(k+1, l, cf, B);
    lincomb(symsize(Mo), k+1-l, FA, allFa+l, cf);
    lincomb(symsize(Mo), k+1-l, FB, allFb+l, cf);

    E1 = E1_eq3(Mo, H, Da, Db, Fa, Fb);
    E2 = E2_eq5(Mo, Da, Db, F2a, F2b);
    E  = E0+E1+E2;

    double dD = Ddiff(Mo, Da, Db, oldD);
    double dE = E-oldE;
    double Rmax = vecabsmax(2*Mo*Mo, allRa[k]); // Rb[k] is right after Ra[k]
    if(fo) iter_diis_print(k, E, dE, dD, Rmax, fo);
    if(dD < dDmax){
      if(fo) fprintf(fo, "converged\n");
      break;
    }

    diagF(Mo, 1, FA, Ca, Va);
    diagF(Mo, 1, FB, Cb, Vb);

    k++;
  }

  if(fo) E012_print(E, E0, E1, E2, fo);
  vecsum(Mo*Mv, Dmp, dEdFa, dEdFb);

  free(Fs);
  free(Rs);
  free(allFa);
  free(B);

  free(Fa);
  free(Fb);
  free(oldD);
  free(Fmpa);
  free(Fmpb);
  free(Xa);
  free(Xb);
  free(FaXa);
  free(FbXb);
  free(sa);
  free(sb);
  free(F2a);
  free(F2b);
  free(FA);
  free(FB);
  free(dEdFa);
  free(dEdFb);
  return E;
}

double scf_diis_r(int N, double E0,
    double * C, double * V, double * D, double * Dmp,
    int maxit, int memit, double dDmax, int * alo, int * alv,
    double * H, double * Hmp, double * mmmm, double * pmmm,
    basis * bo, basis * bv, mol * m, qmdata * qmd, FILE * fo){

  int Mo = bo->M;
  int Mv = bv->M;

  double * F    = malloc(sizeof(double)*symsize(Mo));
  double * oldD = malloc(sizeof(double)*symsize(Mo));
  double * Fmp  = malloc(sizeof(double)*Mo*Mv);
  double * X    = malloc(sizeof(double)*Mo*Mv);
  double * FX   = malloc(sizeof(double)*Mo*Mo);
  double * s    = malloc(sizeof(double)*Mo);
  double * F2   = malloc(sizeof(double)*symsize(Mo));
  double * Feff = malloc(sizeof(double)*symsize(Mo));

  double * B  = calloc(sizeof(double)*(symsize(maxit+1)+(maxit+1)),1);
  double * cf = B + symsize(maxit+1);
  double ** allF = malloc(2*maxit*sizeof(double *));
  double ** allR = allF + maxit;
  int K1 = MIN(maxit,memit);
  double * Fs = malloc(sizeof(double)*K1*symsize(Mo));
  double * Rs = malloc(sizeof(double)*K1*Mo*Mo);
  for(int i=0; i<K1; i++){
    allF[i] = Fs + i * symsize(Mo);
    allR[i] = Rs + i * Mo*Mo;
  }

  double E1 = 0.0;
  double E2 = 0.0;
  double E  = E0+E1+E2;
  vecset(symsize(Mo), oldD, 0.0);
  int k = 0, l = 0;

  while(k < maxit){

    double oldE = E;
    D_eq9(N/2, Mo, C, D);

    if(k>=K1){
      allF[k] = allF[k-K1];
      allR[k] = allR[k-K1];
      if(l==k-K1){
        l++;
      }
    }

    F_eq4_r(D, H, F, alo, mmmm, bo, m, qmd);
    F2_8_7_14_15_6_r(D, Hmp, Fmp, X, FX, s, F2, alo, alv, pmmm, bo, bv, m, qmd);
    dEdF_r(D, X, FX, s, Fmp, Dmp, alo, bo, bv, qmd);
    Heff_r(Dmp, F, F2, allF[k], alo, alv, pmmm, bo, bv, m, qmd);

    errormx(Mo, allR[k], allF[k], D);
    addB(k, B, Mo, allR);
    l = fcoef(k+1, l, cf, B);
    lincomb(symsize(Mo), k+1-l, Feff, allF+l, cf);

    E1 = E1_eq3_r(Mo, H, D, F);
    E2 = E2_eq5_r(Mo, D, F2);
    E  = E0+E1+E2;

    double dD = Ddiff(Mo, D, D, oldD);
    double dE = E-oldE;
    double Rmax = vecabsmax(Mo*Mo, allR[k]);
    if(fo) iter_diis_print(k, E, dE, dD, Rmax, fo);
    if(dD < dDmax){
      if(fo) fprintf(fo, "converged\n");
      break;
    }

    diagF(Mo, 1, Feff, C, V);

    k++;
  }

  if(fo) E012_print(E, E0, E1, E2, fo);

  free(Fs);
  free(Rs);
  free(allF);
  free(B);

  free(oldD);
  free(s);
  free(F);
  free(Fmp);
  free(F2);
  free(Feff);
  free(FX);
  free(X);
  return E;
}

