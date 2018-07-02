#include "qm.h"
#include "eq.h"
#include "matrix.h"

static double Ddiff(int M, double * Da, double * Db, double * oldD){
  double dD = 0.0;
  for(int i=0; i<symsize(M); i++){
    double ab = Da[i]+Db[i];
    double d  = oldD[i] - ab;
    dD += d*d;
    oldD[i] = ab;
  }
  return dD;
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
    if(fo){
      if(!k){
        fprintf(fo, " it %3d     E = % 17.10lf\n", k+1, E);
      }
      else{
        fprintf(fo, " it %3d     E = % 17.10lf    dE = % 17.10lf    dD = % 5.2e\n", k+1, E, dE, dD);
      }
    }
    if(dD < dDmax){
      if(fo) fprintf(fo, "converged\n");
      break;
    }

    mx_id(Mo, Ca);
    veccp(symsize(Mo), Fw, FA);
    jacobi(Fw, Ca, Va, Mo, 1e-15, 20, NULL);
    eigensort(Mo, Va, Ca);
    mx_id(Mo, Cb);
    veccp(symsize(Mo), Fw, FB);
    jacobi(Fw, Cb, Vb, Mo, 1e-15, 20, NULL);
    eigensort(Mo, Vb, Cb);
    k++;
  }

  if(fo){
    fprintf(fo, "\n");
    fprintf(fo, " (E0   = %20.10lf)\n", E0);
    fprintf(fo, " (E0+1 = %20.10lf)\n", E0+E1);
    fprintf(fo, " (E2   = %20.10lf)\n", E2);
    fprintf(fo, "  E    = %20.10lf\n",  E);
  }

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

/*---------------------------------------------------------------------------*/

static inline void errormx(int M, double * R, double * F, double * D){
  /* R := FD-DF */
  mx_symmultsymmx(M, R, F, D);
  mx_antisym(M, R);
  return;
}

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

    if(k>=K1){
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
    if(fo){
      if(!k){
        fprintf(fo, " it %3d     E = % 17.10lf\n", k+1, E);
      }
      else{
        fprintf(fo, " it %3d     E = % 17.10lf    dE = % 17.10lf    dD = % 5.2e\n", k+1, E, dE, dD);
      }
    }
    if(dD < dDmax){
      if(fo) fprintf(fo, "converged\n");
      break;
    }

    mx_id(Mo, Ca);
    jacobi(FA, Ca, Va, Mo, 1e-15, 20, NULL);
    eigensort(Mo, Va, Ca);
    mx_id(Mo, Cb);
    jacobi(FB, Cb, Vb, Mo, 1e-15, 20, NULL);
    eigensort(Mo, Vb, Cb);
    k++;
  }

  if(fo){
    fprintf(fo, "\n");
    fprintf(fo, " (E0   = %20.10lf)\n", E0);
    fprintf(fo, " (E0+1 = %20.10lf)\n", E0+E1);
    fprintf(fo, " (E2   = %20.10lf)\n", E2);
    fprintf(fo, "  E    = %20.10lf\n",  E);
  }

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

