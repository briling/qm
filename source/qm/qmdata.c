#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "matrix.h"
#include "qm.h"

#define END(X) ( qmd->X + (X##_size)/sizeof(*(qmd->X)) )
#define CA_L_PRINT(F, Q, PAR, L) ca_l_print(Q, qmd->PAR, #PAR, qmd->L, qmd->n ## L, F);
#define D_L_PRINT(F, Q, PAR, L) d_l_print(Q, qmd->PAR, #PAR, qmd->L, qmd->n ## L, F);

static qmdata * qmdata_alloc(int Qmax,
                             int nLo, int nLv,
                             int nLp, int Nga,
                             int Nqa, int Nq1a,
                             int Nfb, int Nf1b,
                             int Nfo,
                             int Nub, int Nu1b,
                             int Nvb){

  int N = Qmax + 1;

  size_t used_size = N * sizeof(int);
  size_t Lo_size   = N * sizeof(int);
  size_t Lv_size   = N * sizeof(int);
  size_t Lp_size   = N * sizeof(int);
  size_t p_size    = N * sizeof(int) * nLo;

  size_t ea_size      = N * sizeof(double)        ;
  size_t fa_size      = N * sizeof(double) * nLo  ;
  size_t f2a_size     = N * sizeof(double) * nLv  ;
  size_t r0_size      = N * sizeof(double)        ;
  size_t cf_size      = N * sizeof(double) * nLp  ;
  size_t bub_size     = N * sizeof(ca_t)          ;
  size_t bu1b_size    = N * sizeof(ca_t)          ;
  size_t aub_size     = N * sizeof(ca_t)   * nLo  ;
  size_t au0b_size    = N * sizeof(ca_t)   * nLo  ;
  size_t au1b_size    = N * sizeof(ca_t)   * nLv  ;
  size_t afb_size     = N * sizeof(ca_t)   * nLo  ;
  size_t af0b_size    = N * sizeof(ca_t)   * nLo  ;
  size_t af1b_size    = N * sizeof(ca_t)   * nLv  ;
  size_t afo_size     = N * sizeof(ca_t)   * nLo  ;
  size_t aof_size     = N * sizeof(ca_t)   * nLp  ;
  size_t agb_size     = N * sizeof(double) * nLo  ;
  size_t ag0b_size    = N * sizeof(double) * nLo  ;
  size_t ag1b_size    = N * sizeof(double) * nLv  ;
  size_t bvb_size     = N * sizeof(ca_t)          ;
  size_t avb_size     = N * sizeof(ca_t)   * nLo  ;
  size_t ag6_size     = N * sizeof(ca_t)   * nLo  ;
  size_t q_list_size  = N * sizeof(ql_t)          ;
  size_t ga_size      = Nga  * sizeof(ga_t)       ;
  size_t qa_size      = Nqa  * sizeof(qa_t)       ;
  size_t q1a_size     = Nq1a * sizeof(qa_t)       ;
  size_t qq_list_size = symsize(N) * sizeof(qql_t);
  size_t fb_size      = Nfb   * sizeof(cabm_t)    ;
  size_t f1b_size     = Nf1b  * sizeof(cabm_t)    ;
  size_t fo_size      = Nfo   * sizeof(cabm_t)    ;
  size_t ub_size      = Nub   * sizeof(cabl_t)    ;
  size_t u1b_size     = Nu1b  * sizeof(cabl_t)    ;
  size_t vb_size      = Nvb   * sizeof(cablm_t)   ;

  size_t size = sizeof(qmdata)  + p_size +
    used_size    + Lo_size      + Lv_size      + Lp_size      +
    ea_size      + fa_size      + f2a_size     + r0_size      +
    cf_size      + bub_size     + bu1b_size    + aub_size     +
    au0b_size    + au1b_size    + afb_size     + af0b_size    +
    af1b_size    + afo_size     + aof_size     + agb_size     +
    ag0b_size    + ag1b_size    + bvb_size     + avb_size     +
    ag6_size     + q_list_size  + ga_size      + qa_size      +
    q1a_size     + qq_list_size + fb_size      + f1b_size     +
    fo_size      + ub_size      + u1b_size     + vb_size      ;

  qmdata * qmd  = calloc(1, size);
  if(!qmd){
    return NULL;
  }
  qmd->Qmax    = Qmax;
  qmd->nLo     = nLo;
  qmd->nLv     = nLv;
  qmd->nLp     = nLp;

  qmd->ea      = (double *)(qmd + 1);
  qmd->fa      = END(ea  );
  qmd->f2a     = END(fa  );
  qmd->r0      = END(f2a );
  qmd->cf      = END(r0  );
  qmd->agb     = END(cf  );
  qmd->ag0b    = END(agb );
  qmd->ag1b    = END(ag0b);

  qmd->bub     = (ca_t *)END(ag1b);
  qmd->bu1b    = END(bub );
  qmd->aub     = END(bu1b);
  qmd->au0b    = END(aub );
  qmd->au1b    = END(au0b);
  qmd->afb     = END(au1b);
  qmd->af0b    = END(afb );
  qmd->af1b    = END(af0b);
  qmd->afo     = END(af1b);
  qmd->aof     = END(afo );
  qmd->bvb     = END(aof );
  qmd->avb     = END(bvb );
  qmd->ag6     = END(avb );
  qmd->fb      = (cabm_t *)END(ag6);
  qmd->f1b     = END(fb );
  qmd->fo      = END(f1b);
  qmd->ub      = (cabl_t *)END(fo);
  qmd->u1b     = END(ub);
  qmd->vb      = (cablm_t *)END(u1b);
  qmd->ga      = (ga_t *)END(vb);
  qmd->qa      = (qa_t *)END(ga);
  qmd->q1a     = END(qa);

  qmd->qq_list = (qql_t *)END(q1a);
  qmd->q_list  = (ql_t *)END(qq_list);
  qmd->used    = (int *)END(q_list);
  qmd->Lo      = END(used);
  qmd->Lv      = END(Lo);
  qmd->Lp      = END(Lv);
  qmd->p       = END(Lp);

  return qmd;
}

qmdata * qmdata_read(FILE * f){
  char   s[256];
  char * s1;
  const char type_of_model[] = "m 3 K 0 fo 0\n";
  qmdata * qmd;
  int Q, Lo,Lv,Lp, nLo,nLv,nLp;

  do{
    while(fscanf(f, " %255[^$ ]", s) == 1);
    if(fscanf(f, "%255[$A-z]", s) != 1){
      return NULL;
    }
  } while(strcmp(s, "$qmdata"));

  if(!fgets (s, sizeof(s), f))  return NULL;
  if(!fgets (s, sizeof(s), f))  return NULL;
  s1 = s;
  while(*s1 == ' ') s1++;
  if(strcmp(s1, type_of_model)) return NULL;

  long p     = ftell(f);
  int  Qmax  = 0;
  int  LO    = 0;
  int  LV    = 0;
  int  LP    = 0;
  int  Nga   = 0;
  int  Nqa   = 0;
  int  Nq1a  = 0;
  int  Nfb   = 0;
  int  Nf1b  = 0;
  int  Nfo   = 0;
  int  Nub   = 0;
  int  Nu1b  = 0;
  int  Ngb   = 0;
  int  Ng10b = 0;
  int  Nvb   = 0;

  while ((fgets(s, sizeof(s), f)) && (!strchr(s, '$'))){
    int i1,i2,i3,i4,i5,i6,i7,i8,i9;
    double t1,t2,t3,t4,t5;
    if(sscanf(s, " n %d%d%d%d", &Q, &Lo, &Lv, &Lp) == 4){
      if(Q>Qmax) Qmax = Q;
      if(Lo>LO)  LO   = Lo;
      if(Lv>LV)  LV   = Lv;
      if(Lp>LP)  LP   = Lp;
    }
    else if(sscanf(s, " ga %d%d%d%d%d%d%lf", &i1, &i2, &i3, &i4, &i5, &i6, &t1) == 7){
      Nga++;
    }
    else if(sscanf(s, " qa %d%d%d%d%lf", &i1, &i2, &i3, &i4, &t1) == 5){
      Nqa++;
    }
    else if(sscanf(s, " q1a %d%d%d%d%lf", &i1, &i2, &i3, &i4, &t1) == 5){
      Nq1a++;
    }
    else if(sscanf(s, " fb %d%d%d%d%d a** %lf%lf%lf%lf", &i1, &i2, &i3, &i4, &i5, &t1, &t2, &t3, &t4) == 9){
      Nfb++;
    }
    else if(sscanf(s, " f1b %d%d%d%d%d b* %lf%lf%lf%lf%lf", &i1, &i2, &i3, &i4, &i5, &t1, &t2, &t3, &t4, &t5) == 10){
      Nf1b++;
    }
    else if(sscanf(s, " fo %d%d%d%d%d b* %lf%lf%lf%lf%lf", &i1, &i2, &i3, &i4, &i5, &t1, &t2, &t3, &t4, &t5) == 10){
      Nfo++;
    }
    else if(sscanf(s, " ub %d%d%d%d%d b* %lf%lf%lf%lf%lf", &i1, &i2, &i3, &i4, &i5, &t1, &t2, &t3, &t4, &t5) == 10){
      Nub++;
    }
    else if(sscanf(s, " u1b %d%d%d%d%d b* %lf%lf%lf%lf%lf", &i1, &i2, &i3, &i4, &i5, &t1, &t2, &t3, &t4, &t5) == 10){
      Nu1b++;
    }
    else if(sscanf(s, " gb %d%d%d%d%d%d%d%d%d a* %lf%lf%lf%lf", &i1, &i2, &i3, &i4, &i5, &i6, &i7, &i8, &i9, &t1, &t2, &t3, &t4) == 13){
      Ngb++;
    }
    else if(sscanf(s, " g10b %d%d%d%d%d%d%d%d%d a* %lf%lf%lf%lf", &i1, &i2, &i3, &i4, &i5, &i6, &i7, &i8, &i9, &t1, &t2, &t3, &t4) == 13){
      Ng10b++;
    }
    else if(sscanf(s, " vb %d%d%d%d%d%d%d%d%d b* %lf%lf%lf%lf%lf", &i1, &i2, &i3, &i4, &i5, &i6, &i7, &i8, &i9, &t1, &t2, &t3, &t4, &t5) == 14){
      Nvb++;
    }
  }
#if 0
  printf("Qmax  = %d\n", Qmax );
  printf("LO    = %d\n", LO   );
  printf("LV    = %d\n", LV   );
  printf("LP    = %d\n", LP   );
  printf("Nga   = %d\n", Nga  );
  printf("Nqa   = %d\n", Nqa  );
  printf("Nq1a  = %d\n", Nq1a );
  printf("Nfb   = %d\n", Nfb  );
  printf("Nf1b  = %d\n", Nf1b );
  printf("Nfo   = %d\n", Nfo  );
  printf("Nub   = %d\n", Nub  );
  printf("Nu1b  = %d\n", Nu1b );
  printf("Ngb   = %d\n", Ngb  );
  printf("Ng10b = %d\n", Ng10b);
  printf("Nvb   = %d\n", Nvb  );
#endif
  nLo = LO + 1;
  nLv = LV + 1;
  nLp = LP + 1;
  qmd = qmdata_alloc(Qmax, nLo, nLv, nLp, Nga, Nqa, Nq1a, Nfb, Nf1b, Nfo, Nub, Nu1b, Nvb);

  if(!qmd){
    return NULL;
  }

  fseek(f, p, SEEK_SET);

  while(fscanf(f, " n %d%d%d%d", &Q, &Lo, &Lv, &Lp) == 4){
    qmd->used[Q] = 1;
    qmd->Lo[Q]   = Lo;
    qmd->Lv[Q]   = Lv;
    qmd->Lp[Q]   = Lp;
    for(int j=0; j<=Lo; j++){
      if(fscanf(f, "%d", qmd->p+Q*nLo+j) != 1){
        goto hell;
      }
    }
  }

  Nga   = 0;
  Nqa   = 0;
  Nq1a  = 0;
  Nfb   = 0;
  Nf1b  = 0;
  Nfo   = 0;
  Nub   = 0;
  Nu1b  = 0;
  Nvb   = 0;
  while ((fgets(s, sizeof(s), f)) && (!strchr(s, '$'))){
    double t1,t2,t3,t4,t5;
    int i1,i2,i3,i4,i5,i6,i7,i8,i9;
    if(sscanf(s, " ea %d%lf", &Q, &t1) == 2){
      qmd->ea[Q] = t1;
    }
    else if(sscanf(s, " fa %d%d%lf", &Q, &Lo, &t1) == 3){
      qmd->fa[Q*nLo+Lo] = t1;
    }
    else if(sscanf(s, " f2a %d%d%lf", &Q, &Lv, &t1) == 3){
      qmd->f2a[Q*nLv+Lv] = t1;
    }
    else if(sscanf(s, " r0 %d%lf", &Q, &t1) == 2){
      qmd->r0[Q] = t1;
    }
    else if(sscanf(s, " cf %d%d%lf", &Q, &Lp, &t1) == 3){
      qmd->cf[Q*nLp+Lp] = t1;
    }
    else if(sscanf(s, " bub %d%lf%lf", &Q, &t1, &t2) == 3){
      qmd->bub[Q].c = t1;
      qmd->bub[Q].a = t2;
    }
    else if(sscanf(s, " bu1b %d%lf%lf", &Q, &t1, &t2) == 3){
      qmd->bu1b[Q].c = t1;
      qmd->bu1b[Q].a = t2;
    }
    else if(sscanf(s, " aub %d%d%lf%lf", &Q, &Lo, &t1, &t2) == 4){
      qmd->aub[Q*nLo+Lo].c = t1;
      qmd->aub[Q*nLo+Lo].a = t2;
    }
    else if(sscanf(s, " au0b %d%d%lf%lf", &Q, &Lo, &t1, &t2) == 4){
      qmd->au0b[Q*nLo+Lo].c = t1;
      qmd->au0b[Q*nLo+Lo].a = t2;
    }
    else if(sscanf(s, " au1b %d%d%lf%lf", &Q, &Lv, &t1, &t2) == 4){
      qmd->au1b[Q*nLv+Lv].c = t1;
      qmd->au1b[Q*nLv+Lv].a = t2;
    }
    else if(sscanf(s, " afb %d%d%lf%lf", &Q, &Lo, &t1, &t2) == 4){
      qmd->afb[Q*nLo+Lo].c = t1;
      qmd->afb[Q*nLo+Lo].a = t2;
    }
    else if(sscanf(s, " af0b %d%d%lf%lf", &Q, &Lo, &t1, &t2) == 4){
      qmd->af0b[Q*nLo+Lo].c = t1;
      qmd->af0b[Q*nLo+Lo].a = t2;
    }
    else if(sscanf(s, " af1b %d%d%lf%lf", &Q, &Lv, &t1, &t2) == 4){
      qmd->af1b[Q*nLv+Lv].c = t1;
      qmd->af1b[Q*nLv+Lv].a = t2;
    }
    else if(sscanf(s, " afo %d%d%lf%lf", &Q, &Lo, &t1, &t2) == 4){
      qmd->afo[Q*nLo+Lo].c = t1;
      qmd->afo[Q*nLo+Lo].a = t2;
    }
    else if(sscanf(s, " aof %d%d%lf%lf", &Q, &Lp, &t1, &t2) == 4){
      qmd->aof[Q*nLp+Lp].c = t1;
      qmd->aof[Q*nLp+Lp].a = t2;
    }
    else if(sscanf(s, " agb %d%d%lf", &Q, &Lo, &t1) == 3){
      qmd->agb[Q*nLo+Lo] = t1;
    }
    else if(sscanf(s, " ag0b %d%d%lf", &Q, &Lo, &t1) == 3){
      qmd->ag0b[Q*nLo+Lo] = t1;
    }
    else if(sscanf(s, " ag1b %d%d%lf", &Q, &Lv, &t1) == 3){
      qmd->ag1b[Q*nLv+Lv] = t1;
    }
    else if(sscanf(s, " bvb %d%lf%lf", &Q, &t1, &t2) == 3){
      qmd->bvb[Q].c = t1;
      qmd->bvb[Q].a = t2;
    }
    else if(sscanf(s, " avb %d%d%lf%lf", &Q, &Lo, &t1, &t2) == 4){
      qmd->avb[Q*nLo+Lo].c = t1;
      qmd->avb[Q*nLo+Lo].a = t2;
    }
    else if(sscanf(s, " ag6 %d%d%lf%lf", &Q, &Lo, &t1, &t2) == 4){
      qmd->ag6[Q*nLo+Lo].c = t1;
      qmd->ag6[Q*nLo+Lo].a = t2;
    }
    else if(sscanf(s, " ga %d%d%d%d%d%d%lf", &Q, &i1, &i2, &i3, &i4, &i5, &t1) == 7){
      qmd->ga[Nga].q   =  Q;
      qmd->ga[Nga].lu  =  i1;
      qmd->ga[Nga].lv  =  i2;
      qmd->ga[Nga].lu1 =  i3;
      qmd->ga[Nga].lv1 =  i4;
      qmd->ga[Nga].l   =  i5;
      qmd->ga[Nga].v   =  t1;
      Nga++;
    }
    else if(sscanf(s, " qa %d%d%d%d%lf", &Q, &i1, &i2, &i3, &t1) == 5){
      qmd->qa[Nqa].q  = Q;
      qmd->qa[Nqa].lu = i1;
      qmd->qa[Nqa].lv = i2;
      qmd->qa[Nqa].l  = i3;
      qmd->qa[Nqa].qa = t1;
      Nqa++;
    }
    else if(sscanf(s, " q1a %d%d%d%d%lf", &Q, &i1, &i2, &i3, &t1) == 5){
      qmd->q1a[Nq1a].q  = Q;
      qmd->q1a[Nq1a].lu = i1;
      qmd->q1a[Nq1a].lv = i2;
      qmd->q1a[Nq1a].l  = i3;
      qmd->q1a[Nq1a].qa = t1;
      Nq1a++;
    }
    else if(sscanf(s, " fb %d%d%d%d%d a** %lf%lf%lf%lf", &i1, &i2, &i3, &i4, &i5, &t1, &t2, &t3, &t4) == 9){
      qmd->fb[Nfb].qu = i1;
      qmd->fb[Nfb].qv = i2;
      qmd->fb[Nfb].lu = i3;
      qmd->fb[Nfb].lv = i4;
      qmd->fb[Nfb].m  = i5;
      qmd->fb[Nfb].c0 = t1;
      qmd->fb[Nfb].a  = t4;
      qmd->fb[Nfb].b  = 0.0;
      Nfb++;
    }
    else if(sscanf(s, " f1b %d%d%d%d%d b* %lf%lf%lf%lf%lf", &i1, &i2, &i3, &i4, &i5, &t1, &t2, &t3, &t4, &t5) == 10){
      qmd->f1b[Nf1b].qu = i1;
      qmd->f1b[Nf1b].qv = i2;
      qmd->f1b[Nf1b].lu = i3;
      qmd->f1b[Nf1b].lv = i4;
      qmd->f1b[Nf1b].m  = i5;
      qmd->f1b[Nf1b].c0 = t1;
      qmd->f1b[Nf1b].a  = t4;
      qmd->f1b[Nf1b].b  = t5;
      Nf1b++;
    }
    else if(sscanf(s, " fo %d%d%d%d%d b* %lf%lf%lf%lf%lf", &i1, &i2, &i3, &i4, &i5, &t1, &t2, &t3, &t4, &t5) == 10){
      qmd->fo[Nfo].qu = i1;
      qmd->fo[Nfo].qv = i2;
      qmd->fo[Nfo].lu = i3;
      qmd->fo[Nfo].lv = i4;
      qmd->fo[Nfo].m  = i5;
      qmd->fo[Nfo].c0 = t1;
      qmd->fo[Nfo].a  = t4;
      qmd->fo[Nfo].b  = t5;
      Nfo++;
    }
    else if(sscanf(s, " ub %d%d%d%d%d b* %lf%lf%lf%lf%lf", &i1, &i2, &i3, &i4, &i5, &t1, &t2, &t3, &t4, &t5) == 10){
      qmd->ub[Nub].qu = i1;
      qmd->ub[Nub].qv = i2;
      qmd->ub[Nub].lu = i3;
      qmd->ub[Nub].lv = i4;
      qmd->ub[Nub].l  = i5;
      qmd->ub[Nub].c0 = t1;
      qmd->ub[Nub].a  = t4;
      qmd->ub[Nub].b  = t5;
      Nub++;
    }
    else if(sscanf(s, " u1b %d%d%d%d%d b* %lf%lf%lf%lf%lf", &i1, &i2, &i3, &i4, &i5, &t1, &t2, &t3, &t4, &t5) == 10){
      qmd->u1b[Nu1b].qu = i1;
      qmd->u1b[Nu1b].qv = i2;
      qmd->u1b[Nu1b].lu = i3;
      qmd->u1b[Nu1b].lv = i4;
      qmd->u1b[Nu1b].l  = i5;
      qmd->u1b[Nu1b].c0 = t1;
      qmd->u1b[Nu1b].a  = t4;
      qmd->u1b[Nu1b].b  = t5;
      Nu1b++;
    }
    else if(sscanf(s, " vb %d%d%d%d%d%d%d%d%d b* %lf%lf%lf%lf%lf", &i1, &i2, &i3, &i4, &i5, &i6, &i7, &i8, &i9, &t1, &t2, &t3, &t4, &t5) == 14){
      qmd->vb[Nvb].qu  = i1;
      qmd->vb[Nvb].qv  = i2;
      qmd->vb[Nvb].lu  = i3;
      qmd->vb[Nvb].lv  = i4;
      qmd->vb[Nvb].lu1 = i5;
      qmd->vb[Nvb].lv1 = i6;
      qmd->vb[Nvb].l   = i7;
      qmd->vb[Nvb].l1  = i8;
      qmd->vb[Nvb].m   = i9;
      qmd->vb[Nvb].c0  = t1;
      qmd->vb[Nvb].a   = t4;
      qmd->vb[Nvb].b   = t5;
      Nvb++;
    }

  }

// assume parameters are in the right order
  for(int i=0; i<Nga; i++){
    Q = qmd->ga[i].q;
    qmd->q_list[Q].ga++;
  }
  for(int i=0; i<Nqa; i++){
    Q = qmd->qa[i].q;
    qmd->q_list[Q].qa++;
  }
  for(int i=0; i<Nq1a; i++){
    Q = qmd->q1a[i].q;
    qmd->q_list[Q].q1a++;
  }
  for(Q=0; Q<Qmax; Q++){
    qmd->q_list[Q+1].ga  += qmd->q_list[Q].ga;
    qmd->q_list[Q+1].qa  += qmd->q_list[Q].qa;
    qmd->q_list[Q+1].q1a += qmd->q_list[Q].q1a;
  }
#if 0
  for(int i=0; i<Qmax+1; i++){
    printf("%d %d\n", i, qmd->q_list[i].ga);
  }
#endif

  for(int i=0; i<Nfb; i++){
    int qu = qmd->fb[i].qu;
    int qv = qmd->fb[i].qv;
    if(qu>qv){
      GOTOHELL;
    }
    int qq = mpos(qu,qv);
    qmd->qq_list[qq].fb++;
  }
  for(int i=0; i<Nf1b; i++){
    int qu = qmd->f1b[i].qu;
    int qv = qmd->f1b[i].qv;
    int qq = qu<=qv?mpos(qu,qv):mpos(qv,qu);
    qmd->qq_list[qq].f1b++;
  }
  for(int i=0; i<Nfo; i++){
    int qu = qmd->fo[i].qu;
    int qv = qmd->fo[i].qv;
    int qq = qu<=qv?mpos(qu,qv):mpos(qv,qu);
    qmd->qq_list[qq].fo++;
  }
  for(int i=0; i<Nub; i++){
    int qu = qmd->ub[i].qu;
    int qv = qmd->ub[i].qv;
    int qq = qu<=qv?mpos(qu,qv):mpos(qv,qu);
    qmd->qq_list[qq].ub++;
  }
  for(int i=0; i<Nu1b; i++){
    int qu = qmd->u1b[i].qu;
    int qv = qmd->u1b[i].qv;
    int qq = qu<=qv?mpos(qu,qv):mpos(qv,qu);
    qmd->qq_list[qq].u1b++;
  }
  for(int i=0; i<Nvb; i++){
    int qu = qmd->vb[i].qu;
    int qv = qmd->vb[i].qv;
    int qq = qu<=qv?mpos(qu,qv):mpos(qv,qu);
    qmd->qq_list[qq].vb++;
  }

  for(int qq=1; qq<symsize(Qmax+1); qq++){
    qmd->qq_list[qq].fb  += qmd->qq_list[qq-1].fb;
    qmd->qq_list[qq].f1b += qmd->qq_list[qq-1].f1b;
    qmd->qq_list[qq].fo  += qmd->qq_list[qq-1].fo;
    qmd->qq_list[qq].u1b += qmd->qq_list[qq-1].u1b;
    qmd->qq_list[qq].ub  += qmd->qq_list[qq-1].ub;
    qmd->qq_list[qq].vb  += qmd->qq_list[qq-1].vb;
  }
#if 0
  for(int j=1; j<=Qmax; j++){
    for(int i=1; i<=j; i++){
      printf("%2d %2d   %3d\n", i, j, qmd->qq_list[mpos(i,j)].fb);
    }
  }
#endif

  return qmd;

hell:
  free(qmd);
  return NULL;
}

static void ca_l_print(int Q, ca_t * par, const char s[], int * Lq, int nL, FILE * f){
  for(int L=0; L<=Lq[Q]; L++){
    fprintf(f, "%s %3d %3d  % 12.7lf % 12.7lf\n",
        s, Q, L, par[Q*nL+L].c, par[Q*nL+L].a);
  }
  return;
}

static void d_l_print(int Q, double * par, const char s[], int * Lq, int nL, FILE * f){
  for(int L=0; L<=Lq[Q]; L++){
    fprintf(f, "%s %3d %3d % 12.7lf\n",
        s, Q, L, par[Q*nL+L]);
  }
  return;
}

void qmdata_print(FILE * f, qmdata * qmd){
  int Q,L;

  for(Q=0; Q<=qmd->Qmax; Q++){
    if(!qmd->used[Q]) continue;
    fprintf(f, "n %3d   %2d %2d %2d   ", Q, qmd->Lo[Q], qmd->Lv[Q], qmd->Lp[Q]);
    for(L=0; L<=qmd->Lo[Q]; L++){
      fprintf(f, "%2d ", qmd->p[Q*(qmd->nLo)+L]);
    }
    fprintf(f, "\n");
  }

  for(Q=0; Q<=qmd->Qmax; Q++){
    if(!qmd->used[Q]) continue;
    fprintf(f, "ea %3d % 14.7lf\n", Q, qmd->ea[Q]);
  }

  for(Q=0; Q<=qmd->Qmax; Q++){
    if(!qmd->used[Q]) continue;
    D_L_PRINT(f, Q, fa,  Lo);
    D_L_PRINT(f, Q, f2a, Lv);
  }

  for(Q=0; Q<=qmd->Qmax; Q++){
    if(!qmd->used[Q]) continue;
    fprintf(f, "r0 %3d %12.7lf\n", Q, qmd->r0[Q]);
  }

  for(Q=0; Q<=qmd->Qmax; Q++){
    if(!qmd->used[Q]) continue;
    D_L_PRINT(f, Q, cf,  Lp);
  }

  for(Q=0; Q<=qmd->Qmax; Q++){
    if(!qmd->used[Q]) continue;
    fprintf(f, "bub %3d  % 12.7lf % 12.7lf\n", Q, qmd->bub[Q].c, qmd->bub[Q].a);
    CA_L_PRINT(f, Q, aub, Lo);
  }

  for(Q=0; Q<=qmd->Qmax; Q++){
    if(!qmd->used[Q]) continue;
    CA_L_PRINT(f, Q, afb, Lo);
  }

  for(Q=0; Q<=qmd->Qmax; Q++){
    if(!qmd->used[Q]) continue;
    fprintf(f, "bu1b %3d  % 12.7lf % 12.7lf\n", Q, qmd->bu1b[Q].c, qmd->bu1b[Q].a);
    CA_L_PRINT(f, Q, au0b, Lo);
    CA_L_PRINT(f, Q, au1b, Lv);
  }

  for(Q=0; Q<=qmd->Qmax; Q++){
    if(!qmd->used[Q]) continue;
    CA_L_PRINT(f, Q, af0b, Lo);
    CA_L_PRINT(f, Q, af1b, Lv);
  }

  for(Q=0; Q<=qmd->Qmax; Q++){
    if(!qmd->used[Q]) continue;
    CA_L_PRINT(f, Q, afo, Lo);
    CA_L_PRINT(f, Q, aof, Lp);
  }

  for(Q=0; Q<=qmd->Qmax; Q++){
    if(!qmd->used[Q]) continue;
    D_L_PRINT(f, Q, agb,  Lo);
    D_L_PRINT(f, Q, ag0b, Lo);
    D_L_PRINT(f, Q, ag1b, Lv);
  }

  for(Q=0; Q<=qmd->Qmax; Q++){
    if(!qmd->used[Q]) continue;
    CA_L_PRINT(f, Q, avb, Lo);
    fprintf(f, "bvb %3d  % 12.7lf % 12.7lf\n", Q, qmd->bvb[Q].c, qmd->bvb[Q].a);
  }

  for(Q=0; Q<=qmd->Qmax; Q++){
    if(!qmd->used[Q]) continue;
    CA_L_PRINT(f, Q, ag6, Lo);
  }

  for(Q=1; Q<=qmd->Qmax; Q++){
    int bra = qmd->q_list[Q-1].ga;
    int ket = qmd->q_list[Q  ].ga;
    for(int i=bra; i<ket; i++){
      fprintf(f, "ga %3d  %d %d %d %d  %d  % 12.7lf\n", Q, qmd->ga[i].lu, qmd->ga[i].lv, qmd->ga[i].lu1, qmd->ga[i].lv1, qmd->ga[i].l, qmd->ga[i].v);
    }
  }

  for(Q=1; Q<=qmd->Qmax; Q++){
    int bra,ket;
    bra = qmd->q_list[Q-1].qa;
    ket = qmd->q_list[Q  ].qa;
    for(int i=bra; i<ket; i++){
      fprintf(f, "qa %3d  %d %d %d  % 12.7lf\n", Q, qmd->qa[i].lu, qmd->qa[i].lv, qmd->qa[i].l, qmd->qa[i].qa);
    }
    bra = qmd->q_list[Q-1].q1a;
    ket = qmd->q_list[Q  ].q1a;
    for(int i=bra; i<ket; i++){
      fprintf(f, "q1a %3d  %d %d %d  % 12.7lf\n", Q, qmd->q1a[i].lu, qmd->q1a[i].lv, qmd->q1a[i].l, qmd->q1a[i].qa);
    }
  }

  for(int Qv=1; Qv<=qmd->Qmax; Qv++){
    for(int Qu=1; Qu<=Qv; Qu++){
      int qq  = mpos(Qu,Qv);
      int bra, ket;
      bra = qmd->qq_list[qq-1].fb;
      ket = qmd->qq_list[qq  ].fb;
      for(int i=bra; i<ket; i++){
        fprintf(f, "fb %3d %3d   %d %d   %d a** % 12.7lf % 12.7lf % 12.7lf % 12.7lf\n",
            Qu, Qv, qmd->fb[i].lu, qmd->fb[i].lv, qmd->fb[i].m,
            qmd->fb[i].c0, 0.0,0.0, qmd->fb[i].a);
      }
      bra = qmd->qq_list[qq-1].f1b;
      ket = qmd->qq_list[qq  ].f1b;
      for(int i=bra; i<ket; i++){
        fprintf(f, "f1b %3d %3d   %d %d   %d b* % 12.7lf % 12.7lf % 12.7lf % 12.7lf % 12.7lf\n",
            qmd->f1b[i].qu, qmd->f1b[i].qv, qmd->f1b[i].lu,
            qmd->f1b[i].lv, qmd->f1b[i].m,  qmd->f1b[i].c0,
            0.0, 0.0, qmd->f1b[i].a, qmd->f1b[i].b);
      }
    }
  }

  for(int Qv=1; Qv<=qmd->Qmax; Qv++){
    for(int Qu=1; Qu<=Qv; Qu++){
      int qq  = mpos(Qu,Qv);
      int bra = qmd->qq_list[qq-1].fo;
      int ket = qmd->qq_list[qq  ].fo;
      for(int i=bra; i<ket; i++){
        fprintf(f, "fo %3d %3d   %d %d   %d b* % 12.7lf % 12.7lf % 12.7lf % 12.7lf % 12.7lf\n",
            qmd->fo[i].qu, qmd->fo[i].qv, qmd->fo[i].lu, qmd->fo[i].lv, qmd->fo[i].m,
            qmd->fo[i].c0, 0.0,0.0, qmd->fo[i].a, qmd->fo[i].b);
      }
    }
  }

  for(int Qv=1; Qv<=qmd->Qmax; Qv++){
    for(int Qu=1; Qu<=Qv; Qu++){
      int qq  = mpos(Qu,Qv);
      int bra,ket;
      bra = qmd->qq_list[qq-1].ub;
      ket = qmd->qq_list[qq  ].ub;
      for(int i=bra; i<ket; i++){
        fprintf(f, "ub %3d %3d   %d %d   %d b* % 12.7lf % 12.7lf % 12.7lf % 12.7lf % 12.7lf\n",
            qmd->ub[i].qu, qmd->ub[i].qv, qmd->ub[i].lu, qmd->ub[i].lv, qmd->ub[i].l, qmd->ub[i].c0, 0.0,0.0, qmd->ub[i].a, qmd->ub[i].b);
      }
      bra = qmd->qq_list[qq-1].u1b;
      ket = qmd->qq_list[qq  ].u1b;
      for(int i=bra; i<ket; i++){
        fprintf(f, "u1b %3d %3d   %d %d   %d b* % 12.7lf % 12.7lf % 12.7lf % 12.7lf % 12.7lf\n",
            qmd->u1b[i].qu, qmd->u1b[i].qv, qmd->u1b[i].lu, qmd->u1b[i].lv, qmd->u1b[i].l, qmd->u1b[i].c0, 0.0,0.0, qmd->u1b[i].a, qmd->u1b[i].b);
      }
    }
  }

  for(int Qv=1; Qv<=qmd->Qmax; Qv++){
    for(int Qu=1; Qu<=Qv; Qu++){
      int qq  = mpos(Qu,Qv);
      int bra = qmd->qq_list[qq-1].vb;
      int ket = qmd->qq_list[qq  ].vb;
      for(int i=bra; i<ket; i++){
        fprintf(f, "vb %3d %3d   %d %d %d %d   %d %d  %d  b* % 12.7lf % 12.7lf % 12.7lf % 12.7lf % 12.7lf\n",
            qmd->vb[i].qu, qmd->vb[i].qv, qmd->vb[i].lu, qmd->vb[i].lv, qmd->vb[i].lu1, qmd->vb[i].lv1,
            qmd->vb[i].l, qmd->vb[i].l1, qmd->vb[i].m, qmd->vb[i].c0, 0.0,0.0, qmd->vb[i].a, qmd->vb[i].b);
      }
    }
  }
  return;
}

