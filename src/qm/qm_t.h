
#ifndef QM_T_H
#define QM_T_H

typedef struct{
  double c;
  double a;
} ca_t;

typedef struct{
  int q;
  int lu;
  int lv;
  int lu1;
  int lv1;
  int l;
  double v;
} ga_t;

typedef struct{
  int q;
  int lu;
  int lv;
  int l;
  double qa;
} qa_t;

typedef struct{
  int qu;
  int qv;
  int lu;
  int lv;
  int m;
  double c0;
  double a;
  double b;
} cabm_t;

typedef struct{
  int qu;
  int qv;
  int lu;
  int lv;
  int l;
  double c0;
  double a;
  double b;
} cabl_t;

typedef struct{
  int qu;
  int qv;
  int lu;
  int lv;
  int lu1;
  int lv1;
  int l;
  int l1;
  int m;
  double c0;
  double a;
  double b;
} cablm_t;

typedef struct{
  int ga;
  int qa;
  int q1a;
} ql_t;

typedef struct{
  int fb;
  int f1b;
  int fo;
  int ub;
  int u1b;
  int vb;
} qql_t;

typedef struct{
  int   Qmax;
  int   nLo;   // Lo_max + 1
  int   nLv;   // Lv_max + 1
  int   nLp;   // Lp_max + 1
  int * used;

  int * Lo;
  int * Lv;
  int * Lp;
  int * p;

  double * ea;                // eq2
  double * fa;                // eq26
  double * f2a;               // eq41
  double * r0;                // -
  double * cf;                // eq31
  ca_t   * bub;               // eq49 mm
  ca_t   * bu1b;              // eq49 mp
  ca_t   * aub;               // eq49 mm
  ca_t   * au0b;              // eq49 mp(m)
  ca_t   * au1b;              // eq49 mp(p)
  ca_t   * afb;               // eq47 mm
  ca_t   * af0b;              // eq48 m
  ca_t   * af1b;              // eq48 p
  ca_t   * afo;               // eq50
  ca_t   * aof;               // eq50
  double * agb;               // eq52
  double * ag0b;              // eq52 mp(m)
  double * ag1b;              // eq52 mp(p)
  ca_t   * bvb;               // eq51
  ca_t   * avb;               // eq51
  ca_t   * ag6;               // eq53

  ql_t   * q_list;            // -
  ga_t   * ga;                // eq35
  qa_t   * qa;                // eq52 mm
  qa_t   * q1a;               // eq52 mp

  qql_t   * qq_list;          // -
  cabm_t  * fb;               // eq47 mm (sym)
  cabm_t  * f1b;              // eq48 mp
  cabm_t  * fo;               // eq50
  cabl_t  * ub;               // eq49 mm
  cabl_t  * u1b;              // eq49 mp
  cablm_t * vb;               // eq51
} qmdata;


typedef struct{
  int   M;
  int * Q;  // charge
  int * m;  // azimuthal
  int * l;  // angular
  int * k;  // center
  int * n;  // N wrt center
} basis;

typedef struct {
  double cos_b ;
  double sin_b ;
  double cos_g ;
  double sin_g ;
  double r     ;
} euler;

#endif

