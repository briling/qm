#include "mol.h"
#include "vecn.h"

#define N 2
#define C 8

#define END(S,X) ( (S)->X + (X##_size)/sizeof(*((S)->X)) )

static mol * mol_alloc(mol * mold, int n, int c){

  int q_size = sizeof(int    ) * n     ;
  int m_size = sizeof(double ) * n     ;
  int s_size = sizeof(styp   ) * n     ;
  int r_size = sizeof(double ) * n * 3 ;
  int l_size = sizeof(int    ) *(n + 1);
  int k_size = sizeof(int    ) * c     ;
  int b_size = sizeof(int    ) * c     ;
  int size   = sizeof(mol) + q_size + m_size + s_size + r_size + l_size + k_size + b_size;

  mol * m = realloc(mold, size);
  if(!m){
    GOTOHELL;
  }

  m->n = n;
  m->c = c;

  m->r = (double *) (m + 1);
  m->m = (double *) END(m,r);
  m->q = (int    *) END(m,m);
  m->k = (int    *) END(m,q);
  m->b = (int    *) END(m,k);
  m->l = (int    *) END(m,b);
  m->s = (styp   *) END(m,l);
  return m;
}

static mol * redm(mol * m){

  int n = m->n;
  int c = m->c;

  int q_size = sizeof(int    ) * n     ;
  int m_size = sizeof(double ) * n     ;
  int s_size = sizeof(styp   ) * n     ;
  int r_size = sizeof(double ) * n * 3 ;
  int l_size = sizeof(int    ) *(n + 1);
  int k_size = sizeof(int    ) * c     ;
  int b_size = sizeof(int    ) * c     ;

  mol * tm = mol_alloc(NULL, n, c);

  memmove(tm->s, m->s, s_size);
  memmove(tm->q, m->q, q_size);
  memmove(tm->m, m->m, m_size);
  memmove(tm->r, m->r, r_size);
  memmove(tm->l, m->l, l_size);
  memmove(tm->k, m->k, k_size);
  memmove(tm->b, m->b, b_size);

  free(m);
  return tm;
}

static mol * expm(mol * m){

  int so = m->n;
  int sn = so ? 2*so : N;
  m->n = sn;

  int zo = m->c;
  int zn = zo ? 2*zo : C;
  m->c = zn;

  m = mol_alloc(m, sn, zn);

  int q_size = sizeof(int    ) * so     ;
  int m_size = sizeof(double ) * so     ;
  int s_size = sizeof(styp   ) * so     ;
  int r_size = sizeof(double ) * so * 3 ;
  int l_size = sizeof(int    ) *(so + 1);
  int k_size = sizeof(int    ) * zo    ;
  int b_size = sizeof(int    ) * zo    ;

  mol tm;
  tm.r = (double *) (m + 1);
  tm.m = (double *) END(&tm,r);
  tm.q = (int    *) END(&tm,m);
  tm.k = (int    *) END(&tm,q);
  tm.b = (int    *) END(&tm,k);
  tm.l = (int    *) END(&tm,b);
  tm.s = (styp   *) END(&tm,l);

  memmove(m->s, tm.s, s_size);
  memmove(m->l, tm.l, l_size);
  memmove(m->b, tm.b, b_size);
  memmove(m->k, tm.k, k_size);
  memmove(m->q, tm.q, q_size);
  memmove(m->m, tm.m, m_size);

  return m;
}

static int chckm(mol * m){
  int i,j, a,b, t;
  for(j=0; j<m->n; j++){
    a = m->l[j];
    b = m->l[j+1];
    for(i=a; i<b; i++){
      t = m->k[i];
      if ( t<0 || t>=m->n || t == j){
        return 1;
      }
    }
  }
  return 0;
}

mol * mol_read(FILE * f){

  const double rd = M_PI/180.0;
  mol  * m;
  double r[3];
  char   s[256];
  char   ch, af, zcf;
  styp   t;
  int    q, n,c, a,b, z,mult;
  double mass;
  af   = 1;
  zcf  = 0;
  z    = 0;
  mult = 1;

  do{
    while (fscanf(f, " %255[^$ ]", s) == 1) {
        #if 0
          fprintf(stderr, "%s\n", s);
        #endif
    }
    if (fscanf(f, "%255[$A-z]", s) != 1){
      return NULL;
    }
  } while(strcmp(s, "$molecule"));
  while (fscanf(f, " %255[^$ \n]", s) == 1) {
#if 0
    fprintf(stderr, "%s\n", s);
#endif
    if (s[0] == 'u' ){
      if ((*(strchr(s, '=')+1))=='b'){
        af = 0;
      }
    }
    else if( ! strncmp(s, "charge", 6 ) ){
      sscanf(s, "charge=%d", &z);
    }
    else if( ! strncmp(s, "mult", 4 ) ){
      sscanf(s, "mult=%d", &mult);
    }
    else if (! strncmp(s, "cartesian", 4) ){
      zcf = 0;
      break;
    }
    else if( ! strncmp(s, "z-matrix", 1 ) ){
      zcf = 1;
      break;
    }
  }

  m = calloc(sizeof(mol), 1);
  m = expm(m);

  n = 0;
  c = 0;

  while (fscanf(f, " $%255s", s) != 1){
    if (zcf == 0){
      if (fscanf(f, "%d%lf%lf%lf", &q, r, r+1, r+2) != 4){
        if ( fscanf(f, " %255[^\n]", s) && strstr(s, "set")){
          continue;
        }
        goto hell;
      }
    }
    else if (zcf == 1){
      char   sc;
      int    a1, a2, a3;
      double ab, ac, az;
      switch(n){
        case 0:
          sc = (fscanf(f, "%d", &q) != 1);
          break;
        case 1:
          sc = (fscanf(f, "%d%d%lf", &q, &a1, &ab) != 3);
          break;
        case 2:
          sc = (fscanf(f, "%d%d%lf%d%lf", &q, &a1, &ab, &a2, &ac) != 5);
          break;
        default:
          sc = (fscanf(f, "%d%d%lf%d%lf%d%lf", &q, &a1, &ab, &a2, &ac, &a3, &az) != 7);
          break;
      }
      if(sc || zmat2cart(n, m->r, r, a1-1, a2-1, a3-1, ab, ac*rd, az*rd)){
        goto hell;
      }
    }

#if 0
    printf("%d %lf %lf %lf\n", n, r[0], r[1], r[2]);
#endif

    if((n+1)>(m->n)){
      m = expm(m);
    }

    mass = -1.0;
    fscanf(f, " mass = %lf", &mass);

    t[0] = 0;
    fscanf(f, " type = %7s", t);

    m->q[n]     = q;
    m->r[3*n]   = r[0];
    m->r[3*n+1] = r[1];
    m->r[3*n+2] = r[2];
    strncpy(m->s[n], t, sizeof(styp));
    m->m[n] = mass;

    m->l[n] = c;
    if (fscanf(f, " k%c",&ch) == 1){
      do{
        a = 0; b = 1;
        fscanf(f, "%d(%d)", &a, &b);
        if ((c+1)>(m->c)){
          m = expm(m);
        }
        m->k[c] = a-1;
        m->b[c] = b;
        c++;
      } while((char)getc(f)==',');
    }

    n++;
  }

  m->l[n] = c;
  if ( (m->n > n) || (m->c > c) ){
    m->n = n;
    m->c = c;
    m = redm(m);
  }

  m->z    = z;
  m->mult = mult;

  if (chckm(m)){
    goto hell;
  }

  if (af == 1){
    vecscal(m->n*3, m->r, AB);
  }

  return m;

  hell:
    free(m);
    return NULL;
}

