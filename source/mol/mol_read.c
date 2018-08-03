#include "mol.h"
#include "vecn.h"

#define N 2
#define C 8

static mol * redm(mol * m){

  int n = m->n;
  int c = m->c;

  int size = sizeof(mol    )         +
             sizeof(int    ) * n     + // q
             sizeof(double ) * n     + // m
             sizeof(styp   ) * n     + // s
             sizeof(double ) * n * 3 + // r
             sizeof(int    ) *(n + 1)+ // l
             sizeof(int    ) * c     + // k
             sizeof(int    ) * c     ; // b

  mol * tm = malloc(size);

  tm->n = n;
  tm->c = c;
  tm->s = (styp   *)(tm    + 1    );
  tm->q = (int    *)(tm->s + n    );
  tm->m = (double *)(tm->q + n    );
  tm->r = (double *)(tm->m + n    );
  tm->l = (int    *)(tm->r + n * 3);
  tm->k = (int    *)(tm->l + n + 1);
  tm->b = (int    *)(tm->k + c    );

  memmove(tm->s, m->s, sizeof(styp  )*n    );
  memmove(tm->q, m->q, sizeof(int   )*n    );
  memmove(tm->m, m->m, sizeof(double)*n    );
  memmove(tm->r, m->r, sizeof(double)*n * 3);
  memmove(tm->l, m->l, sizeof(int   )*(n+1));
  memmove(tm->k, m->k, sizeof(int   )*c    );
  memmove(tm->b, m->b, sizeof(int   )*c    );

  free(m);
  return tm;

}

static mol * expm(mol * m){

  int so = m->n;
  int sn;
  if (so == 0){
    sn = N;
  }
  else{
    sn = 2*so;
  }
  m->n = sn;

  int zo = m->c;
  int zn;
  if (zo == 0){
    zn = C;
  }
  else{
    zn = 2*zo;
  }
  m->c = zn;

  int size = sizeof(mol    )          +
             sizeof(int    ) * sn     + // q
             sizeof(double ) * sn     + // m
             sizeof(styp   ) * sn     + // s
             sizeof(double ) * sn * 3 + // r
             sizeof(int    ) *(sn + 1)+ // l
             sizeof(int    ) * zn     + // k
             sizeof(int    ) * zn     ; // b

  m = realloc (m, size );
  m->s = (styp   *)(m    +  1  );
  m->q = (int    *)(m->s + sn  );
  m->m = (double *)(m->q + sn  );
  m->r = (double *)(m->m + sn  );
  m->l = (int    *)(m->r + sn*3);
  m->k = (int    *)(m->l + sn+1);
  m->b = (int    *)(m->k + zn  );

  mol tm;
  tm.s = (styp   *)(m    +  1  );
  tm.q = (int    *)(tm.s + so  );
  tm.m = (double *)(tm.q + so  );
  tm.r = (double *)(tm.m + so  );
  tm.l = (int    *)(tm.r + so*3);
  tm.k = (int    *)(tm.l + so+1);
  tm.b = (int    *)(tm.k + zo  );

  memmove(m->b, tm.b, sizeof(int   )* zo   );
  memmove(m->k, tm.k, sizeof(int   )* zo   );
  memmove(m->l, tm.l, sizeof(int   )*(so+1));
  memmove(m->r, tm.r, sizeof(double)* so*3 );
  memmove(m->m, tm.m, sizeof(double)* so   );
  memmove(m->q, tm.q, sizeof(int   )* so   );

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

  m = malloc(sizeof(mol));
  m = memset(m, 0, sizeof(mol));
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

