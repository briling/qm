
static inline double r3dot (double * u, double * v){
  return (u[0]*v[0] + u[1]*v[1] + u[2]*v[2]);
}

static inline void r3x (double * w, double * u, double * v){
  w[0] =  u[1]*v[2]-u[2]*v[1];
  w[1] = -u[0]*v[2]+u[2]*v[0];
  w[2] =  u[0]*v[1]-u[1]*v[0];
  return;
}

static inline void r3mx (double * u, double * v, double * m){
  u[0] =  v[0]*m[0] + v[1]*m[1] + v[2]*m[2];
  u[1] =  v[0]*m[3] + v[1]*m[4] + v[2]*m[5];
  u[2] =  v[0]*m[6] + v[1]*m[7] + v[2]*m[8];
  return;
}

static inline void r3scal (double * u, double c){
  u[0] *= c;
  u[1] *= c;
  u[2] *= c;
  return;
}

static inline void r3adds (double * u, double * v, double c){
  u[0] += c*v[0];
  u[1] += c*v[1];
  u[2] += c*v[2];
  return;
}

static inline void r3add (double * u, double * v){
  u[0] += v[0];
  u[1] += v[1];
  u[2] += v[2];
  return;
}

static inline void r3min (double * u, double * v){
  u[0] -= v[0];
  u[1] -= v[1];
  u[2] -= v[2];
  return;
}

static inline void r3diff (double * w, double * u, double * v){
  w[0] = u[0] - v[0];
  w[1] = u[1] - v[1];
  w[2] = u[2] - v[2];
  return;
}

static inline void r3cp (double * u, double * v){
  u[0] = v[0];
  u[1] = v[1];
  u[2] = v[2];
  return;
}

static inline void r3cpsc (double * u, double * v, double c){
  u[0] = v[0]*c;
  u[1] = v[1]*c;
  u[2] = v[2]*c;
  return;
}

