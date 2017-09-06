static inline int nel(int q, qmdata * qmd){
  int n = 0;
  for(int l=0; l<=qmd->Lo[q]; l++){
    n += qmd->p[q*(qmd->nLo)+l];
  }
  return n;
}
