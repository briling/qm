for i in mol/*.in ; do
  echo $i;
  if [[ $i = *'.field.'* ]]; then
    ./qm qm_m.in $i ${i/\.in/.x.out1} print:3 field:0.01,0.02,0.03;
  elif [[ $i = *'.grad.'* ]]; then
    ./qm qm_m.in $i ${i/\.in/.x.out1} print:3 task:grad;
  else
    ./qm qm_m.in $i ${i/\.in/.x.out1} print:3 ;
  fi
done
for i in mol/{CH2_m1,CH2_m3,phenylsilatrane,SudanBlackB}.in ; do
  echo $i;
  ./qm qm_m.in $i ${i/\.in/.nodiis.x.out1} print:3 diis:0 it:64;
done

for i in mol/*.x.out ; do
  echo $i;
  diff --suppress-common-lines --side-by-side --report-identical-files $i ${i}1 ;
done

rm mol/*.out1

