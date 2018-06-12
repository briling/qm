for i in mol/*.in ; do
  echo $i;
  if [[ $i = *'.field.'* ]]; then
    ./qm qm_m.in $i ${i/\.in/.x.out1} print:3 field:0.01,0.02,0.03;
  else
    ./qm qm_m.in $i ${i/\.in/.x.out1} print:3;
  fi
done

for i in mol/*.x.out ; do
  diff --suppress-common-lines --side-by-side --report-identical-files $i ${i}1 ;
done

rm mol/*.x.out1

