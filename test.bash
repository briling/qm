for i in mol/*.in ; do
  echo $i;
  ./qm qm_m.in $i ${i/\.in/.x.out1} ;
done

for i in mol/*.x.out ; do
  diff --suppress-common-lines --side-by-side --report-identical-files $i ${i}1 ;
done

rm mol/*.x.out1

