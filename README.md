## qm
An implementation from scratch of the electronic structure model from \[[1]\].

---

### requirements
* `GNU/Linux` or `Cygwin`
* `make`
* `gcc`

### build
```
make qm
```

### usage
```
./qm qm_m.in <molecule>.{in,out}
```
#### command line options:
* `tol:%lf`   – scf convergence criterion
* `it:%d`     – number of iterations
* `read:%s`   – file name for reading scf vectors 
* `write:%s`  – file name for saving scf vectors 

---

### files

`qm_m.in`  –
the file with the set of parameters
(slightly modified `qm.in` \[[1a]\]).

`mol/*.in`  –
input files with molecular geometries,
the format is described in `README` \[[1a]\].
Our program reads the section `$molecule` only.

`mol/*.x.out`  –
corresponding output files from our program.

`mol/*.p11.out` –
corresponding output files from PRIRODA-11
(we used `bin/p1` \[[1a]\]).

---

### references

<a name="ref1">\[1\]</a>
D. N. Laikov, J. Chem. Phys. **135**, 134120 (2011).
DOI: [10.1063/1.3646498](http://dx.doi.org/10.1063/1.3646498)

<a name="ref1a">\[1a\]</a>
Supplementary material of \[[1]\].

[1]: #ref1
[1a]: #ref1a

