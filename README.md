
## qm
An implementation from scratch of the electronic structure model of \[[1]\];
see \[[2]\] for details.

---

### requirements
* `GNU/Linux` or `Cygwin`
* `make`
* `gcc >= 4.7`

### build
```
make qm
make test
```

### usage
```
./qm qm_m.in <molecule>.{in,out}
```
#### command line options:
* `conv:%lf`  – scf convergence criterion
* `it:%d`     – number of iterations
* `print:%d`  – printing options (1 – default, 2 – print scf vectors, 3 – print atomic charges and bond orders)
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
corresponding output files from Priroda-11
(we used `bin/p1` \[[1a]\]).

---

### references

<a name="ref1">\[1\]</a>
D. N. Laikov,  [J. Chem. Phys.][L2011] **135**, 134120 (2011).

<a name="ref1a">\[1a\]</a>
Supplementary material of \[[1]\].

<a name="ref2">\[2\]</a>
K. R. Briling, [J. Chem. Phys.][B2017] **147**, 157101 (2017).

[1]: #ref1
[1a]: #ref1a
[2]: #ref2
[L2011]:https://doi.org/10.1063/1.3646498
[B2017]:https://doi.org/10.1063/1.5000525

<a href="http://m.maploco.com/details/9af5rrkn"> <img src="http://www.maploco.com/vmap/s/9693527.png" width=1 > </a>

