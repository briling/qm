
## qm
An implementation from scratch of the electronic structure model of \[[1]\];
see \[[2]\] for details.

---

### requirements
* `GNU/Linux` or `Cygwin`
* `gcc >= 4.7`
* `make`

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
* `task:%s`     - calculate energy (`energy`) or gradient (`grad`)
* `conv:%lf`    – scf convergence criterion (rms change in density matrix)
* `it:%d,%d`    – number of iterations / size of diis subspace
* `print:%d`    – printing options (`1` – default, `2` – print scf vectors, `3` – print atomic charges and bond orders)
* `read:%s`     – file name for reading scf vectors
* `write:%s`    – file name for saving scf vectors
* `diis:%d`     – scf algorithm (`0` – straightforward procedure, `1` – diis (default))
* `restrict:%d` – if the system has an even number of electrons, use `0` – spin-unrestricted or `1` – spin-restricted orbitals
* `field:%lf,%lf,%lf` – applied electric field (with the opposite sign), i.e.
<img src="http://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;\nabla\phi\equiv-\vec&space;E" title="\nabla\phi\equiv-\vec E" height="22" />

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

`mol/*.p17.out` –
corresponding output files from Priroda-17.

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

