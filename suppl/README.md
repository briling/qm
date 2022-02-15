This is a tool to prepare the triple products 
<img src="http://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;B_{mm'm''}^{\ell\;\;\ell'\;\;\ell''}" title="B_{mm'm''}^{\ell\ell'\ell''}" height="22" />
of real spherical harmonics:


<img src="http://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;B_{mm'm''}^{\ell\;\;\ell'\;\;\ell''}=\displaystyle2\sqrt{\pi}\int_0^\pi{}{\rm\,d}\phi\int_0^{2\pi}\sin\theta{}{\rm\,d}\theta\:\:Y_{\ell{}m}(\theta,\phi)\:Y_{\ell'm'}(\theta,\phi)\:Y_{\ell''m''}(\theta,\phi)."
     title="B_{mm'm''}^{\ell\ell'\ell''}=2\sqrt\pi\int_0^{\pi}d\phi\int_0^{2\pi}\sin\theta{}d\theta{}Y_{\ell{}m}(\theta,\phi)Y_{\ell'm'}(\theta,\phi)Y_{\ell''m''}(\theta,\phi)"
     height="44" />


#### Options ($x denotes an integer cli argument):
|||
| ----------- | ----------- |
 |   `./B c $m`                         | prints a real spherical harmonic <img src="http://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;Y_{\ell{}m}" title="Y_{\ell{}m}" height="15" /> as a sum of complex spherical harmonics <img src="http://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;Y_\ell^{\pm{}m}" title="Y_\ell^{\pm{}m}" height="22" />|
|    `./B s $j1 $j2 $j3 $m1 $m2 $m3`     |prints Wigner's 3-jm symbol (for complex spherical harmonics)
|    `./B r $j1 $j2 $j3 $m1 $m2 $m3`     |prints its analogue for real spherical harmonics
|    `./B i $j1 $j2 $j3 $m1 $m2 $m3`     |prints the integral <img src="http://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;B_{m_1m_2m_3}^{j_1\;\;j_2\;\;j_3}" title="B_{m_1m_2m_3}^{j_1j_2j_3}" height="22" />
 |   `./B t`                             |prints a table with non-zero integrals (the max j is hard-coded to be 3 but can be easily changed)
|   `./B f`                             |prints fome `C` code
