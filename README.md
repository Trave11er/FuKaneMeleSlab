# Tight-binding model for a topological insulator slab [ARCHIVED]

**DISCLAIMER** This is bad code kept only for reproducibility reasons (see below)

- Microscopic tight-binding model of topological a insulator slab as in:

  - G. Siroki, P. D. Haynes, D. K. K. Lee and V. Giannini Phys. Rev. Mater. 1, 024201 (2017)

- The original tight-binding Hamiltonian is presented here:

  - L. Fu, C. L. Kane, and E. J. Mele, Phys. Rev. Lett. 98, 106803 (2007)

### Getting started
To compile run:a
```
g++ -std=c++11 -o slab slab.cc -llapack -lblas -lgfortran -larmadillo
```
To run execute
```
./slab
```
which will print the k-point and save the corresponding eigenvalues in the file output_bs.txt

### Dependencies

Armadillo, LAPACK, BLAS

## Authors

* **Gleb Siroki**

## License

This project is licensed under the MIT License - see the LICENSE.md file for details
