# cddlib

The C-library cddlib is a C implementation of the *Double Description Method*
of Motzkin et al. for generating all vertices (i.e. extreme points) and extreme
rays of a general convex polyhedron in R<sup>d</sup> given by a system of
linear inequalities:

> P = { x=(x<sub>1</sub>, ..., x<sub>d</sub>)<sup>T</sup> : b - A·x ≥ 0 }

where A is a given m×d real matrix, b is a given m-vector and 0 is the m-vector
of all zeros.

The program can be used for the reverse operation (i.e. convex hull
computation). This means that one can move back and forth between an inequality
representation and a generator (i.e. vertex and ray) representation of a
polyhedron with cddlib. Also, cddlib can solve a linear programming problem,
i.e. a problem of maximizing and minimizing a linear function over P.

## Report an Issue
If you find some problems or bugs, please kindly report them to our [issue
tracker](https://github.com/cddlib/cddlib/issues).

## Documentation
The
[documentation](ftp://ftp.math.ethz.ch/users/fukudak/cdd/cddlibman/cddlibman.html)
is still incomplete but contains descriptions of the most important functions.
Please look at the [examples](src/) for how the library functions can be used
in user's C-codes.

## H-representation & V-polyhedron
One convenient feature of cddlib is the ability to handle essentially any data.
More precisely, it can generate an H-representation of a V-polyhedron which is
not full-dimensional, and it can generate a V-representation of an H-polyhedron
which has no extreme points.

## Numerical Problems
A little caution is in order. Many people have seen numerical problems when the
floating version of cddlib is used. As we all know, floating-point computation
might not give a correct answer, especially when an input data is very
sensitive to a small perturbation. When some strange behavior is observed, it
is always wise to create a rational approximation of the input (for example,
one can replace 0.3333333 with 1/3) and to compute it with cddlib compiled with
[GMP Rational](https://gmplib.org) to see what a correct behavior should be.
Whenever the time is not important, it is safer to use GMP rational arithmetic.

If you need speedy computation with floating-point arithmetic, you might want
to *play with* the constant `dd_almostzero` defined in [cdd.h](lib-src/cdd.h):

    #define dd_almostzero  1.0E-6

This number is used to recognize whether a number is zero: a number whose
absolute value is smaller than `dd_almostzero` is considered zero, and nonzero
otherwise. You might want to change this to modify the behavior of cddlib.
Another thing one can do is scaling. If the values in one column of an input is
of smaller magnitude than those in another column, scale one so that they
become comparable.

## Build the Latest Released Version

Download the most recent tarball from our [Releases
page](https://github.com/cddlib/cddlib/releases) and build cddlib with

```
tar zxf cddlib-*.tar.gz
cd cddlib-*
./configure
make
```

To install cddlib to `/usr/local` type
```
sudo make install
```

To link your own code with cddlib, invoke your compiler with something like
```
gcc -I/usr/local/include -L/usr/local/lib YOUR_CODE.c -lcdd
```

or, if you want to use rational arithmetic provided by
[GMP](https://gmplib.org), with
```
gcc -I/usr/local/include -L/usr/local/lib -DGMPRATIONAL YOUR_CODE.c -lcddgmp -lgmp
```

## Build from the Source Code Repository

Similar to the above, you can build from our [github repository](https://github.com/cddlib/cddlib) with
```
git clone https://github.com/cddlib/cddlib.git
cd cddlib
./bootstrap
./configure
make
```

## Contact the Author
The library cddlib is free software, but if cddlib turns out to be useful,
please kindly send a note or paper to the [author](./AUTHORS) mentioning what
purpose and how cddlib has been used for. The most powerful support for free
software development is user's appreciation.
