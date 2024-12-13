**Fixed:**

* The cddlib.pc file now has correct linkage to cddlib (the non-gmp version).
  To link against the gmp version, you can use the new cddgmp.pc file.
  For instance::

    pkg-config cddlib --cflags --libs

  will fetch all compiler and link flags for the non-gmp version of cddlib,
  and::

    pkg-config cddgmp --cflags --libs

  will get them for the gmp version of cddlib.
