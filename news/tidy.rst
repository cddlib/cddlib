**Removed:**

* Remove the stale `dd_Remove` declaration; use `dd_MatrixRowRemove`.

**Fixed:**

* Make `./bootstrap` and `./configure` work without autoconf-archive (AX_*) macros.
* Fix out-of-tree builds so generated GMP headers are found during compilation.
* Only declare `ddd_*` helpers in `cddmp.h` for backends that implement them.
