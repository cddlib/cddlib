**Fixed:**

* cddlib is now thread-safe. We protect static variables with `_Thread_local` directives in the implementation, as provided by the C11 standard of the C language and supported by all supported recent compilers without adding any extra flags.

