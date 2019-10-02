Converting NCL functions to NComp functions
===========================================
1. Search `$NCLSRC/ni/src/lib/nfp/wrapper.c` for relevant `NclRegisterFunc` entry with desired NCL function name as the third argument, "linint2" for example:
```
    NclRegisterFunc(linint2_W,args,"linint2",nargs);
```

The first argument is the NCL "C wrapper" function that needs to be adapted to work in NComp.

2. Search for `linint2_W` in `$NCLSRC/ni/src/lib/nfp/*W.c` ("linint2W.c" in this case). Copy the function into a new file in NComp corresponding to the function name, `$NCOMP/src/wrappers/linint2.c` in this case. Note that this function will still need to be adapted and it not ready for inclusion in NComp at this point.

3. Determine which Fortran routines are needed, usually wrapped in an NGCALLF macro in a loop towards the end of the `_W` function.

4. Search `$NCLSRC/ni/src/lib/nfpfort/*.f` for necessary Fortran routines, and copy files to `$NCOMP/src/fortran` directory as needed. This Fortran code should not be modified for NComp. Add the new `.f` file to the `libncomp_fort_la_SOURCES` list in `$NCOMP/src/fortran/Makefile.am`.

5. Modify `$NCOMP/src/wrappers/funcname.c`
	a. Remove NCL-specific types
	b. Remove NCL stack parameter retrieval (`NclGetArgValue` calls) -- arrays should be passed to the function as pointers to `ncomp_array` structs.
    c . unbundle all elements of ncomp_array struct into local variables and pointers, keep same names as original `_W` function
	d. Change return value to be any error code set in Fortran (usually there's an int "ier" variable used for this).
	e. Remove NCL error handling (NhlPError), replace with printf and return non-zero.
    f. Translate utility function names from NCL to NComp functions as needed (see `$NCOMP/src/util.c`).

6. Add the new `.c` file to the `libncomp_wrappers_la_SOURCES` list in `$NCOMP/src/wrappers/Makefile.am`.

7. Add the new function signature to `$NCOMP/src/ncomp/wrapper.h`.

