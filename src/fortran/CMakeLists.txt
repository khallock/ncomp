# Create a library called "Hello" which includes the source file "hello.cxx".
# The extension is already found. Any number of sources could be listed here.
add_library (ncomp_fort linint2.f)

# Make sure the compiler can find include files for our Hello library
# when other libraries or executables link to Hello
target_include_directories (ncomp_fort PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# create a Fortran static library
IF(NOT CMAKE_Fortran_CREATE_STATIC_LIBRARY)
  SET(CMAKE_Fortran_CREATE_STATIC_LIBRARY
      "<CMAKE_AR> cr <TARGET> <LINK_FLAGS> <OBJECTS> "
      "<CMAKE_RANLIB> <TARGET> ")
ENDIF(NOT CMAKE_Fortran_CREATE_STATIC_LIBRARY)

# compile a Fortran file into an object file
IF(NOT CMAKE_Fortran_COMPILE_OBJECT)
  SET(CMAKE_Fortran_COMPILE_OBJECT
    "<CMAKE_Fortran_COMPILER> -o <OBJECT> <DEFINES> <FLAGS> -c <SOURCE>")
ENDIF(NOT CMAKE_Fortran_COMPILE_OBJECT)
