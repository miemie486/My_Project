#ifndef STRCHAR
#define STRCHAR

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include "def.h"

void print_int_vector( std::string* desc, size_t n, size_t* a );
void print_matrix( std::string* desc, size_t m, size_t n, double* a, size_t lda );

#endif /* STRCHAR */
