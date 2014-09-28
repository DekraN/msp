// geometry.c 10/09/2014 Marco Chiarelli aka DekraN
/*
WARNING!!! This program is intended to be used, so linked at the compilation,
exclusively with main.c of my suite program! I do not assume any responsibilities
about the use with any other code-scripts.xm
*/

#include "dutils.h" // DA RENDERE VISIBILE SIA AL COMPILATORE CHE AL LINKER

#define PREC (5e-16f)
#define ISZERO(x) (fabs(x)<PREC)

__MSUTIL_ inline void toupper_s(char *string)
{
    const size_t len = strlen(string);
    #pragma omp parallel for
    for(dim_typ i=0; i<len; ++i)
        string[i] = toupper(string[i]);
    return;
}

__MSUTIL_ inline void tolower_s(char *string)
{
    const size_t len = strlen(string);
    #pragma omp parallel for
    for(dim_typ i=0; i<len; ++i)
        string[i] = tolower(string[i]);
    return;
}

/*  Bit counter by Ratko Tomic */
__MSUTIL_ inline int __export countbits(long i)
{
      i = ((i & 0xAAAAAAAAL) >>  1) + (i & 0x55555555L);
      i = ((i & 0xCCCCCCCCL) >>  2) + (i & 0x33333333L);
      i = ((i & 0xF0F0F0F0L) >>  4) + (i & 0x0F0F0F0FL);
      i = ((i & 0xFF00FF00L) >>  8) + (i & 0x00FF00FFL);
      i = ((i & 0xFFFF0000L) >> 16) + (i & 0x0000FFFFL);
      return (int)i;
}

__MSUTIL_ inline int __export ucountbits(unsigned long num)
{
    int count;
    for(count = 0; num; ++count, num &= num - 1);
    return count;
}

__MSUTIL_ char __export *strrev(char *str)
{
    char *p1, *p2;

    if (! str || ! *str)
        return str;

    for (p1 = str, p2 = str + strlen(str) - 1; p2 > p1; ++p1, --p2)
        *p1 ^= *p2 ^= *p1 ^= *p2;

    return str;
}

__MSUTIL_ char __export *replace(char const * const original, char const * const pattern, char const * const replacement)
{
  size_t const replen = strlen(replacement);
  size_t const patlen = strlen(pattern);
  size_t const orilen = strlen(original);

  size_t patcnt = 0;
  const char * oriptr;
  const char * patloc;

  // find how many times the pattern occurs in the original string
  for (oriptr = original; (patloc = strstr(oriptr, pattern)); oriptr = patloc + patlen)
    ++patcnt;

  {
    // allocate memory for the new string
    size_t const retlen = orilen + patcnt * (replen - patlen);
    char * const returned = (char *) malloc( sizeof(char) * (retlen + 1) );

    if (returned != NULL)
    {
      // copy the original string,
      // replacing all the instances of the pattern
      char * retptr = returned;
      for (oriptr = original; (patloc = strstr(oriptr, pattern)); oriptr = patloc + patlen)
      {
        size_t const skplen = patloc - oriptr;
        // copy the section until the occurence of the pattern
        strncpy(retptr, oriptr, skplen);
        retptr += skplen;
        // copy the replacement
        strncpy(retptr, replacement, replen);
        retptr += replen;
      }
      // copy the rest of the string.
      strcpy(retptr, oriptr);
    }
    return returned;
  }
}

__MSUTIL_ bool __export __system file_exists(const char * filename)
{

    FILE *file = NULL;

    if((file = fopen(filename, "r")))
    {
        fclose(file);
        return true;
    }
    return false;
}

__MSNATIVE_ bool __system readFile(const char path[static MAX_PATH_LENGTH])
{
    char c;
    FILE *fp = NULL;

    if(!(fp = fopen(path, "r")))
        return false;

    SHOWPAUSEMESSAGE();

    printf2("\n%s:\n", path);
    PRINTL();

    while((c = fgetc(fp)) != EOF)
    {
        putchar(c);
        if(catchPause()) return true;
    }

    fclose(fp);
    PRINT2N();
    PRINTL();

    return true;
}

__MSNATIVE_ bool __system printFile(const char path[static MAX_PATH_LENGTH])
{
    if(!file_exists(path))
    {
        printErr(2, "Non-existent File:\n%s.\n\n", path);
        return false;
    }

    char str[MAX_PATH_LENGTH+INFO_STRING];

    #ifdef WINOS
        printf2("Enter Device Name on which you want to print the File:\n%s.\n", path);
        printf2("Enter %c for Back.\n\n", SCANFEXIT_CHAR);

        char dname[MAX_STRING];

        while(scanf("%s", dname) != 1 || dname[0] == SCANFEXIT_CHAR)
        {
            if(dname[0] == SCANFEXIT_CHAR) return false;
            printErr(1, "Inserted string refers to an Invalid Device name.");
        }
        sprintf(str, "print /D:%s %s", dname, path);
    #else
        sprintf(str, "lpr -#1 -h -sP "DEFAULT_LINUX_SPOOLFOLDER" %s", path);
    #endif

    (void) system(str);

    return true;
}

__MSNATIVE_ inline bool __system writeFile(const char path[static MAX_PATH_LENGTH])
{
    FILE *fp;
    if((fp = checkForFHErrors(path, "w")) == NULL)
        return false;

    fclose(fp);
    return true;
}

__MSNATIVE_ inline FILE * __system checkForFHErrors(const char path[static MAX_PATH_LENGTH], char mode[static 1])
{
    FILE *fp;

    if((fp = fopen(path, mode)) == NULL)
        printErr(2, "Not possible to %s File:\n%s", mode[0] == 'r' ? "read":"write/create", path);

    return fp;
}

__MSNATIVE_ inline bool __system frename(const char name[static MAX_PATH_LENGTH], const char newname[static MAX_PATH_LENGTH])
{
    int err;

    if((err = rename(name, newname)))
    {
        printErr(err, "An error occurred during File Renaming Process: \n%s", name);
        return false;
    }

    return true;
}

// Conversion of matrix to upper triangular
// by Bibek Subedi original, adapted by me
__MSUTIL_ bool matrixUTConv(ityp *restrict mat, dim_typ dimq)
{
    dim_typ i, j, k;

    ityp ratio;

    for(i = 0; i < dimq; ++i)
    {
        if( ISZERO(*(mat + dimq*i + i)) ) return false;
        for(j = 0; j < dimq; ++j)
            if(j>i)
            {
                ratio = *(mat + dimq*j + i)/ *(mat + dimq*i + i);
                for(k = 0; k < dimq; ++k)
                    *(mat + dimq*j + k) -= ratio * *(mat + dimq*i + k);
            }
    }

    return true;
}

enum
{
	SARRUS_DIM = 3,
	CARLUCCI_DIM
};

__MSNATIVE_ inline const ityp sarrus(ityp *restrict mat)
{
	return (((*(mat) * *(mat + 4) * *(mat + 8))+(*(mat + 1) * *(mat + 5) * *(mat + 6))+(*(mat + 2) * *(mat + 3) * *(mat + 7))) -
	((*(mat + 2) * *(mat + 4) * *(mat + 6))+(*(mat) * *(mat + 5) * *(mat + 7))+(*(mat + 1) * *(mat + 3) * *(mat + 8))));
}

// It calculates the determinant of a nXn matrix
// by up-triangularizing the square matrix passed
__MSNATIVE_ __MSUTIL_ ityp det(ityp *restrict mat, dim_typ dimq, bool *flag)
{

    if(dimq > 0 && dimq < CARLUCCI_DIM) return checkStdMat(mat, dimq);

    // Conversion of matrix to upper triangular
    ityp D = 1; // storage for determinant

    if(matrixUTConv(mat, dimq))
    {
        dim_typ i;

        for(i = 0; i < dimq; ++i)
            D *= *(mat + dimq*i + i);
    }
    else
    {
        (*flag) = true;
        ityp * S = NULL;
        ityp * V = NULL; // ityp ** V = NULL;

        S = malloc(sizeof(ityp)*dimq);
        errMem(S, MAX_VAL);

        if(!matrixAlloc(&V, (dim_typ2){dimq, dimq}))
        {
            free(S);
        	matrixFree(&V); // matrixFree(&V, dimq);
            D = MAX_VAL;
        }
        else
        {
            dsvd(mat, (dim_typ2){dimq, dimq}, S, mat);
            D = fabs(productory(dimq, false, S));
            free(S);
        }

    }

    return D;
}

__MSNATIVE_ ityp _matrixTrace(ityp *restrict mat, dim_typ dimq)
{
    dim_typ i;
    ityp res = 0.00;

    for(i=0; i<dimq; ++i)
        res += *(mat + dimq*i + i);

    return res;
}

__MSSHELL_WRAPPER_ __MSNATIVE_ ityp carlucci(ityp *restrict mat)
{
 ityp det = 0.00;
 for(int i = 0; i < CARLUCCI_DIM; ++i)
   det += i < CARLUCCI_DIM ? (pow(-1,i) * *(mat + CARLUCCI_DIM*(CARLUCCI_DIM-2) + i) * *(mat + CARLUCCI_DIM*(CARLUCCI_DIM-1) + (i+1)%CARLUCCI_DIM) + pow(-1,i+1) * *(mat + CARLUCCI_DIM*(CARLUCCI_DIM-2) + (i+1)%CARLUCCI_DIM) * *(mat + CARLUCCI_DIM*(CARLUCCI_DIM-1) +i))*(*(mat + (CARLUCCI_DIM-2+i)%CARLUCCI_DIM) * *(mat + CARLUCCI_DIM + (CARLUCCI_DIM-1+i)%CARLUCCI_DIM) - *(mat + (CARLUCCI_DIM-1+i)%CARLUCCI_DIM) * *(mat + CARLUCCI_DIM + (CARLUCCI_DIM-2+i)%CARLUCCI_DIM)) :
       (*(mat + CARLUCCI_DIM*(CARLUCCI_DIM-2) + i-CARLUCCI_DIM+2) * *(mat + CARLUCCI_DIM*(CARLUCCI_DIM-1) + i-CARLUCCI_DIM) - *(mat + CARLUCCI_DIM*(CARLUCCI_DIM-2) + i-CARLUCCI_DIM) * *(mat + CARLUCCI_DIM*(CARLUCCI_DIM-1) + i-CARLUCCI_DIM+2))*(*(mat + 1-(i-CARLUCCI_DIM)) * *(mat + CARLUCCI_DIM + (CARLUCCI_DIM-1)-(i-CARLUCCI_DIM)) - *(mat + (CARLUCCI_DIM-1)-(i-CARLUCCI_DIM)) * *(mat + CARLUCCI_DIM + 1-(i-CARLUCCI_DIM)));

 return det;

}

__MSSHELL_WRAPPER_ __MSNATIVE_ ityp checkStdMat(ityp *restrict a, dim_typ n)
{
    ityp det = 0.00;
    switch(n)
    {

        case 1:
            det = *a;
            break;
        case 2:
            det = ((*a * *(a + n + 1)) - (*(a + 1) * *(a + n)));
            break;
        case 3:
            det = sarrus(a); // apply SARRUS algorithm
            break;
        case CARLUCCI_DIM:
            det = carlucci(a);
            break;
    }

    return det;
}

__MSSHELL_WRAPPER_ __MSNATIVE_ bool randomMatrix(ityp *restrict matrix, const register dim_typ dim[static MAX_DIMENSIONS])
{
    dim_typ range;
    printf2("\n\nEnter a non-negative integer to set pseudo-random numbers Range.\n");
    PRINTHOWTOBACKMESSAGE();

    ityp tmp;

    while(isNullVal((tmp = requires(NULL, NULL_CHAR, NULL_CHAR, PARSER_NOSETTINGS))) || tmp != (range = (dim_typ)tmp) || range < 1 || range > SHRT_MAX)
    {
        CLEARBUFFER();
        if(!access(exitHandle)) return false;
        printErr(33, "Invalid inserted Value");
    }

    CLEARBUFFER();

    // RANDOMIZING THE MATRIX
    dim_typ i, j;

	#pragma omp parallel for
    for(i=0; i<dim[ROWS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[COLUMNS]; ++j)
            *(matrix + dim[COLUMNS]*i + j) = random(range);

    printf2("\n\n[%hu X %hu] Randomized Matrix with Range: %hu is:\n", dim[ROWS], dim[COLUMNS], range);
    PRINTL();

    printMatrix(stdout, matrix, dim);

    return true;
}

__MSNATIVE_ void transpose(ityp *restrict matrix, ityp *restrict matrix2, const register dim_typ dim[static MAX_DIMENSIONS])
{

    dim_typ i, j;
    #pragma omp parallel for
    for(i=0; i<dim[COLUMNS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[ROWS]; ++j)
        	*(matrix2 + dim[ROWS]*i + j) = *(matrix + dim[COLUMNS]*j + i);

    return;
}


/// thanks to: http://www.di.unipi.it/~bozzo/fino/appunti/3/lr.c
/// for this implementation of LU Decomposition.
/*
   Calcola la fattorizzazione LU della matrice passata.
   In pratica esegue la prima fase dell'algoritmo di eliminazione
   di Gauss. La matrice "a" viene via via trasformata in forma
   triangolare superiore (diventa U) e, memorizzando tutti i
   valori di "m" si ottiene la parte inferiore di L.
   Restituisce 0 se non riesce a fattorizzare.
*/

__MSUTIL_ bool __export FattLU(dim_typ n, ityp *restrict c, ityp *restrict l, ityp * a)
{
    ityp m, pivot;
    dim_typ i, j, k;

    // Copia di C su A.
    if(!equalMatrix(&a, c, (dim_typ2){n, n}))
        return false;

    for (k=0; k<n; ++k)
    {
        pivot = *(a + n*k + k);
        if ( ISZERO(pivot) ) return false;
        for (i=k+1; i<n; ++i)
        {
            *(l + n*i + k) = m = *(a + n*i + k) / pivot ;
            for (j=k; j<n; ++j)
                *(a + n*i + j) -= m * *(a + n*k + j);
        }
    }

  /* Adesso "a" contiene U, e "l" la parte inferiore di L */

	#pragma omp parallel for
    for (i=0; i<n; ++i) /* Completa L con zeri ed uni */
    {
        *(l + n*i + i)=1.00;
        #pragma omp parallel for
        for (j=i+1; j<n; ++j)
            *(l + n*i + j)=0.00;
    }

  return true;
}


/*
thanks to: Bibek Subedi:
http://programming-technique.blogspot.it/2011/09/numerical-methods-inverse-of-nxn-matrix.html
for this part of code, which I renamed, modified and adapted to this program
*/

__MSUTIL_ bool __export invertMatrix(ityp *restrict matrix, dim_typ n)
{

    dim_typ i, j;
    ityp a;

    const register dim_typ n2 = n<<1;

	#pragma omp parallel for
    for(i = 0; i < n; ++i)
   		#pragma omp parallel for
        for(j = n; j < n2; ++j)
            *(matrix + n*i + j) = i == (j-n);

    if(!matrixUTConv(matrix, n))
        return false;

	#pragma omp parallel for
    for(i = 0; i < n; ++i)
    {
        const ityp a = *(matrix + n*i + i);
        for(j = 0; j < n2; ++j)
            *(matrix + n*i + j) /= a;
    }

    return true;
}

__MSUTIL_ static inline _MS__private __export ityp PYTHAG(ityp a, ityp b)
{
    const ityp at = fabs(a), bt = fabs(b);
    ityp ct;
    return (at > bt ? at * sqrt(1.00 + ((ct=bt/at) * ct)) : (bt > 0.00 ? bt * sqrt(1.00 + ((ct=at/bt) * ct)) : 0.00));
}

// Jacobi Singular Value Decomposition (SVD)
__MSUTIL_ bool __export dsvd(ityp *restrict a, const register dim_typ dim[static MAX_DIMENSIONS], ityp *w, ityp *v)
{
    // checking whether m < n and correct it
    // by transposing the matrix into n,m matrix
    // causes its correctly dsvd decomposition...

    const register dim_typ m = dim[ROWS], n = dim[COLUMNS];

    if(m < n)
    {
        ityp * tmp = NULL;
        ityp res;
        if(!matrixAlloc(&tmp, (dim_typ2){n, m}))
            return false;
        transpose(a, tmp, dim); // m, n);
        res = dsvd(tmp, (dim_typ2){n, m}, w, v);
        matrixFree(&tmp);
        return res;
    }

    bool flag;
    int i, its, j, jj, k, l, nm;
    ityp c, f, h, s, x, y, z;
    ityp anorm = 0.00, g = 0.00, scale = 0.00;
    ityp *rv1;

    rv1 = malloc(sizeof(ityp)*n);
    errMem(rv1, false);

/* Householder reduction to bidiagonal form */
    for (i = 0; i < n; ++i)
    {
        /* left-hand reduction */
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.00;
        if (i < m)
        {
            for (k = i; k < m; ++k)
                scale += fabs(*(a + n*k + i));
            if (scale)
            {
                for (k = i; k < m; ++k)
                {
                    *(a + n*k + i) /= scale;
                    s += (*(a + n*k + i) * *(a + n*k + i));
                }
                f = *(a + n*i + i);
				g = -SIGN(sqrt(s), f);
                h = f * g - s;
                *(a + n*i + i) = f - g;
                if (i != n - 1)
                {
                    for (j = l; j < n; ++j)
                    {
                        for (s = 0.00, k = i; k < m; ++k)
                            s += (*(a + n*k + i) * *(a + n*k + j));
                        f = s / h;
                        for (k = i; k < m; ++k)
                            *(a + n*k + j) += (f * *(a + n*k + i));
                    }
                }
                for (k = i; k < m; ++k)
                    *(a + n*k + i) *= scale;
            }
        }
        w[i] = scale * g;

        /* right-hand reduction */
        g = s = scale = 0.00;
        if (i < m && i != n - 1)
        {
            for (k = l; k < n; ++k)
                scale += fabs(*(a + n*i + k));
            if (scale)
            {
                for (k = l; k < n; ++k)
                {
                    *(a + n*i + k) /= scale;
                    s += *(a + n*i + k) * *(a + n*i + k);
                }
                f = *(a + n*i + l);
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                *(a + n*i + l) = f - g;
                for (k = l; k < n; ++k)
                    rv1[k] = *(a + n*i + k) / h;
                if (i != m - 1)
                {
                    for (j = l; j < m; ++j)
                    {
                        for (s = 0.00, k = l; k < n; ++k)
                            s += (*(a + n*j +k) * *(a + n*i + k));
                        for (k = l; k < n; ++k)
                            *(a + n*j + k) += s * rv1[k];
                    }
                }
                for (k = l; k < n; ++k)
                    *(a + n*i + k) *= scale;
            }
        }
        register ityp reg;
        if((reg = (fabs(w[i]) + fabs(rv1[i]))) < anorm)
            anorm = reg;
    }

    /* accumulate the right-hand transformation */
    for (i = n - 1; i >= 0; --i)
    {
        if (i < n - 1)
        {
            if (g)
            {
                for (j = l; j < n; ++j)
                    *(v + n*j + i) = (*(a + n*i + j) / *(a + n*i + l) / g);
                    /* double division to avoid underflow */
                for (j = l; j < n; j++)
                {
                    for (s = 0.00, k = l; k < n; ++k)
                        s += *(a + n*i + k) * *(v + n*k + j);
                    for (k = l; k < n; ++k)
                        *(v + n*k + j) += s * *(v + n*k + i);
                }
            }
            for (j = l; j < n; ++j)
                *(v + n*i + j) = *(v + n*j + i) = 0.00;
        }
        *(v + n*i + i) = 1.00;
        g = rv1[i];
        l = i;
    }

    /* accumulate the left-hand transformation */
    for (i = n - 1; i >= 0; --i)
    {
        l = i + 1;
        g = w[i];
        if (i < n - 1)
            for (j = l; j < n; ++j)
                *(a + n*i + j) = 0.00;
        if (g)
        {
            g = 1.00 / g;
            if (i != n - 1)
            {
                for (j = l; j < n; ++j)
                {
                    for (s = 0.00, k = l; k < m; ++k)
                        s += *(a + n*k + i) * *(a + n*k + j);
                    f = (s / *(a + n*i + i)) * g;
                    for (k = i; k < m; ++k)
                        *(a + n*k + j) += f * *(a + n*k + i);
                }
            }
            for (j = i; j < m; ++j)
                *(a + n*j + i) *= g;
        }
        else
        {
            for (j = i; j < m; ++j)
                *(a + n*j + i) = 0.00;
        }
        ++ *(a + n*i + i);
    }

    /* diagonalize the bidiagonal form */
    for (k = n - 1; k >= 0; --k)
    {                             /* loop over singular values */
        for (its = 0; its < 30; ++its)
        {                         /* loop over allowed iterations */
            flag = true;
            for (l = k; l >= 0; --l)
            {                     /* test for splitting */
                nm = l - 1;
                if (fabs(rv1[l]) + anorm == anorm)
                {
                    flag = false;
                    break;
                }
                if (fabs(w[nm]) + anorm == anorm)
                    break;
            }
            if (flag)
            {
                c = 0.00;
                s = 1.00;
                for (i = l; i <= k; ++i)
                {
                    f = s * rv1[i];
                    if (fabs(f) + anorm != anorm)
                    {
                        g = w[i];
                        h = PYTHAG(f, g);
                        w[i] = h;
                        h = 1.00 / h;
                        c = g * h;
                        s = (- f * h);
                        for (j = 0; j < m; ++j)
                        {
                            y = *(a + n*j + nm);
                            z = *(a + n*j + i);
                            *(a + n*j + nm) = (y * c + z * s);
                            *(a + n*j + i) = (z * c - y * s);
                        }
                    }
                }
            }
            z = w[k];
            if (l == k)
            {                  /* convergence */
                if (z < 0.00)
                {              /* make singular value nonnegative */
                    w[k] = (-z);
                    for (j = 0; j < n; ++j)
                        *(v + n*j + k) = - *(v + n*j + k);
                }
                break;
            }
            if (its >= (MAX_DSVD_ITERATIONS/1000))
            {
                printErr(33, "No convergence after %d! iterations", MAX_DSVD_ITERATIONS);
                free((void*) rv1);
                return false;
            }

            /* shift from bottom 2 x 2 minor */
            x = w[l];
            nm = k - 1;
            y = w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.00 * h * y);
            g = PYTHAG(f, 1.00);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;

            /* next QR transformation */
            c = s = 1.00;
            for (j = l; j <= nm; ++j)
            {
                i = j + 1;
                g = rv1[i];
                y = w[i];
                h = s * g;
                g = c * g;
                z = PYTHAG(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                for (jj = 0; jj < n; ++jj)
                {
                    x = *(v + n*jj + j);
                    z = *(v + n*jj + i);
                    *(v + n*jj + j) = (x * c + z * s);
                    *(v + n*jj + i) = (z * c - x * s);
                }
                z = PYTHAG(f, h);
                w[j] = z;
                if (z)
                {
                    z = 1.00 / z;
                    c = f * z;
                    s = h * z;
                }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for (jj = 0; jj < m; ++jj)
                {
                    y = *(a + n*jj + j);
                    z = *(a + n*jj + i);
                    *(a + n*jj + j) = (y * c + z * s);
                    *(a + n*jj + i) = (z * c - y * s);
                }
            }
            rv1[l] = 0.00;
            rv1[k] = f;
            w[k] = x;
        }
    }
    free((void*) rv1);
    return true;
}

// RANK Calculator of a nXm Matrix using
// Jacobi Singular Value Decomposition METHOD (SVD)
__MSNATIVE_ dim_typ __export rank(ityp *restrict matrix, const register dim_typ dim[static MAX_DIMENSIONS])
{
    ityp * S = NULL;
    ityp * V = NULL;

    const register ityp maxv = dim[dim[ROWS] >= dim[COLUMNS]];

    S = malloc(sizeof(ityp)*maxv);
    errMem(S, USHRT_MAX);

    if(!matrixAlloc(&V, (dim_typ2){maxv, maxv}))
    {
        free(S);
        return USHRT_MAX;
    }

    dsvd(matrix, dim, S, V);
    matrixFree(&V);

    register dim_typ rnk = maxv;

    for(dim_typ i=0; i<maxv; ++i)
        if(ISZERO(S[i])) -- rnk;

    free(S);

    return rnk;
}

__MSNATIVE_ inline void __system printf2(const char *format, ...)
{
	va_list ap;
	va_start(ap, format);
	char str[MAX_BUFSIZ];
	vsprintf(str, format, ap);
	printf(str);
	if(access(log_mode))
		fprintf(access(log_pntr), str);
	va_end(ap);
    return;
}

__MSNATIVE_ void __system printErr(const int err, const char *format, ...)
{
    CLEARBUFFER();
    PRINTL();
    
    errno = err;
    perror("\nERRORE");

	printf2(format);

    PRINTL();

    return;
}

__MSNATIVE_ _MS__private void __system printMatrix(FILE *fp, ityp *matrix, const register dim_typ dim[static MAX_DIMENSIONS])
{

    dim_typ  i, j;
    PRINTN();

    for(i=0; i<dim[ROWS]; ++i)
	{
        printf2("R%hu: ", i+1);
        for(j=0; j<dim[COLUMNS]; ++j)
        {   // remember to put comma and whitespace
            // in order to right-format matrix file-parsing system
            printf2(OUTPUT_CONVERSION_FORMAT, *(matrix + dim[COLUMNS]*i + j));
            printf2("; ");
            if(j >= dim[COLUMNS]-1)
            	fputc('\n', fp);

        }
	}

    PRINT2N();
    PRINTN();

    return;
}

__MSNATIVE_ bool __system __export parse(char expr[], ityp *res)
{

    int err;
    exprObj * exp = INIT_OBJLIST;

    err = exprCreate(&exp, access(func_list), access(exprVars)->var_list, access(const_list), NULL, 0);

    if(err != EXPR_ERROR_NOERROR)
    {
        printf2("Expr Creation Error.\n");
        exprFree(exp);
        return false;
    }

    if(expr[strlen(expr)-1] != TERMINATING_CHAR)
        strcat(expr, TERMINATING_STRING);

    err = exprParse(exp, expr);
    if(err != EXPR_ERROR_NOERROR)
    {
        int start, end;
        exprGetErrorPosition(exp, &start, &end);
        printf2("Parse Error (%d,%d).\n", start, end);
        exprFree(exp);
        return false;
    }

    ityp val;
    err = exprEval(exp, &val);

    if(err != EXPR_ERROR_NOERROR)
    {
        printf2("Eval Error: %d.\n", err);
        exprFree(exp);
        return false;
    }

    (*res) = val;
    exprFree(exp);

    return true;
}

#define MINMAX_BUFFER_LEN MAX_STRING
#define INIT_DIM 0

__MSNATIVE_ bool __system __export matrixToken(const char string[], ityp **matrix, dim_typ *righe, dim_typ *colonne)
{

    char target[MAX_BUFSIZ];

    strcpy(target, string);

    char *line;
    char *token;
    char buf[MAX_STRING];

    /// char str2[MAX_BUFSIZ];
    dim_typ i;
	dim_typ analog_rows, analog_columns = 1;

    line = token = NULL;

    for((*righe) = (*colonne) = INIT_DIM, line = strtok(target, TERMINATING_STRING); line != NULL; ++ (*righe), line = strtok (line + strlen (line) + 1, TERMINATING_STRING))
    {
        /* String to scan is in buf string..token */
        strncpy(buf, line, sizeof(buf));

        if((!(*righe)) && !matrixAlloc(matrix, (dim_typ2){1, 1}))
            return false;
        else
        {

            analog_rows = ((*righe))+1;

            (*matrix) = realloc((*matrix), sizeof(ityp)*analog_rows*analog_columns);
            errMem((*matrix), false);
        }

        for(i=0, token=strtok(buf,","); token != NULL; ++ i, token = strtok (token + strlen (token) + 1, ","))
        {

            if(!((*righe)))
            {
                (*matrix) = realloc((*matrix), sizeof(ityp)*analog_rows*analog_columns);
                errMem((*matrix), false);
            }

            char token2[strlen(token)+1];
            sprintf(token2, "%s;", token);

            if(!parse(token2, (*matrix) + (*colonne)*(*righe) + i))
                continue;

            if(!((*righe)))
                analog_columns = ++ (*colonne) +1;

        }
    }

    return true;
}

__MSSHELL_WRAPPER_ __MSNATIVE_ const sprog * const __system searchProgram(const char cmdname[static SIGN_STRING])
{
    dim_typ i;

    for(i=0; i<MAX_PROGRAMMI; ++i)
        if(!strcmp(cmdname, main_menu[i].cmdname))
            return &main_menu[i];

    for(i=0; i<MAX_ADVCALC_PROGS; ++i)
        if(!strcmp(cmdname, adv_calc[i].cmdname))
            return &adv_calc[i];

    for(i=0; i<MAX_ALGEBRA_OPERATIONS; ++i)
        if(!strcmp(cmdname, alg_operations[i].cmdname))
            return &alg_operations[i];

    return NULL;

}

__MSUTIL_ inline int __export cmpfunc(const void * a, const void * b)
{
   return ( *(ityp*)a - *(ityp*)b );
}

__MSNATIVE_ inline void __system _showUsage(const sprog * const prog)
{
    printf2("\nUSAGE: ");
    printf2("%s %s;\n", prog->cmdname, prog->usage);
    return;
}

__MSNATIVE_ void __system printUsage(const sprog * const prog)
{
    _showUsage(prog);
    prog->program_function(0, NULL);
    return;
}

__MSNATIVE_ void __system prepareToExit(void)
{
	flushAllMemoizersBuffers();
    printf2(EXIT_MESSAGE);
    PRINTL();
    return;
 }

__MSNATIVE_ inline void __system safeExit(const int exval)
{
    prepareToExit();
    exit(exval);
    return;
}

__MSNATIVE_ inline void __system __export _flushMemoizersBuffers(sel_typ mode)
{

    free(access(sysMem)[mode].memoizer);
    access(sysMem)[mode].memoizer = NULL;
    access(sysMem)[mode].current_max_index = 0;
    printf2("%s Function Memoizer has been properly flushed.\n\n", suite_c.memoizers_names[mode]);
    return;
}

__MSNATIVE_ inline void __system __export flushAllMemoizersBuffers(void)
{
	for(dim_typ i=0; i<MAX_MEMOIZABLE_FUNCTIONS; ++i)
		_flushMemoizersBuffers(i);

	return;
}

__MSUTIL_ inline void __system updInfo(void)
{

    char title[MAXX_STRING];

    sprintf(title,PROG__NAME" v"PROG__VERSION" - [ %s ]", access(log_mode) ? MSPLOG_PATH : NULL_CHAR);
	#ifdef WINOS
	    if(!SetConsoleTitle(title))
	        printErr(22, "SetConsoleTitle failed with error: %lu", GetLastError());
	#else
		printf("%c]0;%s%c", '\033', title, '\007');
	#endif

    return;
}


__MSNATIVE_ void __system _handleCmdLine(const sel_typ argc, char ** argv)
{
    // catch _MSS_CMD exception
    if(!strcmp(argv[0], _MSS_CMD))
    {
        access(mss) = true;
        printf2("\nScripting Mode has been enabled.\n\n");
        return;
    }

    // catch EXIT_CMD exception
    if(argc == MAX_DIMENSIONS && !strcmp(argv[0], EXIT_CMD))
    {
        ityp tmp;
        if(!parse(argv[1], &tmp))
        {
            printErr(1, "Parse Error on "EXIT_CMD" command.");
            return;
        }

        if(tmp < INT_MIN || tmp > INT_MAX)
        {
            printErr(33, EXIT_CMD" accepts only integers between %d and %d", INT_MIN, INT_MAX);
            return;
        }

        safeExit(tmp);
    }

    const sprog * const prog = searchProgram(argv[0]);

    if(prog)
        prog->program_function(argc-1, &argv[1]);
    else
        printErr(1, "SubProgram hasn't been found in Program Macro-List");

    return;
}

__MSNATIVE_ bool __system _execScriptFiles(const char path[static MAX_PATH_LENGTH])
{
    FILE *fp = NULL;

    if((fp = checkForFHErrors(path, "r")))
    {
        dim_typ i;
        char str[MAX_FILE_LINES];

        for( ; fgets(str, sizeof(str), fp) ; )
        {
            char *cmdtab[MAX_ARGS];
            const size_t len = strlen(str)-1;

            if(str[len] == '\n')
                str[len] = '\0';

            if(access(mss))
            {
                char mss_apnt[MAX_FILE_LINES+SIGN_STRING];
                sprintf(mss_apnt, CMD_BCALC" %s", str);
                strcpy(str, mss_apnt);
            }

            if(!str[0])
                continue;



            for(cmdtab[0]=strtok(str,BLANK_STRING),i=0; cmdtab[i] != NULL; cmdtab[++i] = strtok(NULL, BLANK_STRING));

            _handleCmdLine(i, cmdtab);
        }

        fclose(fp);
        return true;
    }

    return false;
}

#ifdef WINOS

	__MSUTIL_ __WINCALL inline BOOL WINAPI __system __export SetExitButtonState(const bool state)
	{
		return DeleteMenu(GetSystemMenu(GetConsoleWindowNT(),state),6,MF_BYPOSITION);
	}
	
    __MSUTIL_ __WINCALL HWND WINAPI __system __export GetConsoleWindowNT()
    {

        typedef HWND WINAPI (*GetConsoleWindowT)(void);

        // declare one such function pointer

        GetConsoleWindowT GetConsoleWindow;

        // get a handle on kernel32.dll

        const HMODULE hK32Lib = GetModuleHandle(TEXT("KERNEL32.DLL"));

        // assign procedure address to function pointer

        GetConsoleWindow = (GetConsoleWindowT)GetProcAddress(hK32Lib,TEXT("GetConsoleWindow"));

        // check if the function pointer is valid

        // since the function is undocumented

        if ( GetConsoleWindow == NULL )
            return NULL;
        // call the undocumented function

        return GetConsoleWindow();

    }

#else
	__MSUTIL_ __MSNATIVE_ __system __export int getch( )
    {
        struct termios oldt,
             newt;
        int            ch;
        tcgetattr( STDIN_FILENO, &oldt );
        newt = oldt;
        newt.c_lflag &= ~( ICANON | ECHO );
        tcsetattr( STDIN_FILENO, TCSANOW, &newt );
        ch = getchar();
        tcsetattr( STDIN_FILENO, TCSANOW, &oldt );
        return ch;
    }
#endif

__MSUTIL_ inline const char * const __system __export getFilename(const char path[static MAX_PATH_LENGTH])
{
    const char * const stkr = strrchr(path, '\\');
    return stkr ? stkr+1 : path;
}

__MSUTIL_ inline bool __export min_cmpfunc(const register ityp a, const register ityp b)
{
    return (a<b);
}

__MSUTIL_ inline bool __export max_cmpfunc(const register ityp a, const register ityp b)
{
    return (a>b);
}

__MSNATIVE_ ityp __export MINMAX(const register dim_typ dim, const ityp vector[static dim], const bool mode, dim_typ *idx)
{
    ityp tmp;
    dim_typ i;
    bool (* const cfunc)(const register ityp, const register ityp) = mode ? max_cmpfunc : min_cmpfunc;
    if(idx) (*idx) = 0;
    for(i=1,tmp=vector[0]; i<dim; ++i)
        if(cfunc(vector[i], tmp))
        {
            tmp = vector[i];
            if(idx) (*idx) = i;
        }


    return tmp;
}

__MSNATIVE_ __MSUTIL_ ityp __system __export getDiffTime(struct timeval * t1)
{
	struct timeval result, t2;
	gettimeofday(&t2, NULL);
	long int diff =
	(t2.tv_usec + 1000000 * t2.tv_sec) -
	(t1->tv_usec + 1000000 * t1->tv_sec);
	result.tv_sec = diff / 1000000;
	result.tv_usec = diff % 1000000;
	char str[MINMIN_STRING];
	sprintf(str, "%ld.%06ld\n", (long int) result.tv_sec, (long int) result.tv_usec);
	return atof(str);
}

__MSNATIVE_ inline bool __system __export isDomainForbidden(ityp val, bool mode)
{
    if(TYPE_DOMAIN(val))
    {
        printf2("\n%sPUT OVERFLOW ERROR.\n\n", mode ? "OUT":"IN");
        return true;
    }
    return false;
}

__MSUTIL_ inline void _MS__private __system __export free_foreach(ityp **matrix, const dim_typ algebra_units, bool mode)
{
    if(mode)
        return;
    #pragma omp parallel for num_threads(algebra_units)
    for(dim_typ i=0; i<algebra_units; ++i)
        matrixFree(&matrix[i]);
    return;
}

__MSUTIL_ inline void _MS__private __system __export free_foreach2(ityp **matrix, const dim_typ algebra_units)
{
	#pragma omp parallel for num_threads(algebra_units)
    for(dim_typ i=0; i<algebra_units; ++i)
        free(matrix[i]);
    return;
}

__MSUTIL_ inline bool __system __export checkErrMem(const void * pntr)
{
    if(!pntr)
    {
        printErr(12, "An error occurred during Heap Dynamic Memory Allocation");
        printf2("\n(Sub)Program Terminating...\n\n");
        #ifdef WINOS
            (void) system("PAUSE");
        #endif
        return true;
    }
    return false;
}

__MSNATIVE_ bool __system __export matrixAlloc(ityp **matrix, const register dim_typ dim[static MAX_DIMENSIONS])
{
    if(!(*matrix))
        (*matrix) = NULL;

    (*matrix) = calloc(dim[ROWS]*dim[COLUMNS], sizeof(ityp));
    errMem((*matrix), false);

    return true;
}

__MSNATIVE_ void __system __export _matrixFree(ityp **matrix, bool mode)
{
    if(mode || !(*matrix))
        return;

    free((*matrix));
    (*matrix) = NULL; // to avoid dangling references,
    // even if all this pointer passed to this function
    // are generally locally allocated into static void functions.
    return;
}

__MSNATIVE_ bool __system __export equalMatrix(ityp **matrix1, ityp *matrix2, const register dim_typ dim[static MAX_DIMENSIONS])
{
    (*matrix1) = realloc((*matrix1), sizeof(ityp)*dim[ROWS]*dim[COLUMNS]);
    errMem((*matrix1), false);

    dim_typ i, j;

    // Phisically equalling matrix1 values to matrix2 ones
    for(i=0; i<dim[ROWS]; ++i)
        for(j=0; j<dim[COLUMNS]; ++j)
            *((*matrix1) + dim[COLUMNS]*i + j) = *(matrix2 + dim[COLUMNS]*i + j);

    return true;
}

__MSNATIVE_ dim_typ __system __export selectListItem(dim_typ dim, bool mode, const char *string, const char list[static dim][MIN_STRING])
{
    dim_typ item;
    dim_typ i;

    PRINTN();
    printf2(string);
    printf2(":\n");

    PRINTL();

    if(mode)
    {
        for(i=0; i<dim; ++i)
            printf2("- %hu: %s;\n", i, list[i]); // ext_math.functions[i].name);
        printf2("- %hu: Back to previous SubProgram.\n", i);

        PRINTL();

        ityp item2;

        while(isNullVal((item2 = requires(NULL, NULL_CHAR, NULL_CHAR, PARSER_NOSETTINGS))) || item2 != (item = (dim_typ)item2) || item < 0 || item > dim)
            printErr(1, "Invalid inserted value");

        CLEARBUFFER();
    }
    else
    {
        for(i=0; i<dim; ++i)
            printf2("- %c: %s;\n", i+'A', list[i]);
        printf2("- %c: Back to previous SubProgram.\n", i+'A');

        PRINTL();

        sel_typ item2;

        do item2 = toupper(getch());
        while(item2 < 'A' || item2 > dim+'A');

        item = item2-'A';
    }

    return item;
}

__MSNATIVE_ inline void __system sigexit(void)
{
    access(sigresult) = true;
    return;
}

__MSNATIVE_ inline bool __system __export catchPause()
{
	#ifdef WINOS
    	signal(SIGINT, (__p_sig_fn_t) sigexit);
    #else
    	signal(SIGINT, (__sighandler_t) sigexit);
    #endif
    if(access(sigresult))
    {
        printf2("\nPress any key to continue Program Action.\n");
        printf2("or %c to stop Execution.\n", EXIT_CHAR);
        access(sigresult) = false;
        if(getch() == EXIT_CHAR) return true;
    }

    return false;
}


__MSNATIVE_ void __system getVarList(void)
{
    void *cookie = NULL;
    EXPRTYPE vval;
    char *vname;

    printf2("\nVariable list items:\n\n");

    exprValList * const tmp = access(exprVars)->var_list;
    
    cookie = exprValListGetNext(tmp, &vname, &vval, NULL, NULL);
    // cookie = exprValListGetNext(vlist, &vname, &vval, NULL, NULL);

    do
    {
        printf2("%s = ", vname);
		printf2(OUTPUT_CONVERSION_FORMAT, vval);
		printf2(";\n");
        cookie = exprValListGetNext(tmp/*suite.exprVars->var_list*/, &vname, &vval, NULL, cookie);
    }
    while(cookie);

    return;
}

__MSNATIVE_ ityp __system __export requires(const char *cmd_string, const char *string, const char *res_string, const unsigned options)
{

    char buf[MAX_BUFSIZ];
    exprObj *e = INIT_OBJLIST;
    int err;
    ityp val;
	struct timeval tvBegin;
    jmp_buf jumper;
    int start, end;


    /* Set error buffer */
    err = setjmp(jumper);
    if(err && e)
        exprFree(e);


    /* Gets an expression */
    if(cmd_string)
        strcpy(buf, cmd_string);
    else
    {
        printf2(string);
        PRINT2N();
        (void) gets(buf);
    }

    access(exitHandle) = INVALID_EXITHANDLE;

    // To check whether User wants to exit from requiringINPUT or Not.
    if((!strcmp(buf, MATRIXBACK_COMMAND)) || !buf[0])
    {

        if(!strcmp(buf, MATRIXBACK_COMMAND))
            access(exitHandle) = EXITHANDLE_BACKCMD;

        if(!buf[0])
            access(exitHandle) = EXITHANDLE_EXIT;

        access(sigresult) = false;
        return NULL_VAL;
    }

    if(buf[strlen(buf)-1] != TERMINATING_CHAR)
        strcat(buf, TERMINATING_STRING);

    PRINTN();

    // Creates expression
    err = exprCreate(&e, access(func_list), access(exprVars)->var_list, access(const_list), NULL, 0);
    if(err != EXPR_ERROR_NOERROR)
        {
        printf2("Expr Creation Error.\n");
        longjmp(jumper, err);
        }

    /* Parse expr */
    err = exprParse(e, buf);
    if(err != EXPR_ERROR_NOERROR)
        {
        exprGetErrorPosition(e, &start, &end);
        printf2("Parse Error (%d,%d).\n", start, end);

        longjmp(jumper, err);
        }

    /* Eval expression */
    const bool assert = (options & PARSER_SHOWDIFFTIME) == PARSER_SHOWDIFFTIME;

    if(assert)
    	gettimeofday(&tvBegin, NULL);

    err = exprEval(e, &val);

    if(err != EXPR_ERROR_NOERROR)
        printf2("Eval Error: %d.\n", err);

    if((options & PARSER_SHOWRESULT) == PARSER_SHOWRESULT)
    {

        PRINTN();
        printf2(res_string);
        printf2(": ");
        printf2(OUTPUT_CONVERSION_FORMAT, val);
        printf2(".\n\n");

    }

    if(assert)
    {
        PRINTL();
        printf2("Average Time: %.*f;\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
        PRINTL();
    }

    if((options & PARSER_SAVERESULT) == PARSER_SAVERESULT)
        *(access(exprVars)->e_ANS) = val;

	if((options & PARSER_SHOWVARLIST) == PARSER_SHOWVARLIST)
    	getVarList();
    	
    exprFree(e);
    return val;
}

__MSNATIVE_ bool __system insertDims(dim_typ *righe, dim_typ *colonne)
{
    const dim_typ old_dims[MAX_DIMENSIONS] =
    {
        (*righe),
        (*colonne)
    };

    printf2("\nEnter ROWS and COLUMNS as expected:\n\[ROWS]\n[COLUMNS].\n");

    PRINTHOWTOBACKMESSAGE();

    ityp tmp, tmp2;
    tmp = tmp2 = 0.00;

    // uint64_t tmp3[MAX_DIMENSIONS];

    while(isNullVal((tmp = requires(NULL, NULL_CHAR, "Inserted ROWS", PARSER_SHOWRESULT))) || tmp != ((*righe) = (dim_typ)tmp) || (*righe) < 1 || (*righe) > MAX_RIGHE ||
        (isNullVal((tmp = requires(NULL, NULL_CHAR, "Inserted COLUMNS", PARSER_SHOWRESULT))) || tmp != ((*colonne) = (dim_typ)tmp) || (*colonne) < 1 || (*colonne) > MAX_COLONNE))
    {
        CLEARBUFFER();
        if(!access(exitHandle)) // if(tmp3[ROWS] == NULL_VAL || tmp3[COLUMNS] == NULL_VAL)
        {
            (*righe) = old_dims[ROWS];
            (*colonne) = old_dims[COLUMNS];
            return false;
        }
        printErr(33, "Invalid [ROWS COLUMNS] format.\nYou have to insert non-negative ROWS and COLUMNS,\n\
and must be respectively less than: %hu and %hu", MAX_RIGHE, MAX_COLONNE);
    }
    return true;
}

__MSNATIVE_ bool __system insertDim(dim_typ *dim, bool mode)
{
    const dim_typ old_dim = (*dim);
    dim_typ max_dim;

    PRINTN();

    if(mode > COLUMNS)
    {
        max_dim = MAX_RIGHE_PER_COLONNE;
        printf2("Enter Quad Matrix DIMENSION.");
    }
    else
    {
        max_dim = mode ? MAX_COLONNE : MAX_RIGHE;
        printf2("Enter Matrix %s.", mode ? "COLUMNS":"ROWS");
    }

    PRINTHOWTOBACKMESSAGE();

    ityp tmp;

    while(isNullVal((tmp = requires(NULL, NULL_CHAR, "Inserted DIMENSION is", PARSER_SHOWRESULT))) || tmp != ((*dim) = (dim_typ)tmp) || (*dim) < 1 || (*dim) > max_dim)
    {
        CLEARBUFFER();
        if(!access(exitHandle))
        {
            (*dim) = old_dim;
            return false;
        }
        printErr(33, "Invalid inserted Value.\nMust be a non-negative integer less than: %hu", max_dim);
    }

    // (*dim) = tmp2;



    return true;
}

__MSSHELL_WRAPPER_ __MSNATIVE_ void showNewtonDifferenceTable(dim_typ n, ityp x[static n], ityp y[MAX_NEWTON_DIFFTABLES_DIM][MAX_NEWTON_DIFFTABLES_DIM], bool mode)
{

    printf2("\n***********%s Difference Table ***********\n", mode ? "FORWARD" : "BACKWARD");

    //display Difference Table

    dim_typ i, j;

    for(i=0;i<n;++i)

    {

        PRINTT();
        printf2(OUTPUT_CONVERSION_FORMAT, x[i]);
        printf2("\t%.*f", DEFAULT_PRECISION, x[i]);

        for(j=0; (mode ? (j<(n-i)):(j<=i)) ;++j)
        {
            PRINTT();
            printf2(OUTPUT_CONVERSION_FORMAT, y[i][j]);
        }

        PRINTN();

    }

    return;
}

__MSNATIVE_ bool __system insertNMMatrix(ityp **matrix, const register dim_typ dim[static MAX_DIMENSIONS])
{
    if(!matrixAlloc(matrix, dim))
        return false;

    volatile char tmp=0;
    dim_typ i, j;
    dim_typ start_col_index;


    // we must seek for backtracking...
    // BY Inserting somewhere a 'return false' statement.
    // (FINAL BACKTRACKING RESPONSE)
    for(i=start_col_index=0; i<dim[ROWS]; ++i)
        for(j=start_col_index; j<dim[COLUMNS]; ++j)
        {
            tmp = insertElement((*matrix), (dim_typ2){i, j}, dim[COLUMNS], false);

            if(tmp == -1)
            {
                if(j > 0)
                    j -= 2;
                else if(i > 0)
                {
                    i -= 2;
                    start_col_index = dim[COLUMNS]-1;
                    break;
                }
                else
                {
                    matrixFree(matrix);
                    return false;
                }
            }

        }

    return true;
}

__MSNATIVE_ volatile char __system insertElement(ityp *restrict matrix, const register dim_typ dim[static MAX_DIMENSIONS], const register dim_typ columns, bool square)
{

    if(square && (dim[ROWS] == MAX_RIGHE_PER_COLONNE || dim[COLUMNS] == MAX_RIGHE_PER_COLONNE))
    {
        printErr(33, "MAX ROWS per COLUMNS Reached");
        return 0;
    }

    if(dim[ROWS] == MAX_RIGHE)
    {
        printErr(33, "MAX ROWS Reached");
        return 0;
    }

    if(dim[COLUMNS] == MAX_COLONNE)
    {
        printErr(33, "MAX COLUMNS Reached");
        return 0;
    }

    // PRINTN();
    printf2("\nEnter [%hu,%hu] Matrix Element.\n", dim[ROWS]+1, dim[COLUMNS]+1);

    sel_typ tmp;
    tmp = 1;
    // tmp = true;

    do
    {
        CLEARBUFFER();

        char str[MIN_STRING];

        sprintf(str, "[%hu,%hu] Matrix Element correctly inserted", dim[ROWS]+1, dim[COLUMNS]+1);


        tmp = !((*(matrix + columns*dim[ROWS] + dim[COLUMNS]) = requires(NULL, NULL_CHAR, str, PARSER_SHOWRESULT)) == NULL_VAL);

        CLEARBUFFER();

        if(access(exitHandle) == EXITHANDLE_BACKCMD)
            return -1;

        if(access(exitHandle) == EXITHANDLE_EXIT && (!dim[ROWS]) && (!dim[COLUMNS]))
            return -2;

    }
    while(tmp != 1 && (!dim[ROWS]) && !dim[COLUMNS]);

    CLEARBUFFER();

    return tmp;
}

__MSNATIVE_ inline volatile bool __system checkBackTracking(volatile char tmp, dim_typ *colonna)
{
    if(tmp < -1)
    {
        if((*colonna) > 0)
            (*colonna) -= 2;
        else
            return false;
    }
    return true;
}

__MSNATIVE_ inline volatile __system sel_typ checkBackTracking2(volatile char tmp, dim_typ *riga, dim_typ *colonna, dim_typ *start_col_index, dim_typ colonne)
{
    if(tmp == -1)
    {
        if((*colonna) > 0)
        {
            (*colonna) -= 2;
            return 0;
        }
        else if((*riga) > 1)
        {

            (*riga) -= 2;
            (*start_col_index) = colonne-1;
            return 1;
        }
        return 2;
    }

    return 3;
}

__MSNATIVE_ bool __system __export enterMatrix(ityp **matrix, dim_typ *righe, dim_typ *colonne, bool square, bool view)
{

    printf2("\n\nEnter%s Matrix", square ? " Quad" : "");
    printf2(" by following procedure you'll see.\n");

    volatile char tmp;
    dim_typ start_col_index;

    tmp = 1;
    start_col_index = 0;

    if(square)
    {
        if(!insertDim(righe, MAX_DIMENSIONS))
        {
            matrixFree(matrix);
            return false;
        }
        (*colonne) = (*righe);
    }

    else if(!insertDims(righe, colonne)) return false;

    CLEARBUFFER();

    if(!matrixAlloc(matrix, (dim_typ2){*righe, *colonne}))
    {
        matrixFree(matrix);
        return false;
    }

    dim_typ i;

    for(i=0; i<(*colonne); ++i)
    {
        tmp = insertElement((*matrix), (dim_typ2){0, i}, *colonne, square);
        
        if(!checkBackTracking(tmp, &i))
        {
            matrixFree(matrix);
            return false;
        }
    }

    dim_typ j;
    volatile sel_typ back_tracking;

    for(i = 1; i<(*righe) && __pmode__ != ALGOPS_DOTPRODUCT; ++i)
        for(j=start_col_index; j<(*colonne); ++j)
        {
            tmp = insertElement((*matrix), (dim_typ2){i, j}, *colonne, square);

            if((back_tracking = checkBackTracking2(tmp, &i, &j, &start_col_index, (*colonne))) == 1)
                break;

            else if(back_tracking == 2)
            {
                matrixFree(matrix);
                return false;
            }



            if(tmp == -1)
            {
                if(j > 0)
                    j -= 2;
                else if(i > 1)
                {

                    i -= 2;
                    start_col_index = (*colonne)-1;
                    break;
                }
                else
                {
                    matrixFree(matrix);
                    return false;
                }
            }

        }


    if(view)
    {
        printf2("\nInserted [%hu X %hu] Matrix is:\n\n", *righe, *colonne);
        printMatrix(stdout, (*matrix), (dim_typ2){*righe, *colonne});
    }

    return true;
}

