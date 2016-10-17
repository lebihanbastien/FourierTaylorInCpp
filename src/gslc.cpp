#include "gslc.h"

//----------------------------------------------------------------------------------------
//
// Vector <--> matrix
//
//----------------------------------------------------------------------------------------
/**
 * \brief Transform a vector into a matrix with a given shift in the initial vector
 **/
void gslc_vectorToMatrix(gsl_matrix *m, const double y[], int rows, int columns, int shift)
{
    int i,j,k;
    for(i=1; i<=rows ; i++)
        for(j=1; j<=columns; j++)
        {
            k = rows*(i-1)+j;
            gsl_matrix_set(m, i-1, j-1, y[shift+k-1]);
        }
}

/**
 * \brief Transform a matrix into a vector with a given shift in the final vector
 **/
void gslc_matrixToVector(double y[], const gsl_matrix *m, int rows, int columns, int shift)
{
    int i,j,k;
    for(i=1; i<=rows ; i++)
        for(j=1; j<=columns; j++)
        {
            k = rows*(i-1)+j;
            y[shift+k-1] = gsl_matrix_get(m, i-1, j-1);
        }
}


//----------------------------------------------------------------------------------------
// Printing a real matrix
//----------------------------------------------------------------------------------------
/**
 * \brief Print a real matrix
 **/
void gslc_matrix_printf(const gsl_matrix *M)
{
    int im = M->size1;
    int jm = M->size2;

    for(int i = 0; i < im ; i++)
    {
        for(int j=0; j <jm; j++) printf("%+5.15e ", gsl_matrix_get(M, i, j));
        printf("\n");
    }
}

/**
 * \brief Print an approximation of a real matrix: the coefficients under 1e-5 are displayed as null.
 **/
void gslc_matrix_printf_approx(const gsl_matrix *M)
{
    int im = M->size1;
    int jm = M->size2;
    double c;
    for(int i = 0; i < im ; i++)
    {
        for(int j=0; j <jm; j++)
        {
            c = gsl_matrix_get(M, i, j);
            if(fabs(c) > 1e-5) printf("%+1.3e ", c);
            else printf("%+1.3e ", 0.0);
        }
        printf("\n");
    }
}

/**
 * \brief Print a real matrix into a txt file.
 **/
void gslc_matrix_fprintf(const gsl_matrix *M, char* fileName)
{

    FILE *f;
    f = fopen(fileName, "w");
    int im = M->size1;
    int jm = M->size2;

    for(int i = 0; i < im ; i++)
    {
        for(int j=0; j <jm; j++) fprintf(f, "%+5.16e ", gsl_matrix_get(M, i, j));
        fprintf(f, "\n");
    }

    fclose(f);
}


//----------------------------------------------------------------------------------------
// Complex numbers
//----------------------------------------------------------------------------------------
/**
 * \brief Create a gsl_complex from 2 doubles
 **/
gsl_complex gslc_complex(double real, double imag)
{
    gsl_complex one_c;
    GSL_SET_COMPLEX(&one_c, real, imag);
    return one_c;
}

/**
 * \brief Create a gsl_complex from a complex number
 **/
gsl_complex gslc_complex(cdouble x)
{
    gsl_complex one_c;
    GSL_SET_COMPLEX(&one_c, creal(x), cimag(x));
    return one_c;
}

/**
 * \brief Create a complex number from a gsl_complex
 **/
cdouble gslc_complex(gsl_complex c)
{
    return GSL_REAL(c) + I*GSL_IMAG(c);
}

//----------------------------------------------------------------------------------------
// Printing a complex matrix
//----------------------------------------------------------------------------------------
/**
 * \brief Print a complex matrix.
 **/
void gslc_matrix_complex_printf(const gsl_matrix_complex *M)
{

    int im = M->size1;
    int jm = M->size2;

    for(int i = 0; i < im ; i++)
    {
        for(int j=0; j <jm; j++)
        {
            /*if( GSL_IMAG(gsl_matrix_complex_get(M, i, j)) != 0.0)*/ printf("%+1.10e%+1.10ei ", GSL_REAL(gsl_matrix_complex_get(M, i, j)), GSL_IMAG(gsl_matrix_complex_get(M, i, j)));
            //else printf("%+1.0e ", GSL_REAL(gsl_matrix_complex_get(M, i, j)));

        }
        printf("\n");
    }
}

/**
 * \brief Print the real part of a complex matrix into a txt file.
 **/
void gslc_matrix_complex_printf_real(const gsl_matrix_complex *M)
{

    int im = M->size1;
    int jm = M->size2;

    for(int i = 0; i < im ; i++)
    {
        for(int j=0; j <jm; j++) printf("%+5.5e ", GSL_REAL(gsl_matrix_complex_get(M, i, j)));
        printf("\n");
    }
}

/**
 * \brief Print the imag part of a complex matrix into a txt file.
 **/
void gslc_matrix_complex_printf_imag(const gsl_matrix_complex *M)
{

    int im = M->size1;
    int jm = M->size2;

    for(int i = 0; i < im ; i++)
    {
        for(int j=0; j <jm; j++) printf("%+5.5e ", GSL_IMAG(gsl_matrix_complex_get(M, i, j)));
        printf("\n");
    }
}

/**
 * \brief Print an approximation of a complex matrix: the coefficients under 1e-5 are displayed as null.
 **/
void gslc_matrix_complex_printf_approx(const gsl_matrix_complex *M)
{
    int im = M->size1;
    int jm = M->size2;
    double c;

    printf("real(M) = \n");
    for(int i = 0; i < im ; i++)
    {
        for(int j=0; j <jm; j++)
        {
            c  = GSL_REAL(gsl_matrix_complex_get(M, i, j));
            if(fabs(c) > 1e-5) printf("%+1.3e ", c);
            else printf("%+1.3e ", 0.0);
        }
        printf("\n");
    }

    printf("imag(M) = \n");
    for(int i = 0; i < im ; i++)
    {
        for(int j=0; j <jm; j++)
        {
            c = GSL_IMAG(gsl_matrix_complex_get(M, i, j));
            if(fabs(c) > 1e-5) printf("%+1.3e ", c);
            else printf("%+1.3e ", 0.0);
        }
        printf("\n");
    }
}

/**
 * \brief Print a complex matrix into a txt file.
 **/
void gslc_matrix_complex_fprintf(const gsl_matrix_complex *M, char* fileName)
{
    FILE *f;

    f = fopen(fileName, "w");
    int im = M->size1;
    int jm = M->size2;

    for(int i = 0; i < im ; i++)
    {
        for(int j=0; j <jm; j++) fprintf(f, "%+5.16e %+5.16e ", GSL_REAL(gsl_matrix_complex_get(M, i, j)), GSL_IMAG(gsl_matrix_complex_get(M, i, j)));
        fprintf(f, "\n");
    }

    fclose(f);
}

/**
 * \brief Print a complex matrix into a txt file, separating the real and imaginary part, in order to be easier to read.
 **/
void gslc_matrix_complex_fprintf_pretty(const gsl_matrix_complex *M, char* fileName)
{
    FILE *f;
    f = fopen(fileName, "w");
    int im = M->size1;
    int jm = M->size2;

    fprintf(f,"real part:\n");
    for(int i = 0; i < im ; i++)
    {
        for(int j=0; j <jm; j++) fprintf(f, "%+5.16e ", GSL_REAL(gsl_matrix_complex_get(M, i, j)));
        fprintf(f, "\n");
    }

    fprintf(f,"imag part:\n");
    for(int i = 0; i < im ; i++)
    {
        for(int j=0; j <jm; j++) fprintf(f, "%+5.16e ", GSL_IMAG(gsl_matrix_complex_get(M, i, j)));
        fprintf(f, "\n");
    }
    fclose(f);
}

/**
 * \brief Print the real part of a complex matrix into a txt file.
 **/
void gslc_matrix_complex_fprintf_real(const gsl_matrix_complex *M, char* fileName)
{
    FILE *f;

    f = fopen(fileName, "w");
    int im = M->size1;
    int jm = M->size2;

    for(int i = 0; i < im ; i++)
    {
        for(int j=0; j <jm; j++) fprintf(f, "%+2.0e ", GSL_REAL(gsl_matrix_complex_get(M, i, j)));
        fprintf(f, "\n");
    }

    fclose(f);
}


//----------------------------------------------------------------------------------------
// Complex to real matrices
//----------------------------------------------------------------------------------------
void gslc_matrix_complex_to_matrix(const gsl_matrix_complex *Mc, gsl_matrix *M)
{
    int im = M->size1;
    int jm = M->size2;

    for(int i = 0; i < im ; i++)
    {
        for(int j=0; j <jm; j++) gsl_matrix_set(M, i, j, GSL_REAL(gsl_matrix_complex_get(Mc, i, j)));
    }
}


//----------------------------------------------------------------------------------------
// Allocation
//----------------------------------------------------------------------------------------
/**
 * \brief Initialize and return an array of GSL matrices
 **/
gsl_matrix** gslc_matrix_array_calloc(int size1, int size2, int M)
{
    gsl_matrix **DAT = (gsl_matrix**) malloc((M)*sizeof(gsl_matrix*));
    for(int i = 0; i< M; i++) DAT[i] = gsl_matrix_calloc(size1,size2);
    return DAT;
}

/**
 * \brief Free an array of GSL matrices
 **/
void gslc_matrix_array_free(gsl_matrix ** P, int M)
{
    for(int i = 0; i< M; i++) gsl_matrix_free(P[i]);
}

/**
 * \brief Initialize and return an array of GSL vectors
 **/
gsl_vector** gslc_vector_array_calloc(int size1, int M)
{
    gsl_vector **DAT = (gsl_vector**) malloc((M)*sizeof(gsl_vector*));
    for(int i = 0; i< M; i++)  DAT[i] = gsl_vector_calloc(size1);
    return DAT;
}

/**
 * \brief Free an array of GSL vectors
 **/
void gslc_vector_array_free(gsl_vector ** P, int M)
{
    for(int i = 0; i< M; i++) gsl_vector_free(P[i]);
}


//----------------------------------------------------------------------------------------
// Norm
//----------------------------------------------------------------------------------------
/**
 *  \brief L2 norm of a complex matrix: res = ||Mc||_2
 **/
double gslc_matrix_complex_L2(const gsl_matrix_complex *Mc)
{
    //Init
    double res = 0.0;

    //Loop
    for(unsigned int i = 0; i < Mc->size1; i++)
    {
        for(unsigned int j=0; j < Mc->size2; j++) res += gsl_complex_abs2((gsl_matrix_complex_get(Mc, i, j)));
    }

    return sqrt(res);
}

/**
 *  \brief L2 norm of a difference of complex matrices: res = ||M1c - M2c ||_2
 **/
double gslc_matrix_complex_diff_L2(const gsl_matrix_complex *M1c, const gsl_matrix_complex *M2c)
{
    //------------------------------------------------------------------------------------
    //Check
    //------------------------------------------------------------------------------------
    if(M1c->size1 != M2c->size1 ||  M1c->size2 != M2c->size2)
    {
        printf("gslc_matrix_complex_diff_L2. Dimension mismatch. 0 is returned\n");
        return 0;
    }

    //------------------------------------------------------------------------------------
    //Compute the difference
    //------------------------------------------------------------------------------------
    //Init Mt
    gsl_matrix_complex *Mt = gsl_matrix_complex_alloc(M1c->size1, M1c->size2);
    //Mt = M1c
    gsl_matrix_complex_memcpy(Mt, M1c);
    //Mt = M1c - M2c
    gsl_matrix_complex_sub(Mt, M2c);
    //res = L2(Mt)
    double res = gslc_matrix_complex_L2(Mt);

    //------------------------------------------------------------------------------------
    //Free and return
    //------------------------------------------------------------------------------------
    gsl_matrix_complex_free(Mt);
    return res;
}
