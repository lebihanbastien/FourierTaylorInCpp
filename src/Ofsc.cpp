#include "Ofsc.h"

int OFS_ORDER;

/**
 * \file ofs.cpp
 * \brief Fourier series template class (src)
 * \author BLB
 * \date May 2015
 * \version 1.0
 */


//----------------------------------------------------------------------------------------
//Create
//----------------------------------------------------------------------------------------
/**
 *  \brief Default constructor of the class Ofsc.
 */
Ofsc::Ofsc()
{
    order = OFS_ORDER;
    coef = new cdouble[2*order+1];
    this->zero(); //all coefficients to zero
}

/**
 *  \brief Constructor with a given order.
 */
Ofsc::Ofsc(const int newOrder)
{
    order = newOrder;
    coef = new cdouble[2*order+1];
    this->zero(); //all coefficients to zero
}

/**
 *  \brief Constructor from a given Ofsc object.
 */
Ofsc::Ofsc(Ofsc const& b)
{
    order = b.order;
    coef = new cdouble[2*order+1];
    for(int i = -order ; i<= order; i++) coef[i+order] = b.coef[i+order];
}

//----------------------------------------------------------------------------------------
//Copy
//----------------------------------------------------------------------------------------
/**
 *  \brief  Copy from a given Ofs object (only the coefficients).
 */
Ofsc& Ofsc::ccopy(Ofsc const& b)
{
    if(order != b.order)
    {
        cout << "Erreur in Ofsc::ccopy: orders do not match. Nothing is done." << endl;
        return *this;
    }
    else
    {
        for(int i = -order ; i<= order; i++) coef[i+order] = b.coef[i+order];
        return *this;
    }
}


//----------------------------------------------------------------------------------------
//Delete
//----------------------------------------------------------------------------------------
/**
 *  \brief Default destructor of the class Ofsc.
 */
Ofsc::~Ofsc()
{
    if(coef != NULL) delete[] coef;
    coef = 0;
}

//----------------------------------------------------------------------------------------
//Setters
//----------------------------------------------------------------------------------------
/**
 *  \brief Sets a coefficient at a given position in the serie.
 */
void Ofsc::setCoef(cdouble const& value, int const& pos)
{
    if(fabs(pos) <= order) coef[pos+order] = value;
    else cout << "Error in Ofsc::setCoef: position is out of scope. No coefficient is set." << endl;
}


/**
 *  \brief Adds a coefficient at a given position in the serie.
 */
void Ofsc::addCoef(cdouble const& value, int const& pos)
{
    if(fabs(pos) <= order) coef[pos+order] += value;
    else cout << "Error in Ofsc::addCoef: position is out of scope\n" << endl;
}


//----------------------------------------------------------------------------------------
//Operators (+=, -=, ...)
//----------------------------------------------------------------------------------------
/**
 *  \brief  An operator. Constructor from a given Ofsc object (only the coefficients).
 */
Ofsc& Ofsc::operator = (Ofsc const& b)
{
    if(this != &b)
    {
        if(order != b.order)
        {
            order = b.order;
            if(coef != NULL) delete coef;
            coef = new cdouble[2*order+1];
            for(int i = -order ; i<= order; i++) coef[i+order] = b.coef[i+order]; //this->setCoef(b.getCoef(i), i);
        }
        else
        {
            for(int i = -order ; i<= order; i++) coef[i+order] = b.coef[i+order];
        }
    }
    return *this; //same object if returned
}

/**
 *  \brief  An operator. Adds all coefficients term by term  from a given Ofs object.
 *
 *   Allows b.order != order.
 */
Ofsc& Ofsc::operator += (Ofsc const& b)
{
    if(b.order > order) //if b.order > order, a new array of coefficients must be set
    {
        //Copy coef into temporary array
        cdouble temp[2*order+1];
        for(int i = 0 ; i< 2*order + 1; i++) temp[i] = coef[i];
        //Recreate a good array
        delete coef;
        coef = new cdouble[2*b.order+1];
        //Store the coefficients again
        for(int i = -order ; i<= order; i++) coef[i+b.order] = temp[i+order];
        order = b.order;
    }

    //Adding the coefficients
    for(int i = -b.order ; i <= b.order ; i++) coef[i+order] += b.coef[i+b.order];
    return *this;
}

/**
 *  \brief  An operator. Subtracts all coefficients term by term  from a given Ofs object.
 *
 *   Allows b.order != order.
 */
Ofsc& Ofsc::operator -= (Ofsc const& b)
{
    if(b.order > order) //if b.order > order, a new array of coefficients must be set
    {
        //Copy coef into temporary array
        cdouble temp[2*order+1];
        for(int i = 0 ; i< 2*order + 1; i++) temp[i] = coef[i];
        //Recreate a good array
        delete coef;
        coef = new cdouble[2*b.order+1];
        //Store the coefficients again
        for(int i = -order ; i<= order; i++) coef[i+b.order] = temp[i+order];
        order = b.order;
    }

    //Adding the coefficients
    for(int i = -b.order ; i <= b.order ; i++) coef[i+order] -= b.coef[i+b.order];
    return *this;
}

/**
 *  \brief  An operator. Multiplies all coefficients by a given \c cdouble coefficient.
 */
Ofsc& Ofsc::operator *= (cdouble const& c)
{
    for(int i=0; i<2*order+1; i++) coef[i] *= c;
    return *this;
}

/**
 *  \brief  An operator. Divides all coefficients by a given \c cdouble coefficient.
 */
Ofsc& Ofsc::operator /= (cdouble const& c)
{
    for(int i=0; i<2*order+1; i++) coef[i] /= c;
    return *this;
}

//----------------------------------------------------------------------------------------
//Getters
//----------------------------------------------------------------------------------------
/**
 *  \brief  Gets the order of the serie.
 */
int Ofsc::getOrder() const
{
    return order;
}

/**
 *  \brief  Gets the pointer address of the Ofsc object
 */
Ofsc* Ofsc::getAddress() const
{
    return (Ofsc*) this;
}

/**
 *  \brief  Gets the coefficient at a given position. cdouble case
 */
cdouble Ofsc::getCoef(int pos) const
{
    if(fabs(pos) <= order)  return coef[pos+order];
    else
    {
        cout << "Warning in Ofsc::getCoef: position " << pos << " is out of scope. 0.0 is returned by default." << endl;
        return  0.0+0.0*I;
    }
}


//----------------------------------------------------------------------------------------
//Zeroing
//----------------------------------------------------------------------------------------

/**
 *  \brief  Sets all coefficients to zero. cdouble case
 */
void Ofsc::zero()
{
    for(int i = -order ; i<= order; i++) coef[i+order] = 0.0+0.0*I;
}

/**
 *  \brief  Is the Ofsc object equal to zero, at order ofs_order?
 */
bool Ofsc::isnull(const int ofs_order) const
{
    for(int i = -min(ofs_order, order); i <= min(ofs_order, order); i++)
    {
        if(cabs(getCoef(i)) != 0.0) return false;
    }
    return true;
}


//----------------------------------------------------------------------------------------
// Functions (evaluate)
//----------------------------------------------------------------------------------------
/**
 *  \brief  Evaluates the Ofsc object at angle theta (theta = nt) and at a certain order eff_order.
 *
 *  The sum is given in the form: \f$ c_0 + \sum \limits_{k = 1}^{J} \cos(k\theta) (c_k + c_{-k}) + i\sin(k\theta) (c_k - c_{-k}) \f$
 *  in order to avoid any additionnal roundoff errors by using cos and sin functions as little as possible.
 */
cdouble Ofsc::fevaluate(double const cR[], double const sR[], int eff_order) const
{
    cdouble result = 0+0.0*I;
    //Order -k and +k are calculated together to avoid any additionnal roundoff errors
    for(int k=eff_order; k>=1; k--) result+= cR[k-1]*(coef[k+order] + coef[-k+order]) + I*sR[k-1]*(coef[k+order] - coef[-k+order]);
    //Order 0
    result+= coef[0+order];
    //Return result
    return result;
}

/**
 *  \brief  Evaluates the Ofsc object at angle theta (theta = nt) and at a certain order eff_order.
 *
 *  The sum is given in the form: \f$ c_0 + \sum \limits_{k = 1}^{J} \cos(k\theta) (c_k + c_{-k}) + i\sin(k\theta) (c_k - c_{-k}) \f$
 *  in order to avoid any additionnal roundoff errors by using cos and sin functions as little as possible.
 */
cdouble Ofsc::evaluate(double const& theta, int eff_order) const
{
    cdouble result = 0+0.0*I;
    double cR[eff_order];
    double sR[eff_order];

    cR[0] = cos(theta);
    sR[0] = sin(theta);
    for(int i = 1; i< eff_order; i++)
    {
        cR[i] =  cR[i-1]*cR[0] - sR[i-1]*sR[0];
        sR[i] =  sR[i-1]*cR[0] + cR[i-1]*sR[0];
    }

    //Order -k and +k are calculated together to avoid any additionnal roundoff errors
    for(int k=eff_order; k>=1; k--) result+= cR[k-1]*(coef[k+order] + coef[-k+order]) + I*sR[k-1]*(coef[k+order] - coef[-k+order]);
    //Order 0
    result+= coef[0+order];//getCoef(0);


    return result;
}

/**
 *  \brief  Evaluates the Ofsc object at angle theta (theta = nt).
 *
 *  The sum is given in the form: \f$ c_0 + \sum \limits_{k = 1}^{J} \cos(k\theta) (c_k + c_{-k}) + i\sin(k\theta) (c_k - c_{-k}) \f$
 *  in order to avoid any additionnal roundoff errors by using cos and sin functions as little as possible.
 */
cdouble Ofsc::evaluate(double const& theta) const
{
    cdouble result = 0+0.0*I;
    double cR[order];
    double sR[order];

    cR[0] = cos(theta);
    sR[0] = sin(theta);
    for(int i = 1; i< order; i++)
    {
        cR[i] =  cR[i-1]*cR[0] - sR[i-1]*sR[0];
        sR[i] =  sR[i-1]*cR[0] + cR[i-1]*sR[0];
    }

    //Order -k and +k are calculated together to avoid any additionnal roundoff errors
    for(int k=order; k>=1; k--) result+= cR[k-1]*(coef[k+order] + coef[-k+order]) + I*sR[k-1]*(coef[k+order] - coef[-k+order]);
    //Order 0
    result+= coef[0+order];//getCoef(0);

    return result;
}

//----------------------------------------------------------------------------------------
// Functions (operations)
//----------------------------------------------------------------------------------------
/**
 *  \brief  An operation. Adds the product: \c this \f$  += m a \f$ at a certain order eff_order.
 *
 *  Note 1: can be used in place.
 *  Note 2: For speed, the use of addCoef and getCoef has been discarded,
 *          which means no test at all in this routine!
 */
void Ofsc::ofs_smult(Ofsc const& a, cdouble c, int eff_order)
{
    //Sum
    for(int i = -eff_order; i <= eff_order; i++)
    {
        //addCoef(c*a.getCoef(i), i);
        coef[i+order] += c*a.coef[i+order];
    }
}

/**
 *  \brief  An operation. Computes the product: \c this \f$  = this a \f$ at a certain order eff_order.
 */
void Ofsc::ofs_mult_inline(cdouble c, int eff_order)
{
    //Sum
    for(int i = -eff_order; i <= eff_order; i++)
    {
        coef[i+order] = c*coef[i+order];
    }
}

/**
 *  \brief  An operation. Adds the product: \c this \f$ += a \times b \f$.
 *           Works fine when a.order = b.order which is the default case.
 */
void Ofsc::ofs_sprod(Ofsc const& a, Ofsc const& b)
{
    int psup, pinf;
    //Product
    for(int n=-order ; n<= order; n++)
    {
        psup = min(n+order,  order);
        pinf = max(n-order, -order);
        for(int p = pinf; p <= psup; p++) coef[n+order] += a.coef[p+order]*b.coef[n-p+order];
    }
}

/**
 *  \brief  An operation. Sets the sum-product: \c this \f$  = m_a a + m_b b \f$.
 *
 *  Note: can be used in place.
 */
void Ofsc::ofs_fsum(Ofsc const& a, cdouble const& ma, Ofsc const& b, cdouble const& mb)
{
    if(order != a.order ||  order != b.order)
    {
        cout << "Error using fsum: the order does not match. Initial Ofsc is returned" << endl;
    }
    else
    {
        for(int i=-order; i<=order; i++)
        {
            coef[i+order] = ma*a.coef[i+order] + mb*b.coef[i+order];//setCoef(ma*a.getCoef(i)+mb*b.getCoef(i), i);
        }
    }
}

//----------------------------------------------------------------------------------------
// Functions (operations) involving the time domain
//----------------------------------------------------------------------------------------
/**
 *  \brief An operation. Performs the expansion this \f$ = a^{\alpha} \f$ using inverse FFT, the power function in time domain, and finally direct FFT.
 **/
void Ofsc::ofs_pows(Ofsc const& a, cdouble const& alpha)
{
    this->tfs_from_ofs(a);
    this->tfs_pows(alpha);
    this->tfs_to_ofs_inline();
}

/**
 *  \brief An operation. Performs the division this \f$ = a / b \f$ using inverse FFT, the division in time domain, and finally direct FFT.
 **/
void Ofsc::ofs_div(Ofsc const& a, Ofsc const& b, Ofsc& temp)
{
        //this = a in tfs
        this->tfs_from_ofs(a);

        //temp = b in tfs
        temp.tfs_from_ofs(b);

        //this = a/b in tfs
        this->tfs_div_inline(temp);

        //back in frequency domain (ofs)
        this->tfs_to_ofs_inline();
}

//----------------------------------------------------------------------------------------
// Frequency domain <--> Time domain
//----------------------------------------------------------------------------------------
/**
 *  \brief  From Frequency domain to time domain.
 */
void Ofsc::tfs_from_ofs(Ofsc const& a)
{
    //---------------------
    //if wrong orders, stop
    //---------------------
    if(order < a.getOrder())
    {
        cout << "tfs_from_ofs. Wrong orders. Stop." << endl;
        return;
    }

    int N = 2*order+1;
    //---------------------
    //Copy the coefficients in coef
    //---------------------
    for(int i = 0; i < N; i++) coef[i] = ((Ofsc)a).evaluate(i*2*M_PI/((double)N));
}

/**
 *  \brief  Inline from Frequency domain to time domain.
 */
 void Ofsc::tfs_from_ofs_inline(Ofsc& temp)
{
    //---------------------
    //if wrong orders, stop
    //---------------------
    if(order != temp.getOrder())
    {
        cout << "tfs_from_ofs. Wrong orders. Stop." << endl;
        return;
    }

    //---------------------
    //Copy in temp
    //---------------------
    temp.ccopy(*this);

    int N = 2*order+1;
    //---------------------
    //Copy the coefficients in coef
    //---------------------
    for(int i = 0; i < N; i++) coef[i] = temp.evaluate(i*2*M_PI/((double)N));
}

/**
 *  \brief  From Time domain to Frequency domain.
 */
 void Ofsc::tfs_to_ofs(Ofsc const& a)
{
    //---------------------
    // If wrong orders, stop
    //---------------------
    if(order > a.getOrder())
    {
        cout << "tfs_to_ofs. Wrong orders. Stop." << endl;
        return;
    }

    //---------------------
    // FFT structures
    //---------------------
    int N = 2*a.getOrder()+1;
    gsl_vector_complex *data_fft          = gsl_vector_complex_calloc(N);
    gsl_fft_complex_wavetable *wavetable  = gsl_fft_complex_wavetable_alloc (N);
    gsl_fft_complex_workspace *workspace  = gsl_fft_complex_workspace_alloc (N);

    //---------------------
    //Copy a in data_fft
    //---------------------
    for(int i = 0; i < N; i++) gsl_vector_complex_set(data_fft, i, gslc_complex(creal(a.coef[i]), cimag(a.coef[i])));

    //---------------------
    //FFT
    //---------------------
    gsl_fft_complex_forward (data_fft->data, 1, N, wavetable, workspace);

    //---------------------
    //Order 0
    //---------------------
    this->setCoef(+GSL_REAL(gsl_vector_complex_get(data_fft, 0))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, 0))/(double)N,  0);
    //Version without setCoef
    //this->coef[order] = +GSL_REAL(gsl_vector_complex_get(data_fft, 0))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, 0))/(double)N;

    //---------------------
    //Order n
    //---------------------
    for(int i = 1; i<= order; i++)
    {
        //Negative frequecies
        this->setCoef(+GSL_REAL(gsl_vector_complex_get(data_fft, N-i))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, N-i))/(double)N, -i);
        //Positive frequencies
        this->setCoef(+GSL_REAL(gsl_vector_complex_get(data_fft, i))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, i))/(double)N,  i);
        //Version without setCoef
        //this->coef[order-i] = +GSL_REAL(gsl_vector_complex_get(data_fft, N-i))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, N-i))/(double)N;
        //this->coef[order+i] = +GSL_REAL(gsl_vector_complex_get(data_fft, i))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, i))/(double)N;
    }

    //------------------------
    // Free
    //------------------------
    gsl_vector_complex_free(data_fft);
    gsl_fft_complex_wavetable_free (wavetable);
    gsl_fft_complex_workspace_free (workspace);
}

/**
 *  \brief  Inline from Time domain to Frequency domain.
 */
 void Ofsc::tfs_to_ofs_inline()
{
    //---------------------
    // FFT structures
    //---------------------
    int N = 2*order+1;
    gsl_vector_complex *data_fft          = gsl_vector_complex_calloc(N);
    gsl_fft_complex_wavetable *wavetable  = gsl_fft_complex_wavetable_alloc (N);
    gsl_fft_complex_workspace *workspace  = gsl_fft_complex_workspace_alloc (N);

    //---------------------
    //Copy a in data_fft
    //---------------------
    for(int i = 0; i < N; i++) gsl_vector_complex_set(data_fft, i, gslc_complex(creal(this->coef[i]), cimag(this->coef[i])));

    //---------------------
    //FFT
    //---------------------
    gsl_fft_complex_forward (data_fft->data, 1, N, wavetable, workspace);

    //---------------------
    //Order 0
    //---------------------
    this->setCoef(+GSL_REAL(gsl_vector_complex_get(data_fft, 0))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, 0))/(double)N,  0);
    //Version without setCoef
    //this->coef[order] = +GSL_REAL(gsl_vector_complex_get(data_fft, 0))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, 0))/(double)N;

    //---------------------
    //Order n
    //---------------------
    for(int i = 1; i<= order; i++)
    {
        //Negative frequecies
        this->setCoef(+GSL_REAL(gsl_vector_complex_get(data_fft, N-i))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, N-i))/(double)N, -i);
        //Positive frequencies
        this->setCoef(+GSL_REAL(gsl_vector_complex_get(data_fft, i))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, i))/(double)N,  i);
        //Version without setCoef
        //this->coef[order-i] = +GSL_REAL(gsl_vector_complex_get(data_fft, N-i))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, N-i))/(double)N;
        //this->coef[order+i] = +GSL_REAL(gsl_vector_complex_get(data_fft, i))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, i))/(double)N;
    }

    //------------------------
    // Free
    //------------------------
    gsl_vector_complex_free(data_fft);
    gsl_fft_complex_wavetable_free (wavetable);
    gsl_fft_complex_workspace_free (workspace);
}


//----------------------------------------------------------------------------------------
// TFS operations
//----------------------------------------------------------------------------------------

// pows
//----------------------------------------------------------------------------------------
/**
 *  \brief  An operation. Performs the power this = a^alpha in time domain.
 */
void Ofsc::tfs_pows(Ofsc const& a, cdouble const& alpha)
{
    for(int i = 0; i < 2*order+1; i++) coef[i] = cpow(a.coef[i], alpha);
}

/**
 *  \brief  An operation. Performs the power this = this^alpha in time domain.
 */
void Ofsc::tfs_pows(cdouble const& alpha)
{
    for(int i = 0; i < 2*order+1; i++) coef[i] = cpow(coef[i], alpha);
}

// sprod
//----------------------------------------------------------------------------------------

/**
 *  \brief  An operation. Computes the product: \c this \f$ = this \times b \f$ in time domain.
 */
void Ofsc::tfs_prod_inline(Ofsc const& b)
{
    for(int k=0 ; k< 2*order+1; k++) coef[k] = coef[k]*b.coef[k];
}


/**
 *  \brief  An operation. Adds the product: \c this \f$ += a \times b \f$ in time domain.
 */
void Ofsc::tfs_sprod(Ofsc const& a, Ofsc const& b)
{
    for(int k=0 ; k< 2*order+1; k++) coef[k] += a.coef[k]*b.coef[k];
}


// sdiv
//----------------------------------------------------------------------------------------
/**
 *  \brief  An operation. Computes the division: \c this \f$ = this / b \f$ in time domain.
 */
void Ofsc::tfs_div_inline(Ofsc const& b)
{
    for(int k=0 ; k< 2*order+1; k++) coef[k] = coef[k]/b.coef[k];
}

/**
 *  \brief  An operation. Adds the division: \c this \f$ += a / b \f$ in time domain.
 */
void Ofsc::tfs_sdiv(Ofsc const& a, Ofsc const& b)
{
    for(int k=0 ; k< 2*order+1; k++) coef[k] += a.coef[k]/b.coef[k];
}



//----------------------------------------------------------------------------------------
//Print
//----------------------------------------------------------------------------------------
/**
 *  \brief  A stream operator
 */
std::ostream& operator << (std::ostream& stream, Ofsc const& ofs)
{
    //Coefficients
    for(int i = -min(ofs.order,ofs.order) ; i<= min(ofs.order,ofs.order) ; i++)
    {
        stream << setw(3) << setiosflags(ios::right) << std::showpos << i << "   ";
        stream <<  setiosflags(ios::scientific) << setprecision(15);
        stream << creal(ofs.coef[i+ofs.order]) << "  " << cimag(ofs.coef[i+ofs.order]) << endl;

    }
    return stream;
}


//----------------------------------------------------------------------------------------
// Reading an OFS from a text file
//----------------------------------------------------------------------------------------
/**
 * \fn void inline readOFS_txt(Ofsc& xFFT, string filename, int fftN)
 * \brief Reading an Ofsc object from a text file.
 */
void readOFS_txt(Ofsc& xFFT, string filename)
{
    //Init
    ifstream readStream;
    double ct, cr, ci;
    int fftN = xFFT.getOrder();

    //Reading
    readStream.open((filename+".txt").c_str());
    for(int i = -fftN; i<=fftN; i++)
    {
        readStream >> ct;  //current order
        readStream >> cr;  //real part
        readStream >> ci;  //imag part
        xFFT.setCoef(cr+I*ci, i);
    }
    readStream.close();
}

//----------------------------------------------------------------------------------------
//Derivation
//----------------------------------------------------------------------------------------
/**
 *  \brief  An operation. Set the time derivative of object \c a with pulsation \f$ \omega = n \f$, so that \c this \f$ = \dot{a} \f$.
 */
void Ofsc::dot(Ofsc const& a, double const& n)
{
    //d(a_k)/dt = k*n*I*a_k
    for(int k=-order; k<= order; k++) coef[k+order] = k*n*I*a.coef[k+order]; //this->setCoef(k*n*I*a.getCoef(k), k);
}
