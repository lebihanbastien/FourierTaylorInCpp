#ifndef OFTS_H
#define OFTS_H

#include <fstream>
#include <vector>
#include <sstream>

#include "Ofsc.h"
#include "FTA.h"


class Oftsc
{
private:

    int nv;            //number of variables
    int order;         //order
    int cnv;           //number of variables of the coefficients
    int corder;        //order of the coefficients
    Ofsc **coefs;      //matrix of coefficients

public:

    //------------------
    //Create
    //------------------
    /**
     *  \brief Default constructor of the class Oftsc.
     */
    //Oftsc();
    /**
     *  \brief Constructor with given order/number of variables both for the Fourier-Taylor series and the coefficients.
     *  \param newNv: number of variables of the serie
     *  \param newOrder: order of the serie
     *  \param newCnv: number of variables of the coefficients
     *  \param newOrder: order of the coefficients
     */
    Oftsc(int newNv, int newOrder, int newCnv, int newCorder);
    /**
     *  \brief Constructor from a given Oftsc object (without any link).
     *  \param b:  a reference to the Oftsc object to copy in the new object
     */
    Oftsc(Oftsc const& b);


    //------------------
    //Getters
    //------------------
    /**
     *  \brief  Gets the order of the serie.
     *  \return the order of the serie as an \c int
     */
    int getOrder() const;

    /**
     *  \brief  Gets the number of variables of serie.
     *  \return the number of variables of the serie as an \c int
     */
    int getNV() const;

    /**
     *  \brief  Gets the adress of the coefficient at order \c ord and position \c pos
     *  \return a pointer to the desired coefficient
     */
    Ofsc* getCA(int const& ord, int const& pos) const;

    /**
     *  \brief  Gets the pointer address of the Oftsc object
     */
    Oftsc* getAddress() const;

    //------------------
    //Delete
    //------------------
    /**
     *  \brief Default destructor of the class Oftsc. WARNING: potential memory leak here, through the terms of type Oftsh.
     */
    ~Oftsc();

    //------------------
    //Zeroing
    //------------------
    /**
     *  \brief  Sets all coefficients to zero.
     */
    void zero();

    //------------------
    //Evaluate
    //------------------
    void evaluate(cdouble const  X[], Ofsc& z, int const& m, int const& ofs_order) const;
    void sevaluate(cdouble const  X[], Ofsc& z, int const& m, int const& ofs_order) const;
    void evaluate(cdouble const  X[], Ofsc* z, int const& m) const;

    //------------------
    //Friendly streaming
    //------------------
    friend std::ostream& operator << (std::ostream& stream, Oftsc const& ofts);
    void printNonNullCoeffs();
};


//----------------------------------------------------------------------------------------
//
//          Reading & writing
//
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
// Text format, read
//----------------------------------------------------------------------------------------
/**
 * \brief Reads a given \c Ofsc  object within a \c Oftsc, in txt format.
 **/
void readOFS_txt(Ofsc &xFFT, ifstream &readStream);
/**
 * \brief Reads a given \c Oftsc  object, in txt format.
 **/
int  readOFTS_txt(Oftsc &x, string filename);
/**
 * \brief Reads a given vector W of type \c Oftsc  in a txt files of the form "filename+i.txt", with i = 0, length(W)-1.
 **/
void readVOFTS_txt(vector<Oftsc >  &W, string filename);


//----------------------------------------------------------------------------------------
// Binary format, read
//----------------------------------------------------------------------------------------
/**
 * \brief Reads a given \c Ofsc  object within a \c Oftsc, in bin format.
 **/
void readOFS_bin(Ofsc  &xFFT, fstream &myfile);
/**
 * \brief Reads a given \c Oftsc  object, in bin format.
 **/
int readOFTS_bin(Oftsc &W, string filename);
/**
 * \brief Reads a given \c Oftsc  object, in bin format, at order n
 **/
int readOFTS_bin(Oftsc &W, string filename, int n);
/**
 * \brief Reads a given vector W of type \c Oftsc  in binary files of the form "filename+i.bin", with i = 0, length(W)-1.
 **/
void readVOFTS_bin(vector<Oftsc >  &W, string filename);

//----------------------------------------------------------------------------------------
// Binary format, copy in lesser dimension series
//----------------------------------------------------------------------------------------
/**
 * \brief Reads a given \c Oftsc  object, in bin format, and transform it into another Oftsc object, of lesser dimensions
 *        More precisely:
 *        We suppose that W is a very sparse Fourier-Taylor series, with only non-null coefficients along the dimension dim.
 *        Hence, it is possible to entirely store W into a simpler series W1, with only one dimension.
 **/
int fromOFTStoOFTS_bin(Oftsc &W, Oftsc &W1, int dim);

/**
 * \brief Reads some Oftsc  objects, in bin format,
 *        and transform it into other Oftsc objects, of lesser dimensions.
 *
 *        Indeed, we know that, if the center-stable or center-unstable manifolds are
 *        computed using the parameterization method, the series in the dimensions
 *        0, and 3 of Wh (parameterization in TFC coordinates) are of the form:
 *
 *          Wh[0] = sum_(k=0)^n c_k s_5^k
 *
 *        I.e. they are FT series of only one variable, namely the variable s5,
 *        last variable of the reduced variables (s1, s2, s3, s4, s5).
 *
 *        To limit the amount of memory needed to stored these series, we can use
 *        one-dimensional FT series. This is done via the present routine.
 *
 *        This routine makes use of fromOFTStoOFTS_bin to read and store into the less dimensions objects.
 *        Then, the results are stored in binary files.
 **/
void fromVOFTStoVOFTS_bin(Oftsc &W, Oftsc &W1, string filein, string fileout);
//----------------------------------------------------------------------------------------
// Text format, write
//----------------------------------------------------------------------------------------
/**
 * \brief Writes a given vector W of type \c Oftsc  in a txt files of the form "filename+i.txt", with i = 0, length(W)-1.
 **/
void  writeVOFTS_txt(vector<Oftsc > &W, string filename);

//----------------------------------------------------------------------------------------
// Binary format, write
//----------------------------------------------------------------------------------------
/**
 * \brief Writes a given \c Ofsc  object within a \c Oftsc, in bin format.
 **/
void  writeOFS_bin(Ofsc const &xFFT, fstream &myfile);
/**
 * \brief Writes a given \c Oftsc  object, in bin format.
 **/
void  writeOFTS_bin(Oftsc  const &W, string filename);


#endif // OFTS_H
