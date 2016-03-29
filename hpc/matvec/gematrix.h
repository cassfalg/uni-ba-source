#ifndef HPC_MATVEC_GEMATRIX_H
#define HPC_MATVEC_GEMATRIX_H 1

#include <cassert>
#include <cstdlib>
#include <hpc/matvec/densevector.h>

namespace hpc { namespace matvec {

enum StorageOrder {
    ColMajor,
    RowMajor
};

template <typename T, typename I>
    struct GeMatrixView;

template <typename T, typename I>
    struct GeMatrixConstView;

template <typename T, typename I=std::size_t>
struct GeMatrix
{
    typedef T                              ElementType;
    typedef I                              Index;

    typedef GeMatrix<T,Index>              NoView;
    typedef GeMatrixConstView<T,Index>     ConstView;
    typedef GeMatrixView<T,Index>          View;

    typedef DenseVectorConstView<T,Index>  ConstVectorView;
    typedef DenseVectorView<T,Index>       VectorView;


    GeMatrix(Index numRows, Index numCols,
             StorageOrder order=ColMajor)
        : numRows(numRows), numCols(numCols),
          incRow(order==ColMajor ? 1: numCols),
          incCol(order==RowMajor ? 1: numRows),
          data(new T[numRows*numCols])
    {
    }

    ~GeMatrix()
    {
        delete[] data;
    }

    const ElementType &
    operator()(Index i, Index j) const
    {
        assert(i<numRows && j<numCols);
        return data[i*incRow + j*incCol];
    }

    ElementType &
    operator()(Index i, Index j)
    {
        assert(i<numRows && j<numCols);
        return data[i*incRow + j*incCol];
    }

    ConstView
    operator()(Index i, Index j, Index m, Index n) const
    {
        assert(i+m<=numRows);
        assert(j+n<=numCols);
        return ConstView(m, n, &(operator()(i,j)), incRow, incCol);
    }

    View
    operator()(Index i, Index j, Index m, Index n)
    {
        assert(i+m<=numRows);
        assert(j+n<=numCols);
        return View(m, n, &(operator()(i,j)), incRow, incCol);
    }

    ConstVectorView
    row(Index i) const
    {
        return ConstVectorView(numCols, &(operator()(i,0)), incCol);
    }

    VectorView
    row(Index i)
    {
        return VectorView(numCols, &(operator()(i,0)), incCol);
    }

    ConstVectorView
    col(Index j) const
    {
        return ConstVectorView(numRows, &(operator()(0,j)), incRow);
    }

    VectorView
    col(Index j)
    {
        return VectorView(numRows, &(operator()(0,j)), incRow);
    }

    const Index     numRows, numCols, incRow, incCol;
    ElementType*    data;
};

template <typename T, typename I=std::size_t>
struct GeMatrixView
{
    typedef T                               ElementType;
    typedef I                               Index;
    typedef GeMatrix<T,Index>               NoView;
    typedef GeMatrixConstView<T,Index>      ConstView;
    typedef GeMatrixView<T,Index>           View;

    typedef DenseVectorConstView<T,Index>   ConstVectorView;
    typedef DenseVectorView<T,Index>        VectorView;

    GeMatrixView(Index numRows, Index numCols,
                 T *data,
                 Index incRow, Index incCol)
        : numRows(numRows), numCols(numCols),
          incRow(incRow), incCol(incCol), data(data)
    {
    }

    GeMatrixView(const GeMatrixView &rhs)
        : numRows(rhs.numRows), numCols(rhs.numCols),
          incRow(rhs.incRow), incCol(rhs.incCol),
          data(rhs.data)
    {
    }

    const ElementType &
    operator()(Index i, Index j) const
    {
        assert(i<numRows && j<numCols);
        return data[i*incRow + j*incCol];
    }

    ElementType &
    operator()(Index i, Index j)
    {
        assert(i<numRows && j<numCols);
        return data[i*incRow + j*incCol];
    }

    ConstView
    operator()(Index i, Index j, Index m, Index n) const
    {
        assert(i+m<=numRows);
        assert(j+n<=numCols);
        return ConstView(m, n, &(operator()(i,j)), incRow, incCol);
    }


    View
    operator()(Index i, Index j, Index m, Index n)
    {
        assert(i+m<=numRows);
        assert(j+n<=numCols);
        return View(m, n, &(operator()(i,j)), incRow, incCol);
    }

    ConstVectorView
    row(Index i) const
    {
        return ConstVectorView(numCols, &(operator()(i,0)), incCol);
    }

    VectorView
    row(Index i)
    {
        return VectorView(numCols, &(operator()(i,0)), incCol);
    }

    ConstVectorView
    col(Index j) const
    {
        return ConstVectorView(numRows, &(operator()(0,j)), incRow);
    }

    VectorView
    col(Index j)
    {
        return VectorView(numRows, &(operator()(0,j)), incRow);
    }


    const Index     numRows, numCols, incRow, incCol;
    ElementType*    data;
};

template <typename T, typename I=std::size_t>
struct GeMatrixConstView
{
    typedef T                               ElementType;
    typedef I                               Index;
    typedef GeMatrix<T,Index>               NoView;
    typedef GeMatrixConstView<T,Index>      ConstView;
    typedef GeMatrixView<T,Index>           View;

    typedef DenseVectorConstView<T,Index>   ConstVectorView;
    typedef DenseVectorView<T,Index>        VectorView;

    GeMatrixConstView(Index numRows, Index numCols,
                      const T *data,
                      Index incRow, Index incCol)
        : numRows(numRows), numCols(numCols),
          incRow(incRow), incCol(incCol), data(data)
    {
    }

    GeMatrixConstView(const GeMatrixConstView &rhs)
        : numRows(rhs.numRows), numCols(rhs.numCols),
          incRow(rhs.incRow), incCol(rhs.incCol),
          data(rhs.data)
    {
    }

    const ElementType &
    operator()(Index i, Index j) const
    {
        assert(i<numRows && j<numCols);
        return data[i*incRow + j*incCol];
    }

    ConstView
    operator()(Index i, Index j, Index m, Index n) const
    {
        assert(i+m<=numRows);
        assert(j+n<=numCols);
        return ConstView(m, n, &(operator()(i,j)), incRow, incCol);
    }

    ConstVectorView
    row(Index i) const
    {
        return ConstVectorView(numCols, &(operator()(i,0)), incCol);
    }

    ConstVectorView
    col(Index j) const
    {
        return ConstVectorView(numRows, &(operator()(0,j)), incRow);
    }

    const Index         numRows, numCols, incRow, incCol;
    const ElementType*  data;
};



} } // namespace matvec, hpc

#endif // HPC_MATVEC_GEMATRIX_H
