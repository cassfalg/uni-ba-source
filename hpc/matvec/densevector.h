#ifndef HPC_MATVEC_DENSEVECTOR_H
#define HPC_MATVEC_DENSEVECTOR_H 1

#include <cassert>
#include <cstdlib>

namespace hpc { namespace matvec {

template <typename T, typename I>
    struct DenseVectorConstView;

template <typename T, typename I>
    struct DenseVectorView;

template <typename T, typename I=std::size_t>
struct DenseVector
{
    typedef T                               ElementType;
    typedef I                               Index;

    typedef DenseVector<T,Index>            NoView;
    typedef DenseVectorConstView<T,Index>   ConstView;
    typedef DenseVectorView<T,Index>        View;

    DenseVector(Index length)
        : length(length), inc(1), data(new T[length])
    {
    }

    ~DenseVector()
    {
        delete[] data;
    }

    const ElementType &
    operator()(Index i) const
    {
        assert(i<length);
        return data[i];
    }

    ElementType &
    operator()(Index i)
    {
        assert(i<length);
        return data[i];
    }

    ConstView
    operator()(Index i, Index num, Index stride=1) const
    {
        assert(i+stride*num<=length);
        return ConstView(num, &(operator()(i)), stride);
    }

    View
    operator()(Index i, Index num, Index stride=1)
    {
        assert(i+stride*num<=length);
        return View(num, &(operator()(i)), stride);
    }

    const Index     length, inc;
    ElementType*    data;

};

template <typename T, typename I=std::size_t>
struct DenseVectorView
{
    typedef T                           ElementType;
    typedef I                           Index;

    typedef DenseVector<T,Index>            NoView;
    typedef DenseVectorConstView<T,Index>   ConstView;
    typedef DenseVectorView<T,Index>        View;

    DenseVectorView(Index length, ElementType *data, Index inc)
        : length(length), inc(inc), data(data)
    {
    }

    const ElementType &
    operator()(Index i) const
    {
        assert(i<length);
        return data[i*inc];
    }

    ElementType &
    operator()(Index i)
    {
        assert(i<length);
        return data[i*inc];
    }

    ConstView
    operator()(Index i, Index num, Index stride=1) const
    {
        assert(i+stride*num<=length);
        return ConstView(num, &(operator()(i)), inc*stride);
    }

    View
    operator()(Index i, Index num, Index stride=1)
    {
        assert(i+stride*num<=length);
        return View(num, &(operator()(i)), inc*stride);
    }

    const Index     length, inc;
    ElementType*    data;
};

template <typename T, typename I=std::size_t>
struct DenseVectorConstView
{
    typedef T                           ElementType;
    typedef I                           Index;

    typedef DenseVector<T,Index>            NoView;
    typedef DenseVectorConstView<T,Index>   ConstView;
    typedef DenseVectorView<T,Index>        View;

    DenseVectorConstView(Index length, const ElementType *data, Index inc)
        : length(length), inc(inc), data(data)
    {
    }

    const ElementType &
    operator()(Index i) const
    {
        assert(i<length);
        return data[i*inc];
    }

    ConstView
    operator()(Index i, Index num, Index stride=1) const
    {
        assert(i+stride*num<=length);
        return ConstView(num, &(operator()(i)), inc*stride);
    }

    const Index         length, inc;
    const ElementType*  data;
};


} } // namespace matvec, hpc

#endif // HPC_MATVEC_DENSEVECTOR_H
