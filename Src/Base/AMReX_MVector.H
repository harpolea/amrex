#ifndef AMREX_MVECTOR_H_
#define AMREX_MVECTOR_H_

#include <algorithm>
#include <vector>
#include <memory>
#include <AMReX_BLassert.H>
#ifdef AMREX_SPACEDIM
#include <AMReX_Array.H>
#include <AMReX_TypeTraits.H>
#endif
#ifdef AMREX_USE_CUDA
#include <AMReX_CudaAllocators.H>
#endif

namespace amrex {
/**
* \brief This class is a thin wrapper around std::vector.  Unlike vector,
* Vector::operator[] provides bound checking when compiled with
* DEBUG=TRUE.
*/
#ifdef AMREX_USE_CUDA
template <class T, class Allocator=CudaManagedAllocator<T> >
#else
template <class T, class Allocator=std::allocator<T> >
#endif
class MVector
    :
        Vector<T, Allocator>
{
public:

    using std::vector<T, Allocator>::vector;

    T& operator[] (std::size_t i)
    {
	BL_ASSERT(static_cast<long>(i) < this->size());
	return this->std::vector<T, Allocator>::operator[](i);
    }

    const T& operator[] (std::size_t i) const
    {
	BL_ASSERT(static_cast<long>(i) < this->size());
	return this->std::vector<T, Allocator>::operator[](i);
    }

    //! get access to the underlying data pointer
    T* dataPtr () { return this->data(); }
    //! get access to the underlying data pointer
    const T* dataPtr () const { return this->data(); }

    long size() const {return static_cast<long>(std::vector<T, Allocator>::size());}

};

}
#endif
