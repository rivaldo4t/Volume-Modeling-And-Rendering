
#pragma once
#include "Vector.hpp"
#include "Matrix.hpp"
//#include "Forms.h"
#include "Color.hpp"
#include <vector>
#include <map>
#include <string>
//#include <boost/shared_ptr.hpp>
#include <memory>
#include <cstring>

//#include "VolumeGrid.h"
//#include "SparseGrid.h"

namespace lux
{


//-----------------------------------------------------------------------------
// Setting up logic to be able to determine the data type of the gradient 
template <typename U>
struct GradType
{
   typedef int GType;
};

template<>
struct GradType<float>
{
   typedef Vector GType;
};

template<>
struct GradType<Vector>
{
   typedef Matrix GType;
};



//-----------------------------------------------------------------------------

template< typename U >
class Volume 
{
  public:

    Volume(){}

   virtual ~Volume(){}

   typedef U volumeDataType;
   typedef typename GradType<U>::GType volumeGradType;

   virtual const volumeDataType eval( const Vector& P ) const { volumeDataType base; return base; }
   virtual const volumeGradType grad( const Vector& P ) const { volumeGradType base; return base; }


};

typedef Volume<float>* VolumeFloatPtr;
typedef Volume<Color>* VolumeColorPtr;
typedef Volume<Vector>* VolumeVectorPtr;
typedef Volume<Matrix>* VolumeMatrixPtr;



}