//*******************************************************************
//
//   Vector.h
//
// 3D vector class in the namespace lux
//
//
//
//*******************************************************************

#pragma once
#include <cmath>
#define _USE_MATH_DEFINES
#include <math.h>

namespace lux
{
	//! Vector is a 3D vector class
	class Vector
	{
	  public:

	   Vector(){ xyz[0] = xyz[1] = xyz[2] = 0; }

	   Vector(const Vector& v)
	   { 
		  xyz[0] = v.xyz[0];
		  xyz[1] = v.xyz[1];
		  xyz[2] = v.xyz[2]; 
	   }
   
	   Vector(const float a, const float b, const float c)
	   {
		  xyz[0] = a;
		  xyz[1] = b;
		  xyz[2] = c; 
	   }

	   ~Vector(){}

	   //!  Set all three components
	   void set( const float vx, const float vy, const float vz )
	   {
		  xyz[0] = vx;
		  xyz[1] = vy;
		  xyz[2] = vz;
	   }

	   //! Add two vectors together
	   const Vector operator+        (const Vector& v) const 
	   { 
		  return Vector(xyz[0]+v.xyz[0], xyz[1]+v.xyz[1], xyz[2]+v.xyz[2]); 
	   }
  
	   //! Subtract one vector from another
	   const Vector operator-        (const Vector& v) const
	   { 
		  return Vector(xyz[0]-v.xyz[0], xyz[1]-v.xyz[1], xyz[2]-v.xyz[2]); 
	   }

	   //! Unary minus
	   friend const Vector operator- (const Vector& v)
	   { return Vector(-v.xyz[0],-v.xyz[1],-v.xyz[2]); }

	   //! Multiplication of a constant with a vector
	   friend const Vector operator* (const float w, const Vector& v)
	   { return v*w; }
	  
	   //! Multiplication of a vector with a constant
	   const Vector operator*        (const float v) const
	   { return Vector(xyz[0]*v, xyz[1]*v, xyz[2]*v); }

	   const Vector operator/        (const float v) const
	   { return Vector(xyz[0]/v, xyz[1]/v, xyz[2]/v); }

	   //! Inner product
	   const float operator*        (const Vector& v) const  
	   { return (xyz[0]*v.xyz[0] + xyz[1]*v.xyz[1] + xyz[2]*v.xyz[2]); }
  
	   //! cross product
	   const Vector operator^        (const Vector& v) const 
	   { return Vector(xyz[1]*v.xyz[2] - xyz[2]*v.xyz[1], 
			   xyz[2]*v.xyz[0] - xyz[0]*v.xyz[2], 
			   xyz[0]*v.xyz[1] - xyz[1]*v.xyz[0]); }

	   Vector& operator=       (const Vector& v)
	   { xyz[0] = v.xyz[0]; xyz[1] = v.xyz[1]; xyz[2] = v.xyz[2]; return *this; }
  
	   Vector& operator+=      (const Vector& v)
	   { xyz[0] += v.xyz[0]; xyz[1] += v.xyz[1]; xyz[2] += v.xyz[2]; return *this; }
  
	   Vector& operator-=      (const Vector& v)
	   { xyz[0] -= v.xyz[0]; xyz[1] -= v.xyz[1]; xyz[2] -= v.xyz[2]; return *this; }
  
	   Vector& operator*=      (const float v)
	   { xyz[0] *= v; xyz[1] *= v; xyz[2] *= v; return *this; }
  
	   Vector& operator/=      (const float v)
	   { xyz[0] /= v; xyz[1] /= v; xyz[2] /= v; return *this; }
  

	   const float& operator[] (const int v) const { return xyz[v]; }
			 float& operator[] (const int v)       { return xyz[v]; }
	   const float& operator() (const int v) const { return xyz[v]; }

	   const float X() const { return xyz[0]; }
	   const float Y() const { return xyz[1]; }
	   const float Z() const { return xyz[2]; }

	   const float magnitude() const 
	   { return sqrt( xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2] ); }

	   const float magnitudeSquared() const
	   { return xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2]; }
   
	   const Vector unitvector() const { return *this/magnitude(); }

	   void normalize() 
	   { float mag = magnitude(); xyz[0] /= mag; xyz[1] /= mag; xyz[2] /= mag; }

	//  Comparisons

	   const bool operator==         (const Vector& v) const
		   { return ( xyz[0]==v.xyz[0] && xyz[1]==v.xyz[1] && xyz[2]==v.xyz[2] ); }
  
	   const bool operator!=         (const Vector& v) const
		   { return ( xyz[0]!=v.xyz[0] || xyz[1]!=v.xyz[1] || xyz[2]!=v.xyz[2] ); }
  
	   const bool operator<          (const Vector& v) const
		   { return ( magnitude() < v.magnitude() ); }
  
	   const bool operator<=         (const Vector& v) const
		   { return ( magnitude() <= v.magnitude() ); }
  
	   const bool operator>          (const Vector& v) const
		   { return ( magnitude() > v.magnitude() ); }
  
	   const bool operator>=         (const Vector& v) const
		   { return ( magnitude() >= v.magnitude() ); }

	   // Test for parallel
	   const bool operator||         (const Vector& v) const
		   { return (  fabs((*this)*v) == v.magnitude()*((*this).magnitude()) ); }
   
	   Vector rodriguesRotation(const lux::Vector& axis, float angle)
	   {
		   float theta = angle;// *3.14f / 180;
		   auto v = *this;
		   return v * cos(theta) + (axis ^ v) * sin(theta) + axis * (axis * v) * (1 - cos(theta));
	   }

	  private:
	  float xyz[3];
	};
}