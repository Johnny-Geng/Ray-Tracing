// A simple library that contains some basic math constants and functions, as well as a definition of 3D vector.
// When I constructed this library, I partially referred to some of my coursework in ITP 380 class.

#pragma once

#include <cmath>
#include <limits>
#include <float.h>

namespace Math
{
	const double Pi = 3.1415926535;
	const double MAX = DBL_MAX;
	const double EPI = 0.000001;
	template <typename T>
	[[nodiscard]] T Max(const T& a, const T& b)
	{
		return (a < b ? b : a);
	}

	template <typename T>
	[[nodiscard]] T Min(const T& a, const T& b)
	{
		return (a < b ? a : b);
	}

	template <typename T>
	[[nodiscard]] T Clamp(const T& value, const T& lower, const T& upper)
	{
		return Min(upper, Max(lower, value));
	}

	[[nodiscard]] inline double Sqrt(double value)
	{
		return sqrt(value);
	}

	[[nodiscard]] inline double Tan(double angle)
	{
		return tan(angle);
	}

	[[nodiscard]] inline bool NearZero(double val, double epsilon = 0.00000001)
	{
		if (fabs(val) <= epsilon)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

}

// 3D Vector
class Vector3
{
public:
	double x;
	double y;
	double z;

	Vector3()
		:x(0.0)
		,y(0.0)
		,z(0.0)
	{}

	explicit Vector3(double inX, double inY, double inZ)
		:x(inX)
		,y(inY)
		,z(inZ)
	{}

	// Set all three components in one line
	void Set(double inX, double inY, double inZ)
	{
		x = inX;
		y = inY;
		z = inZ;
	}

	// Vector addition (a + b)
	[[nodiscard]] friend Vector3 operator+(const Vector3& a, const Vector3& b)
	{
		return Vector3(a.x + b.x, a.y + b.y, a.z + b.z);
	}

	// Vector subtraction (a - b)
	[[nodiscard]] friend Vector3 operator-(const Vector3& a, const Vector3& b)
	{
		return Vector3(a.x - b.x, a.y - b.y, a.z - b.z);
	}

	// Component-wise multiplication
	[[nodiscard]] friend Vector3 operator*(const Vector3& left, const Vector3& right)
	{
		return Vector3(left.x * right.x, left.y * right.y, left.z * right.z);
	}

	// Scalar multiplication
	[[nodiscard]] friend Vector3 operator*(const Vector3& vec, double scalar)
	{
		return Vector3(vec.x * scalar, vec.y * scalar, vec.z * scalar);
	}

	// Scalar multiplication
	[[nodiscard]] friend Vector3 operator*(double scalar, const Vector3& vec)
	{
		return Vector3(vec.x * scalar, vec.y * scalar, vec.z * scalar);
	}

	// Scalar *=
	Vector3& operator*=(double scalar)
	{
		x *= scalar;
		y *= scalar;
		z *= scalar;
		return *this;
	}

	// Vector +=
	Vector3& operator+=(const Vector3& right)
	{
		x += right.x;
		y += right.y;
		z += right.z;
		return *this;
	}

	// Vector -=
	Vector3& operator-=(const Vector3& right)
	{
		x -= right.x;
		y -= right.y;
		z -= right.z;
		return *this;
	}

	// Length squared of vector
	[[nodiscard]] double LengthSq() const
	{
		return (x*x + y*y + z*z);
	}

	// Length of vector
	[[nodiscard]] double Length() const
	{
		return (Math::Sqrt(LengthSq()));
	}

	// Normalize this vector
	void Normalize()
	{
		double length = Length();
		x /= length;
		y /= length;
		z /= length;
	}

	// Normalize the provided vector
	[[nodiscard]] static Vector3 Normalize(const Vector3& vec)
	{
		Vector3 temp = vec;
		temp.Normalize();
		return temp;
	}

	// Reflect V about (normalized) N
	[[nodiscard]] static Vector3 Reflect(const Vector3& v, const Vector3& n)
	{
		return 2.0 * Vector3::Dot(v, n) * n - v;
	}

	// Dot product between two vectors (a dot b)
	[[nodiscard]] static double Dot(const Vector3& a, const Vector3& b)
	{
		return (a.x * b.x + a.y * b.y + a.z * b.z);
	}

	// Cross product between two vectors (a cross b)
	[[nodiscard]] static Vector3 Cross(const Vector3& a, const Vector3& b)
	{
		Vector3 temp;
		temp.x = a.y * b.z - a.z * b.y;
		temp.y = a.z * b.x - a.x * b.z;
		temp.z = a.x * b.y - a.y * b.x;
		return temp;
	}

	// Get distance between two points
	[[nodiscard]] static double Distance(const Vector3& a, const Vector3& b)
	{
		return (a - b).Length();
	}
};
