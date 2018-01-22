#pragma once



#ifndef PCG_MATH_DBL_H
#define PCG_MATH_DBL_H




//#define _USE_MATH_DEFINES

#include <cmath>
//#include <cfloat>
#include <limits>
//#include <array>

//#include <amp_math.h>
//#include "Constant.h"
//#include "Fraction.h"
//#include "UtilityTypefinder.h"
//#include "PCGMath.h"

using namespace std;

namespace PCG
{
	template <class Real> struct Math;
	//template<> struct Math<double>;
	//template<> struct Math<float>;


	template<> struct Math<double>
	{
	//public:

		typedef double	BaseType;
		//typedef double	RealType;
		//typedef double	ResultType;
		//typedef double&	RefercenceType;
		//typedef double*	PointerType;
		//typedef typename PCG::Utility::TypeFinder<double>::ConstArgType	ArgType;
		//typedef typename PCG::Utility::TypeFinder<double>::ConstArgType	GradeType;
		//typedef typename PCG::Utility::TypeFinder<double>::ConstArgType	SignedGradeType;

		//static const ::numeric_limits<double> limits_real;


		static const double PCG_PI;
		static const double PCG_PI_DIV_2;
		static const double PCG_1_DIV_PI;
		static const double PCG_1_DIV_2PI;
		static const double PCG_1_DIV_4PI;
		static const double PCG_4PI;
		static const double PCG_4_DIV_3PI;
		static const double PCG_2PI;

		static const double PCG_PHI;
		static const double PCG_PHI_MAIOR;
		static const double PCG_PHI_MINOR;
		static const double PCG_E;
		static const double PCG_SQRT_E;

		//static const double GAUSS_THRESHOLD(0.056375536236835322586065; /*(exp ( -pow (1.0 / (1.0 - exp (-0.5)), 2.0) * 0.5 ))*/ 
		//static const double GAUSS_MAX		( return 2.39821599079104275109315; /*(1.0 / (1.0 - exp (-0.5)))*/ 	
		//static const double MAP_GAUSS		( return 2.8757199692429309);
		//static const double MAP_INV_GAUSS	( return 0.34773900473461694);

		static const double PCG_SQRT2;
		static const double PCG_SQRT5;
		static const double PCG_SQRT_1_DIV_3;
		static const double PCG_SQRT_1_DIV_2;
		static const double PCG_1_DIV_SQRT5;
		static const double PCG_LOG_1_DIV_2;
		static const double PCG_LOG10;

		static const double PCG_1_DIV_LOG_1_DIV_2;
		static const double PCG_1_DIV_LOG_PHI;
		static const double PCG_RAD_TO_DEG;
		static const double PCG_DEG_TO_RAD;

		static const double PCG_IND;
		static const double PCG_QNAN;
		static const double PCG_NINF;
		static const double PCG_PINF;
		static const double PCG_DIGITS;
		static const double PCG_DIGITS10;
		static const double PCG_MAX;
		static const double PCG_ROOT2_MAX;
		static const double PCG_ROOT3_MAX;
		static const double PCG_ROOT4_MAX;
		static const double PCG_PRECISE_MAX;
		static const double PCG_PRECISE_BASE;
		static const double PCG_MIN;
		static const double PCG_DENORM;
		static const double PCG_EPSILON;
		static const double PCG_NEAR1;
		//static const double ZERO_TOLERANCE;
		

		static double Increment(double value) // 9
		{

			//static const numeric_limits<double> limit_double;
			//return ::nextafter(value, limit_double.infinity());
			//return ::nexttoward(value, PCG::Math<double>::PCG_PINF);
			switch (::fpclassify(value))
			{
			case FP_INFINITE:	return ::signbit(value) ? -PCG::Math<double>::PCG_MAX : PCG::Math<double>::PCG_PINF;
			case FP_SUBNORMAL:	return ::signbit(value) ? -0.0 : PCG::Math<double>::PCG_MIN;
			case FP_NAN:		return PCG::Math<double>::PCG_MIN;
			case FP_ZERO:		return PCG::Math<double>::PCG_MIN;
			default:			return value == -PCG::Math<double>::PCG_MIN ? -0.0 : ::_nextafter(value, PCG::Math<double>::PCG_PINF);
			}
		}

		static double Decrement(double value) // 12
		{
			//static const numeric_limits<double> limit_double;
			//return ::_nextafter(value, PCG::Math<double>::PCG_NINF);

			switch (::fpclassify(value))
			{
			case FP_INFINITE:	return ::signbit(value) ? PCG::Math<double>::PCG_NINF : PCG::Math<double>::PCG_MAX;
			case FP_SUBNORMAL:	return ::signbit(value) ? -PCG::Math<double>::PCG_MIN : 0.0;
			case FP_NAN:		return -PCG::Math<double>::PCG_MIN;
			case FP_ZERO:		return -PCG::Math<double>::PCG_MIN;
			default:			return value == PCG::Math<double>::PCG_MIN ? 0.0 : ::_nextafter(value, PCG::Math<double>::PCG_NINF);
			}
		}

		static double Nextafter(double value, double target) 
		{
			double target_(PCG::Math<double>::Filter(target, 0.0, PCG::Math<double>::PCG_NINF, PCG::Math<double>::PCG_PINF));
			double value_(PCG::Math<double>::Filter(value, 0.0, PCG::Math<double>::PCG_NINF, PCG::Math<double>::PCG_PINF));

			if (target_ < value_) { return Decrement(value_); }
			else if (target_ > value_) { return Increment(value_); }
			else { return value_; }
			//else { return PCG::Math<double>::Filter(value, 0.0, PCG::Math<double>::PCG_NINF, PCG::Math<double>::PCG_PINF); }
		}

		// ******** filtering anomalies ********

		// replace Zero
		static double FilterZero(double x)	// 1
		{
			return	::fpclassify(x) == FP_ZERO || ::isnan(x) ?
				(::signbit(x) ? -PCG_MIN : PCG_MIN) :
				x;
		}

		//static double FilterZero(double x, double replace_nzero, double replace_zero)	
		//{
		//	return	::fpclassify(x) == FP_ZERO ?
		//		(::signbit(x) ? -replace_nzero : replace_zero) :
		//		x;
		//}

		// replace NAN values with 0
		static double FilterNAN(double x)	// #FilterNAN, #MakeValid // 34
		{
			return	::_isnan(x) ? 0.0 : x;
		}

		// replace denormalized values with defined safe values
		static double FilterDNORM(double x, double replace_NDNORM = -0.0, double replace_DNORM = 0.0) // 10
		{
			return	::fpclassify(x) == FP_SUBNORMAL ?
				(::signbit(x) ? replace_NDNORM : replace_DNORM) :
				x;

			//if (::fpclassify(x) == FP_SUBNORMAL)
			//{
			//	if (::signbit(x))
			//	{
			//		return replace_NDNORM;
			//	}
			//	else
			//	{
			//		return replace_DNORM;
			//	}
			//}
			//else
			//{
			//	return x;
			//}

		}

		// replace NAN values with defined safe value
		static double FilterNAN(double x, double replace_NAN) // 20
		{
			return	::_isnan(x) ? replace_NAN : x;
		}

		// replace infinite values with signed maximum values
		static double FilterINF(double x) // 	#FilterINF, #MakeFinite // 1
		{
			return	::_finite(x) ? x : ::_copysign(PCG::Math<double>::PCG_MAX, x);
		}

		// replace infinite values with defined safe values // 1
		static double FilterINF(double x, double replace_PINF, double replace_NINF)
		{
			if (::_finite(x)) return x;
			else
			{
				if (::signbit(x))
				{
					return replace_NINF;
				}
				return replace_PINF;
			}

		}

		// replace all anomalies with defined safe values
		static double Filter(double x, double replace_NAN = 0.0, double replace_NINF = -PCG::Math<double>::PCG_MAX, double replace_PINF = PCG::Math<double>::PCG_MAX, double replace_NDNORM = -0.0, double replace_DNORM = 0.0, double replace_NZERO = -0.0, double replace_ZERO = 0.0)	// #Filter, #FilterAnomaly // 67
		{
			//if (::_isnan(x))
			//{
			//	return replace_NAN;
			//}
			//else if (::_finite(x)) return x;
			//else
			//{
			//	if (x < 0)
			//	{
			//		return replace_NINF;
			//	}
			//	return replace_PINF;
			//}
			switch (::fpclassify(x))
			{
			case FP_NAN:		return replace_NAN;
			case FP_INFINITE:	return ::signbit(x) ? replace_NINF : replace_PINF;
			case FP_SUBNORMAL:	return ::signbit(x) ? replace_NDNORM : replace_DNORM;
			case FP_ZERO:		return ::signbit(x) ? replace_NZERO : replace_ZERO;
			default:			return x;
			}
		}

		// ******** Basic Operations ********


		// Plus(a, b) = a + b  // for compatibility with higher order fuzzy numbers:
		static inline double Plus(double a, double b) // #Plus
		{
			return	a + b;
		}

		// Minus(a, b)  = a - b // for compatibility with higher order fuzzy numbers:
		static inline double Minus(double a, double b) // #Minus
		{
			return	a - b;
		}

		// Multiply(a, b)  = a * b // for compatibility with higher order fuzzy numbers:
		static inline double Multiply(double a, double b) // #Multiply
		{
			return	a * b;
		}

		// Divide(a, b)  = a / b // for compatibility with higher order fuzzy numbers:
		static inline double Divide(double a, double b) // #Divide
		{
			return	a / b;
		}

		// logical complement of x
		//static double LogicComplement(double x) // #ComplementLogic, #LogicNot, #LogicComplement
		//{
		//	return	1.0 - x;
		//}

		// complement of x
		static double Complement(double x) // #Complement // 4
		{
			return	-x;
		}

		// multiplicative inverse of x
		static double Inverse(double x)  // #Inverse // 18
		{
			return	1.0 / x;
		}
		
		// [-∞, ∞]->[0, ∞]	// absolute value of x
		static double Abs(double x)	// #Abs // 52
		{
			//return	x < 0.0 ? -x : x;

			return	::fabs(x);
		}

		// [-∞, ∞]->[-∞, ∞]	// sigmoid continous exponent function (symmetric)
		static double Pow(double x, double exponent) // #PowSigmoid, #PowContinuous, #PlotPow, #Pow // 30
		{
			return	::_copysign(::pow(::fabs(x), exponent), x);
		}

		static double Sqrt(double x)	// [-∞, ∞]->[-∞, ∞]	// sigmoid continous exponent function
		{
			return	::_copysign(::sqrt(::fabs(x)), x);
		}

		static double Cbrt(double x)	// [-∞, ∞]->[-∞, ∞]	// sigmoid continous exponent function: #PowSigmoid, #PowContinuous, #PlotPow, #Pow
		{
			return	::_copysign(::cbrt(::fabs(x)), x);
		}

		// safe natural logarithmic function (symmetric)
		static double Log(double x)	// #Log
		{
			// ALGORITHM 1:
			//return	::log(x < 0.0 ? 0.0 : x);
			//return	::log(::fabs(x));

			// ALGORITHM 2: fast / symmetric
			if (::signbit(x))
			{
				return	-::log(::fabs(x));
			}
			else
			{
				return	::log(x);
			}
		}

		static double Log(double x, double numeric_base)
		{
			// ALGORITHM (1):
			//return	::log2 (x) / ::log2 (numeric_base); // TODO: safe? 
			//return	::log2 (x < 0.0 ? 0.0 : x) / ::log2 (numeric_base < 0.0 ? 0.0 : numeric_base);
			//return	::log2(x < 0.0 ? 0.0 : x) / ::log2(::fabs(numeric_base));

			// ALGORITHM 2: symmetric
			if (::signbit(x))
			{
				if (::signbit(numeric_base))
				{
					return	-::log2(::fabs(x)) / -::log2(::fabs(numeric_base));
				}
				else
				{
					return	-::log2(::fabs(x)) / ::log2(numeric_base);
				}
			}
			else
			{
				if (::signbit(numeric_base))
				{
					return	::log2(x) / -::log2(::fabs(numeric_base));
				}
				else
				{
					return	::log2(x) / ::log2(numeric_base);
				}
			}



			// ALGORITHM (2): safe
			//return	FilterNAN 
			//		(
			//			::log2(x < 0.0 ? 0.0 : x) / ::log2(::fabs(numeric_base))
			//		);
		}


		static double Log2(double x)	// #Log2
		{
			// ALGORITHM 1:
			//return	::log2(x < 0.0 ? 0.0 : x);

			// ALGORITHM 2: fast / symmetric
			if (::signbit(x))
			{
				return	-::log2(::fabs(x));
			}
			else
			{
				return	::log2(x);
			}
		}

		static double Log10(double x)	// #Log10
		{
			// ALGORITHM 1:
			//return	::log10(x < 0.0 ? 0.0 : x);

			// ALGORITHM 2: fast / symmetric
			if (::signbit(x))
			{
				return	-::log10(::fabs(x));
			}
			else
			{
				return	::log10(x);
			}
		}



		static double CopySign(double x, double y)	// #CopySign, #CopySgn
		{
			// ALGORITHM (1):
			return	::_copysign(x, y);

			// ALGORITHM (2):
			/*
			return	y ? x * (y / :::fabs (y)) : x;
			*/
		}

		static double Sign(double x)	// #Sign, #Sgn
		{
			// Alogrithm 1:
			//return	::_copysign (x != 0.0, x);

			// Alogrithm 2:
			return (x > 0.0) - (x < 0.0);
		}



		

		// ******** sorting ********

		// 2 numbers: Min Max
		// 3 numbers: Min Mid Max
		// 4 numbers: Min MidMin MidMax Max

		static double Min(double a, double b)	// #Min
		{
			// ALGORITHM 1:
			return	a < b ? a : b;

			// ALGORITHM 2: slow?
			//return	Concurrency::fast_math::fmin(a, b);
		}

		static double Min(double a, double b, double c)
		{
			return	Min(Min(a, b), c);
		}

		static double Min(double a, double b, double c, double d)
		{
			return	Min(Min(Min(a, b), c), d);
		}

		static double Max(double a, double b)	// #Max
		{
			// ALGORITHM 1:
			return	a < b ? b : a;

			// ALGORITHM 2: slow?
			//return	Concurrency::fast_math::fmax(a, b);
		}

		static double Max(double a, double b, double c)
		{
			return	Max(Max(a, b), c);
		}

		static double Max(double a, double b, double c, double d)
		{
			return	Max(Max(Max(a, b), c), d);
		}

		static double Mid(double a, double b, double c)	// #Mid // #Core
		{
			return	Max
			(
				Max
				(
					Min(a, b),
					Min(b, c)
				),
				Min(a, c)
			);
		}

		static double MidMax(double a, double b, double c, double d)	// #MidMax, #Mid2, #Mid3of4, #CoreMax
		{
			return	Max
			(
				Max
				(
					Max
					(
						Max
						(
							Max
							(
								Min(a, b),
								Min(a, c)
							),
							Min(a, d)
						),
						Min(b, c)
					),
					Min(b, d)
				),
				Min(c, d)
			);
		}

		static double MidMin(double a, double b, double c, double d)	// #MidMin, #Mid1, #Mid2of4, #CoreMin
		{
			return	Min
			(
				Min
				(
					Min
					(
						Min
						(
							Min
							(
								Max(a, b),
								Max(a, c)
							),
							Max(a, d)
						),
						Max(b, c)
					),
					Max(b, d)
				),
				Max(c, d)
			);
		}


		static double MinAbs(double a, double b)	// #MinAbs
		{
			return	::fabs(a) < ::fabs(b) ? a : b;
		}

		static double MinAbs(double a, double b, double c)
		{
			return	MinAbs(MinAbs(a, b), c);
		}

		static double MinAbs(double a, double b, double c, double d)
		{
			return	MinAbs(MinAbs(MinAbs(a, b), c), d);
		}

		static double MaxAbs(double a, double b)	// #MaxAbs
		{
			return	::fabs(a) < ::fabs(b) ? b : a;
		}

		static double MaxAbs(double a, double b, double c)
		{
			return	MaxAbs(MaxAbs(a, b), c);
		}

		static double MaxAbs(double a, double b, double c, double d)
		{
			return	MaxAbs(MaxAbs(MaxAbs(a, b), c), d);
		}

		static double MidAbs(double a, double b, double c)	// #MidAbs
		{
			return	MaxAbs
			(
				MaxAbs
				(
					MinAbs(a, b),
					MinAbs(b, c)
				),
				MinAbs(a, c)
			);
		}

		static double MidMaxAbs(double a, double b, double c, double d)	//  #MidMaxAbs, #Mid2Abs
		{
			return	MaxAbs
			(
				MaxAbs
				(
					MaxAbs
					(
						MaxAbs
						(
							MaxAbs
							(
								MinAbs(a, b),
								MinAbs(a, c)
							),
							MinAbs(a, d)
						),
						MinAbs(b, c)
					),
					MinAbs(b, d)
				),
				MinAbs(c, d)
			);
		}

		static double MidMinAbs(double a, double b, double c, double d)	//  #MidMinAbs, #Mid1Abs
		{
			return	MinAbs
			(
				MinAbs
				(
					MinAbs
					(
						MinAbs
						(
							MinAbs
							(
								MaxAbs(a, b),
								MaxAbs(a, c)
							),
							MaxAbs(a, d)
						),
						MaxAbs(b, c)
					),
					MaxAbs(b, d)
				),
				MaxAbs(c, d)
			);
		}

		// ******** Clamp ********

		// Clamp(x, min, max) = min <= x <= max
		static double Clamp(double x, double min, double max)	// #Clamp
		{
			// ALGORITHM (1): Does NOT work with Pulse!
			//return (x < min) ? min : ((max < x) ? max : x);	

			// ALGORITHM (2): nested
			return	Min(Max(x, min), max);
		}

		// Clamp(x)) = 0 <= x <= 1
		static double Clamp(double x)
		{
			return	Min(Max(x, 0.0), 1.0);
		}

		// Clamp(x, max)) = 0 <= x <= max
		static double Clamp(double x, double max)
		{
			return	Min(Max(x, 0.0), max);
		}

		//int Clamp (int x, int a, int b)
		//{
		//	// ALGORITHM (1):
		//	return (x < a) ? a : ((b < x) ? b : x);

		//	// ALGORITHM (2): nested
		//	/*
		//	return	Min (Max (x, a), b);
		//	*/
		//}





		// ******** linear interpolation ********

		// Linear Interpolation:
		// y = y0 + x * (y1 - y0) 
		// x[0, 1] - >y[y0, y1]
		static double Lerp(double x, double y0, double y1)	// #Interpolate, #Lerp
		{
			return	y0 + x * (y1 - y0);

			/*
			Input (x, y0, y1)	| Output (y)			| Description
			x < 0				  y < y0
			0 <= x <= 1			  y0 <= y <= y1
			x > 1				  y > y1
			y0 == y1			  y0
			*/
		}

		// Inverse Linear Interpolation:
		// x = (y - y0) / (y1 - y0) 
		// y[y0, y1] -> x[0, 1]
		static double LerpInv(double y, double y0, double y1)	// #LerpInv
		{
			// ALGORITHM (1):
			return	(y - y0) / (y1 - y0);

			// ALGORITHM (2): Does NOT work with Pulse or Ramp!
			//return	FilterNAN ((y - y0) / (y1 - y0));

			/*
			Input (y, y0, y1)	| Output (x)			| Description
			y < y0				x < 0
			y0 <= y <= y1		0 <= x <= 1
			y > 1					x > 1
			y0 == y1				INF
			*/
		}

		// Signed Linear Interpolation:
		// y = y0 + (x + 1.0) * 0.5 * (y1 - y0) 
		// x[-1, 1] -> y[y0, y1]
		static double LerpSgn(double x, double y0, double y1)	// level0	#LerpSgn, #LerpSigned
		{
			return	y0 + (x + 1.0) * 0.5 * (y1 - y0);

			/*
			Input (x, y0, y1)	| Output (y)			| Description
			x < -1				y < y0
			-1 <= x <= 1		y0 <= y <= y1
			x > 1					y > y1
			y0 == y1				y0
			*/
		}

		// Inverse Signed Linear Interpolation: 
		// x = (2.0 * (y - y0) / (y1 - y0)) - 1.0 
		// y[y0, y1] -> x[-1, 1]
		static double LerpSgnInv(double y, double y0, double y1)
		{
			return	(2.0 * (y - y0) / (y1 - y0)) - 1.0;

			/*
			Input (y, y0, y1)	| Output (x)			| Description
			y < y0				  x < -1
			y0 <= y <= y1         -1 <= x <= 1
			y > 1				  x > 1
			y0 == y1			  INF
			*/
		}

		// Generalized Linear Interpolation: 
		// y = y0 + (x - x0) * (y1 - y0) / (x1 - x0)
		// x[x0, x1] -> y[y0, y1]
		static double LerpGen(double x, double x0, double x1, double y0, double y1) // #Transform, #Translate, #Interpolate, #Lerp
		{
			// ALGORITHM 1: fast
			if (x0 == x1 && x == x0)
			{
				return	(y0 + y1)*0.5;
			}
			else
			{
				return	y0 + (x - x0) * (y1 - y0) / (x1 - x0);
			}

			// ALGORITHM 2:
			//return	 x0 == x1 && x == x0 ? (y0+y1)*0.5 : y0 + (x - x0) * (y1 - y0) / (x1 - x0);

			//return	 y0 + (x - x0) * (y1 - y0) / (x1 - x0);	// ERROR: (x == x0 == x1) -> PCG_IND
			//return	FilterNAN (y0 + (x - x0) * (y1 - y0) / (x1 - x0), (y0+y1)*0.5);
		}

		// Inverse Generalized Linear Interpolation: 
		// x = x0 + (y - y0) * (x1 - x0) / (y1 - y0)
		// y[y0, y1] -> x[x0, x1]
		static double LerpGenInv(double y, double x0, double x1, double y0, double y1)
		{
			// ALGORITHM 1: fast
			if (y0 == y1 && y == y0)
			{
				return	(x0 + x1)*0.5;
			}
			else
			{
				return	x0 + (y - y0) * (x1 - x0) / (y1 - y0);
			}

			// ALGORITHM 2:
			//return	 y0 == y1 && y == y0 ? (x0+x1)*0.5 : x0 + (y - y0) * (x1 - x0) / (y1 - y0);

			//return	x0 + (y - y0) * (x1 - x0) / (y1 - y0);	// ERROR: (y == y0 == y1) -> PCG_IND
			//return	FilterNAN (x0 + (y - y0) * (x1 - x0) / (y1 - y0), (x0+x1)*0.5);
		}

		// ******** allometry ********

		// y = a * pow (x, b) + c;
		// a = (y - c) / pow (x, b);
		// b = log ((y - c) / a) / log (x);
		// c = y - a * pow (x, b);
		// x = pow ((y - c) / a, 1 / b);

		// Allometric Interpolation:
		// y = (max_val - min_val) * pow (x, exponent) + min_val;
		static double AlloLerp(double x, double y0, double y1, double exponent) // #NLerp ?, #AlloLerp
		{
			// y = (y1 - y0) * pow (x, exponent) + y0;
			//return	FilterNAN ((y1 - y0) * ::pow (x, exponent) + y0);
			//return	(y1 - y0) * ::copysign (::pow (::fabs(x), exponent), x) + y0;

			return	(y1 - y0) * PCG::Math<double>::Pow(x, exponent) + y0;
		}

		// Inverse Allometric Interpolation:
		// x = pow ((y - min_val) / (max_val - min_val), 1 / exponent); 
		static double AlloLerpInv (double y, double y0, double y1, double exponent)
		{
			// x = pow ((y - y0) / (y1 - y0), 1 / exponent);
			//return	FilterNAN (::pow ((y - y0) / (y1 - y0), 1.0 / exponent));
			//double y_ ((y - y0) / (y1 - y0));
			//return	::copysign (::pow (::fabs(y_), 1.0 / exponent), y_);

			return	PCG::Math<double>::Pow((y - y0) / (y1 - y0), 1 / exponent);
		}

		static double AlloLerpInvE (double x, double y0, double y1, double y)
		{
			// y = (y1 - y0) * pow (x, exponent) + y0;
			// (y - y0) = (y1 - y0) * pow (x, exponent);
			// (y - y0) / (y1 - y0) =  pow (x, exponent);
			// log ((y - y0) / (y1 - y0)) / log (x) = exponent
			//return  log ((y - y0) / (y1 - y0)) / log (x);		
			//double y_ ((y - y0) / (y1 - y0));
			//return	::copysign (log (::fabs (y_)), y_) / ::copysign (log (::fabs(x)), x);

			return PCG::Math<double>::FilterNAN(PCG::Math<double>::Log((y - y0) / (y1 - y0), x), 1.0);
		}

		static double Allometric(double x, double a, double b, double c)	// #Allometric
		{
			//return	FilterNAN (a * ::pow (x, b) + c);
			//return	FilterNAN (a * ::copysign (::pow (::fabs(x), b), x) + c);
			return	a * PCG::Math<double>::Pow(x, b) + c;
		}

		static double AllometricA(double x, double y, double b, double c)	// #AllometricA, #AllometricInvA, #AllometricScale
		{
			//return	FilterNAN ((y - c) / ::pow (x, b));
			//return	PCG::Math<double>::SolveDivide(PCG::Math<double>::SolveMinus(y, c), PCG::Math<double>::Pow(x, b));
			return	(y - c) / PCG::Math<double>::Pow(x, b);
		}

		static double AllometricB(double x, double y, double a, double c)	// #AllometricB, #AllometricInvB, #AllometricPower
		{
			//return	FilterNAN (::log ((y - c) / a) / ::log (x));
			//return	PCG::Math<double>::Log((y - c) / a) / PCG::Math<double>::Log(x);
			//return	PCG::Math<double>::SolveLog(PCG::Math<double>::SolveDivide(PCG::Math<double>::SolveMinus(y, c), a), x);
			return	PCG::Math<double>::FilterNAN(PCG::Math<double>::Log((y - c) / a, x), 1.0);
		}

		static double AllometricC(double x, double y, double a, double b)	// #AllometricC, #AllometricInvC, #AllometricDisplacement
		{
			//return	FilterNAN (y - a * ::pow (x, b));
			//return	a * PCG::Math<double>::Pow(x, b) + c;
			//return	PCG::Math<double>::SolveMinus(y, a * PCG::Math<double>::Pow(x, b)) =  c;
			return	y - a * PCG::Math<double>::Pow(x, b);
		}

		static double AllometricX(double y, double a, double b, double c)	// #AllometricX, #AllometricInv
		{
			//return	FilterNAN (::pow ((y - c) / a, 1 / b));
			//double y_ ((y - c) / a);
			//return	FilterNAN (::copysign (::pow (::fabs(y_), 1.0 / b), y_));
			//return	PCG::Math<double>::SolvePow(PCG::Math<double>::SolveDivide(PCG::Math<double>::SolveMinus(y, c), a), PCG::Math<double>::Inverse(b));
			return	PCG::Math<double>::Pow((y - c) / a, 1 / b);
		}


		// ******** Pattern Generation ********

		static double Step(double x, double a)	// #Step
		{
			return	x >= a;
		}

		static double HStep(double x, double a)	// #HStep
		{
			if (x < a)
			{
				return  0.0;
			}
			else if (x > a)
			{
				return  1.0;
			}
			else
			{
				return  0.5;
			}
		}

		static double Ramp(double x, double a, double b)	// #Ramp
		{
			// ALGORITHM (1): fast
			if (::_finite(a))
			{
				if (a < b)
				{
					return	x < a ? 0.0 : x < b ? (x - a) / (b - a) : 1.0;
				}
				else if (b < a)
				{
					return	x < b ? 1.0 : x < a ? (x - a) / (b - a) : 0.0;
				}
				else // if (a == b)
				{
					return	x < a ? 0.0 : 1.0;
				}
			}
			else
			{
				return	a != b || (a < 0.0 || x == b);
			}

			// ALGORITHM (2): nested
			//return	Clamp 
			//		(
			//			LerpInv (x, a, b), 
			//			0.0, 
			//			1.0
			//		);
		}

		//double Ramp (double x, double a, double b, double c, double d)
		//{
		//	return	PCG::Math<double>::Ramp(x, a, b)*0.5 + PCG::Math<double>::Ramp(x, c, d)*0.5;
		//}

		static double HRamp(double x, double a, double b)	// #HRamp
		{
			// ALGORITHM (1):
			if (::_finite(a))
			{
				if (a < b)
				{
					return	x < a ? 0.0 : x < b ? (x - a) / (b - a) : 1.0;
				}
				else if (b < a)
				{
					return	x < b ? 1.0 : x < a ? (x - a) / (b - a) : 0.0;
				}
				else // if (a == b)
				{
					return	x < a ? 0.0 : x == a ? 0.5 : 1.0;
				}
			}
			else
			{
				if (a != b)
				{
					return	1.0;
				}
				else
				{
					if (x == b)
					{
						return	0.5;
					}
					else
					{
						return	a < 0.0;
					}
				}
			}

			// ALGORITHM (2): nested version
			/*
			return	Clamp
			(
			FilterNAN (LerpInv (x, a, b), 0.5),
			0.0,
			1.0
			);
			*/
		}


		//double HRamp (double x, double a, double b, double c, double d)
		//{
		//	return	PCG::Math<double>::HRamp(x, a, b)*0.5 + PCG::Math<double>::HRamp(x, c, d)*0.5;
		//}

		// ******** Wave Pattern Generation ********

		// WaveSin(origin, length, origin) = 0.0;
		static double WaveSin(double x, double length, double origin)	// #WaveSin
		{
			static const double PCG_2PI(2.0 * ::acos(-1.0));

			return	::sin(PCG_2PI * (x - origin) / length);
		}

		static double WaveSinInv(double y, double length, double origin)
		{
			static const double TWO_PI_INV(1.0 / (2.0 * ::acos(-1.0)));

			return	origin + ::asin(y) * length  * TWO_PI_INV;
		}

		// WaveCos(origin, length, origin) = 1.0;
		static double WaveCos(double x, double length, double origin)	// #WaveCos
		{
			static const double PCG_2PI(2.0 * ::acos(-1.0));

			return	::cos(PCG_2PI * (x - origin) / length);
		}

		static double WaveCosInv(double y, double length, double origin)
		{
			static const double TWO_PI_INV(1.0 / (2.0 * ::acos(-1.0)));

			return	origin + ::acos(y) * length  * TWO_PI_INV;
		}

		static double UWaveSin(double x, double length, double origin)	// #UWaveSin, #WaveSinUnsigned
		{
			return	LerpInv(WaveSin(x, length, origin), -1, 1);
		}

		static double UWaveSinInv(double y, double length, double origin)
		{
			return	WaveSinInv(Lerp(y, -1, 1), length, origin);
		}

		static double UWaveCos(double x, double length, double origin)	// #UWaveCos, #WaveCosUnsigned
		{
			return	LerpInv(WaveCos(x, length, origin), -1, 1);
		}

		static double UWaveCosInv(double y, double length, double origin)
		{
			return	WaveCosInv(Lerp(y, -1, 1), length, origin);
		}

		static double WavePulse(double x, double length, double origin)	// #WavePulse, #PulseWave ?
		{
			return	0.5 - WaveCos(x, length, origin) * 0.5;
		}

		static double WavePulseInv(double y, double length, double origin)
		{
			return	WaveCosInv(2.0 * (0.5 - y), length, origin);
		}

		// ******** Rectangular Pattern Generation ********

		static double Pulse2(double x, double a, double b)	// #Pulse2
		{
			// ALGORITHM (1):

			return	x >= a && x <= b;

			// ALGORITHM (2): nested version
			/*
			return	Min // Or
			(
			Step (x, a), // x >= a
			Step (b, x)  // b >= x
			);
			*/

			// ALGORITHM (3): safe and nested
			/*
			return	Min // Or
			(
			Step  // x >= a
			(
			x,
			Min (a, b)
			),
			Step   // b >= x
			(
			Max (a, b),
			x
			)
			);
			*/
		}

		static double HPulse2(double x, double a, double b)	// #HPulse2
		{

			if (x < a || x > b)
			{
				return  0.0;
			}
			else if (x > a && x < b)
			{
				return  1.0;
			}
			else
			{
				return  0.5;
			}
		}

		static double Pulse2Density(double a, double b)	// #Pulse2Density, #Pulse2Density
		{
			return	1.0 / (1.0 + (b - a));
		}


		// ******** Triangular Pattern Generation ********

		static double Pulse3(double x, double a, double b, double c)	// #Pulse3
		{
			return	Clamp
			(
				Min
				(
					LerpInv(x, a, b),
					LogicNot(LerpInv(x, b, c))
				),
				0.0,
				1.0
			);
		}

		static double HPulse3(double x, double a, double b, double c)	// #HPulse3
		{
			// ALGORITHM (1):
			if ((a != c) && ((a == b) && (x == a)) || ((b == c) && (x == b)))
			{
				return  0.5;
			}
			else
			{
				return	Pulse3(x, a, b, c);
			}

			// ALGORITHM (2):
			//return	Clamp
			//		(
			//			Min
			//			(
			//				FilterNAN (LerpInv (x, a, b),  0.5),
			//				FilterNAN (LogicNot (LerpInv (x, b, c)),  0.5)
			//			),
			//			0.0, 
			//			1.0
			//		);
		}

		// center of gravity for the triangle[a,b,c]
		static double Pulse3Mean(double a, double b, double c)	// #Pulse3Center	// #Pulse3Mean // #Pulse3COG
		{
			return	(a+b+c) / 3.0;
		}

		static double Pulse3Median(double a, double b, double c)	// #Pulse3Median
		{
			double A1((b - a) * 0.5);
			double A2((c - b) * 0.5);
			double Am((A1 + A2) * 0.5);
			double s;
			double m;

			if (A1 < A2)
			{
				s = ::sqrt(Am / A2);
				m = Lerp(s, c, b);
				return	FilterNAN(m, b);
			}
			else
			{
				s = ::sqrt(Am / A1);
				m = Lerp(s, a, b);
				return	FilterNAN(m, b);
			}
		}



		static double Pulse3Fractile(double f, double a, double b, double c)	// #Pulse3Quantile // #Pulse3Fractile
		{
			double A1((b - a) * 0.5);
			double A2((c - b) * 0.5);
			double Am((A1 + A2) * f);
			double s;
			double m;

			if (Am < A1)
			{
				s = ::sqrt(Am / A1);
				m = Lerp(s, a, b);
				return	FilterNAN(m, a);
			}
			else
			{
				s = ::sqrt((A1 + A2 - Am) / A2);
				m = Lerp(s, c, b);
				return	FilterNAN(m, c);
			}
		}

		static double Pulse3Density(double a, double b, double c)	// #Pulse3Density, #Pulse3Density
		{
			//return	2.0 / (c - a);

			return	2.0 / (2.0 + (c - a)); // #Crispness
		}


		//double Pulse3Distribute (double y, double a, double b, double c)	// #Pulse3Distribute
		//{
		//	double A1 ((b-a) * 0.5);
		//	double A2 ((c-b) * 0.5);
		//	double A12 (A1+A2);
		//	double y_ (y * A12);
		//
		//	if (y_ < A1)
		//	{
		//		return	Lerp (::sqrt (LerpInv (y_, 0.0, A1)), a, b);
		//	}
		//	else // if (y_ < A12)
		//	{
		//		return	Lerp (::sqrt (LerpInv (y_, A12, A1)), c, b);
		//	}
		//}
		//
		//double Pulse3Distribute (double y, double a, double b, double c, double exponent)
		//{
		//	double inv_power (1.0 / exponent);
		//	double A1 ((b-a) * inv_power);
		//	double A2 ((c-b) * inv_power);
		//	double A12 (A1+A2);
		//	double y_ (y * A12);
		//
		//	if (y_ < A1)
		//	{
		//		return	Lerp (::pow (LerpInv (y_, 0.0, A1),  inv_power), a, b);
		//	}
		//	else if (y_ < A12)
		//	{
		//		return	Lerp (::pow (LerpInv (y_, A12, A1), inv_power), c, b);
		//	}
		//	else // if (exponent == 0)
		//	{
		//		return	y < 0.5 ? a : c;	
		//	}
		//}

		// ******** Trapezoidal Pattern Generation ********

		// membership-function / height for the trapezoid[a,b,c,d] at x
		static double Pulse4(double x, double a, double b, double c, double d)	// #Pulse4
		{
			// ALGORITHM (1): fast
			if (x < a || d < x)
			{
				return	0.0;
			}
			else if (x <= b)
			{
				double r1((x - a) / (b - a));
				return	r1 < 1.0 ? r1 : 1.0;
				//return	PCG::Math<double>::Min (PCG::Math<double>::LerpInv (x, a, b), 1.0);
				//return	PCG::Math<double>::Min ((x - a) / (b - a), 1.0);
			}
			else if (x <= d)
			{
				double r2(1.0 - (x - c) / (d - c));
				return	r2 < 1.0 ? r2 : 1.0;
				//return	PCG::Math<double>::Min (1.0 - PCG::Math<double>::LerpInv (x, c, d), 1.0);
				//return	PCG::Math<double>::Min (1.0 - (x - c) / (d - c), 1.0);
			}
			else return 0.0;


			// ALGORITHM (2): nested version
			//return	/*FilterNAN
			//			(*/
			//				PCG::Math<double>::Clamp
			//				(
			//					PCG::Math<double>::Min
			//					(
			//						PCG::Math<double>::LerpInv (x, a, b), 
			//						1.0 - PCG::Math<double>::FilterNAN(PCG::Math<double>::LerpInv (x, c, d))
			//						//PCG::Math<double>::LerpInv (c-x, d, c)
			//					),
			//					0.0, 
			//					1.0
			//				);
			//			//);

			//double ab (/*PCG::Math<double>::FilterNAN*/(PCG::Math<double>::LerpInv (x, a, b)));
			//double cd ((1- PCG::Math<double>::FilterNAN(PCG::Math<double>::LerpInv (x, c, d))));
			//double min_ab_cd (PCG::Math<double>::Min (ab, cd));
			//double clamp01 (PCG::Math<double>::Clamp (min_ab_cd,0,1));
			//return clamp01;

		}

		static double HPulse4(double x, double a, double b, double c, double d)	// #HPulse4
		{
			// ALGORITHM (1):

			if ((a != d) && (((x == a) && (a == b)) || ((x == c) && (c == d))))
			{
				return  0.5;
			}
			else
			{
				return	Pulse4(x, a, b, c, d);
			}

			// ALGORITHM (2): nested version
			/*
			return	Clamp
			(
			Min
			(
			FilterNAN (LerpInv (x, a, b),  0.5),
			FilterNAN (1.0 - LerpInv (x, c, d),  0.5)
			),
			0.0,
			1.0
			);
			*/
		}

		// area for the trapezoid[a,b,c,d]
		static double Pulse4Area(double a, double b, double c, double d)	// #Pulse4Area
		{
			return	(c - b + d - a) * 0.5 ;
		}

		// area left of x for the trapezoid[a,b,c,d]
		static double Pulse4Low(double x, double a, double b, double c, double d)	// #Pulse4Low // #Pulse4Left
		{
			double x0(PCG::Math<double>::Clamp(x, a, b));
			double y0(PCG::Math<double>::Ramp(x0, a, b));
			double x1(PCG::Math<double>::Clamp(x, b, c));
			double x2(PCG::Math<double>::Clamp(x, c, d));
			double y2(PCG::Math<double>::Ramp(x2, d, c));
			double A0(0.5 * y0 * (x0 - a));	// PCG::Math<double>::SegmentArea (a, x0, 0.0, y0)
			double A1(x1 - b);	//PCG::Math<double>::SegmentArea (b, x1, 1.0, 1.0)
			double A2(0.5 * (1.0 + y2) * (x2 - c));	//PCG::Math<double>::SegmentArea (c, x2, 1.0, y2)

													//return	(::isnan (A0) ? 0.0 : A0) + (::isnan (A1) ? 0.0 : A1) + (::isnan (A2) ? 0.0 : A2); 
			return	A0 + A1 + A2;
		}

		// area right of x for the trapezoid[a,b,c,d]
		static double Pulse4High(double x, double a, double b, double c, double d)	// #Pulse4High // #Pulse4Right
		{
			double x0(PCG::Math<double>::Clamp(x, a, b));
			double y0(PCG::Math<double>::Ramp(x0, a, b));
			double x1(PCG::Math<double>::Clamp(x, b, c));
			double x2(PCG::Math<double>::Clamp(x, c, d));
			double y2(PCG::Math<double>::Ramp(x2, d, c));
			double A0(0.5 * (1.0 + y0) * (b - x0));	// PCG::Math<double>::SegmentArea (x0, b, y0, 1.0)
			double A1(c - x1);	//PCG::Math<double>::SegmentArea (x1, c, 1.0, 1.0)
			double A2(0.5 * y2 * (d - x2));	//PCG::Math<double>::SegmentArea (x2, d, y2, 0.0)

			//return	(::isnan (A0) ? 0.0 : A0) + (::isnan (A1) ? 0.0 : A1) + (::isnan (A2) ? 0.0 : A2); 
			return	A0 + A1 + A2;
		}

		static double Pulse4LowInv(double area, double a, double b, double c, double d)	// #Pulse4LowInv
		{
			double A1((b - a) * 0.5);
			double A2(c - b);
			double A3((d - c) * 0.5);
			double s;
			double x;

			if (area <= 0.0)
			{
				return a;
			}
			if (area < A1)
			{
				s = ::sqrt(area / A1);
				x = a + s * (b - a); // m = Lerp (s, a, b); // y0 + x * (y1 - y0)
				return	::isnan(x) ? a : x; // FilterNAN (x, a);
			}
			else if ((area - A1) < A2)
			{
				return b + area - A1;
			}
			else if ((area - A1 - A2)< A3)
			{
				s = ::sqrt(((A1 + A2 + A3) - area) / A3);
				x = d + s * (c - d); // x = Lerp (s, d, c); // y0 + x * (y1 - y0)
				return	::isnan(x) ? d : x; // FilterNAN (x, d);
			}
			else
			{
				return	d;
			}
		}

		static double Pulse4HighInv(double area, double a, double b, double c, double d)	// #Pulse4HighInv
		{
			double A1((b - a) * 0.5);
			double A2(c - b);
			double A3((d - c) * 0.5);
			double s;
			double x;

			if (area <= 0.0)
			{
				return d;
			}
			if (area < A3)
			{
				s = ::sqrt(area / A3);
				x = d + s * (c - d); // x = Lerp (s, d, c); // y0 + x * (y1 - y0)
				return	::isnan(x) ? d : x; // FilterNAN (x, d);
			}
			else if ((area - A3) < A2)
			{
				return c - (area - A3);
			}
			else if ((area - A3 - A2)< A1)
			{
				s = ::sqrt(((A1 + A2 + A3) - area) / A1);
				x = a + s * (b - a); // m = Lerp (s, a, b); // y0 + x * (y1 - y0)
				return	::isnan(x) ? a : x; // FilterNAN (x, a);
			}
			else
			{
				return	a;
			}
		}

		// degree of triangluar <-> rectangular shape for the trapezoid[a,b,c,d]
		static double Pulse4Shape(double a, double b, double c, double d)	// #Pulse4CentroidY, #Pulse4TriRect, #Pulse4Shape
		{
			// ALGORITHM 1: [0.3333, 0.5] centroid height = y 
			//return	a < d ? (1.0 + (b - c) / (a + b - c - d)) / 3.0 : 0.5;

			// ALGORITHM 2: [0, 1] normalized -> triangle vs rectangle
			return	 PCG::Math<double>::LerpInv(a < d ? (1.0 + (b - c) / (a + b - c - d)) / 3.0 : 0.5, 1.0/3.0, 1.0/2.0); 
		}

		static double Pulse4Density(double a, double b, double c, double d)	// 	// #Pulse4Height, #Pulse4Density ?
		{
			//return	Inverse (Pulse4Area (a, b, c, d)); // #Density
			//return	2.0 / ((d - a) + (c - b));
			//return	Inverse (1.0 + (c - b + d - a) * 0.5);
			return	2.0 / (2.0 + (c - b + d - a)); // #Crispness
		}
		
		// lateral excentricity for the trapezoid[a,b,c,d]
		static double Pulse4ExentricityLateral(double a, double b, double c, double d)	// #Pulse4ExentricityLateral
		{
			return	RatioToBias(b-a, d-c);
		}

		// central excentricity for the trapezoid[a,b,c,d]
		static double Pulse4ExentricityCentral(double a, double b, double c, double d)	// #Pulse4ExentricityCentral
		{
			return	RatioToBias(c-b, (b-a) * 0.5 + (d-c) * 0.5);
		}

		// length of membership-function for the trapezoid[a,b,c,d]
		static double Pulse4Length(double a, double b, double c, double d)	// #Pulse4Length
		{
			return	::sqrt(::pow(b - a, 2.0) + 1.0) +
				c - b +
				::sqrt(::pow(d - c, 2.0) + 1.0);
		}

		// perimeter of the trapezoid[a,b,c,d]
		static double Pulse4Perimeter(double a, double b, double c, double d)	// #Pulse4Perimeter
		{
			return	::sqrt(::pow(b - a, 2.0) + 1.0) +
				c - b +
				::sqrt(::pow(d - c, 2.0) + 1.0) +
				d - a;
		}

		// center of gravity for the trapezoid[a,b,c,d]
		static double Pulse4Mean(double a, double b, double c, double d)	// #Pulse4Center	// #Pulse4Mean // #Pulse4COG
		{
			return	(c * c - a * a - a * b - b * b + c * d + d * d) / (3.0 * (c - b + d - a));
		}

		// calculate point(x) that divides trapezoid[a,b,c,d] into equal-sized areas
		static double Pulse4Median(double a, double b, double c, double d)	// #Pulse4Median
		{
			static const double q(0.5);

			double A1((b - a) * 0.5);
			double A2(c - b);
			double A3((d - c) * 0.5);
			double A123(A1 + A2 + A3);
			double Am(A123 * q);
			double s;
			double m;

			if (Am < A1)
			{
				s = ::sqrt(Am / A1);
				m = a + s * (b - a); // m = Lerp (s, a, b); // y0 + x * (y1 - y0)
				return	::isnan(m) ? a : m; // FilterNAN (m, a);
			}
			else if (Am < (A1 + A2))
			{
				s = (Am - A1) / A2;
				m = b + s * (c - b); // m = Lerp (s, b, c); // y0 + x * (y1 - y0)
				return	m; // FilterNAN (m, (b+c)*0.5);
			}
			else
			{
				s = ::sqrt((A123 - Am) / A3);
				m = d + s * (c - d); // m = Lerp (s, d, c); // y0 + x * (y1 - y0)
				return	::isnan(m) ? d : m; // FilterNAN (m, d);
			}
		}

		// calculate intersection-point(x) of slopes ([a,b] and [d,c]) for the trapezoid[a,b,c,d]
		static double Pulse4Mode(double a, double b, double c, double d)	// #Pulse4Peak	// #Pulse4Mode
		{

			if (a < b)
			{
				if (c < d)
				{
					//double s = PCG::Geo2::RayIntersectionToScale(a, 0.0, b, 1.0, c, 1.0, d, 0.0);
					return PCG::Math<double>::Lerp((a - d) / (c + a - d - b), a, b);
				}
				else
				{
					return c;
				}
			}
			else
			{
				if (c < d)
				{
					return b;
				}
				else
				{
					return PCG::Math<double>::FilterNAN((b + c) / 2.0);
				}
			}
		}

		// calculate point(x) that divides trapezoid[a,b,c,d] into a fraction(f) and (1-f) for it's area
		static double Pulse4Fractile(double f, double a, double b, double c, double d)	// #Pulse4Distribute, #Pulse4Quantile, #Pulse4Percentile, #Pulse4Lerp, #Pulse4AreaLerp, #Pulse4Fractile
		{
			double A1((b - a) * 0.5);
			double A2(c - b);
			double A3((d - c) * 0.5);
			double A123(A1 + A2 + A3);
			double Ax(A123 * f);
			double s;
			double x;

			if (Ax < A1)
			{
				s = ::sqrt(Ax / A1);
				x = a + s * (b - a); // x = Lerp (s, a, b); // y0 + x * (y1 - y0)
				return	::isnan(x) ? a : x; // FilterNAN (x, a);
			}
			else if (Ax < (A1 + A2))
			{
				s = (Ax - A1) / A2;
				x = b + s * (c - b); // x = Lerp (s, b, c); // y0 + x * (y1 - y0)
				return	x; // FilterNAN (x, (b+c)*0.5);
			}
			else
			{
				s = ::sqrt((A123 - Ax) / A3);
				x = d + s * (c - d); // x = Lerp (s, d, c); // y0 + x * (y1 - y0)
				return	::isnan(x) ? d : x; // FilterNAN (x, d);
			}
		}

		static double Pulse4DistributePow(double f, double a, double b, double c, double d, double power_inv) // cumulative distribution function, CDF? Quantile?
		{
			//double power_inv(1.0 / exponent);
			double A1((b - a) * power_inv);	// Area of the first (left) triangle
			double A2(c - b); // Area of the central rectangle (plateau, core)
			double A3((d - c) * power_inv);	// Area of the second (right) triangle
			double A123(A1 + A2 + A3);	// Area of the whole trapezoid
			double Ax(f * A123); // Area left of x
			double f_;
			double x;

			if (Ax < A1) // left slope
			{
				f_ = Ax / A1; // LerpInv (Ax, 0, A1)
				f_ = ::pow(f_, power_inv); // adjust for slope
				x = Lerp(f_, a, b); // find x
									//return	FilterNAN (Lerp (::pow (LerpInv (Ax, 0.0, A1), power_inv), a, b), a);
									//return	FilterNAN (Lerp (::pow (Ax / A1, power_inv), a, b), a);
				return	FilterNAN(x, a);
			}
			else if (Ax < (A1 + A2)) // central plateau
			{
				f_ = (Ax - A1) / A2; // LerpInv (Ax, A1, A1+A2)
				x = Lerp(f_, b, c);
				//return	Lerp (LerpInv (Ax, A1, A12), b, c);
				//return	Lerp ((Ax - A1) / A2, b, c);
				return	x;
			}
			else /*if (Ax < A123)*/ // right slope
			{
				f_ = (A123 - Ax) / A3; // LerpInv (Ax, A1+A2+A3, A1+A2)
				f_ = ::pow(f_, power_inv); // adjust for slope
				x = Lerp(f_, d, c); // find x
									//return	FilterNAN (Lerp (::pow (LerpInv (Ax, A123, A12), power_inv), d, c), d);
									//return	FilterNAN (Lerp (::pow ((A123-Ax) / A3, power_inv), d, c), d);
				return	FilterNAN(x, d);
			}

			// LerpInv(Ax, A123, A12) = (Ax - A123) / (A12 - A123);
			// LerpInv(Ax, (A1 + A2 + A3), (A1 + A2)) = ((A1 + A2 + A3) - Ax) / A3;

			// LerpInv(Ax, A1, A12) = (Ax - A1) / (A12 - A1);
			// LerpInv(Ax, A1, (A1 + A2)) = (Ax - A1) / A2;

			// LerpInv(Ax, 0, A1) = (Ax - 0) / (A1 - 0);
			// LerpInv(Ax, 0, A1) = Ax / A1;
		}

		static double Pulse4DistributeXY(double x, double y, double a, double b, double c, double d)	// #Pulse4DistributeXY
		{
			static const double half_pi(0.5 * ::acos(-1.0));
			y = LerpInv(::asin(y), 0.0, half_pi); // Linear Distribution

			return	Lerp(x, Lerp(y, a, b), Lerp(y, d, c));
		}

		// ******** Inclined Rectangular Pattern Generation ********

		static double SegmentArea(double x1, double x2, double y1, double y2)	// #SegmentArea
		{
			return	0.5 * (y1 + y2) * (x2 - x1);
		}

		static double SegmentCOG(double x1, double x2, double y1, double y2)	// #SegmentCOG, #SegmentMean
		{
			//if (y1 < y2)
			//{
			//	double area_rectangle((x2 - x1) * y1);
			//	double area_triangle(((x2 - x1) * (y2 - y1)) / 2.0);
			//	double centroid_rectangle(x1 + (x2 - x1) / 2.0);
			//	double centroid_triangle(x1 + (x2 - x1) * 2.0 / 3.0);
			//	double centroid((centroid_rectangle * area_rectangle + centroid_triangle * area_triangle) / (area_rectangle + area_triangle));
			//	//double centroid(((x1 + (x2 - x1) / 2.0) * ((x2 - x1) * y1) + (x1 + (x2 - x1) * 2.0 / 3.0) * (((x2 - x1) * (y2 - y1)) / 2.0)) / (((x2 - x1) * y1) + (((x2 - x1) * (y2 - y1)) / 2.0)));
			//	return centroid; // (x1 (2 y1 + y2) + x2 (y1 + 2 y2))/(3 (y1 + y2))
			//}
			//else 
			//{
			//	double area_rectangle((x2 - x1) * y2);
			//	double area_triangle(((x2 - x1) * (y1 - y2)) / 2.0);
			//	double centroid_rectangle(x1 + (x2 - x1) / 2.0);
			//	double centroid_triangle(x1 + (x2 - x1) * 1.0 / 3.0);
			//	//double centroid(((x1 + (x2 - x1) / 2) * ((x2 - x1) * y2) + (x1 + (x2 - x1) * 1 / 3) * (((x2 - x1) * (y1 - y2)) / 2)) / (((x2 - x1) * y2) + (((x2 - x1) * (y1 - y2)) / 2)));
			//	double centroid((centroid_rectangle * area_rectangle + centroid_triangle * area_triangle) / (area_rectangle + area_triangle));
			//	return centroid; // (x1 (2 y1 + y2) + x2 (y1 + 2 y2))/(3 (y1 + y2))
			//}
			return	(x1 * (2.0 * y1 + y2) + x2 * (y1 + 2.0 * y2)) / (3.0 * (y1 + y2));
		}

		static double SegmentAreaAbs(double x1, double x2, double y1, double y2)	// #SegmentAreaAbs
		{
			double x3(Clamp(LerpGen(0.0, y1, y2, x1, x2), x1, x2));
			double y3(LerpGen(x3, x1, x2, y1, y2));

			return	::fabs(SegmentArea(x1, x3, y1, y3)) + ::fabs(SegmentArea(x3, x2, y3, y2));
		}

		static double SegmentCOGAbs(double x1, double x2, double y1, double y2)	// #SegmentCOGAbs, #SegmentMeanAbs
		{
			double x3(Clamp(LerpGen(0.0, y1, y2, x1, x2), x1, x2));
			double y3(LerpGen(x3, x1, x2, y1, y2));
			double a1(::fabs(SegmentArea(x1, x3, y1, y3)));
			double a2(::fabs(SegmentArea(x3, x2, y3, y2)));

			return	(SegmentCOG(x1, x3, y1, y3) * a1 + SegmentCOG(x3, x2, y3, y2) * a2) / (a1 + a2);
		}

		static double SegmentFractile(double f, double x1, double x2, double y1, double y2) // #SegmentFractile
		{
			double As(SegmentArea(x1, x2, y1, y2));
			double Aq(f * As);
			double xq(LerpGenInv(0, x1, x2, y1, y2));
			double A0;
			double At;
			double s;
			double m;
			double dy(::fabs(y2 - y1));

			if (dy < PCG::Math<double>::PCG_EPSILON || !::_finite(dy))
			{
				m = x1 + f * (x2 - x1); // m = Lerp (f, x1, x2); y0 + x * (y1 - y0)
				return	m;
			}
			else if (y1 < y2)
			{
				A0 = (x1 - xq)*y1*0.5;
				At = A0 + As;
				s = ::sqrt((A0 + Aq) / At);
				m = xq + s * (x2 - xq); // m = Lerp (s, xq, x2);
				return	::isnan(m) ? x1 : m; // FilterNAN (m, x1);
			}
			else //if (y2 < y1)
			{
				A0 = (xq - x2)*y2*0.5;
				At = A0 + As;
				s = ::sqrt((At - Aq) / At);
				m = xq + s * (x1 - xq); // m = Lerp (s, xq, x1);
				return	::isnan(m) ? x2 : m; // FilterNAN (m, x2);
			}
		}

		// ********  ********

		static double PulseInterval(double a0, double a1, double b0, double b1)	// #PulseInterval
		{
			// ALGORITHM (1):

			return	(a0 >= b0 && a1 <= b1) || (a0 <= b1 && a1 >= b0);

			// ALGORITHM (2): nested version
			/*
			return	Max // Or
			(
			Min // And
			(
			Step (a0, b0),	// a0 >= b0
			Step (b1, a1)	// b1 >= a1
			),
			Min // And
			(
			Step (b1, a0),	// b1 >= a0
			Step (a1, b0)	// a1 >= b0
			)
			);
			*/
		}

		static double HPulseInterval(double a0, double a1, double b0, double b1)	// #HPulseInterval
		{
			return	Max // Or
			(
				Min // And
				(
					HStep(a0, b0), // a0 >= b0
					HStep(b1, a1)	// b1 >= a1
				),
				Min // And
				(
					HStep(b1, a0), // b1 >= a0
					HStep(a1, b0)	// a1 >= b0
				)
			);
		}

		// ******** convert scale ********

		// convert decimal digits to scale:
		static double DigitsToScale(double digits)	// #DigitsToScale
		{
			// ALGORITHM (1):
			//return	::exp (digits * PCG_LOG10);

			// ALGORITHM (2):
			return ::pow(10.0, -digits);

		}

		// convert scale to decimal digits:
		static double ScaleToDigits(double scale)	// #DigitsToScale
		{
			return -PCG::Math<double>::Log10(scale);
		}

		// convert numeric digits to scale:
		static double DigitsToScale(double digits, double numeric_base)
		{
			// ALGORITHM (1):
			//return	::exp (digits * ::log (numeric_base));

			// ALGORITHM (2):
			//return ::pow (numeric_base, -digits);
			//return PCG::Math<double>::Pow(numeric_base, -digits);
			return	::_copysign(::pow(::fabs(numeric_base), -digits), numeric_base);
		}

		// convert scale to numeric digits:
		static double ScaleToDigits(double scale, double numeric_base)	// #DigitsToScale
		{
			return -PCG::Math<double>::Log(scale, numeric_base);
		}

		// ******** Nearest Integer Operations ********

		static double Floor(double x)	// #Floor
		{
			return	::floor(x);
		}

		static double Floor(double x, double scale)
		{
			return	scale ? ::floor(x / scale) * scale : x;
		}

		static double Floor(double x, double a, double b)
		{
			// ALGORITHM (1): Fast
			return	a + Floor(x - a, b - a);

			// ALGORITHM (2): Nested
			//return	Lerp (Floor (LerpInv (x, a, b)), a, b);
		}

		static double Ceiling(double x)	// #Ceiling, #Ceil
		{
			return	::ceil(x);
		}

		static double Ceiling(double x, double scale)
		{
			return	scale ? ::ceil(x / scale) * scale : x;
		}

		static double Ceiling(double x, double a, double b)
		{
			// ALGORITHM (1): Fast
			return	a + Ceiling(x - a, b - a);

			// ALGORITHM (2): Nested
			//return	Lerp (Ceiling (LerpInv (x, a, b)), a, b);
		}

		static double Round(double x)	// #Round
		{
			// ALGORITHM (1): fast
			return	::_copysign
			(
				::floor(::fabs(x) + 0.5),
				x
			);

			// ALGORITHM (2): slow
			//return ::round(x);
		}

		static double Round(double x, double scale)
		{
			return	scale ? Round(x / scale) * scale : x;
		}

		static double Round(double x, double a, double b)
		{
			// ALGORITHM (1): Fast
			return	a + Round(x - a, b - a);

			// ALGORITHM (2): Nested
			//return	Lerp (Round (LerpInv (x, a, b)), a, b);
		}

		static double Trunc(double x)	// #Trunc
		{
			// ALGORITHM (1):
			return	::_copysign
			(
				::floor(::fabs(x)),
				x
			);

			// ALGORITHM (2): 
			//return Concurrency::precise_math::trunc(x);
		}

		static double Trunc(double x, double scale)
		{
			return	scale ? Trunc(x / scale) * scale : x;
		}

		static double Trunc(double x, double a, double b)
		{
			// ALGORITHM (1): Fast
			return	a + Trunc(x - a, b - a);

			// ALGORITHM (2): Nested
			//return	Lerp (Trunc (LerpInv (x, a, b)), a, b);
		}

		static double Nint(double x)	// #Nint
		{
			return ::lrint(x);
		}

		static double Nint(double x, double scale)
		{
			return	scale ? Nint(x / scale) * scale : x;
		}

		static double Nint(double x, double a, double b)
		{
			// ALGORITHM (1): Fast
			return	a + Nint(x - a, b - a);

			// ALGORITHM (2): Nested
			//return	Lerp (Nint (LerpInv (x, a, b)), a, b);
		}

		static double Group(double x)	// #Group
		{
			return	::floor(x + 0.5);
		}

		static double Group(double x, double width)
		{
			return	width ? ::floor((x + 0.5 * width) / width) * width : x;
		}

		static double Group(double x, double width, double offset)
		{
			return	offset + Group(x - offset, width);
		}

		// ******** fractional part ********

		static double Mod(double x)	// #FracPart, #FractionalPart , #Mod, #ModSigned
		{
			/*

			Mod (fractional part):

			/  /  /
			/  /  /
			-2 -1  0  1  2
			/  /  /
			/  /  /

			*/

			// ALGORITHM (1): fast
			//double intPart;
			//return	::modf (x, &intPart);	
			//return	::modf(x, &intPart);

			// ALGORITHM (2): fast
			return	::copysign
			(
				::fabs(x) - ::floor(::fabs(x)),
				x
			);

			// ALGORITHM (3): slow
			//return	::fmod (x, 1.0);
		}

		static double Mod(double x, double scale)
		{
			return	scale ? Mod(x / scale) * scale : scale;
		}

		static double Mod(double x, double a, double b)
		{
			return	a + Mod(x - a, b - a);
		}

		static double ModComplement(double x, double a, double b) // #WrapComplement, #ComplementWrap, #WrapComp, #CompWrap ?
		{
			return	Mod(x + (b - a) * 0.5, a, b);
		}

		// ******** remainder ********

		static double Remainder(double x)	// #Remainder
		{
			return	::remainder(x, 1.0);
		}

		static double Remainder(double x, double scale)
		{
			//return	scale ? PCG::Math<double>::Remainder(x / scale) * scale : scale;
			return	scale ? ::remainder(x, scale) : scale;
		}

		static double Remainder(double x, double a, double b)
		{
			return	a + PCG::Math<double>::Remainder(x - a, b - a);
		}

		// ******** fractional value ********

		static double Wrap(double x) // 	// #FractValue, #FractionalValue, #Wrap, #Mod, #Modulo, #Cyclic, #CyclicWrap, #ModWrap
		{
			/*

			Wrap (fractional value):

			  /  /  /  /  /
			 /  /  /  /  /
			-2 -1  0  1  2

			*/

			return x - ::floor(x);
		}

		static double Wrap(double x, double scale)
		{
			// Wrap (x, scale) = x - ::floor (x / scale) * scale


			// ALGORITHM (1): fast
			return	scale ? Wrap(x / scale) * scale : 0.0;

			// ALGORITHM (2): safe
			//return	FilterNAN (Wrap (x / scale) * scale);
		}

		static double Wrap(double x, double a, double b)
		{
			// ALGORITHM (1): fast

			return	a + Wrap(x - a, b - a);

			// ALGORITHM (2): nested
			/*
			return	Lerp
			(
			Wrap
			(
			LerpInv (x, a, b)
			),
			a,
			b
			);
			*/
		}

		static double WrapComplement(double x, double a, double b) // #WrapComplement, #ComplementWrap, #WrapComp, #CompWrap ?
		{
			return	Wrap(x + (b - a) * 0.5, a, b);
		}

		// ******** fractional complement ********

		static double Reflect(double x)	// #Reflect, #CyclicReflect, #ModReflect
		{
			/*

			Reflect (fractional complement):

			\   / \   / \   /
 			 \ /   \ /   \ /
			 -2 -1  0  1  2

			*/

			//ALGORITHM 1:
			//double integral;
			//double fractional (::modf (::fabs(x), &integral));
			//return ::fmod (integral, 2.0) ? 1.0 - fractional : fractional;


			//ALGORITHM 2:		
			double floor_x(::floor(x));
			double wrap_x(x - floor_x);
			return ::fmod(floor_x, 2.0) ? 1.0 - wrap_x : wrap_x;


		}

		static double Reflect(double x, double scale)
		{
			// ALGORITHM (1): fast
			return	scale ? Reflect(x / scale) * scale : 0.0;

			// ALGORITHM (2): safe
			//return	FilterNAN (Reflect (x / scale) * scale);		
		}

		static double Reflect(double x, double a, double b)
		{
			// ALGORITHM (1): Fast
			return	a + Reflect(x - a, b - a);

			// ALGORITHM (2): Nested
			//return	Lerp (Reflect (LerpInv (x, a, b)), a, b);
		}

		static double ReflectComplement(double x, double a, double b) // #WrapComplement, #ComplementWrap, #WrapComp, #CompWrap ?
		{
			return	Reflect(x + (b - a) * 0.5, a, b);
		}

		// ********  ********
		static double SignedToUnsigned(double x)	// #Unsigned, #SignedToUnsigned
		{
			return	0.5 + x * 0.5;
		}

		static double UnsignedToSigned(double x)	// #Signed, #UnsignedToSigned, #UnsignedInv ?
		{
			return	2.0 * (x - 0.5);
		}

		// ******** sigmoid normalization ********

		// x[-∞, ∞]->y[-1, 1]
		static double Sigmoid(double x)		// #Sigmoid, #Normalize, #Canonize, #Unify
		{
			if (::_finite(x))
			{
				return	x / (1.0 + ::fabs(x));
				//return	::copysign (1 - 1 / (1 + ::fabs(x)), x)
			}
			else
			{
				return	::copysign(1.0, x);
			}
		}

		static double SigmoidInv(double y)
		{
			return	y / (1.0 - ::fabs(y));
		}

		static double Sigmoid(double x, double exponent)
		{
			//return Sigmoid (::pow (x, exponent), ::pow (scale, exponent));

			double x_(::pow(::fabs(x), exponent));

			if (::_finite(x_))
			{
				return	::copysign(x_ / (1.0 + x_), x);
			}
			else
			{
				return	::copysign(1.0, x);
			}
		}

		static double SigmoidInv(double y, double exponent)
		{
			//return	::pow (SigmoidInv (y, ::pow (scale, exponent)), 1.0 / exponent);

			double a(y / (1.0 - ::fabs(y)));

			//double b (_copysign(::pow (::fabs(a), 1.0 / exponent), a));

			return	_copysign(::pow(::fabs(a), 1.0 / exponent), a);
		}

		static double Sigmoid(double x, double exponent, double scale)
		{
			//return Sigmoid (::pow (x, exponent), ::pow (scale, exponent));

			double x_(::pow(::fabs(x), exponent));

			if (::_finite(x_))
			{
				return	::copysign(x_ / (scale + x_), x);
			}
			else
			{
				return	::copysign(1.0, x);
			}
		}

		static double SigmoidInv(double y, double exponent, double scale)
		{
			//return	::pow (SigmoidInv (y, ::pow (scale, exponent)), 1.0 / exponent);

			double a(scale * y / (1.0 - ::fabs(y)));

			//double b (_copysign(::pow (::fabs(a), 1.0 / exponent), a));

			return	_copysign(::pow(::fabs(a), 1.0 / exponent), a);
		}

		// ********  ********

		// x[-∞, ∞]->y[0, 1]
		static double USigmoid(double x)		// #USigmoid, #SigmoidUnsigned, #SigmoidSgn, #USigmoid
		{
			return	0.5 + Sigmoid(x) * 0.5;
		}

		// y[0, 1]->x[-∞, ∞]
		static double USigmoidInv(double y)
		{
			return	SigmoidInv(2.0 * (y - 0.5));
		}

		static double USigmoid(double x, double exponent, double scale)
		{
			return	0.5 + Sigmoid(x, exponent, scale) * 0.5;
		}

		static double USigmoidInv(double y, double exponent, double scale)
		{
			return	SigmoidInv(2.0 * (y - 0.5), exponent, scale);
		}

		// ********  ********

		// [-∞, ∞]->[-1, 1]
		static double Tanh(double x)	
		{
			// Tanh (x) = SigmoidLog (2*x) = SigmoidLog (2*x, numeric_base = exp(1.0)) = SigmoidLog (x, scale = 0.5, numeric_base = exp(1.0))

			// ALGORITHM (1):
			return	::tanh(x);

			// ALGORITHM (2):
			//double e2x (::exp (2.0 * x));
			//return	(e2x - 1.0) / (e2x + 1.0);
		}

		// y[0, 1]->x[-∞, ∞]
		static double TanhInv(double y)
		{
			//return	0.5 * ::log((1.0 + y) / (1.0 - y));
			return	::atanh(y);			
		}

		// [-∞, ∞]->[-1, 1]
		static double SigmoidLog(double x)		// #SigmoidExp, #Fuzzify, #SigmoidLog
		{
			// Tanh (x) = SigmoidLog (2*x) = SigmoidLog (2*x, numeric_base = exp(1.0)) = SigmoidLog (x, scale = 0.5, numeric_base = exp(1.0))

			double exp_x(::exp(x));
			//return	_finite (exp_x) ? (exp_x - 1.0) / (exp_x + 1.0) : _copysign (1.0, x);

			if (_finite(exp_x))
			{
				return	(exp_x - 1.0) / (exp_x + 1.0);
			}
			else
			{
				return	_copysign(1.0, x);
			}
		}

		// y[0, 1]->x[-∞, ∞]
		static double SigmoidLogInv(double y)	// #SigmoidLogInv
		{
			// ALGORITHM (1): fast
			return	 ::log((1.0 + y) / (1.0 - y));

			// ALGORITHM (2): safe
			//double y_ (Clamp (y, -1.0, 1.0));
			//return	 ::log ((1.0 + y_) / (1.0 - y_));

		}

		// [-∞, ∞]->[-1, 1]
		static double SigmoidLog(double x, double scale, double numeric_base)		// #SigmoidExp //#SigmoidLog
		{
			// Tanh (x) = SigmoidLog (2*x) = SigmoidLog (2*x, numeric_base = exp(1.0)) = SigmoidLog (x, scale = 0.5, numeric_base = exp(1.0))

			//return	::FilterNAN ((::pow (numeric_base, x / scale) - 1) / (::pow (numeric_base, x / scale) + 1), ::copysign (1, x));

			double exp_x(::pow(numeric_base, x / scale));

			if (::_finite(exp_x))
			{
				return	(exp_x - 1.0) / (exp_x + 1.0);
			}
			else
			{
				return	::_copysign(1.0, x);
			}
		}

		static double SigmoidLogInv(double y, double scale, double numeric_base)
		{
			// ALGORITHM (1): fast
			return	::log2((1.0 + y) / (1.0 - y)) / ::log2(numeric_base) * scale;

			// ALGORITHM (2): safe
			//return	PCG::Math<double>::Log ((1.0 + y) / (1.0 - y), numeric_base) * scale;

		}

		// [-∞, ∞]->[-1, 1]
		static double SigmoidLog2(double x)		// #SigmoidExp, #Fuzzify, #SigmoidLog
		{
			// Tanh (x) = SigmoidLog (2*x) = SigmoidLog (2*x, numeric_base = exp(1.0)) = SigmoidLog (x, scale = 0.5, numeric_base = exp(1.0))

			double exp_x(::exp2(x));
			//return	_finite (exp_x) ? (exp_x - 1.0) / (exp_x + 1.0) : _copysign (1.0, x);

			if (_finite(exp_x))
			{
				return	(exp_x - 1.0) / (exp_x + 1.0);
			}
			else
			{
				return	_copysign(1.0, x);
			}
		}

		// y[0, 1]->x[-∞, ∞]
		static double SigmoidLog2Inv(double y)	// #SigmoidLogInv
		{
			// ALGORITHM (1): fast
			return	 ::log2((1.0 + y) / (1.0 - y));

			// ALGORITHM (2): safe
			//double y_ (Clamp (y, -1.0, 1.0));
			//return	 ::log ((1.0 + y_) / (1.0 - y_));

		}

		// [-∞, ∞]->[-1, 1]
		static double SigmoidLog2(double x, double scale)		// #SigmoidExp //#SigmoidLog
		{
			// Tanh (x) = SigmoidLog (2*x) = SigmoidLog (2*x, numeric_base = exp(1.0)) = SigmoidLog (x, scale = 0.5, numeric_base = exp(1.0))

			//return	::FilterNAN ((::pow (numeric_base, x / scale) - 1) / (::pow (numeric_base, x / scale) + 1), ::copysign (1, x));

			double exp_x(::exp2(x / scale));

			if (::_finite(exp_x))
			{
				return	(exp_x - 1.0) / (exp_x + 1.0);
			}
			else
			{
				return	::_copysign(1.0, x);
			}
		}

		static double SigmoidLog2Inv(double y, double scale)
		{
			// ALGORITHM (1): fast
			return	::log2((1.0 + y) / (1.0 - y)) * scale;

			// ALGORITHM (2): safe
			//return	PCG::Math<double>::Log ((1.0 + y) / (1.0 - y), numeric_base) * scale;

		}

		// ********  ********

		// [-∞, ∞]->[0, 1]
		static double USigmoidLog(double x)	// USigmoidLog, #SigmoidLogUnsigned, #SigmoidLogSgn // x[-∞, ∞]->y[0, 1]
		{
			return	1.0 / (1.0 + ::exp(-x));
		}

		static double USigmoidLogInv(double y)
		{
			//return	::log(y / (1.0 - y)); // Logit

			// ALGORITHM (1): fast
			return	-::log(1.0 / y - 1.0);

			// ALGORITHM (2): safe
			//return	-::log (1.0 / Clamp (y, 0.0, 1.0) - 1.0);

		}

		static double USigmoidLog(double x, double scale, double numeric_base)
		{
			return	1.0 / (1.0 + ::pow(numeric_base, -x / scale));
		}

		static double USigmoidLogInv(double y, double scale, double numeric_base)
		{
			//return	-Log (1.0 / y - 1.0, numeric_base) * scale;
			//return	-::log (1.0 / y - 1.0) / (::log (numeric_base) / scale);
			return	scale * -::log2(1.0 / y - 1.0) / ::log2(numeric_base);
			//return	-::log (1.0 / Clamp (y, 0.0, 1.0) - 1.0) / (::log (numeric_base) / scale);
		}

		// ********  ********

		static double SigmoidGamma(double x) 
		{
			double x_(::abs(x));
			return ::copysign(::sqrt(1.0 - 1.0 / (x_*x_ + 2 * x_ + 1.0)), x);
		}

		static double SigmoidGammaInv(double y)
		{
			double y_(::abs(y));
			return ::copysign(1.0 / ::sqrt(1.0 - y_*y_) - 1.0, y);
		}

		static double SigmoidGamma(double x, double exponent)
		{
			double x_(::abs(x));
			return ::copysign(::pow(1.0 - ::pow(1.0 / (x_ + 1.0), exponent), 1.0 / exponent), x);
		}

		static double SigmoidGammaInv(double y, double exponent)
		{
			double y_(::abs(y));
			return ::copysign(1.0 / ::pow(1.0 - ::pow(y_,exponent), 1.0/exponent) - 1.0, y);
		}

		static double SigmoidGamma(double x, double exponent, double scale) // #SigmoidRel, #SigmoidGamma // Relativistic Sigmoid Function
		{
			double x_(::abs(x / scale) );
			return ::copysign(::pow(1.0 - ::pow(1.0 / (x_ + 1.0), exponent), 1.0 / exponent), x);
		}

		static double SigmoidGammaInv(double y, double exponent, double scale)
		{
			double y_(::abs(y));
			return ::copysign((1.0 / ::pow(1.0 - ::pow(y_, exponent), 1.0 / exponent) - 1.0) * scale, y) ;
		}

		// ******** continuous sigmoid plotting ********

		static double PlotPow(double x, double exponent)	// [-∞, ∞]->[-∞, ∞]	// sigmoid continous exponent function: #PowSigmoid, #PowContinuous, #PlotPow
		{
			return	::copysign(::pow(::fabs(x), exponent), x);
		}

		static double PlotPowInv(double y, double exponent)
		{
			return	::copysign(::pow(::fabs(y), 1.0 / exponent), y);
		}

		static double PlotLog(double x)	// [-∞, ∞]->[-∞, ∞]	// sigmoid continuous log function: #LogSigmoid, #SymLog, #LogContinuous, #PlotLog
		{
			//return	::copysign(::log(::fabs(x) + 1.0), x);

			return	::copysign(::log1p(::fabs(x)), x);
		}

		static double PlotLogInv(double y)
		{
			//return	::copysign(::exp(::fabs(y)) - 1.0, y);

			return	::copysign(::expm1(::fabs(y)), y);
		}

		static double PlotLog(double x, double numeric_base)
		{
			// ALGORITHM (1): fast
			return	::copysign(::log2(1.0 + ::abs(x)), x) / ::log2(numeric_base);
			//return	::copysign(::log1p(::fabs(x)), x) / ::log(::fabs(numeric_base));

			// ALGORITHM (2): safe
			//return	FilterNAN (::copysign (::log (1.0 + ::abs (x)), x) / ::log (numeric_base));

			// ALGORITHM (3): symmetric numeric_base
			//double y_ (numeric_base - 1.0);
			//return	FilterNAN (::copysign (::log (1.0 + ::abs (x)), x) / ::copysign (::log (1.0 + ::abs (y_)), y_));
		}

		static double PlotLogInv(double y, double numeric_base)
		{
			// ALGORITHM (1): fast
			return	::copysign(::pow(numeric_base, ::fabs(y)) - 1.0, y);

			// ALGORITHM (2): symmetric numeric_base
			//double y_ (numeric_base - 1.0);
			//return	::copysign (::pow (::abs(y_) + 1.0, ::fabs(y)) - 1.0, y * y_);
		}

		static double PlotLog2(double x)	// [-∞, ∞]->[-∞, ∞]	// sigmoid continuous log function: #LogSigmoid, #SymLog, #LogContinuous, #PlotLog
		{
			return	::copysign(::log2(::fabs(x) + 1.0), x);
		}

		static double PlotLog2Inv(double y)
		{
			return	::copysign(::exp2(::fabs(y)) - 1.0, y);
		}

		static double PlotLog10(double x)	// [-∞, ∞]->[-∞, ∞]	// sigmoid continuous log function: #LogSigmoid, #SymLog, #LogContinuous, #PlotLog
		{
			return	::copysign(::log10(::fabs(x) + 1.0), x);
		}

		static double PlotLog10Inv(double y)
		{
			return	::copysign(::pow(10.0, ::fabs(y)) - 1.0, y);
		}

		// ******** Fuzzy Comparison ********

		static double FuzzyEqual(double x, double y)	// #FuzzyEqual
		{
			//return	::exp (	-::pow (y - x, 2.0));
			return ::exp(-::fabs(x - y));
			//return ::pow (10, -::fabs (x - y));
		}

		static double FuzzyEqual(double x, double y, double scale, double exponent, double numeric_base)
		{
			//return	::exp (	-::pow (y - x, 2.0) * s);
			//return ::exp (-::fabs ((x - y) * s));

			return	::pow(numeric_base, -::pow(::fabs(y - x) / scale, exponent)* scale);
		}

		static double FuzzyLess(double x, double y)	// #FuzzyLess
		{
			//return	LogicNot
			//			(
			//				::exp 
			//				(
			//					-::pow 
			//					(
			//						Max					
			//						(
			//							y - x, 
			//							0.0
			//						),
			//						2.0
			//					)
			//				)
			//			);

			return	LogicNot
			(
				::exp
				(
					-Max
					(
						y - x,
						0.0
					)
				)
			);
		}

		static double FuzzyLess(double x, double y, double scale, double exponent, double numeric_base)
		{
			return	LogicNot
			(
				::pow
				(
					numeric_base,
					-::pow
					(
						Max
						(
							y - x,
							0.0
						) / scale,
						exponent
					) * scale
				)
			);

		}

		static double FuzzyGreater(double x, double y)	// #FuzzyGreater
		{
			//return	LogicNot
			//			(
			//				::exp 
			//				(
			//					-::pow 
			//					(
			//						Max					
			//						(
			//							x - y, 
			//							0.0
			//						),
			//						2.0
			//					)
			//				)
			//			);

			return	LogicNot
			(
				::exp
				(
					-Max
					(
						x - y,
						0.0
					)
				)
			);

			//return 1.0 - ::exp (Min (y - x, 0.0));
		}

		static double FuzzyGreater(double x, double y, double scale, double exponent, double numeric_base)
		{
			return	LogicNot
			(
				::pow
				(
					numeric_base,
					-::pow
					(
						Max
						(
							x - y,
							0.0
						) / scale,
						exponent
					) * scale
				)
			);

			//return 1.0 - ::exp (Min ((y - x) * s, 0.0));
		}

		// ******** SLOPE ********

		// ******** bias ********

		static double SlopeBias(double x, double bias)	// #SlopeBias, #MapBias
		{
			static const double LOG_HALF_INV(1.0 / ::log(0.5));

			return	::pow
			(
				x,
				::log(bias) * LOG_HALF_INV
			);
		}

		static double SlopeBiasInv(double y, double bias)
		{
			static const double PCG_LOG_1_DIV_2(::log(0.5));

			return	::pow
			(
				y,
				PCG_LOG_1_DIV_2 / ::log(bias)
			);
		}

		static double RatioToBias(double weight_a, double weight_b) // #RatioToBias,  #WeightToBias, #BiasFromRatio
		{
			// ALGORITHM (1): safe for INF
			double a_(::_finite(weight_a) ? (_finite(weight_b) ? ::fabs(weight_a) : 0.0) : 1.0);
			double b_(::_finite(weight_b) ? (_finite(weight_a) ? ::fabs(weight_b) : 0.0) : 1.0);

			return	(a_ + b_) ? b_ / (a_ + b_) : 0.5;

			// ALGORITHM (2): fast
			//return	(weight_a + weight_b) ? weight_b / (weight_a + weight_b) : 0.5;
		}

		static double PowerToBias(double exponent)	// #PowerToBias, #BiasFromPower
		{
			// ::exp (exponent * PCG_LOG_1_DIV_2) = bias;

			static const double PCG_LOG_1_DIV_2(::log(0.5));

			return	::exp(exponent * PCG_LOG_1_DIV_2);
		}

		static double BiasToPower(double bias)	// #BiasToPower
		{
			// exponent = ::log (bias) * LOG_HALF_INV

			static const double LOG_HALF_INV(1.0 / ::log(0.5));

			return	::log(bias) * LOG_HALF_INV;
		}

		//double inv_bias (SlopeBiasInv (0.5, bias)); // Invert Bias!!!
		static double InvertBias(double bias)	// #InvertBias
		{
			return	SlopeBiasInv(0.5, bias);
		}

		static double TranslateBias(double observed, double expected) // #TranslateBias, #SlopeBiasTranslate, #SlopeBiasTransform, #SlopeProbabilityToBias, #BiasTranslate
		{
			//return	SlopeBias (expected, SlopeBias (observed, SlopeBiasInv (0.5, observed))); // Transform observed
			//return	SlopeBiasInv (observed, SlopeBias (expected, SlopeBias (observed, SlopeBiasInv (0.5, observed)))); // Get inv_bias

			//return	SlopeBias (expected, SlopeBiasInv (0.5, observed)); // Get bias
			return	SlopeBiasInv(expected, observed); // Get bias
		}

		static double TranslateBias(double observed, double expected, double gain) // #TranslateBias, #BiasTranslate
		{
			//return	SlopeBias (expected, SlopeBias (observed, SlopeBiasInv (0.5, observed))); // Transform observed
			//return	SlopeBiasInv (observed, SlopeBias (expected, SlopeBias (observed, SlopeBiasInv (0.5, observed)))); // Get inv_bias

			return	SlopeBias(expected, SlopeGain(SlopeBiasInv(0.5, observed), gain)); // Get bias
		}

		// ******** gain ********

		static double SlopeGain(double x, double gain)	// #SlopeGain, #MapGain
		{
			// SlopeSmooth(x) :=  SlopeGain (x, ::exp (::log(0.5) / 2.0)) = SlopeGain (x, 0.707106781186548)
			// SlopeSmoothStrong(x) := SlopeGain (x, ::exp (::log(0.5) / 3.0)) = SlopeGain (x, 0.793700525984100)

			if (x < 0.5)
			{
				return	0.5 * SlopeBias(x * 2.0, 1.0 - gain);
			}
			else
			{
				return	0.5 *
					(
						2.0 -
						SlopeBias
						(
							2.0 - x * 2.0,
							1.0 - gain
						)
						);
			}

		}

		static double SlopeGain(double x, double  gain, double inflection)
		{

			return	SlopeBias
			(
				SlopeGain
				(
					SlopeBiasInv
					(
						x,
						inflection
					),
					gain
				),
				inflection
			);
		}

		//double inv_gain (LogicNot (SlopeBiasInv (0.5, LogicNot (gain)))); // Invert Gain!!!
		static double InvertGain(double gain) // #InvertGain
		{
			return	LogicNot(SlopeBiasInv(0.5, LogicNot(gain)));
		}

		static double SlopeGainInv(double x, double gain)
		{
			double comp_gain(LogicNot(gain));
			double x2(x * 2.0);

			if (x < 0.5)
			{
				return	0.5 * SlopeBiasInv(x2, comp_gain);
			}
			else
			{
				return	0.5 *
					(
						2.0 -
						SlopeBiasInv
						(
							2.0 - x2,
							comp_gain
						)
						);
			}
		}

		static double SlopeGainInv(double x, double  gain, double inflection)
		{
			return	SlopeGain(x, InvertGain(gain), inflection);
		}

		// ******** generalized gain function ********

		static double SlopeTransform(double x, double bias_in, double gain, double inflection, double bias_out)	// #SlopeTransform, #MapTransform
		{
			return	PCG::Math<double>::SlopeBias(PCG::Math<double>::SlopeGain(PCG::Math<double>::SlopeBias(x, bias_in), gain, inflection), bias_out);
		}

		static double SlopeTransformInv(double y, double bias_in, double gain, double inflection, double bias_out)
		{
			return	PCG::Math<double>::SlopeBiasInv(PCG::Math<double>::SlopeGainInv(PCG::Math<double>::SlopeBiasInv(y, bias_out), gain, inflection), bias_in);
		}

		// ********  ********
		static double SlopeArc(double x)	// #MapArc, #SlopeArc
		{
			return	::sqrt(1.0 - ::pow(1.0 - x, 2.0));
		}

		static double SlopeArcInv(double y)
		{
			return	1.0 - ::sqrt(1.0 - ::pow(y, 2.0));
		}

		static double SlopeArc(double x, double bias)	// #SlopeArc, #SlopeArcBias, #SlopeBiasArc, #MapArc
		{
			double exponent(BiasToPower(bias));
			return	::pow(1.0 - ::pow(1.0 - x, 1.0 / exponent), exponent);
		}

		static double SlopeArcInv(double y, double bias)
		{
			double exponent(BiasToPower(bias));
			return	1.0 - ::pow(1.0 - ::pow(y, 1.0 / exponent), exponent);
		}

		static double SlopeSmooth(double x)	// #SlopeSmooth, #MapSmooth, 
		{
			//return	x * x * (3.0 - 2.0 * x);	// not symmetric and not invertible!
			// Use instead: SlopeGain (x, ::exp (::log(0.5) * 0.5)) = 0.707106781186548

			static const double gain(::exp(::log(0.5) / 2.0));
			return	SlopeGain(x, gain);
		}

		static double SlopeSmoothInv(double y)
		{
			static const double gain(::exp(::log(0.5) / 2.0));
			return	SlopeGainInv(y, gain);
		}

		static double SlopeSmoothStrong(double x)	// #SlopeSmoothStrong, #MapSmooth2, 
		{
			//return	6.0*::pow(x, 5.0) - 15.0*::pow(x, 4.0) + 10.0*::pow(x, 3.0);	// not symmetric and not invertible!
			// Use instead: SlopeGain (x, ::exp (::log(0.5) / 3.0)) = 0.793700525984100

			static const double gain(::exp(::log(0.5) / 3.0));
			return	SlopeGain(x, gain);
		}

		static double SlopeSmoothStrongInv(double y)
		{
			static const double gain(::exp(::log(0.5) / 3.0));
			return	SlopeGainInv(y, gain);
		}

		static double SlopeSin(double x)	// #SlopeSin, #MapSin
		{
			static const double half_pi(0.5 * ::acos(-1.0));
			return	::sin(Lerp(x, 0.0, half_pi));
		}

		static double SlopeSinInv(double y)
		{
			static const double half_pi(0.5 * ::acos(-1.0));
			return	LerpInv(::asin(y), 0.0, half_pi);
		}

		static double SlopeSmoothSin(double x)	// #SlopeSmoothSin, #MapSmoothSin
		{
			static const double half_pi(0.5 * ::acos(-1.0));
			return	0.5 * (1.0 + ::sin(Lerp(x, -half_pi, half_pi)));
		}

		static double SlopeSmoothSinInv(double y)
		{
			static const double half_pi(0.5 * ::acos(-1.0));
			return	LerpInv(::asin(y * 2.0 - 1.0), -half_pi, half_pi);
		}

		static double SlopeSplit(double x)	// #SlopeSplit, #MapSplit, SlopeSolarize, SlopeLiminality ?
		{
			return	::fabs
			(
				x - 0.5
			) * 2.0;
		}

		// ******** spline ********

		static double SplineCatRom(double x, double s0, double s1, double s2, double s3)	// #SplineCatRom 
		{
			//		| -1  2 -1  0 |
			//		|  3 -5  0  2 | * 0.5 = Catmull-Rom Spline
			//		| -3  4  1  0 |
			//		|  1 -1  0  0 |

			//= 0.5 *(	(2 * s1) + 
			//			(-s0 + s2) * x +
			//			(2*s0 - 5*s1 + 4*s2 - s3) * x*x +
			//			(-s0 + 3*s1- 3*s2 + s3) * x*x*x
			//		);

			return	0.5 *
				(
				(2.0 * s1) +
					(-s0 + s2) * x +
					(2.0 * s0 - 5.0 * s1 + 4.0 * s2 - s3) * x * x +
					(-s0 + 3.0 * s1 - 3.0 * s2 + s3) * x * x * x
					);
		}

		static double SplineCatRomInv(double y, double s0, double s1, double s2, double s3)	// #SplineCatRomInv 
		{
			//	x = (sqrt(s1^2 - 8 s1 s2 - 2 s1 s3 + 8 s1 y + 16 s2^2 - 8 s2 s3 - 16 s2 y + s3^2 + 8 s3 y) + 3 s1 - 4 s2 + s3)/(2 (s1 - 2 s2 + s3))


			return (::sqrt(::pow(s1, 2.0) - 8.0 *s1 *s2 - 2.0 *s1 *s3 + 8.0 *s1 *y + 16 * ::pow(s2, 2.0) - 8.0* s2 *s3 - 16.0* s2 *y + ::pow(s3, 2.0) + 8.0* s3* y) + 3.0 *s1 - 4.0 *s2 + s3) / (2.0* (s1 - 2.0* s2 + s3));
		}

		// ********  ********

		static double Mean(double a, double b)	// #Mean
		{
			return	0.5 * (a + b);
		}

		static double Mean(double a, double b, double c)
		{
			static const double inv3(1.0 / 3.0);
			return	inv3 * (a + b + c);
		}

		static double Mean(double a, double b, double c, double d)
		{
			return	0.25 * (a + b + c + d);
		}


		static double MeanWeighted(double a, double weight_a, double b, double weight_b)	// #MeanWeighted
		{
			double wxy(::fabs(weight_a) + ::fabs(weight_b));

			if (wxy)
			{
				return	(a * weight_a + b * weight_b) / wxy;
			}
			else
			{
				return	Mean(a, b);
			}
		}

		static double MeanWeighted(double a, double weight_a, double b, double weight_b, double c, double weight_c)
		{
			double wxyz(::fabs(weight_a) + ::fabs(weight_b) + ::fabs(weight_c));

			if (wxyz)
			{
				return (a * weight_a + b * weight_b + c * weight_c) / wxyz;
			}
			else
			{
				return	Mean(a, b, c);
			}
		}

		static double MeanWeighted(double a, double weight_a, double b, double weight_b, double c, double weight_c, double d, double weight_d)
		{
			double wxyzw(::fabs(weight_a) + ::fabs(weight_b) + ::fabs(weight_c) + ::fabs(weight_d));

			if (wxyzw)
			{
				return (a * weight_a + b * weight_b + c * weight_c + d * weight_d) / wxyzw;
			}
			else
			{
				return	Mean(a, b, c, d);
			}
		}

		// ********  ********

		static double CyclicMean(double a, double b, double m1, double m2)	// #CyclicMean, #CyclicMean
		{
			double wrap_min1(Wrap(Min(a, b), m1, m2));
			double wrap_min2(wrap_min1 + m2 - m1);
			double wrap_max(Wrap(Max(a, b), m1, m2));

			if (::fabs(wrap_max - wrap_min1) < ::fabs(wrap_max - wrap_min2))
			{
				return	Wrap(Mean(wrap_min1, wrap_max), m1, m2);
			}
			else
			{
				return	Wrap(Mean(wrap_min2, wrap_max), m1, m2);
			}
		}

		static double CyclicMeanWeighted(double a, double weight_a, double b, double weight_b, double m1, double m2)	// #CyclicMeanWeighted
		{
			double wrap_min1(Wrap(Min(a, b), m1, m2));
			double wrap_min2(wrap_min1 + m2 - m1);
			double wrap_max(Wrap(Max(a, b), m1, m2));
			double weight_min(a < b ? weight_a : weight_b);
			double weight_max(a < b ? weight_b : weight_a);

			if (::fabs(wrap_max - wrap_min1) < ::fabs(wrap_max - wrap_min2))
			{
				return	Wrap(MeanWeighted(wrap_min1, weight_min, wrap_max, weight_max), m1, m2);
			}
			else
			{
				return	Wrap(MeanWeighted(wrap_min2, weight_min, wrap_max, weight_max), m1, m2);
			}
		}

		static double CyclicLerpAscend(double x, double y0, double y1, double a, double b)	// #CyclicLerpAscend
		{
			double d(::fabs(b - a));
			double y0_(Wrap(y0, a, b));
			double y1_(Wrap(y1, a, b));

			if (y0_ < y1_)
			{
				return	Wrap(Lerp(x, y0_, y1_), a, b);
			}
			else
			{
				return	Wrap(Lerp(x, y0_, y1_ + d), a, b);
			}
		}

		static double CyclicLerpDescend(double x, double y0, double y1, double a, double b)	// #CyclicLerpDescend
		{
			double d(::fabs(b - a));
			double y0_(Wrap(y0, a, b));
			double y1_(Wrap(y1, a, b));

			if (y1_ < y0_)
			{
				return	Wrap(Lerp(x, y0_, y1_), a, b);
			}
			else
			{
				return	Wrap(Lerp(x, y0_ + d, y1_), a, b);
			}
		}

		static double CyclicLerpShort(double x, double y0, double y1, double a, double b)	// #CyclicLerpShort
		{
			double d(::fabs(b - a));
			double y0_(Wrap(y0, a, b));
			double y1_(Wrap(y1, a, b));
			double dy(::fabs(y1_ - y0_));

			if (dy < (d*0.5))
			{
				return	Wrap(Lerp(x, y0_, y1_), a, b);
			}
			else if (y0_ < y1_)
			{
				return	Wrap(Lerp(x, y0_ + d, y1_), a, b);
			}
			else
			{
				return	Wrap(Lerp(x, y0_, y1_ + d), a, b);
			}
		}

		static double CyclicLerpLong(double x, double y0, double y1, double a, double b)	// #CyclicLerpLong
		{
			double d(::fabs(b - a));
			double y0_(Wrap(y0, a, b));
			double y1_(Wrap(y1, a, b));
			double dy(::fabs(y1_ - y0_));

			if ((d*0.5) < dy)
			{
				return	Wrap(Lerp(x, y0_, y1_), a, b);
			}
			else if (y0_ < y1_)
			{
				return	Wrap(Lerp(x, y0_ + d, y1_), a, b);
			}
			else
			{
				return	Wrap(Lerp(x, y0_, y1_ + d), a, b);
			}
		}

		static double CyclicLerp(double x, double y0, double y1, double a, double b)	// #CyclicLerp
		{
			double y0_(Wrap(y0, a, b));
			double y1_(Wrap(y1, a, b));

			return	Wrap(Lerp(x, y0_, y1_), a, b);
		}

		// ******** Logic Operations ********

		/*************************************************************/
		//	binary operators:
		//	1100	= a
		//	1010	= b
		//	xxxx	= operator (a, b)
		/*************************************************************/
		//	0000	Contradiction / False (a,b) = (a && !a) && (b && !b)
		//	0000	Contradiction / FalseWeak (a,b) = (a && !a) || (b && !b)
		//	0001	Nor (a,b) = !a && !b
		//	0010	ConverseNonimplication / InhibitB (a,b) = !a && b
		//	0011	NEGATION / COMPLEMENT / NotA (a,b) = !a && (b || !b)
		//	0100	MaterialNonimplication / InhibitA (a,b) = a && !b
		//	0101	NEGATION / COMPLEMENT / NotB (a,b) = (a || !a) && !b
		//	0110	ExclusiveDisjunction / EXCLUSION / XOr (a,b) = (!a && b) || (a && !b)
		//	0111	NAnd (a,b) = (!a) || (!b)
		//	1000	CONJUNCTION / INTERSECTION / And (a,b)	= a && b
		//	1001	Biconditional / XNOR / EQ / Equivalence (a,b)	= (a && b) || (!a && !b)
		//	1010	Projection / IdentityB (a,b) = (a || !a) && b
		//	1011	MaterialImplication / Implication (a,b) = !a || b
		//	1100	Projection / IdentityA (a,b) = a && (b || !b)
		//	1101	ConverseImplication / Replication (a,b) = a || !b
		//	1110	DISJUNCTION / UNION / Or (a,b) = a || b
		//	1111	Tautology / True (a,b) = (a || !a) || (b || !b)	
		//	1111	Tautology / TrueWeak (a,b) = (a || !a) && (b || !b)
		//			Liminality(a) = Maybe (a) = Pulse3 (a, 0.0, 0.5, 1.0)
		//			Matrix = (TT && (a && b)) || (TF && (a && !b)) || (FT && (!a && b)) || (FF && (!a && !b))
		//			FuzzyMatrix = MeanWeighted (TT, (a && b), TF, (a && !b), FT, (!a && b), FF, (!a && !b))
		//			OrAnd = ((a || b) && !cruciality) || ((a && b) && cruciality)
		//			FuzzyOrAnd = Lerp (cruciality, (a || b), (a && b));
		//			Switch = (!if_a & else_c) || (if_a & then_b)
		//			FuzzySwitch = Lerp (if_a, else_c, then_b)	
		//			Blend = Switch (alpha, Matrix (a, b, TT, TF, FT, FF), a)	
		//			FuzzyBlend = FuzzySwitch (alpha, FuzzyMatrix (a, b, TT, TF, FT, FF), a)	
		//			Not(a) = LogicNot (a)
		/*************************************************************/

		// logical complement of a:
		static double LogicNot(double a)	// #LogicNot, #LogicNot
		{
			return	1.0 - a;
		}

		static double LogicLiminal(double a)	// #LogicLiminal
		{
			return	Pulse3(a, 0.0, 0.5, 1.0);
		}

		static double LogicLiminal(double a, double threshold)
		{
			return	Pulse3(a, 0.0, threshold, 1.0);
		}

		//	1110	DISJUNCTION / UNION / Or (a,b)	= a || b
		static double LogicProbOr(double a, double b)	// #LogicProbOr	
		{
			return	a + b - a * b;
		}

		//	1000	CONJUNCTION / INTERSECTION / And (a,b)	= a && b
		static double LogicProbAnd(double a, double b)	// #LogicProbAnd	
		{
			return	a * b;
		}

		//	0110	EXCLUSION / XOr (a,b)	= (!a && b) || (a && !b)
		static double LogicProbXOr(double a, double b)	// #LogicProbXOr	
		{
			return	PCG::Math<double>::LogicProbOr
			(
				PCG::Math<double>::LogicProbAnd
				(
					PCG::Math<double>::LogicNot(a),
					b
				),
				PCG::Math<double>::LogicProbAnd
				(
					a,
					PCG::Math<double>::LogicNot(b)
				)
			);
		}

		static double LogicProbOrAnd(double a, double b, double cruciality)	// #LogicProbFuzzyOrAnd
		{
			// OrAnd = ((a || b) && !cruciality) || ((a && b) && cruciality)
			return	PCG::Math<double>::LogicProbOr
			(
				PCG::Math<double>::LogicProbAnd
				(
					PCG::Math<double>::LogicProbOr(a, b),
					PCG::Math<double>::LogicNot(cruciality)
				),
				PCG::Math<double>::LogicProbAnd
				(
					PCG::Math<double>::LogicProbAnd(a, b),
					cruciality
				)
			);
		}

		static double LogicProbFuzzyOrAnd(double a, double b, double cruciality)	// #LogicProbFuzzyOrAnd
		{
			return	PCG::Math<double>::Lerp(cruciality, PCG::Math<double>::LogicProbOr(a, b), PCG::Math<double>::LogicProbAnd(a, b));
		}

		static double LogicProbMatrix(double a, double b, double TT, double TF, double FT, double FF)	// #LogicProbMatrix
		{
			// Matrix = (TT && (a && b)) || (TF && (a && !b)) || (FT && (!a && b)) || (FF && (!a && !b))
			return	PCG::Math<double>::LogicProbOr
			(
				PCG::Math<double>::LogicProbOr
				(
					PCG::Math<double>::LogicProbOr
					(
						PCG::Math<double>::LogicProbAnd(TT, PCG::Math<double>::LogicProbAnd(a, b)),
						PCG::Math<double>::LogicProbAnd(TF, PCG::Math<double>::LogicProbAnd(a, PCG::Math<double>::LogicNot(b)))
					),
					PCG::Math<double>::LogicProbAnd(FT, PCG::Math<double>::LogicProbAnd(PCG::Math<double>::LogicNot(a), b))
				),
				PCG::Math<double>::LogicProbAnd(FF, PCG::Math<double>::LogicProbAnd(PCG::Math<double>::LogicNot(a), PCG::Math<double>::LogicNot(b)))
			);

		}

		static double LogicProbFuzzyMatrix(double a, double T, double M, double F, double threshold = 0.5)	// #LogicProbFuzzyMatrix
		{
			return	PCG::Math<double>::MeanWeighted
			(
				T, PCG::Math<double>::Ramp(a, threshold, 1.0),
				M, PCG::Math<double>::Pulse3(a, 0.0, threshold, 1.0),
				F, PCG::Math<double>::Ramp(a, threshold, 0.0)
			);
		}

		static double LogicProbFuzzyMatrix(double a, double b, double TT, double TF, double FT, double FF)	// #LogicProbFuzzyMatrix
		{
			return	PCG::Math<double>::MeanWeighted
			(
				TT, PCG::Math<double>::LogicProbAnd(a, b),
				TF, PCG::Math<double>::LogicProbAnd(a, PCG::Math<double>::LogicNot(b)),
				FT, PCG::Math<double>::LogicProbAnd(PCG::Math<double>::LogicNot(a), b),
				FF, PCG::Math<double>::LogicProbAnd(PCG::Math<double>::LogicNot(a), PCG::Math<double>::LogicNot(b))
			);
		}

		static double LogicProbFuzzyMatrix(double a, double b, double c, double TTT, double TTF, double TFT, double TFF, double FTT, double FTF, double FFT, double FFF)	// #LogicProbFuzzyMatrix
		{
			double wTTT(PCG::Math<double>::LogicProbAnd(PCG::Math<double>::LogicProbAnd(a, b), c));
			double wTTF(PCG::Math<double>::LogicProbAnd(PCG::Math<double>::LogicProbAnd(a, b), 1.0 - c));
			double wTFT(PCG::Math<double>::LogicProbAnd(PCG::Math<double>::LogicProbAnd(a, 1.0 - b), c));
			double wTFF(PCG::Math<double>::LogicProbAnd(PCG::Math<double>::LogicProbAnd(a, 1.0 - b), 1.0 - c));
			double wFTT(PCG::Math<double>::LogicProbAnd(PCG::Math<double>::LogicProbAnd(1.0 - a, b), c));
			double wFTF(PCG::Math<double>::LogicProbAnd(PCG::Math<double>::LogicProbAnd(1.0 - a, b), 1.0 - c));
			double wFFT(PCG::Math<double>::LogicProbAnd(PCG::Math<double>::LogicProbAnd(1.0 - a, 1.0 - b), c));
			double wFFF(PCG::Math<double>::LogicProbAnd(PCG::Math<double>::LogicProbAnd(1.0 - a, 1.0 - b), 1.0 - c));

			return	(TTT*wTTT + TTF*wTTF + TFT*wTFT + TFF*wTFF + FTT*wFTT + FTF*wFTF + FFT*wFFT + FFF*wFFF) /
				(wTTT + wTTF + wTFT + wTFF + wFTT + wFTF + wFFT + wFFF);
		}

		static double LogicProbFuzzyMatrix(double a, double b, double TT, double TF, double FT, double FF, double TM, double MT, double FM, double MF, double MM, double threshold = 0.5)
		{
			double wTT(PCG::Math<double>::LogicProbAnd(PCG::Math<double>::Ramp(a, threshold, 1.0), PCG::Math<double>::Ramp(b, threshold, 1.0)));
			double wTF(PCG::Math<double>::LogicProbAnd(PCG::Math<double>::Ramp(a, threshold, 1.0), PCG::Math<double>::Ramp(b, threshold, 0.0)));
			double wFT(PCG::Math<double>::LogicProbAnd(PCG::Math<double>::Ramp(a, threshold, 0.0), PCG::Math<double>::Ramp(b, threshold, 1.0)));
			double wFF(PCG::Math<double>::LogicProbAnd(PCG::Math<double>::Ramp(a, threshold, 0.0), PCG::Math<double>::Ramp(b, threshold, 0.0)));
			double wTM(PCG::Math<double>::LogicProbAnd(PCG::Math<double>::Ramp(a, threshold, 1.0), PCG::Math<double>::Pulse3(b, 0.0, threshold, 1.0)));
			double wMT(PCG::Math<double>::LogicProbAnd(PCG::Math<double>::Pulse3(a, 0.0, threshold, 1.0), PCG::Math<double>::Ramp(b, threshold, 1.0)));
			double wFM(PCG::Math<double>::LogicProbAnd(PCG::Math<double>::Ramp(a, threshold, 0.0), PCG::Math<double>::Pulse3(b, 0.0, threshold, 1.0)));
			double wMF(PCG::Math<double>::LogicProbAnd(PCG::Math<double>::Pulse3(a, 0.0, threshold, 1.0), PCG::Math<double>::Ramp(b, threshold, 0.0)));
			double wMM(PCG::Math<double>::LogicProbAnd(PCG::Math<double>::Pulse3(a, 0.0, threshold, 1.0), PCG::Math<double>::Pulse3(b, 0.0, threshold, 1.0)));

			return	(TT * wTT + TF * wTF + FT * wFT + FF * wFF + TM * wTM + MT * wMT + FM * wFM + MF * wMF + MM * wMM) /
				(wTT + wTF + wFT + wFF + wTM + wMT + wFM + wMF + wMM);
		}

		static double LogicProbSwitch(double if_a, double then_b, double else_c)	// #LogicProbSwitch
		{
			// Switch = (!if_a & else_c) || (if_a & then_b)
			return	PCG::Math<double>::LogicProbOr
			(
				PCG::Math<double>::LogicProbAnd(PCG::Math<double>::LogicNot(if_a), else_c),
				PCG::Math<double>::LogicProbAnd(if_a, then_b)
			);

		}

		static double LogicProbFuzzySwitch(double if_a, double then_b, double else_c)	// #LogicProbFuzzySwitch
		{
			return	PCG::Math<double>::Lerp(if_a, else_c, then_b);
		}

		/*************************************************************/

		static double LogicMinMaxOr(double a, double b)	// #LogicMinMaxOr, #LogicStandardOr, #LogicWeakOr
		{
			return	PCG::Math<double>::Max(a, b);
		}

		static double LogicMinMaxAnd(double a, double b)	// #LogicMinMaxAnd, #LogicStandardAnd
		{
			return	PCG::Math<double>::Min(a, b);
		}

		static double LogicMinMaxXOr(double a, double b)	// #LogicMinMaxXOr
		{
			return	LogicMinMaxOr
			(
				LogicMinMaxAnd
				(
					LogicNot(a),
					b
				),
				LogicMinMaxAnd
				(
					a,
					LogicNot(b)
				)
			);
		}

		/*************************************************************/

		static double LogicStrongOr(double a, double b)	// #LogicStrongOr, #LogicBoundedOr
		{
			return	PCG::Math<double>::Min(1.0, a + b);
		}

		static double LogicStrongAnd(double a, double b)	// #LogicStrongAnd, #LogicBoundedAnd
		{
			return	PCG::Math<double>::Max(0.0, a + b - 1.0);
		}

		/*************************************************************/

		static double LogicEinsteinOr(double a, double b)	// #LogicEinsteinOr
		{
			return	(a + b) / (1.0 + a * b);
		}

		static double LogicEinsteinAnd(double a, double b)	// #LogicEinsteinAnd
		{
			return	(a || b) ? (a * b) / (1.0 - (1.0 - a) * (1.0 - b)) : 0.0;
		}

		/*************************************************************/

		static double LogicHamacherOr(double a, double b)	// #LogicHamacherOr
		{
			return	(1.0 - a * b) ? (a + b - 2.0 * a * b) / (1.0 - a * b) : 1.0;
		}

		static double LogicHamacherAnd(double a, double b)	// #LogicHamacherAnd
		{
			return	(a || b) ? (a * b) / (a + b - a * b) : 0.0;
		}

		/*************************************************************/

		static double LogicDrasticOr(double a, double b)	// #LogicDrasticOr
		{
			return	PCG::Math<double>::Min
			(
				PCG::Math<double>::Max(PCG::Math<double>::LogicNot(PCG::Math<double>::Step(0.0, a)), b),
				PCG::Math<double>::Max(PCG::Math<double>::LogicNot(PCG::Math<double>::Step(0.0, b)), a)
			);
		}

		static double LogicDrasticAnd(double a, double b)	// #LogicDrasticAnd
		{
			return	PCG::Math<double>::Max
			(
				PCG::Math<double>::Min(PCG::Math<double>::Step(a, 1.0), b),
				PCG::Math<double>::Min(PCG::Math<double>::Step(b, 1.0), a)
			);
		}

		// ********  ********

		static double DiffAbs(double a, double b)	// #DiffAbs
		{
			return	::fabs(a - b);
		}

		static double RelativeDifference(double a, double b)	// #RelativeDifference
		{
			return	::fabs((b - a) / Max(a, b));
		}

		static double RelativeDifferenceAbs(double a, double b)	// #RelativeDifference
		{
			//	return	::fabs ((b - a) / MaxAbs (a, b));
			return	::fabs((b - a) / MaxAbs(a, b));
		}

		// ********  ********

		static void IntervalUnite(double a_low, double a_high, double b_low, double b_high, double& result_low, double& result_high)
		{
			result_low	= PCG::Math<double>::Min (a_low, b_low);
			result_high	= PCG::Math<double>::Max (a_high, b_high);
		}

		static void IntervalIntersect(double a_low, double a_high, double b_low, double b_high, double& result_low, double& result_high)
		{
			result_low = PCG::Math<double>::Max(a_low, b_low);
			result_high = PCG::Math<double>::Min(a_high, b_high);
		}

		static void IntervalPlus(double a_low, double a_high, double b_low, double b_high, double& result_low, double& result_high) // #IntervalPlus
		{
			//double r_low (a_low + b_low);
			//double r_high (a_high + b_high);
			//result_low	= PCG::Math<double>::Min (r_low, r_high);
			//result_high	= PCG::Math<double>::Max (r_low, r_high);

			result_low = a_low + b_low;
			result_high = a_high + b_high;
		}

		static void IntervalPlusSolve(double a_low, double a_high, double b_low, double b_high, double& result_low, double& result_high) // #IntervalPlusSolve
		{
			//double r_low (a_low + b_high);
			//double r_high (a_high + b_low);
			//result_low	= PCG::Math<double>::Min (r_low, r_high);
			//result_high	= PCG::Math<double>::Max (r_low, r_high);

			result_low = a_low + b_high;
			result_high = a_high + b_low;
		}

		static void IntervalMinus(double a_low, double a_high, double b_low, double b_high, double& result_low, double& result_high) // #IntervalMinus
		{
			//double r_low (a_low - b_high);
			//double r_high (a_high - b_low);
			//result_low	= PCG::Math<double>::Min (r_low, r_high);
			//result_high	= PCG::Math<double>::Max (r_low, r_high);

			result_low = a_low - b_high;
			result_high = a_high - b_low;
		}

		static void IntervalMinusSolve(double a_low, double a_high, double b_low, double b_high, double& result_low, double& result_high) // #IntervalMinusSolve
		{
			//double r_low (a_low - b_low);
			//double r_high (a_high - b_high);
			//result_low	= PCG::Math<double>::Min (r_low, r_high);
			//result_high	= PCG::Math<double>::Max (r_low, r_high);

			result_low = a_low - b_low;
			result_high = a_high - b_high;
		}

		static void IntervalInverse(double a_low, double a_high, double& result_low, double& result_high) // #IntervalInverse
		{
			// ALGORITHM 1:
			double r1(1.0 / a_low);
			double r2(1.0 / a_high);

			if (0.0 < a_low || 0.0 > a_high)
			{
				result_low = PCG::Math<double>::Min(r1, r2);
				result_high = PCG::Math<double>::Max(r1, r2);
			}
			else
			{
				if (a_low == 0.0 && a_high > 0.0)
				{
					result_low = r2;
					result_high = PCG::Math<double>::PCG_PINF;
				}
				else if (a_low < 0.0 && a_high == 0.0)
				{
					result_low = PCG::Math<double>::PCG_NINF;
					result_high = r1;
				}
				else
					result_low = PCG::Math<double>::PCG_NINF;
				result_high = PCG::Math<double>::PCG_PINF;
			}

		}

		static void IntervalMultiply(double a_low, double a_high, double b_low, double b_high, double& result_low, double& result_high) // #IntervalMultiply
		{
			double r1(a_low * b_low);
			double r2(a_low * b_high);
			double r3(a_high * b_low);
			double r4(a_high * b_high);
			result_low = PCG::Math<double>::Min(r1, r2, r3, r4);
			result_high = PCG::Math<double>::Max(r1, r2, r3, r4);
		}

		static void IntervalMultiplySolve(double a_low, double a_high, double b_low, double b_high, double& result_low, double& result_high) // #IntervalMultiplySolve
		{
			double r1(a_low * b_low);
			double r2(a_low * b_high);
			double r3(a_high * b_low);
			double r4(a_high * b_high);
			result_low = PCG::Math<double>::MidMin(r1, r2, r3, r4);
			result_high = PCG::Math<double>::MidMax(r1, r2, r3, r4);
		}

		static void IntervalDivide(double a_low, double a_high, double b_low, double b_high, double& result_low, double& result_high) // #IntervalDivide
		{
			double b_low_inv;
			double b_high_inv;
			IntervalInverse(b_low, b_high, b_low_inv, b_high_inv);
			IntervalMultiply(a_low, a_high, b_low_inv, b_high_inv, result_low, result_high);
		}

		static void IntervalDivideSolve(double a_low, double a_high, double b_low, double b_high, double& result_low, double& result_high) // #IntervalDivideSolve
		{
			double b_low_inv;
			double b_high_inv;
			IntervalInverse(b_low, b_high, b_low_inv, b_high_inv);
			IntervalMultiplySolve(a_low, a_high, b_low_inv, b_high_inv, result_low, result_high);
		}

		static void IntervalAbs(double a_low, double a_high, double& result_low, double& result_high) // #IntervalAbs
		{
			double r1 = (PCG::Math<double>::Abs(a_low));
			double r2 = (PCG::Math<double>::Abs(a_high));

			if (0.0 < a_low || 0.0 > a_high)
			{
				result_low = PCG::Math<double>::Min(r1, r2);
				result_high = PCG::Math<double>::Max(r1, r2);
			}
			else
			{
				result_low = 0.0;
				result_high = PCG::Math<double>::Max(r1, r2);
			}

			//double r1 (PCG::Math<double>::Abs (a_low));
			//double r2 (PCG::Math<double>::Abs (a_high));

			//if (PCG::Math<double>::Pulse2 (0.0, a_low, a_high))
			//{
			//	result_low	= 0.0;
			//	result_high = PCG::Math<double>::Max (r1, r2);
			//}
			//else
			//{		
			//	result_low	= PCG::Math<double>::Min (r1, r2);
			//	result_high	= PCG::Math<double>::Max (r1, r2);
			//}
		}

		static void IntervalPow(double a_low, double a_high, double b_low, double b_high, double& result_low, double& result_high) // #IntervalPow
		{
			double r1(::_copysign(::pow(::fabs(a_low), b_low), a_low));
			double r2(::_copysign(::pow(::fabs(a_low), b_high), a_low));
			double r3(::_copysign(::pow(::fabs(a_high), b_low), a_high));
			double r4(::_copysign(::pow(::fabs(a_high), b_high), a_high));

			if (0.0 <= a_low || 0.0 >= a_high)
			{
				result_low = PCG::Math<double>::Min(r1, r2, r3, r4);
				result_high = PCG::Math<double>::Max(r1, r2, r3, r4);
			}
			else // if (PCG.Math.Pulse2 (0.0, a_low, a_high) != 0)
			{
				double r5(::pow(0.0, b_low));
				double r6(::pow(0.0, b_high));

				result_low = PCG::Math<double>::Min(PCG::Math<double>::Min(PCG::Math<double>::Min(r1, r2, r3, r4), r5), r6);
				result_high = PCG::Math<double>::Max(PCG::Math<double>::Max(PCG::Math<double>::Max(r1, r2, r3, r4), r5), r6);
			}

			//double r1 (::pow (a_low, b_low));
			//double r2 (::pow (a_low, b_high));
			//double r3 (::pow (a_high, b_low));
			//double r4 (::pow (a_high, b_high));

			//if (PCG::Math<double>::Pulse2 (0.0, a_low, a_high))
			//{
			//	double r5 (::pow (0.0, b_low));
			//	double r6 (::pow (0.0, b_high));

			//	result_low	= PCG::Math<double>::Min (PCG::Math<double>::Min (PCG::Math<double>::Min (r1, r2,  r3, r4), r5), r6);
			//	result_high = PCG::Math<double>::Max (PCG::Math<double>::Max (PCG::Math<double>::Max (r1, r2,  r3, r4), r5), r6);
			//}
			//else
			//{		
			//	result_low	= PCG::Math<double>::Min (r1, r2, r3, r4);
			//	result_high	= PCG::Math<double>::Max (r1, r2, r3, r4);
			//}
		}

		static void IntervalPowSolve(double a_low, double a_high, double b_low, double b_high, double& result_low, double& result_high) // #IntervalPow
		{
			double r1(::_copysign(::pow(::fabs(a_low), b_low), a_low));
			double r2(::_copysign(::pow(::fabs(a_low), b_high), a_low));
			double r3(::_copysign(::pow(::fabs(a_high), b_low), a_high));
			double r4(::_copysign(::pow(::fabs(a_high), b_high), a_high));

			if (PCG::Math<double>::Pulse2(0.0, a_low, a_high))
			{
				double r5(::pow(0.0, b_low));
				double r6(::pow(0.0, b_high));

				result_low = PCG::Math<double>::Min(PCG::Math<double>::Min(PCG::Math<double>::MidMin(r1, r2, r3, r4), r5), r6);
				result_high = PCG::Math<double>::Max(PCG::Math<double>::Max(PCG::Math<double>::MidMax(r1, r2, r3, r4), r5), r6);
			}
			else
			{
				result_low = PCG::Math<double>::MidMin(r1, r2, r3, r4);
				result_high = PCG::Math<double>::MidMax(r1, r2, r3, r4);
			}
		}

		static void IntervalMod(double x_low, double x_high, double& result_low, double& result_high) // #IntervalMod
		{
			//double intPart;
			//return	::modf (x, &intPart);
			double low = PCG::Math<double>::Ceiling(x_low);
			double high = PCG::Math<double>::Floor(x_high);
			double a(PCG::Math<double>::Mod(x_low));
			double b(PCG::Math<double>::Mod(PCG::Math<double>::Clamp(PCG::Math<double>::Decrement(low), x_low, x_high)));
			double c(PCG::Math<double>::Mod(PCG::Math<double>::Min(PCG::Math<double>::Increment(low), x_high)));
			double d(PCG::Math<double>::Mod(PCG::Math<double>::Max(PCG::Math<double>::Decrement(high), x_low)));
			double e(PCG::Math<double>::Mod(PCG::Math<double>::Clamp(PCG::Math<double>::Increment(high), x_low, x_high)));
			double f(PCG::Math<double>::Mod(x_high));
			result_low = PCG::Math<double>::Min(PCG::Math<double>::Min(PCG::Math<double>::Min(PCG::Math<double>::Min(PCG::Math<double>::Min(a, b), c), d), e), f);
			result_high = PCG::Math<double>::Max(PCG::Math<double>::Max(PCG::Math<double>::Max(PCG::Math<double>::Max(PCG::Math<double>::Max(a, b), c), d), e), f);
		}

		static void IntervalWrap(double x_low, double x_high, double& result_low, double& result_high) // #IntervalWrap
		{
			// ALGORITHM 1:
			//double r1 (PCG::Math<double>::Wrap (x_low));
			//double r2 (r1 + (x_high - x_low));

			//result_low	= PCG::Math<double>::Min (r1, r2);
			//result_high	= PCG::Math<double>::Max (r1, r2);

			// ALGORITHM 2: safe
			double a(PCG::Math<double>::Wrap(x_low));
			double b(PCG::Math<double>::Wrap(PCG::Math<double>::Clamp(PCG::Math<double>::Decrement(PCG::Math<double>::Ceiling(x_low)), x_low, x_high)));
			double c(PCG::Math<double>::Wrap(PCG::Math<double>::Max(PCG::Math<double>::Floor(x_high), x_low)));
			double d(PCG::Math<double>::Wrap(x_high));
			result_low = PCG::Math<double>::Min(PCG::Math<double>::Min(PCG::Math<double>::Min(a, b), c), d);
			result_high = PCG::Math<double>::Max(PCG::Math<double>::Max(PCG::Math<double>::Max(a, b), c), d);
		}

		static void IntervalReflect(double x_low, double x_high, double& result_low, double& result_high) // #IntervalWrap
		{
			double a = (PCG::Math<double>::Reflect(x_low));
			double b = (PCG::Math<double>::Reflect(PCG::Math<double>::Min(PCG::Math<double>::Ceiling(x_low), x_high)));
			double c = (PCG::Math<double>::Reflect(PCG::Math<double>::Max(PCG::Math<double>::Floor(x_high), x_low)));
			double d = (PCG::Math<double>::Reflect(x_high));
			result_low = PCG::Math<double>::Min(PCG::Math<double>::Min(PCG::Math<double>::Min(a, b), c), d);
			result_high = PCG::Math<double>::Max(PCG::Math<double>::Max(PCG::Math<double>::Max(a, b), c), d);
		}

		// ********  ********

		static double ForceEquilibrium(double x1, double x2, double eq_a, double eq_b, double eq_c, double eq_d)
		{
			//double	distance (::sqrt((x2-x1)*(x2-x1)));
			double	distance(::abs(x2 - x1));
			double  magnitude(1.0 - PCG::Math<double>::Pulse4(distance, eq_a, eq_b, eq_c, eq_d));
			double  direction(distance < (0.5*(eq_b + eq_c)) ? -1.0 : 1.0);
			return  direction * magnitude;
		}

		static double Equilibrium(double x, double a, double b, double c, double d) // #Equilibrium Force
		{
			// Alogrithm 1: fast
			if (x <= a)
			{
				return	1.0;
			}
			else if (x < b)
			{
				return	(x - b) / (a - b);
			}
			else if (x <= c)
			{
				return 0.0;
			}
			else if (x < d)
			{
				return	-(x - c) / (d - c);
			}
			else
			{
				return	-1.0;
			}

			// Alogrithm 2:
			//double sign (x < (b+c)*0.5 ? (d < a ?  -1.0 : 1.0) : (d < a ?  1.0 : -1.0));
			//return	::copysign (1.0 - PCG::Math<double>::Pulse4(x, a, b, c, d), sign);
		}

		static double Equilibrium2(double x1, double x2, double a, double b, double c, double d) // #Equilibrium Energy
		{
			double F1(PCG::Math<double>::Pulse4High(x1, PCG::Math<double>::PCG_NINF, PCG::Math<double>::PCG_NINF, a, b) - PCG::Math<double>::Pulse4High(x2, PCG::Math<double>::PCG_NINF, PCG::Math<double>::PCG_NINF, a, b));
			double F2(PCG::Math<double>::Pulse4Low(x2, c, d, PCG::Math<double>::PCG_PINF, PCG::Math<double>::PCG_PINF) - PCG::Math<double>::Pulse4Low(x1, c, d, PCG::Math<double>::PCG_PINF, PCG::Math<double>::PCG_PINF));

			return	F1 - F2;
		}

		static double EquilibriumLow(double x, double a, double b) // #EquilibriumLow
		{

			double a_(PCG::Math<double>::Max(0.0, a));
			double b_(PCG::Math<double>::Max(0.0, b));
			double x_(PCG::Math<double>::Clamp(x, a_, b_));
			double ya(PCG::Math<double>::Ramp(a_, b, a));
			double yx(PCG::Math<double>::Ramp(x_, b, a));
			double A0(PCG::Math<double>::Max(0.0, PCG::Math<double>::Min(x, a)));
			double A1(PCG::Math<double>::SegmentArea(a_, x_, ya, yx));

			return	A0 + A1;

		}

		static double EquilibriumHigh(double x, double c, double d) // #EquilibriumHigh
		{
			double c_(PCG::Math<double>::Max(0.0, c));
			double d_(PCG::Math<double>::Max(0.0, d));
			double x_(PCG::Math<double>::Clamp(x, c_, d_));
			double yc(PCG::Math<double>::Ramp(c_, c, d));
			double yx(PCG::Math<double>::Ramp(x_, c, d));
			double A0(PCG::Math<double>::SegmentArea(c_, x_, yc, yx));
			double A1(PCG::Math<double>::Max(0.0, x - d_));

			return	A0 + A1;
		}

		static double EquilibriumSin(double x, double a, double b, double c, double d) // #EquilibriumSin
		{
			double sign(x < (b + c)*0.5 ? (d < a ? 1.0 : -1.0) : (d < a ? -1.0 : 1.0));
			return	::copysign(1.0 - PCG::Math<double>::SlopeSmoothSin(PCG::Math<double>::Pulse4(x, a, b, c, d)), sign);

			//static double Equilibrium (Equilibrium (x, a, b, c, d));
			//return	::copysign (PCG::Math<double>::SlopeSmoothSin (::fabs (equilibrium)), equilibrium);

		}

		// ********  ********

		static double Fibonacci(double nu)	// #Fibonacci
		{
			static const double PCG_PHI((::sqrt(5.0) + 1.0) / 2.0);
			static const double INV_SQRT_5(1.0 / ::sqrt(5.0));
			static const double PCG_PI(::acos(-1.0));

			return	(::pow(PCG_PHI, nu) - ::cos(PCG_PI * nu) * ::pow(PCG_PHI, -nu)) * INV_SQRT_5;
		}


		static double Factorial(double n)	// #Factorial
		{
			//return Gamma (n+1.0);
			return ::tgamma(n + 1.0);
		}

		static double Poisson(double k, double l)	// #Poisson
		{
			//return ::pow (l, k) * ::exp (-l) / Gamma (k + 1.0);
			return ::pow(l, k) * ::exp(-l) / ::tgamma(k + 1.0);
		}


		// ********  ********

		static struct Fuzzy // for compatibility with higher order fuzzy numbers:
		{

			// membership // for compatibility with higher order fuzzy numbers:
			static inline double Membership(double v, double x)
			{
				return x == v;
			}

			static inline double Core(double v) // #
			{
				return	v;
			}

			static inline double Support(double v) // #
			{
				return	v;
			}

			static inline double Section(double v, double height) // #
			{
				return	v;
			}

			static inline double Infimum(double v) // #
			{
				return	v;
			}

			static inline double CoreInfimum(double v) // #
			{
				return	v;
			}

			static inline double CoreSupremum(double v) // #
			{
				return	v;
			}

			static inline double Supremum(double v) // #
			{
				return	v;
			}

			static inline double Mean(double v) // #
			{
				return	v;
			}

			static inline double Median(double v) // #
			{
				return	v;
			}

			static inline double Mode(double v) // #
			{
				return	v;
			}

			static inline double Area(double v) // #
			{
				return	0.0;
			}

			static inline double Width(double v, double height) // #
			{
				return	0.0;
			}

			static inline double Density(double v) // #
			{
				return	1.0;
			}

			static inline double Length(double v) // #
			{
				return	2.0;
			}

			static inline double Perimeter(double v) // #
			{
				return	2.0;
			}

			static inline double DistributeXY(double v, double fraction, double height) // #
			{
				return	v;
			}

			static inline double Fractile(double v, double fraction) // #
			{
				return	v;
			}

			static inline double SolvePlus(double x, double y) // #SolvePlus
			{
				return	x + y;
			}

			static inline double SolveMinus(double x, double y) // #SolveMinus
			{
				return	x - y;
			}

			static inline double SolveMultiply(double x, double y) // #SolveMultiply
			{
				return	x * y;
			}

			static inline double SolveDivide(double x, double y) // #SolveDivide
			{
				return	x / y;
			}

			static inline double SolvePow(double x, double pow) // #SolveDivide
			{
				return	Pow(x, pow);
			}

			static inline double SolveLog(double x, double numeric_base) // #SolveDivide
			{
				return	Log(x, numeric_base);
			}
			
			static inline double SolveLerp(double x, double y0, double y1)
			{
				return Lerp(x, y0, y1);
			}

			static inline double SolveLerpInv(double y, double y0, double y1)
			{
				return LerpInv(y, y0, y1);
			}

			static inline double SolveLerpGen(double x, double x0, double x1, double y0, double y1)
			{
				return LerpGen(x, x0, x1, y0, y1);
			}

			static inline double SolveLerpGenInv(double y, double x0, double x1, double y0, double y1)
			{
				return LerpGen(y, x0, x1, y0, y1);
			}

		} Fuzzy;

	}; // template <> class Math<double>

	//inline double PCG::Math<double>::TEST(double arg)
	//{
	//	return 0.0;
	//}

	const double PCG::Math<double>::PCG_PI(::acos(-1.0));
	const double PCG::Math<double>::PCG_PI_DIV_2(0.5 * ::acos(-1.0));
	const double PCG::Math<double>::PCG_1_DIV_PI(1.0 / ::acos(-1.0));
	const double PCG::Math<double>::PCG_1_DIV_2PI(1.0 / (2.0 * ::acos(-1.0)));
	const double PCG::Math<double>::PCG_1_DIV_4PI(1.0 / (4.0 * ::acos(-1.0)));
	const double PCG::Math<double>::PCG_4PI(4.0 * ::acos(-1.0));
	const double PCG::Math<double>::PCG_4_DIV_3PI(4.0 / 3.0 * ::acos(-1.0));
	const double PCG::Math<double>::PCG_2PI(2.0 * ::acos(-1.0));
	const double PCG::Math<double>::PCG_PHI((::sqrt(5.0) + 1.0) / 2.0);
	const double PCG::Math<double>::PCG_PHI_MAIOR((::sqrt(5.0) - 1.0) / 2.0);
	const double PCG::Math<double>::PCG_PHI_MINOR(1.0 - (::sqrt(5.0) - 1.0) / 2.0);
	const double PCG::Math<double>::PCG_E(::exp(1.0));
	const double PCG::Math<double>::PCG_SQRT_E(::sqrt(::exp(1.0)));

	//const double PCG::Math<double>::GAUSS_THRESHOLD(0.056375536236835322586065; /*(exp ( -pow (1.0 / (1.0 - exp (-0.5)), 2.0) * 0.5 ))*/ }
	//const double PCG::Math<double>::GAUSS_MAX()		{ return 2.39821599079104275109315; /*(1.0 / (1.0 - exp (-0.5)))*/ }		
	//const double PCG::Math<double>::MAP_GAUSS()		{ return 2.8757199692429309);
	//const double PCG::Math<double>::MAP_INV_GAUSS()	{ return 0.34773900473461694);

	const double PCG::Math<double>::PCG_SQRT2(::sqrt(2.0));
	const double PCG::Math<double>::PCG_SQRT5(::sqrt(5.0));
	const double PCG::Math<double>::PCG_SQRT_1_DIV_3(::sqrt(1.0 / 3.0));
	const double PCG::Math<double>::PCG_SQRT_1_DIV_2(::sqrt(0.5));
	const double PCG::Math<double>::PCG_1_DIV_SQRT5(1.0 / ::sqrt(5.0));
	const double PCG::Math<double>::PCG_LOG_1_DIV_2(::log(0.5));
	const double PCG::Math<double>::PCG_LOG10(::log(10.0));

	const double PCG::Math<double>::PCG_1_DIV_LOG_1_DIV_2(1.0 / ::log(0.5));
	const double PCG::Math<double>::PCG_1_DIV_LOG_PHI(1.0 / ::log((::sqrt(5.0) + 1.0) / 2.0));
	const double PCG::Math<double>::PCG_RAD_TO_DEG(180.0 / ::acos(-1.0));
	const double PCG::Math<double>::PCG_DEG_TO_RAD(::acos(-1.0) / 180.0);
	const double PCG::Math<double>::PCG_IND(::sqrt(-1.0));
	const double PCG::Math<double>::PCG_QNAN(-::sqrt(-1.0));
	const double PCG::Math<double>::PCG_NINF(-::std::numeric_limits<double>::infinity());
	const double PCG::Math<double>::PCG_PINF(::std::numeric_limits<double>::infinity());

	const double PCG::Math<double>::PCG_DIGITS(::std::numeric_limits<double>::digits);
	const double PCG::Math<double>::PCG_DIGITS10(DBL_DIG);
	const double PCG::Math<double>::PCG_MAX(DBL_MAX);
	const double PCG::Math<double>::PCG_ROOT2_MAX(::pow(DBL_MAX, 1.0 / 2.0));
	const double PCG::Math<double>::PCG_ROOT3_MAX(::pow(DBL_MAX, 1.0 / 3.0));
	const double PCG::Math<double>::PCG_ROOT4_MAX(::pow(DBL_MAX, 1.0 / 4.0));
	const double PCG::Math<double>::PCG_PRECISE_MAX(::exp2(PCG_DIGITS));
	const double PCG::Math<double>::PCG_PRECISE_BASE(::exp(log1p(PCG_MAX) / PCG_PRECISE_MAX));
	const double PCG::Math<double>::PCG_MIN(DBL_MIN);
	const double PCG::Math<double>::PCG_DENORM(::std::numeric_limits<double>::denorm_min());
	const double PCG::Math<double>::PCG_EPSILON(DBL_EPSILON);
	const double PCG::Math<double>::PCG_NEAR1(::_nextafter(1.0, 0.0));

	//const double PCG::Math<double>::ZERO_TOLERANCE(::pow(10.0, 2.0 - DBL_DIG));

	//const double PCG::Math<double>::REAL_LOG_TEN(::log(10.0));
	//const double PCG::Math<double>::REAL_LOG_HALF(::log(0.5));
	//const double PCG::Math<double>::REAL_LOG_HALF_INV(1.0 / ::log(0.5));
	//const double PCG::Math<double>::REAL_E(::exp(1.0));
	//const double PCG::Math<double>::REAL_PI(::acos(-1.0));
	//const double PCG::Math<double>::REAL_MIN(DBL_MIN);
	//const double PCG::Math<double>::REAL_MAX(DBL_MAX);
	//const double PCG::Math<double>::REAL_EPSILON(DBL_EPSILON);
	//const double PCG::Math<double>::REAL_NEAR_ONE(::nextafter(1.0, 0.0));
	//const double PCG::Math<double>::REAL_ZERO(0.0);
	//const double PCG::Math<double>::REAL_PINF(PCG::Math<double>::limits_real.infinity());
	//const double PCG::Math<double>::REAL_NINF(-REAL_PINF);
	//const double PCG::Math<double>::REAL_PHI((::sqrt(5.0) + 1.0) / 2.0);
	//const double PCG::Math<double>::REAL_PHI_MAIOR((::sqrt(5.0) - 1.0) / 2.0);
	//const double PCG::Math<double>::REAL_PHI_MINOR(1.0 - (::sqrt(5.0) - 1.0) / 2.0);
	//const double PCG::Math<double>::REAL_SQRT_TWO(::sqrt(2.0));
	//const double PCG::Math<double>::REAL_SQRT_TWO_INV(1.0 / REAL_SQRT_TWO);
	//const double PCG::Math<double>::REAL_SQRT_E(::sqrt(REAL_E));
	//const double PCG::Math<double>::REAL_SQRT_E_INV(1.0 / REAL_SQRT_E);



	

//#include "PCGMath.inl"
} // namespace PCG

#endif	// PCG_MATH__DBL_H