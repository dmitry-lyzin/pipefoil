#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <clocale>
#include <cassert>
#include <type_traits>

#define CE      constexpr
#define E       explicit
#define F       friend
#define OP      operator
#define ъ       const
#define This    ​ // этот тип
#define Thisъ   This ъ
#define ㄥ       cis

#define FN( f) constexpr friend auto f( const This& r) -> decltype( r.f()) { return r.f(); }

template< typename T
        , typename = decltype( std::declval<T>().print( std::declval< std::ostream&>()) )
        >
std::ostream& operator <<( std::ostream& os, const T& t)
{
        t.print( os);
        return os;
}

#pragma region // template operators

template< typename... >
using void_t = void;

template< typename T, typename U, typename = void> struct exist_oper_plus_eq:  std::false_type {};
template< typename T, typename U, typename = void> struct exist_oper_minus_eq: std::false_type {};
template< typename T, typename U, typename = void> struct exist_oper_mul_eq:   std::false_type {};
template< typename T, typename U, typename = void> struct exist_oper_div_eq:   std::false_type {};

template< typename T, typename U>
struct exist_oper_plus_eq < T, U, void_t< decltype( std::declval< T>() += std::declval< U>())> >
        : std::true_type
{};
template< typename T, typename U>
struct exist_oper_minus_eq< T, U, void_t< decltype( std::declval< T>() -= std::declval< U>())> >
        : std::true_type
{};
template< typename T, typename U>
struct exist_oper_mul_eq  < T, U, void_t< decltype( std::declval< T>() *= std::declval< U>())> >
        : std::true_type
{};
template< typename T, typename U>
struct exist_oper_div_eq  < T, U, void_t< decltype( std::declval< T>() /= std::declval< U>())> >
        : std::true_type
{};

/*
#define OP_B( o, name, U )                                                        \
template< typename T> CE std::enable_if_t<   exist_oper_##name##_eq< U, T>::value    \
, U> operator o ( U u, const T& t) { u o##= t; return u; };
*/

// U o= T -> U o T
#define OP_B( o, name, U )                      \
template< typename T>                           \
CE auto operator o ( U u, const T& t) ->        \
std::remove_reference_t< decltype( u o##= t )>  \
{ u o##= t; return u; };

// коммутативный оператор с каким-то типом
#define OP_C( o, name, U )                                      \
template< typename T>                                           \
CE std::enable_if_t<  !exist_oper_##name##_eq< T, U>::value     \
                   &&  exist_oper_##name##_eq< U, T>::value     \
, U> operator o ( const T& t, const U& u)                       \
{                                                               \
        U u1( u);                                               \
        u1 o##= t;                                              \
        return u1;                                              \
};

// некоммутативный оператор с каким-то типом
#define OP_Ȼ( o, name, U )                                      \
template< typename T>                                           \
CE std::enable_if_t<  !exist_oper_##name##_eq< T, U>::value     \
                   &&  exist_oper_##name##_eq< U, T>::value     \
, U> operator o ( const T& t, const U& u)                       \
{                                                               \
        U t1( t);                                               \
        t1 o##= u;                                              \
        return t1;                                              \
};

// оба варианта коммутативного оператора с каким-то типом
#define OPS_C( o, name, U ) OP_B( o, name, U ) OP_C( o, name, U )

// оба варианта некоммутативного оператора с каким-то типом
#define OPS_Ȼ( o, name, U ) OP_B( o, name, U ) OP_Ȼ( o, name, U )

#pragma endregion


ъ class Ф{} ф; // флаг для конструкторов из сырых данных
ъ class 𝐈{} 𝐢; // самая мнимая единица на свете (настолько мнимая, что даже память не занимает!)

struct ℂ;
struct ℂ₁;

namespace ce
{
        CE inline double abs  ( double a) { return a >= 0 ? a : -a; }

        //https://gist.github.com/alexshtf/eb5128b3e3e143187794
        double CE sqrt_Newton_Raphson( double r, double curr, double prev)
        {
                return curr == prev ? curr
                                    : sqrt_Newton_Raphson( r, 0.5 * (curr + r / curr), curr);
        }

        /*
        * Constexpr version of the square root
        * \return 
        *   - For a finite and non-negative value of "r", returns an approximation for the square root of "r"; 
        *   - Otherwise, returns NaN
        */
        double CE sqrt( double r)
        {
                return r >= 0 && r < std::numeric_limits< double>::infinity()
                        ? sqrt_Newton_Raphson( r, r, 0)
                        : std::numeric_limits< double>::quiet_NaN();
        }
}

CE ъ double π = 3.14159265358979323846;
#pragma warning( disable: 4455)
CE ъ double operator ""π  ( unsigned long long a) { return π * a;                }
CE ъ double operator ""⅟²(        long double a) { return ce::sqrt( double(a)); }
CE ъ double operator ""⅟²( unsigned long long a) { return ce::sqrt( double(a)); }

CE inline double  ²( double a)
{
        return a * a;
}
CE inline void swap( double &a, double &b)
{
        double a0 = a;
        a = b;
        b = a0;
}
CE inline bool eq  ( double a, double b)
{
        if( (a -= b) < 0 )
                a = -a;
        return a < 1e-14; // epsilon ≠ std::numeric_limits< double>::epsilon();
}

#pragma region // Комплексные Числа ℂ{}, ℂ₁{}

// Комплексное Число, но есть нюанс...
struct ℂ
{
protected:
        double r, i;
public:
        using This = ℂ;

CE      This( 𝐈                 ): r( 0 ), i( 1 ) {}
CE E    This( double s          ): r( s ), i( 0 ) {}
CE      This( double r, double i): r( r ), i( i ) {}

CE      double  abs² () ъ { return  r*r + i*i   ; } // абсолютная величина в квадрате
CE      double  abs  () ъ { return ce::sqrt(abs²());} // абсолютная величина
CE      This    OP - () ъ { return { -r, -i }   ; } // унарный минус
CE      This    OP ~ () ъ { return {  r, -i }   ; } // Сопряжённое (conjugate) число
CE      This    conj () ъ { return {  r, -i }   ; } // Сопряжённое (conjugate) число
CE      This    recip() ъ { This z( r, -i); z /= abs²(); return z; } // 1/z - обратная величина (reciprocal, multiplicative inverse)
CE      double  re   () ъ { return r            ; } // действительная часть
CE      double  im   () ъ { return i            ; } // мнимая часть
CE      double  ℜ   () ъ { return r            ; } // действительная часть
CE      double  ℑ    () ъ { return i            ; } // мнимая часть
CE      double  arg  () ъ
        {
                if( r > +0. ) return atan( i / r );
                if( r < -0. ) return atan( i / r ) + π;

                if( i > +0. ) return π/2;
                if( i < -0. ) return 3π/2;

                return NAN;
        }
        FN( abs²); FN( abs); FN( conj); FN( recip); FN( re); FN( im); FN( ℜ); FN( ℑ); FN( arg);

CE      bool    OP ==( Thisъ& z) ъ { return  eq( this->r, z.r ) &&  eq( this->i, z.i ); }
CE      bool    OP !=( Thisъ& z) ъ { return !eq( this->r, z.r ) || !eq( this->i, z.i ); }

CE      This&   OP +=( Thisъ& z) { r+=z.r; i+=z.i; return *this; }
CE      This&   OP -=( Thisъ& z) { r-=z.r; i-=z.i; return *this; }
CE      This&   OP *=( Thisъ& z) // Умножение векторов как комплексных чисел
        {
                This rot90 = {-i, r}; // поворот на 90° *this;
                rot90 *= z.i;
                *this *= z.r;
                *this += rot90;
                return *this;
        }
CE      This&   OP /=( Thisъ& z) { return (*this *= z.recip());  } // Деление КЧ на другое КЧ
CE      This&   OP /=( ℂ₁ ъ& n);                                  // Деление КЧ на единичное КЧ дает КЧ

CE      This&   OP +=( double s) { r += s;         return *this; }
CE      This&   OP -=( double s) { r -= s;         return *this; }
CE      This&   OP *=( double s) { r *= s; i *= s; return *this; } // Умножение на скаляр
CE      This&   OP /=( double s) { r /= s; i /= s; return *this; } // Деление на скаляр

CE      This&   OP *=( 𝐈 ъ& ) { swap( r, i); r = -r; return *this; } // Умножение на мнимую
CE      This&   OP /=( 𝐈 ъ& ) { swap( r, i); i = -i; return *this; } // Деление на мнимую

// на время представим, что наше комплексное число это вектор...
CE      double  OP , ( Thisъ& v) ъ { return  r*v.r + i*v.i; } // Скалярное произведение (Inner product)
CE      double  OP ^ ( Thisъ& v) ъ { return  r*v.i - i*v.r; } // Псевдоскалярное, Векторное, косое произведение (Outer, cross product)

        void print( std::ostream& os) ъ
        {
                static This last = { NAN, NAN };
                if( last != *this )
                {
                        os << std::setw(8) << this->r << std::setw(12) << this->i << '\n';
                        last = *this;
                }
        };
};

OPS_C( +, plus,  ℂ )
OPS_Ȼ( -, minus, ℂ )
OPS_C( *, mul,   ℂ )
OPS_Ȼ( /, div,   ℂ )


// Единичное Комплексное Число (Unit Complex Number). ℂ₁ = {𝑧 ∈ ℂ: |𝑧| = 1}
struct ℂ₁: public ℂ
{
        using This  = ℂ₁;
        using Super = ℂ;

#ifndef NDEBUG
CE      This( Ф, double r, double i): Super( r, i) {}
#endif
protected:
CE      This( Ф, Super ъ& z     ): Super( z   ) {}
CE      This( double r, double i): Super( r, i) {}

public:
CE      This( 𝐈          ): Super( 0., 1.                       ) {}
CE E    This( Super ъ& z ): Super( z / z.abs()                  ) {}
CE E    This( double sin ): Super( ce::sqrt( 1. - sin*sin), sin ) {}        
        // к-тор берет КЧ и делает из него ЕКЧ, попутно деля *p на длину этого КЧ.
CE      This( Super ъ& z, double *p )
        : Super( z)
        {
                assert( p);
                double len = z.abs();
                *(Super *)(this) /= len;
                *p               /= len; 
        }

CE      double  abs² () ъ { return     1.   ; } // абсолютная величина в квадрате = 1
CE      double  abs  () ъ { return     1.   ; } // абсолютная величина = 1, этож единичное число ;)
CE      This    OP - () ъ { return {-r, -i }; } // унарный минус (ЕКЧ = -ЕКЧ)
CE      This    OP ~ () ъ { return { r, -i }; } // Сопряжённое (conjugate) число
CE      This    conj () ъ { return { r, -i }; } // Сопряжённое (conjugate) число
CE      This    recip() ъ { return { r, -i }; } // 1/z - обратная величина (reciprocal, multiplicative inverse)
CE      double  ψarg () ъ { double x1 = r + 1.; return i>=0 ? -x1 : x1;} // псевдоугол (для сравнений)
        FN( abs²); FN( abs); FN( conj); FN( recip); FN( ψarg);

CE F    This    cis  ( double 𝜑 );// { return { cos(𝜑), sin(𝜑)}; } // фазор угла 𝜑 (единичное КЧ с аргументом 𝜑)

CE      bool    OP > ( Thisъ& z) ъ { return ψarg() > z.ψarg();}
CE      bool    OP < ( Thisъ& z) ъ { return ψarg() < z.ψarg();}

CE      This&   OP *=( Thisъ& z) { Super::OP *=( z); return *this; } // ЕКЧ*ЕКЧ = ЕКЧ
CE      This&   OP /=( Thisъ& z) { return (*this *= z.conj());     } // ЕКЧ/ЕКЧ = ЕКЧ, при этом само деление вырождено
CE      ℂ&     OP /=( ℂ ъ& z ) { return ℂ::OP /=( z); }            // КЧ/ЕКЧ = КЧ

CE      This&   OP *=( 𝐈 ъ& ) { swap( r, i); r = -r; return *this; } // Умножение на мнимую
CE      This&   OP /=( 𝐈 ъ& ) { swap( r, i); i = -i; return *this; } // Деление на мнимую

static  Thisъ 𝟏;
static  Thisъ 𝐢;
};

OPS_C( *, mul,   ℂ₁)
OPS_Ȼ( /, div,   ℂ₁)


// фазор угла 𝜑 (единичное КЧ с аргументом 𝜑)
CE ℂ₁ cis( double 𝜑 )
{       return { cos(𝜑), sin(𝜑) }; }

CE ℂ& ℂ::OP /=( ℂ₁ ъ& n) 
{       return (*this *= n.recip()); }

// комплексная единица (единичная, естественно)
CE ъ ℂ₁ ℂ₁::𝟏 = { 1., 0. };
// комплексная мнимая единица
CE ъ ℂ₁ ℂ₁::𝐢 = { 0., 1. };

// единица
CE ъ ℂ₁ 𝟏 = ℂ₁::𝟏;

CE ℂ OP ""𝐢 ( unsigned long long a) { return { 0., double(a)}; }
CE ℂ OP ""𝐢 (        long double a) { return { 0., double(a)}; }

void   ℂ_test()
{
        static_assert( (5 + 𝐢)*(7 - 6𝐢) / (3 + 𝐢)          == (10 - 11𝐢), "");
        static_assert( (4 + 𝐢)*(5 + 3𝐢) + (3 + 𝐢)*(3 - 2𝐢) == (28 + 14𝐢), "");

        // это просто должно скомпилиться
        CE auto a0 = ℂ₁( 0.1 ) * ℂ( 0, 1 );
        CE ℂ₁ u0 = 𝐢;        
        CE ℂ₁ u1 = 𝟏 / 𝐢;
        CE ℂ  u2 = ℂ₁( 0.1 ) *= 𝐢;
        CE ℂ₁ u3 = ℂ₁( 0.1 );
        CE ℂ₁ u4 = u3 * ℂ₁( 0.1 ) / 𝐢;
        CE ℂ₁ u5 = 𝟏/ℂ₁( 0.1 );
        CE ℂ₁ u6 = 𝟏*ℂ₁( 0.1 );
        CE ℂ  u7 = ℂ( 0, 1 ) * ℂ( 0.995, 0.1 ) ;
        CE ℂ₁ u8 = 𝐢 * ℂ₁( 0.1 );
        CE ℂ  u9 = u7 / u8;
        CE ℂ  ua = u8 / u7;
        CE ℂ  ub = u7 + u9;
        CE ℂ  uc = u7 * 9;

        CE ℂ₁ a = ℂ₁( 0.1);
        static_assert( a / ℂ(3, -2) == ℂ( 0.2142278701015276898, 0.1761519134010184617), "");
        static_assert( a      > 𝟏 && !( a      < 𝟏), "");
        static_assert(~a      > 𝟏 && !(~a      < 𝟏), "");
        static_assert(~a      > a && !(~a      < a), "");
        static_assert( a*𝐢*𝐢*𝐢 > a && !( a*𝐢*𝐢*𝐢 < a), "");
}
#pragma endregion

struct Line
{
        double p; // растояние от начало координат до прямой
        ℂ₁    n̂; // единичный (|n̂| = 1) перпедикуляр к прямой

        // Прямая, заданой нормальным уравнением (n̅, r̅) = 𝐶 т.е. 𝐴𝑥ᵣ + 𝐵𝑦ᵣ = 𝐶, n̅ = {𝐴, 𝐵}, |n̅| > 0
        // \param[in] n̅ - вектор, нормальный к прямой
        // \param[in] 𝐶 - скалярный параметр
CE      Line( ℂ  ъ& n̅, double 𝐶 ): p( 𝐶 ), n̂( n̅, &p ) {}
        // Прямая, задана нормированным уравнением (n̂, r̅) = 𝑝, |n̂| = 1, 𝑝 ⩾ 0
        // \param[in] n̂ - единичный вектор, нормальный к прямой
        // \param[in] 𝑝 - растояние от начало координат до прямой
CE      Line( ℂ₁ ъ& n̂, double 𝑝 ): p( 𝑝 ), n̂( n̂     ) {};
        // Прямая через точки a̅ и b̅
CE      Line( ℂ ъ& a̅, ℂ ъ& b̅   ): Line( (a̅ - b̅)*𝐢, (a̅ ^ b̅) ) {};

CE      bool    OP ==( Line ъ& l) ъ { return eq( p, l.p) && n̂ == l.n̂; }
CE      double  dist ( ℂ   ъ& r̅) ъ { return (n̂, r̅) - p             ; } // растояние между прямой и точкой r̅
CE      ℂ      proj ( ℂ   ъ& r̅) ъ { return r̅ - n̂ * dist( r̅)       ; } // проекция точки r̅ на прямую

        // растояние между точкой r̅ и прямой l
CE F    double dist( ℂ ъ& r̅, Line ъ& l) { return l.dist( r̅); }
        // растояние между прямой l и точкой r̅
CE F    double dist( Line ъ& l, ℂ ъ& r̅) { return l.dist( r̅); }
        // проекция точки r̅ на прямую l
CE F    ℂ OP >>( ℂ ъ& r̅, Line ъ& l ) { return l.proj( r̅); }

/*
friend  std::ostream& OP<<( std::ostream &os, Line ъ& obj )
        {
                os << ℂ( 0.0, obj.p/obj.n.i)
                   << ℂ( obj.p/obj.n.r, 0.0);
                return os;
        };
*/
};
struct Vertical: public Line
{
        // Вертикаль пересекающая ось 𝑋 в точке x0
CE      Vertical( double x0 ): Line( ℂ₁::𝟏, x0) {};
};
struct Horizontal: public Line
{
        // Горизонталь пересекающая ось 𝑌 в точке y0
CE      Horizontal( double y0 ): Line( ℂ₁::𝐢, y0) {};
};
void   Line_test()
{
        CE Vertical l_test( 2 );
        static_assert(       l_test.n̂        == 𝟏, "");
        static_assert(       l_test.n̂ * 𝐢    == 𝐢, "");
        static_assert( dist( l_test, {3, 1}) == 1, "");
}

struct Segment: public Line
{
        ℂ p̅1, p̅2;

CE      Segment( ℂ ъ& p1_, ℂ ъ& p2_ )
        : Line( p1_, p2_ )
        , p̅1( p1_), p̅2( p2_)
        {};

        // TODO стремный конст-ор, как бы его спрятать
CE      Segment( Line ъ& l, ℂ ъ& p1_, ℂ ъ& p2_ )
        : Line( l )
        , p̅1( p1_), p̅2( p2_)
        {};

        // обмен концов (рисоваться будет в другую сторону)
CE      Segment OP - () ъ { return Segment( this->p̅2, this->p̅1 ); }

        void print( std::ostream& os) ъ
        {
                os << p̅1 << p̅2;
        };
};

struct Circle
{
        ℂ     o̅; // центр окружности
        double R; // радиус окружности

CE      Circle( ℂ ъ& center, double radius = 0. )
        : R( radius), o̅( center )
        {};

CE      bool   OP ==( Circle ъ& c ) ъ { return eq( R, c.R) && o̅ == c.o̅; }
CE      Circle OP - (             ) ъ
        {
                Circle a = *this;
                a.R = -a.R;
                return a;
        }

CE      ℂ intersect( Line   ъ& l ) ъ
        {
                double h = dist( o̅, l);          // расстояние от центра окружности до прямой
                //ℂ    p̅ = (o̅ - l.n̂ * h);       // точка проекции центра окружности на прямую
                //ℂ₁   v̂ = l.n̂ * 𝐢;             // направляющий вектор прямой l
                //return p̅ + v̂ * ce::sqrt(²(R) - ²(h)) * sign;
                return o̅ - l.n̂ * ℂ( h, ce::sqrt( ²(R) - ²(h)) );
        }
CE      ℂ intersect( Circle ъ& c2) ъ
        {
                Line radical_line( (c2.o̅ - o̅) * 2. 
                                 , abs²(c2.o̅) - abs²(o̅) - ²(c2.R) + ²(R)
                                 );
                return intersect( radical_line );
        }

CE F    Circle tangent( Circle ъ& c1, Circle ъ& c2, double R )
        {
                ℂ center = Circle( c1.o̅, R+c1.R).intersect( Circle( c2.o̅, R+c2.R) );
                return Circle( center, R );
        }

        // нормаль к касательной к окружностям *this и c2
CE      ℂ₁ tangent_norm( Circle ъ& c2 ) ъ
        {
                double sin𝜑 = c2.R - R;         // 𝜑 - угол между линией центров и касательной
                ℂ₁ â( c2.o̅ - o̅, &sin𝜑 );        // направляющий вектор линии центров
                return â * ℂ₁( sin𝜑) / 𝐢;       // повернуть â на 𝜑-90°
        }

        // касательная к окружностям *this и c2
CE      Line tangent( Circle ъ& c2 ) ъ
        {
                ℂ₁ n̂ = tangent_norm( c2 );
                return Line( n̂, (n̂, o̅) - R );
        }

        // точка касания касательной к окружности исходящей из точки p
CE      ℂ tangent_point( ℂ ъ& p ) ъ
        {
                ℂ₁ n̂ = tangent_norm( {p, 0} );
                return o̅ - n̂*R;
        };

        void print( std::ostream &os) ъ
        {
                CE ъ int segs = 40;

                //CE ъ ℂ₁ m̂¹⁰ = cis( 2π/segs);
                CE ъ ℂ₁ m̂( 2π/segs/10 ); // типа 𝛼 ≈ sin 𝛼
                CE ъ ℂ₁ m̂⁵  = m̂*m̂*m̂*m̂*m̂;
                CE ъ ℂ₁ m̂¹⁰ = m̂⁵*m̂⁵;

                ℂ₁ n̂ = 𝟏;
                for( int i = segs; i --> 0; )
                {
                        os << (o̅ + R*n̂);
                        n̂ *= m̂¹⁰;
                }
        };

        void print00( std::ostream &os, double α1, double α2 ) ъ
        {
                double Δα = α2 - α1; // угол поворота
                // кол. сегментов, на круг - 160 сегментов, примерно
                int segments = static_cast< int>( round( abs( Δα / 2π * 160)));
                ℂ₁ m̂ = cis( Δα / segments);   // ед. век. повернутый на угол Δα/segments

                ℂ  r̅ = abs( R) * cis( α1);    // радиус от центра окружности
                ℂ  s̅ = m̂*r̅ - r̅;               // сегментик, которым рисуем окружность
                r̅ += o̅;
                for( ; segments >= 0; --segments )
                {
                        os << r̅;
                        r̅ += ( s̅ *= m̂ );
                }
        };
};
void   Circle_test()
{
        CE Circle c_test( {0, 0}, 3 );
        CE ℂ o1( 5, 0);
        CE ℂ o2( 5, 5);

        //CE auto b = Circle( o1, -2).tangent( Circle( o2, -2));
        //CE auto a = Line( ℂ( 0.6,  0.8), 5);

        static_assert( Circle( o1,  2).tangent( Circle( o2,  2)) == Line( ℂ( 1.0,   .0), 3), "");
        static_assert( Circle( o1, -2).tangent( Circle( o2,  2)) == Line( ℂ( 0.6,  0.8), 5), "");
        static_assert( Circle( o1,  2).tangent( Circle( o2, -2)) == Line( ℂ( 0.6, -0.8), 1), "");
        static_assert( Circle( o1, -2).tangent( Circle( o2, -2)) == Line( ℂ( 1.0,   .0), 7), "");

        static_assert( c_test.intersect( Line( {3, 0}, {0, 3})) == ℂ( 3,    0), "");
        static_assert( c_test.intersect( Vertical( 3)         ) == ℂ( 3,    0), "");
        static_assert( c_test.intersect( Circle( {4, 0}, 3)   ) == ℂ( 2,-5⅟²), "");

        static_assert( c_test.tangent_point( ℂ( 5, 0)       ) == ℂ(1.8, 2.4), "");
        //static_assert( c_test.tangent_point( ℂ( 5, 0), minus) == ℂ(1.8,-2.4), "");
        /*
        static_assert( c_test.tangent( Circle( {   2, 0}, 1)       ) == Line( ℂ₁::𝟏             ,  3), "");
        static_assert( c_test.tangent( Circle( {2⅟², 0}, 2)       ) == Line( ℂ( 0.5⅟²,  0.5⅟²), -3), "");
        static_assert( c_test.tangent( Circle( {2⅟², 0}, 2), minus) == Line( ℂ( 0.5⅟², -0.5⅟²), -3), "");
        static_assert( c_test.tangent( Circle( {   9, 0},-3)       ) == Line( ℂ(   2./3,  5⅟²/3), -3), "");
        */
        static_assert( tangent( c_test, Circle( {  6, 0}, 3),  3 ) == Circle( {  3,-3*3⅟²},  3), "");
        static_assert( tangent( c_test, Circle( {  2, 0}, 1), 10 ) == Circle( { 13,      0}, 10), "");
}

struct Arc: public Circle
{
        ℂ₁ n̂₁; // единич. вектор из центра дуги на первый конец
        ℂ₁ n̂₂; // единич. вектор из центра дуги на второй конец

CE      Arc( Circle ъ& c, ℂ₁ ъ& _n̂₁, ℂ₁ ъ& _n̂₂)
        : Circle( c), n̂₁(_n̂₁), n̂₂(_n̂₂)
        {}
CE      Arc( Circle ъ& c, ℂ ъ& start_point, ℂ ъ& next_point )
        : Circle( c)
        , n̂₁( start_point- c.o̅ )
        , n̂₂( next_point - c.o̅ )
        {}
CE      Arc( Circle ъ& c, Arc ъ& prev, ℂ ъ& next_point )
        : Arc( c, prev.o̅ + ce::abs(prev.R)*prev.n̂₂, next_point )
        {
                // ставим знак вектору n̂₂ таким, чтоб конечная точка дуги была максимально близко к next_point
                //if( ((next_point - o̅), n̂₂) < 0)
                //        n̂₂ = -n̂₂;
                if( prev.R * R * (prev.n̂₂, n̂₁) < 0 )
                        R = -R;
        }

CE      bool OP ==( Arc ъ& a) ъ
        {
                return  Circle::OP==( a)
                        && n̂₁ == a.n̂₁
                        && n̂₂ == a.n̂₂
                        ;
        }

        void print( std::ostream &os) ъ
        {
                CE ъ int segs = 160;

                //CE ъ ℂ₁ m̂¹⁰ = cis( 2π/segs);
                CE ъ ℂ₁ m̂( 2π/segs/10 ); // типа 𝛼 ≈ sin 𝛼
                CE ъ ℂ₁ m̂⁵ = m̂*m̂*m̂*m̂*m̂;

                ℂ₁ m̂¹⁰ = m̂⁵*m̂⁵;
                bool dir = (R >= 0); // направление рисования
                if( !dir )
                        m̂¹⁰ = ~m̂¹⁰;

                ℂ r̅ = abs(R) * n̂₁;
                os << (o̅ + r̅);
                double a₂₁ = ψarg( n̂₂/n̂₁);

                for( ℂ₁ n̂ = m̂¹⁰; (ψarg(n̂) < a₂₁) == dir; n̂ *= m̂¹⁰ )
                        os << (o̅ + r̅*n̂);
        };
};

// зазор для соблюдения постулата Жуковского-Чаплыгина (Kutta condition)
CE ъ ℂ TE1( 1.,  0.00001); // задняя кромка верх
CE ъ ℂ TE2( 1., -0.00001); // задняя кромка низ

#ifndef NDEBUG
/*
int test1()
{
        // профиль со скруглением сверху носка
        CE double b_abs   =  40.0; // хорда профиля
        CE double D_abs   = 110.0; // диаметр трубы
        CE double s_abs   =   2.7; // толщина стенки трубы
        CE double lr1_abs =   1.0; // радиус скругл. носка 1
        CE double lr2_abs =  18.0; // радиус скругл. носка 2

        CE double R   = D_abs/2./ b_abs; // относительный внешний радиус трубы
        CE double s   = s_abs   / b_abs; // относительная толщина стенки трубы
        CE double ½s  = s / 2.         ; // относительная полутолщина стенки трубы
        CE double lr1 = lr1_abs / b_abs; // относительный радиус скругл. носка 1
        CE double lr2 = lr2_abs / b_abs; // относительный радиус скругл. носка 2

        CE Circle  c_lead({0, ½s}, ½s );
        CE Circle  с_top = tangent( c_lead, TE1, R );
        CE Circle  c_bottom( с_top.o̅, R-s );
        CE Segment s_end = tangent_segment( c_bottom, TE2 );
        CE Circle  c_l1( Vertical( lr1) ^ Circle( с_top.o̅, c_bottom.R+lr1), lr1 );
        CE Circle  c_l2 = tangent( с_top, c_l1, -lr2 );

        CE Arc a_top = chain_arc( TE1      , с_top , c_l2     );
        CE Arc a_l2  = chain_arc( a_top.p̅2 , c_l2  , c_l1     );
        CE Arc a_l1  = chain_arc( a_l2.p̅2  , c_l1  , c_bottom, minus );
        CE Arc a_bottom( c_bottom, a_l1.p̅2, s_end.p̅1 );

        std::cout << "PIPE " << D_abs << 'r' << s_abs << '-' << b_abs << '\n'
                << std::setprecision(5) << std::fixed
                << a_top << a_l2 << -a_l1 << a_bottom << s_end;

        return 0;
}
int test2()
{
        // круглая ПК немного сверху зализана и заостренная ЗК 
        // Cy/Cx = 25 при α ≈ 5° плоская полка α ≈ 0÷5°, срыв при α ≈ 7,5°
        //CE double b_abs   =  50.0; // хорда профиля
        //CE double D_abs   = 160.0; // диаметр трубы
        //CE double s_abs   =   4.7; // толщина стенки трубы
        //CE double lr1_abs =   1.0; // радиус скругл. носка 1
        //CE double lr2_abs =  20.0; // радиус скругл. носка 2

        // просто круглая ПК и заостренная ЗК
        // Cy/Cx = 28 при α ≈ 5,2° очень узкий пик, срыв при α ≈ 5,2°
        CE double b_abs   =  40.0; // хорда профиля
        CE double D_abs   = 110.0; // диаметр трубы
        CE double s_abs   =   2.7; // толщина стенки трубы
        CE double lr1_abs =   1.3; // радиус скругл. носка 1
        CE double lr2_abs =  18.0; // радиус скругл. носка 2

        CE double R   = D_abs/2./ b_abs; // относительный внешний радиус трубы
        CE double s   = s_abs   / b_abs; // относительная толщина стенки трубы
        CE double lr1 = lr1_abs / b_abs; // относительный радиус скругл. носка 1

        CE Circle  c_lead({lr1, lr1}, lr1 );
        CE Circle  с_top = tangent( c_lead, TE1, R );
        CE Circle  c_bottom( с_top.o̅, R-s );
        CE Segment s_start = tangent_segment( c_lead, -c_bottom, minus );
        CE Segment s_end = tangent_segment( c_bottom, TE2 );

        CE Arc a_top = chain_arc( TE1, с_top, c_lead      );
        CE Arc a_lead  ( c_lead  , a_top.p̅2  , s_start.p̅1 );
        CE Arc a_bottom( c_bottom, s_start.p̅2, s_end.p̅1   );

        std::cout << "PIPE " << D_abs << 'x' << s_abs << '-' << b_abs << '\n'
                << std::setprecision(5) << std::fixed
                << a_top << -a_lead << s_start << a_bottom << s_end;

        return 0;
}
int test3()
{
        // ровная нижняя поверхность (предполагался скотч)
        // Cy/Cx = 31 при α = 6,5°

        CE double b_abs   =  50.0; // хорда профиля
        CE double D_abs   = 160.0; // диаметр трубы
        CE double s_abs   =   4.7; // толщина стенки трубы
        CE double lr1_abs =   1.5; // радиус скругл. носка 1
        CE double lr2_abs =  20.0; // радиус скругл. носка 2

        CE double R   = D_abs/2./ b_abs; // относительный внешний радиус трубы
        CE double s   = s_abs   / b_abs; // относительная толщина стенки трубы
        CE double lr1 = lr1_abs / b_abs; // относительный радиус скругл. носка 1
        CE double lr2 = lr2_abs / b_abs; // относительный радиус скругл. носка 1

        CE Circle  c_l1   ( {lr1, lr1}, lr1 );
        CE Circle  c_trail( TE1       , -s  );
        CE Circle  c_bottom = tangent( -c_l1, c_trail, R-s );
        CE Circle  с_top( c_bottom.o̅, R );
        CE Circle  c_l2  = tangent( c_l1, с_top, -lr2, minus );
        CE Segment s_end = tangent_segment( c_l1, TE2, minus );

        CE Arc a_top = chain_arc( TE1     , с_top, c_l2 );
        CE Arc a_l2  = chain_arc( a_top.p̅2, c_l2 , c_l1 );
        CE Arc a_l1( c_l1, a_l2.p̅2, s_end.p̅1 );

        std::cout << "PIPE " << D_abs << 'x' << s_abs << '-' << b_abs << '\n'
                << std::setprecision(5) << std::fixed
                << a_top << a_l2 << a_l1 << s_end;

        return 0;
}
int test4()
{
        // ровная нижняя поверхность (предполагался скотч, но он оторвался)
        // тупой пик, Cy/Cx > 30 α ≈ 5÷9°, Cy/Cx ≈ 34 α ≈ 6.6°, срыв α ≈ 15°
        CE double b_abs   =  50.0; // хорда профиля
        CE double D_abs   = 160.0; // диаметр трубы
        CE double s_abs   =   4.7; // толщина стенки трубы
        CE double lr1_abs =   1.5; // радиус скругл. носка 1
        CE double lr2_abs =  20.0; // радиус скругл. носка 2

        CE double R   = D_abs/2./ b_abs; // относительный внешний радиус трубы
        CE double s   = s_abs   / b_abs; // относительная толщина стенки трубы
        CE double lr1 = lr1_abs / b_abs; // относительный радиус скругл. носка 1
        CE double lr2 = lr2_abs / b_abs; // относительный радиус скругл. носка 1

        CE Circle  c_l1   ( {lr1, lr1}, lr1 );
        CE Circle  c_trail( TE1       , -s  );
        CE Circle  c_bottom = tangent( -c_l1, c_trail, R-s );
        CE Circle  с_top( c_bottom.o̅, R );
        CE Circle  c_l2  = tangent( c_l1, с_top, -lr2, minus );
        CE Segment s_end0 = tangent_segment( c_l1, TE2, minus );

        CE Arc a_top = chain_arc( TE1     , с_top, c_l2 );
        CE Arc a_l2  = chain_arc( a_top.p̅2, c_l2 , c_l1 );
        CE Arc a_l1( c_l1, a_l2.p̅2, s_end0.p̅1 );

        CE Segment s_start( s_end0.p̅1, s_end0 ^ c_bottom );
        CE Segment s_end( s_end0 & c_bottom, TE2 );
        CE Arc a_bottom( c_bottom, s_start.p̅2, s_end.p̅1 );

        std::cout << "PIPE " << D_abs << 'x' << s_abs << '-' << b_abs << '\n'
                << std::setprecision(5) << std::fixed
                << a_top << a_l2 << a_l1 << s_start << a_bottom << s_end;

        return 0;
}*/
int test5()
{
        // круглая ПК немного сверху зализана и заостренная ЗК
        // тупой пик, Cy/Cx > 30 α ≈ 4.5÷9.5°, Cy/Cx ≈ 35 α ≈ 6.6°, срыв α ≈ 15°
        CE double b_abs   =  50.0; // хорда профиля
        CE double D_abs   = 160.0; // диаметр трубы
        CE double s_abs   =   4.7; // толщина стенки трубы
        CE double ler_abs =   1.5; // радиус скругл. носка 1
        CE double lef_abs =  20.0; // радиус скругл. носка 2

        CE double R   = D_abs/2./ b_abs; // относительный внешний радиус трубы
        CE double s   = s_abs   / b_abs; // относительная толщина стенки трубы
        CE double ler = ler_abs / b_abs; // относительный радиус передней кромки
        CE double lef = lef_abs / b_abs; // относительный радиус скругл. передней кромки

        CE Circle c_le     = { {ler, ler}, ler };
        CE Circle c_trail  = { TE1       , s   };
        CE Circle c_bottom = tangent( c_le, c_trail, R-s );
        CE Circle с_top    = { c_bottom.o̅, R };
        CE Circle c_lef    = tangent( с_top, c_le, -lef );
        CE ℂ     p_end    = c_bottom.tangent_point( TE2 );

        CE Arc a_top   ( с_top   , TE1  , c_lef.o̅   );
        CE Arc a_lef   ( c_lef   , a_top, c_le.o̅    );
        CE Arc a_le    ( c_le    , a_lef, c_bottom.o̅);
        CE Arc a_bottom( c_bottom, a_le , p_end     );
        
        // static_asserts
        {
                static_assert(    a_top == Arc{{{0.3656664547982899638, -1.468874254606818974},  1.6  }, {ф, 0.3964584657510687449, 0.9180526591292618166}, {ф,-0.1059312120582786843,  0.9943734601807634466}}, "");
                static_assert(    a_lef == Arc{{{0.238549000328355576, -0.2756261023899031493},  0.4  }, {ф,-0.1059312120582786426, 0.9943734601807632245}, {ф,-0.5636459468333934186,  0.8260164929456842442}}, "");
                static_assert(     a_le == Arc{{{0.03                ,  0.03                 },  0.03 }, {ф,-0.5636459468333931966, 0.8260164929456842442}, {ф, 0.2185328481759700181, -0.9758295928429810973}}, "");
                static_assert( a_bottom == Arc{{{0.3656664547982899638, -1.468874254606818974}, -1.506}, {ф,-0.2185328481759700459, 0.9758295928429812083}, {ф, 0.06316730084589926297, 0.9980029519514677094}}, "");
        }
        
        std::cout << "PIPE " << D_abs << 'x' << s_abs << '-' << b_abs << " r" << ler_abs << " f" << lef_abs << '\n'
                << std::setprecision(5) << std::fixed
                << a_top << a_lef << a_le << a_bottom << p_end << TE2
                ;

        return 0;
}
#endif

int main( unsigned argc, ъ char *argv[])
{
#ifndef NDEBUG
        return test5();
#endif

        static CE char *param_name[] =
        { "диаметр трубы"
        , "толщина стенки трубы"
        , "хорда"
        , "радиус передней кромки"
        , "радиус зализа над передней кромкой"
        };

        static double param[ std::size( param_name)] = { 160.0, 4.7, 50.0, 1.5, 18.0};

        double &D_abs   = param[0];  // диаметр трубы
        double &s_abs   = param[1];  // толщина стенки трубы
        double &b_abs   = param[2];  // хорда профиля
        double &ler_abs = param[3];  // радиус передней кромки
        double &lef_abs = param[4];  // радиус скругл. передней кромки
/*
 #i fn def NDEBUG
        if( argc == 3 )
        {
                std::size_t n = 0;
                if( !strcmp( argv[1], "test") && ( n = atoi( argv[2])) )
                {
                        static CE int (*checks[])() = { &test1, &test2, &test3, &test4, &test5 };
                        if( (n > 0) && (n <= std::size( checks)) )
                                return checks[ n-1 ]();
                }
                setlocale( LC_ALL, "" );
                std::cerr << "Ошибка тестов!\n";
                return 1;
        }
 #e nd if
*/
        if( argc <= std::size( param_name) )
        {
                setlocale( LC_ALL, "" );

                std::cerr << "Примерные параметры запуска:\n  pipefoil";
                for( unsigned i = 0; i < std::size( param); ++i)
                        std::cerr << ' ' << param[ i ];

                std::cerr << " >1.dat\nгде:\n";
                for( unsigned i = 0; i < std::size( param); ++i)
                        std::cerr << "  " << param[ i ] << "\t- " << param_name[ i ] << '\n';

                std::cerr << "  1.dat\t- результат: файл с точками профиля, можно \"продуть\" его в xflr5\n";
                return 1;
        }

        for( unsigned i = 0; i < std::size( param); ++i)
        {
                char *last_char;
                param[i] = strtod( argv[i+1], &last_char);
                if( *last_char)
                {
                        setlocale( LC_ALL, "" );
                        std::cerr << "Ошибка " << (i+1) << "-ого параметра \"" << param_name[i] << "\"\n"
                                "\"" << argv[i+1] << "\" распозналось как \"" << param[i] << "\"\n";
                        return 1;
                }
        }

        double R   = D_abs/2./ b_abs; // относительный внешний радиус трубы
        double s   = s_abs   / b_abs; // относительная толщина стенки трубы
        double ler = ler_abs / b_abs; // относительный радиус передней кромки
        double lef = lef_abs / b_abs; // относительный радиус скругл. передней кромки

        Circle c_le     = { {ler, ler}, ler };
        Circle c_trail  = { TE1       , s   };
        Circle c_bottom = tangent( c_le, c_trail, R-s );
        Circle с_top    = { c_bottom.o̅, R   };
        Circle c_lef    = tangent( с_top, c_le, -lef );
        ℂ     p_end    = c_bottom.tangent_point( TE2 );

        Arc a_top   ( с_top   , TE1  , c_lef.o̅   );
        Arc a_lef   ( c_lef   , a_top, c_le.o̅    );
        Arc a_le    ( c_le    , a_lef, c_bottom.o̅);
        Arc a_bottom( c_bottom, a_le , p_end     );
        
        std::cout << "PIPE " << D_abs << 'x' << s_abs << '-' << b_abs << " r" << ler_abs << " f" << lef_abs << '\n'
                << std::setprecision(5) << std::fixed
                << a_top << a_lef << a_le << a_bottom << p_end << TE2;

        return 0;
}
