#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <clocale>
#include <cassert>

#define CE      constexpr
#define CExplct constexpr explicit
#define CEfrnd  constexpr friend

#define FUNC( f) constexpr friend auto f( const This& x) { return x.f(); }

// оператор со своим типом
#define OPER_THIS( o) constexpr friend This operator o ( This t, const This& a) { return t o##= a; }

// оператор с родительским типом
#define OPER_SUPER( o) \
constexpr friend Super operator o ( This t, const Super& a) { return t o##= a; }; \
constexpr friend Super operator o ( Super t, const This& a) { return t o##= a; }

// коммутативный оператор с каким-то простым типом
#define OPER_COMM( o, Any) \
constexpr friend This operator o ( This t, Any a) { return t o##= a; }; \
constexpr friend This operator o ( Any a, This t) { return t o##= a; }

// некоммутативный оператор с каким-то простым типом
#define OPER_NOCOMM( o, Any) \
constexpr friend This operator o ( This t,        Any a) {             return t o##= a; }; \
constexpr friend This operator o ( Any a, const This& t) { This b( a); return b o##= t; }

const class Ф{} ф; // флаг для конструкторов из сырых данных
const class Ĵ{} ĵ; // самая мнимая единица на свете (умножение на неё поворачивает вектор на 90°)
//CE const Ĵ ⅈ;

class  Vec;
class NVec;

class Sgn
{
        int _;
CE      Sgn( Ф, int a): _( a) {};
public:
CExplct Sgn( bool   a): _( a      ? +1 : -1) {}; // Sgn( false) = -1, Sgn( true) = +1
CExplct Sgn( int    a): _( a >= 0 ? +1 : -1) {};
CExplct Sgn( double a): _( a >= 0 ? +1 : -1) {};

CE      Sgn    operator *=( Sgn s          ) { _ *= s._; return *this; }
CEfrnd  Sgn    operator * ( Sgn s, Sgn    a) { return s *= a;  }
CEfrnd  double operator *=( double a, Sgn s) { return a *= s._;} // Умножение на знак
CEfrnd  double operator * ( Sgn s, double a) { return s._ * a; } // Умножение на знак
CEfrnd  double operator * ( double a, Sgn s) { return s._ * a; } // Умножение на знак
CE      Sgn    operator - () const           { return {ф, -_}; }

friend  class  Vec;
friend  class NVec;
};

namespace ce
{
CE inline double abs  ( double a) { return a >= 0 ? a : -a; }

        double CE sqrt_Newton_Raphson( double x, double curr, double prev)
        {
                return curr == prev ? curr
                                    : sqrt_Newton_Raphson( x, 0.5 * (curr + x / curr), curr);
        }

        /*
        * Constexpr version of the square root
        * \return 
        *   - For a finite and non-negative value of "x", returns an approximation for the square root of "x"; 
        *   - Otherwise, returns NaN
        */
        double CE sqrt( double x)
        {
                return x >= 0 && x < std::numeric_limits< double>::infinity()
                        ? sqrt_Newton_Raphson( x, x, 0)
                        : std::numeric_limits< double>::quiet_NaN();
        }
}

// epsilon ≠ std::numeric_limits< double>::epsilon();
CE const double eps = 1e-14;
CE const double π = 3.14159265358979323846;
#pragma warning( disable: 4455)
CE const double operator ""π  ( unsigned long long a) { return π * a;                }
CE const double operator ""⅟²(        long double a) { return ce::sqrt( double(a)); }
CE const double operator ""⅟²( unsigned long long a) { return ce::sqrt( double(a)); }

CE inline double     ²( double a)
{
        return a * a;
}
CE inline void swap   ( double &a, double &b)
{
        double a0 = a;
        a = b;
        b = a0;
}
CE inline bool eq     ( double a, double b)
{
        if( (a -= b) < 0 )
                a = -a;
        return a < eps;
}

#pragma region // Vec{}, NVec{}
// ℂ
struct Vec
{
protected:
        double x, y;
public:
        using This = Vec;

CE      Vec( Ĵ                   ): x( 0 ), y( 1 ) {}
CExplct Vec( double x_           ): x( x_), y( 0 ) {}
CE      Vec( double x_, double y_): x( x_), y( y_) {}

CE      double  operator , ( const Vec& v) const { return  x*v.x + y*v.y  ; } // Скалярное произведение (Inner product)
CE      double  operator ^ ( const Vec& v) const { return  x*v.y - y*v.x  ; } // Псевдоскалярное, Векторное, косое произведение (Outer, cross product)

CE      double  abs²       () const { return  x*x + y*y   ; } // Длина²
CE      double  abs        () const { return ce::sqrt(abs²());} // Длина
CE      Vec     operator - () const { return { -x, -y }   ; } // поворот на 180°
CE      Vec     operator ~ () const { return {  x, -y }   ; } // Сопряжённое (conjugate) число (зеркально отраженный вектор)
CE      Vec     conj       () const { return {  x, -y }   ; } // Сопряжённое (conjugate) число (зеркально отраженный вектор)
CE      Vec     inv        () const { return conj()/abs²(); } // 1/v - обратная (inverse) величина
CE      double  re         () const { return x            ; } // действительная часть
CE      double  im         () const { return y            ; } // мнимая часть
CE      double  ℜ         () const { return x            ; } // действительная часть
CE      double  ℑ          () const { return y            ; } // мнимая часть
CE      double  arg        () const
{
        if( x > +0. ) return atan( y / x );
        if( x < -0. ) return atan( y / x ) + π;

        if( y > +0. ) return π/2;
        if( y < -0. ) return 3π/2;

        return NAN;
}
        FUNC( abs²); FUNC( abs); FUNC( conj); FUNC( inv); FUNC( re); FUNC( im); FUNC( ℜ); FUNC( ℑ); FUNC( arg);

CE      bool    operator ==( const Vec& v) const { return  eq( this->x, v.x ) &&  eq( this->y, v.y ); }
CE      bool    operator !=( const Vec& v) const { return !eq( this->x, v.x ) || !eq( this->y, v.y ); }

CE      Vec&    operator +=( const Vec& v) {  x+=v.x; y+=v.y; return *this; }
CE      Vec&    operator -=( const Vec& v) {  x-=v.x; y-=v.y; return *this; }
CE      Vec&    operator *=( const Vec& v) // Умножение векторов как комплексных чисел
        {
                Vec rot90 = {-y, x}; // поворот на 90° *this;
                rot90 *= v.y;
                *this *= v.x;
                *this += rot90;
                return *this;
        }
CE      Vec&    operator /=( const Vec& v) {*this *= v.inv(); return *this; } // Деление векторов как комплексных чисел
CE      Vec&    operator /=( const NVec &n);
        OPER_THIS( +); OPER_THIS( -); OPER_THIS( *); OPER_THIS( /);

CE      Vec&    operator +=( double s ) {  x += s;         return *this; }
CE      Vec&    operator -=( double s ) {  x -= s;         return *this; }
CE      Vec&    operator *=( double s ) {  x *= s; y *= s; return *this; } // Умножение на скаляр
CE      Vec&    operator /=( double s ) {  x /= s; y /= s; return *this; } // Деление на скаляр
        OPER_COMM( +, double); OPER_COMM( *, double); OPER_NOCOMM( -, double); OPER_NOCOMM( /, double)

CE      Vec&    operator *=( Ĵ ) { swap( x, y); x = -x; return *this; } // Умножение на мнимую
CE      Vec&    operator /=( Ĵ ) { swap( x, y); y = -y; return *this; } // Деление на мнимую
        OPER_COMM( *, Ĵ)
CEfrnd  Vec     operator / ( Vec v,        Ĵ) { return {  v.y, -v.x}; } //  векторов как комплексных чисел
CEfrnd  Vec     operator / ( Ĵ, const Vec& v) { return {  v.y,  v.x}; } // на мнимую

CE      Vec&    operator *=( Sgn a ) {  x *= a._; y *= a._; return *this; } // Умножение на знак
        OPER_COMM( *, Sgn)

friend  std::ostream& operator<<( std::ostream &os, const Vec& v )
        {
                static Vec last = { NAN, NAN };
                if( last != v )
                {
                        os << std::setw(8) <<v.x << std::setw(12) <<v.y << '\n';
                        last = v;
                }
                return os;
        };
};

// Нормированный вектор (unit vector)
struct NVec: public Vec
{
#ifdef NDEBUG
protected:
#endif
CE      NVec( Ф, const Vec& v      ): Vec( v   ) {}
CE      NVec( Ф, double x, double y): Vec( x, y) {}
public:
        using This = NVec;
        using Super = Vec;

CE      NVec( Ĵ            ): Vec( 0, 1                        ) {}
CExplct NVec( const Vec& v ): Vec( v / v.abs()                 ) {}
CExplct NVec( double sin   ): Vec( ce::sqrt( 1 - sin*sin), sin ) {}        
CE      NVec( const Vec& v, double *p ) // нормализующий всё подряд к-тор: нормализует вектор v и ещё что дадут (*p)
        : Vec( v)
        {
                assert( p);
                double len = v.abs();
                *(Vec *)(this) /= len;
                *p             /= len; 
        }

CE      double  abs²       () const { return  1.         ; } // Длина² = 1, единичный вектор же
CE      double  abs        () const { return  1.         ; } // Длина  = 1, единичный вектор же
CE      NVec    operator - () const { return {ф, -x, -y }; } // поворот на 180°
CE      NVec    operator ~ () const { return {ф,  x, -y }; } // Сопряжённое (conjugate) число (зеркально отраженный вектор)
CE      NVec    conj       () const { return {ф,  x, -y }; } // Сопряжённое (conjugate) число (зеркально отраженный вектор)
CE      NVec    inv        () const { return {ф,  x, -y }; } // 1/v - обратная (inverse) величина
CE      double  ψarg       () const { double x1 = x + 1; return y>=0 ? -x1 : x1;} // псевдоугол (для сравнений)
        FUNC( abs²); FUNC( abs); FUNC( conj); FUNC( inv); FUNC( ψarg);

CEfrnd  NVec    cis        ( double β );// { return {ф, cos(β), sin(β)}; } // Единичный вектор, повернутый на угол φ

CE      bool    operator > ( const NVec& v) const { return ψarg() > v.ψarg();}
CE      bool    operator < ( const NVec& v) const { return ψarg() < v.ψarg();}

CE      NVec&   operator *=( const NVec& v) { Vec::operator*=( v); return *this; } // Умножение единичных векторов как комплексных чисел
CE      NVec&   operator /=( const NVec& v) { *this *= v.conj()  ; return *this; } // Деление векторов как комплексных чисел
        OPER_THIS( *); OPER_THIS( /);

CE      Vec&    operator *=( const  Vec& v) { Vec::operator*=( v); return *this; } // Умножение единичных векторов как комплексных чисел
CE      Vec&    operator /=( const  Vec& v) { *this *= v.conj()  ; return *this; } // Деление векторов как комплексных чисел
        OPER_SUPER( *); OPER_SUPER( /);

CE      NVec&   operator *=( Ĵ ) { swap( x, y); x = -x; return *this; } // Умножение на мнимую
CE      NVec&   operator /=( Ĵ ) { swap( x, y); y = -y; return *this; } // Деление на мнимую
        //OPER_COMM( *, Ĵ); OPER_COMM( /, Ĵ);
CEfrnd  NVec   operator * ( NVec v,        Ĵ) { return {ф, -v.y,  v.x}; } // Умножение на мнимую
CEfrnd  NVec   operator * ( Ĵ, const NVec& v) { return {ф, -v.y,  v.x}; } // Умножение на мнимую
CEfrnd  NVec   operator / ( NVec v,        Ĵ) { return {ф,  v.y, -v.x}; } // Деление на мнимую
CEfrnd  NVec   operator / ( Ĵ, const NVec& v) { return {ф,  v.y,  v.x}; } // Деление мнимой

CE      NVec&    operator *=( Sgn a ) {  x *= a._; y *= a._; return *this; } // Умножение на знак
        OPER_COMM( *, Sgn)

static  const NVec î;
static  const NVec ĵ;
};

CE NVec cis( double β ) { return {ф, cos(β), sin(β) }; }

CE Vec& Vec::operator /=( const NVec& n)
{
        *this *= n.inv();
        return *this;
}

// единичный вектор вдоль оси X
CE const NVec NVec::î = {ф, 1, 0 };
// единичный вектор вдоль оси Y, мнимая единица
CE const NVec NVec::ĵ = {ф, 0, 1 };

// единица
CE const NVec î = NVec::î;

CE Vec operator + ( Ĵ, double a) { return { a, 1}; }
CE Vec operator + ( double a, Ĵ) { return { a, 1}; }
CE Vec operator - ( Ĵ, double a) { return {-a, 1}; }
CE Vec operator - ( double a, Ĵ) { return { a,-1}; }
CE Vec operator * ( Ĵ, double a) { return { 0, a}; }
CE Vec operator * ( double a, Ĵ) { return { 0, a}; }

CE Vec operator ""ĵ ( unsigned long long a) { return { 0, double(a)}; }
CE Vec operator ""ĵ (        long double a) { return { 0, double(a)}; }

void   Vec_test()
{
        static_assert( (5 + ĵ)*(7 - 6ĵ) / (3 + ĵ)          == (10 - 11ĵ), "");
        static_assert( (4 + ĵ)*(5 + 3ĵ) + (3 + ĵ)*(3 - 2ĵ) == (28 + 14ĵ), "");

        // это просто должно скомпилиться
        CE NVec u0 = ĵ;        
        CE NVec u1 = î / ĵ;
        CE NVec u3 = NVec( 0.1 );
        CE NVec u4 = u3 * NVec( 0.1 ) / ĵ;
        CE NVec u5 = î/NVec( 0.1 );
        CE NVec u6 = î*NVec( 0.1 );
        CE  Vec u7 = Vec( 0, 1 ) * Vec( 0.995, 0.1 ) ;
        CE NVec u8 = ĵ * NVec( 0.1 );
        CE  Vec u9 = u7 / u8;

        CE NVec a = NVec( 0.1);
        static_assert( a       > î && !( a       < î), "");
        static_assert(~a       > î && !(~a       < î), "");
        static_assert(~a       > a && !(~a       < a), "");
        static_assert( a*ĵ*ĵ*ĵ > a && !( a*ĵ*ĵ*ĵ < a), "");
}
#pragma endregion

struct Line
{
        double p; // растояние от начало координат до прямой
        NVec   n̂; // единичный (|n̂| = 1) перпедикуляр к прямой

        // Прямая, заданой нормальным уравнением (n̅, r̅) = 𝐶 т.е. 𝐴𝑥ᵣ + 𝐵𝑦ᵣ = 𝐶, n̅ = {𝐴, 𝐵}, |n̅| > 0
        // \param[in] n̅ - вектор, нормальный к прямой
        // \param[in] 𝐶 - скалярный параметр
CE      Line( const  Vec& n̅, double 𝐶 ): p( 𝐶 ), n̂( n̅, &p ) {}
        // Прямая, задана нормированным уравнением (n̂, r̅) = 𝑝, |n̂| = 1, 𝑝 ⩾ 0
        // \param[in] n̂ - единичный вектор, нормальный к прямой
        // \param[in] 𝑝 - растояние от начало координат до прямой
CE      Line( const NVec& n̂, double 𝑝 ): p( 𝑝 ), n̂( n̂     ) {};
        // Прямая через точки a̅ и b̅
CE      Line( const Vec& a̅, const Vec& b̅ ): Line( (a̅ - b̅)*ĵ, (a̅ ^ b̅) ) {};

CE      bool    operator ==( const Line& l) const { return eq( p, l.p) && n̂ == l.n̂; }
CE      double  dist       ( const Vec&  r̅) const { return (n̂, r̅) - p             ; } // растояние между прямой и точкой r̅
CE      Vec     proj       ( const Vec&  r̅) const { return r̅ - n̂ * dist( r̅)       ; } // проекция точки r̅ на прямую

        // растояние между точкой r̅ и прямой l
CEfrnd  double dist( const Vec& r̅, const Line& l) { return l.dist( r̅); }
        // растояние между прямой l и точкой r̅
CEfrnd  double dist( const Line& l, const Vec& r̅) { return l.dist( r̅); }
        // проекция точки r̅ на прямую l
CEfrnd  Vec operator >>( const Vec& r̅, const Line& l) { return l.proj( r̅); }

/*
friend  std::ostream& operator<<( std::ostream &os, const Line& obj )
        {
                os << Vec( 0.0, obj.p/obj.n.y)
                   << Vec( obj.p/obj.n.x, 0.0);
                return os;
        };
*/
};
struct Vertical: public Line
{
        // Вертикаль пересекающая ось 𝑋 в точке x0
CE      Vertical( double x0 ): Line( NVec::î, x0) {};
};
struct Horizontal: public Line
{
        // Горизонталь пересекающая ось 𝑌 в точке y0
CE      Horizontal( double y0 ): Line( NVec::ĵ, y0) {};
};
void   Line_test()
{
        CE Vertical l_test( 2 );
        static_assert(       l_test.n̂        == î, "");
        static_assert(       l_test.n̂ * ĵ    == ĵ, "");
        static_assert( dist( l_test, {3, 1}) == 1, "");
}

struct Segment: public Line
{
        Vec p̅1, p̅2;

CE      Segment( const Vec& p1_, const Vec& p2_ )
        : Line( p1_, p2_ )
        , p̅1( p1_), p̅2( p2_)
        {};

        // TODO стремный конст-ор, как бы его спрятать
CE      Segment( const Line& l, const Vec& p1_, const Vec& p2_ )
        : Line( l )
        , p̅1( p1_), p̅2( p2_)
        {};

        // обмен концов (рисоваться будет в другую сторону)
CE      Segment operator - () const { return Segment( this->p̅2, this->p̅1 ); }

friend  std::ostream& operator<<( std::ostream &os, const Segment& _ )
        {
                os <<_.p̅1 <<_.p̅2;
                return os;
        };
};

struct Circle
{
        Vec    o̅; // центр окружности
        double R; // радиус окружности

CE      Circle( const Vec& center, double radius = 0. )
        : R( radius), o̅( center )
        {};

CE      bool   operator ==( const Circle& c) const { return eq( R, c.R) && o̅ == c.o̅; }
CE      Circle operator - (                ) const
        {
                Circle a = *this;
                a.R = -a.R;
                return a;
        }

CE      Vec intersect( const Line&    l ) const
        {
                double h = dist( o̅, l);               // расстояние от центра окружности до прямой
                //Vec    p̅ = (o̅ - l.n̂ * h);               // точка проекции центра окружности на прямую
                //NVec   v̂ = l.n̂ * ĵ; // направляющий вектор прямой l
                //return p̅ + v̂ * ce::sqrt(²(R) - ²(h)) * sign;
                return o̅ - l.n̂ * Vec( h, ce::sqrt( ²(R) - ²(h)) );
        }
CE      Vec intersect( const Circle& c2 ) const
        {
                Line radical_line( (c2.o̅ - o̅) * 2. 
                                 , abs²(c2.o̅) - abs²(o̅) - ²(c2.R) + ²(R)
                                 );
                return intersect( radical_line );
        }

CEfrnd  Circle tangent( const Circle& c1, const Circle& c2, double R )
        {
                Vec center = Circle( c1.o̅, R+c1.R).intersect( Circle( c2.o̅, R+c2.R) );
                return Circle( center, R );
        }

        // нормаль к касательной к окружностям *this и c2
CE      NVec tangent_norm( const Circle& c2 ) const
        {
                double sinφ = c2.R - R;       // φ - угол между линией центров и касательной
                NVec â( c2.o̅ - o̅, &sinφ );    // направляющий вектор линии центров
                return â * NVec( sinφ) / ĵ;   // повернуть â на φ-90°
        }

        // касательная к окружностям *this и c2
CE      Line tangent( const Circle& c2 ) const
        {
                NVec n̂ = tangent_norm( c2 );
                return Line( n̂, (n̂, o̅) - R );
        }

        // точка касания касательной к окружности исходящей из точки p
CE      Vec tangent_point( const Vec& p ) const
        {
                NVec n̂ = tangent_norm( {p, 0} );
                return o̅ - n̂*R;
        };

        void print( std::ostream &os) const
        {
                CE const int segs = 40;

                //CE const NVec m̂¹⁰ = cis( 2π/segs);
                CE const NVec m̂( 2π/segs/10 ); // типа 𝛼 ≈ sin 𝛼
                CE const NVec m̂⁵  = m̂*m̂*m̂*m̂*m̂;
                CE const NVec m̂¹⁰ = m̂⁵*m̂⁵;

                NVec n̂ = î;
                for( int i = segs; i --> 0; )
                {
                        os << (o̅ + R*n̂);
                        n̂ *= m̂¹⁰;
                }
        };

        void print00( std::ostream &os, double α1, double α2 ) const
        {
                double Δα = α2 - α1; // угол поворота
                // кол. сегментов, на круг - 160 сегментов, примерно
                int segments = static_cast< int>( round( abs( Δα / 2π * 160)));
                NVec m̂ = cis( Δα / segments);   // ед. век. повернутый на угол Δα/segments

                Vec  r̅ = abs( R) * cis( α1);    // радиус от центра окружности
                Vec  s̅ = m̂*r̅ - r̅;               // сегментик, которым рисуем окружность
                r̅ += o̅;
                for( ; segments >= 0; --segments )
                {
                        os << r̅;
                        r̅ += ( s̅ *= m̂ );
                }
        };
friend  std::ostream& operator<<( std::ostream &os, const Circle& _ ) { _.print( os); return os; };
};
void   Circle_test()
{
        CE Circle c_test( {0, 0}, 3 );
        CE Vec o1( 5, 0);
        CE Vec o2( 5, 5);

        //CE auto b = Circle( o1, -2).tangent( Circle( o2, -2));
        //CE auto a = Line( Vec( 0.6,  0.8), 5);

        static_assert( Circle( o1,  2).tangent( Circle( o2,  2)) == Line( Vec( 1.0,   .0), 3), "");
        static_assert( Circle( o1, -2).tangent( Circle( o2,  2)) == Line( Vec( 0.6,  0.8), 5), "");
        static_assert( Circle( o1,  2).tangent( Circle( o2, -2)) == Line( Vec( 0.6, -0.8), 1), "");
        static_assert( Circle( o1, -2).tangent( Circle( o2, -2)) == Line( Vec( 1.0,   .0), 7), "");

        static_assert( c_test.intersect( Line( {3, 0}, {0, 3})) == Vec( 3,    0), "");
        static_assert( c_test.intersect( Vertical( 3)         ) == Vec( 3,    0), "");
        static_assert( c_test.intersect( Circle( {4, 0}, 3)   ) == Vec( 2,-5⅟²), "");

        static_assert( c_test.tangent_point( Vec( 5, 0)       ) == Vec(1.8, 2.4), "");
        //static_assert( c_test.tangent_point( Vec( 5, 0), minus) == Vec(1.8,-2.4), "");
        /*
        static_assert( c_test.tangent( Circle( {   2, 0}, 1)       ) == Line( NVec::î             ,  3), "");
        static_assert( c_test.tangent( Circle( {2⅟², 0}, 2)       ) == Line( Vec( 0.5⅟²,  0.5⅟²), -3), "");
        static_assert( c_test.tangent( Circle( {2⅟², 0}, 2), minus) == Line( Vec( 0.5⅟², -0.5⅟²), -3), "");
        static_assert( c_test.tangent( Circle( {   9, 0},-3)       ) == Line( Vec(   2./3,  5⅟²/3), -3), "");
        */
        static_assert( tangent( c_test, Circle( {  6, 0}, 3),  3 ) == Circle( {  3,-3*3⅟²},  3), "");
        static_assert( tangent( c_test, Circle( {  2, 0}, 1), 10 ) == Circle( { 13,      0}, 10), "");
}

struct Arc: public Circle
{
        NVec n̂₁; // единич. вектор из центра дуги на первый конец
        NVec n̂₂; // единич. вектор из центра дуги на второй конец

CE      Arc( const Circle& c, const NVec&_n̂₁, const NVec&_n̂₂)
        : Circle( c), n̂₁(_n̂₁), n̂₂(_n̂₂)
        {}
CE      Arc( const Circle& c, const Vec& start_point, const Vec& next_point )
        : Circle( c)
        , n̂₁( start_point- c.o̅ )
        , n̂₂( next_point - c.o̅ )
        {}
CE      Arc( const Circle& c, const Arc& prev, const Vec& next_point )
        : Arc( c, prev.o̅ + ce::abs(prev.R)*prev.n̂₂, next_point )
        {
                // ставим знак вектору n̂₂ таким, чтоб конечная точка дуги была максимально близко к next_point
                //if( ((next_point - o̅), n̂₂) < 0)
                //        n̂₂ = -n̂₂;
                if( prev.R * R * (prev.n̂₂, n̂₁) < 0 )
                        R = -R;
        }

CE      bool operator ==( const Arc& a) const
        {
                return  Circle::operator==( a)
                        && n̂₁ == a.n̂₁
                        && n̂₂ == a.n̂₂
                        ;
        }

        void print( std::ostream &os) const
        {
                CE const int segs = 160;

                //CE const NVec m̂¹⁰ = cis( 2π/segs);
                CE const NVec m̂( 2π/segs/10 ); // типа 𝛼 ≈ sin 𝛼
                CE const NVec m̂⁵ = m̂*m̂*m̂*m̂*m̂;

                NVec m̂¹⁰ = m̂⁵*m̂⁵;
                bool dir = (R >= 0);
                if( !dir )
                        m̂¹⁰ = ~m̂¹⁰;

                Vec r̅ = abs(R) * n̂₁;
                os << (o̅ + r̅);
                double a₂₁ = ψarg( n̂₂/n̂₁);

                for( NVec n̂ = m̂¹⁰; (ψarg(n̂) < a₂₁) == dir; n̂ *= m̂¹⁰ )
                        os << (o̅ + r̅*n̂);
        };
friend  std::ostream& operator<<( std::ostream &os, const Arc& _) { _.print( os); return os; };
};

// зазор для соблюдения постулата Жуковского-Чаплыгина (Kutta condition)
CE const Vec TE1( 1.,  0.00001); // задняя кромка верх
CE const Vec TE2( 1., -0.00001); // задняя кромка низ

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

        std::cout << "PIPE " << D_abs << 'x' << s_abs << '-' << b_abs << '\n'
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
        CE Vec    p_end    = c_bottom.tangent_point( TE2 );

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

int main( unsigned argc, const char *argv[])
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
        Vec    p_end    = c_bottom.tangent_point( TE2 );

        Arc a_top   ( с_top   , TE1  , c_lef.o̅   );
        Arc a_lef   ( c_lef   , a_top, c_le.o̅    );
        Arc a_le    ( c_le    , a_lef, c_bottom.o̅);
        Arc a_bottom( c_bottom, a_le , p_end     );
        
        std::cout << "PIPE " << D_abs << 'x' << s_abs << '-' << b_abs << " r" << ler_abs << " f" << lef_abs << '\n'
                << std::setprecision(5) << std::fixed
                << a_top << a_lef << a_le << a_bottom << p_end << TE2;

        return 0;
}
