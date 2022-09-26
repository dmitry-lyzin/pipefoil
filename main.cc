#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <clocale>
#include <cassert>

#define CE      constexpr
#define ₠       constexpr
#define CEstat  constexpr static
#define CExplct constexpr explicit
#define CEfrnd  constexpr friend

// epsilon ≠ std::numeric_limits< double>::epsilon();
CE const double ε = 1e-14;

CE const double π = 3.14159265358979323846;
CE double operator ""π ( unsigned long long a) { return π * a; }

class Sing
{
        int value;
CE      Sing( int x): value(x) {};
public:
CE      operator int   () const { return  value; };
CE      Sing operator- () const { return -value; };
static  const Sing plus;
static  const Sing minus;
};
CE const Sing Sing::plus ( 1 );
CE const Sing Sing::minus(-1 );

CE const Sing plus  = Sing::plus;
CE const Sing minus = Sing::minus;

namespace ce
{
        double CE sqrt_Newton_Raphson( double x, double curr, double prev)
        {
                return curr == prev
                        ? curr
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
                return x >= 0 && x < std::numeric_limits<double>::infinity()
                        ? sqrt_Newton_Raphson(x, x, 0)
                        : std::numeric_limits<double>::quiet_NaN();
        }
}

CE inline double ²( double a)
{
        return a * a;
}

CE inline bool eq( double a, double b)
{
        if( (a -= b) < 0 )
                a = -a;
        return a < ε;
}

class Vec
{
protected:
        double x, y;
public:
CE      Vec(                     ): x(NAN), y(NAN) {}
CE      Vec( double x_, double y_): x( x_), y( y_) {}

CE      double  operator , ( const Vec& v) const { return  x*v.x + y*v.y  ;     } // Скалярное произведение
CE      double  operator ^ ( const Vec& v) const { return  x*v.y - y*v.x  ;     } // Псевдоскалярное (косое) произведение
CEfrnd  double  ²          ( const Vec& v)       { return  (v, v)         ;     } // Длина² (Скалярное произведение самого на себя)
CEfrnd  Vec     L          ( const Vec& v)       { return {  -v.y,   v.x };     } // поворот на 90°
CE      Vec     operator - (             ) const { return {    -x,    -y };     } // поворот на 180°

CE      Vec     operator + ( const Vec& v) const { return { x+v.x, y+v.y };     }
CE      Vec     operator - ( const Vec& v) const { return { x-v.x, y-v.y };     }
CE      Vec     operator * ( double s    ) const { return { x*s  , y*s   };     } // Умножение на скаляр
CE      Vec     operator / ( double s    ) const { return { x/s  , y/s   };     } // Деление на скаляр

CE      Vec&    operator +=( const Vec& v) { x+=v.x, y+=v.y; return *this;  }
CE      Vec&    operator -=( const Vec& v) { x-=v.x, y-=v.y; return *this;  }
CE      Vec&    operator *=( double s    ) { x *= s, y *= s; return *this;  } // Умножение на скаляр
CE      Vec&    operator /=( double s    ) { x /= s, y /= s; return *this;  } // Деление на скаляр

CE      bool    operator ==( const Vec& v) const { return eq( x, v.x) && eq( y, v.y);   }
CE      bool    operator !=( const Vec& v) const { return !(*this == v);                }

friend  std::ostream& operator<<( std::ostream &os, const Vec& _ )
        {
                static Vec last;
                if( last != _ )
                {
                        os << std::setw(8) <<_.x << std::setw(12) <<_.y << '\n';
                        last = _;
                }
                return os;
        };

CEfrnd  double angle( const Vec& v̅ )
        {
                if( v̅.x > +0. ) return atan( v̅.y / v̅.x );
                if( v̅.x < -0. ) return atan( v̅.y / v̅.x ) + π;

                if( v̅.y > +0. ) return π/2;
                if( v̅.y < -0. ) return 3π/2;

                return NAN;
        }
};

class Rot;

// Единичный вектор (unit vector)
class UVec: public Vec
{
        friend class Rot;

CE      UVec( const Vec& v      ): Vec( v    ) {}
CE      UVec( double x, double y): Vec( x, y ) {}

public:                                                  
CE      UVec( double θ                ): Vec( cos(θ), sin(θ)) {}
        // нормализующий всё подряд к-тор: нормализует вектор v̅ и ещё что дадут (*p)
CE      UVec( const Vec& v̅, double *p ): Vec( v̅)
        {
                double l = ce::sqrt((v̅, v̅));
                *this /= l;
                *p    /= l; 
        }

CEfrnd  double  ²          ( const UVec& v)       { return 1.           ;} // Скалярное произведение самого на себя = 1
                                                                           // при поворотах нормализованность сохраняется
CEfrnd  UVec    L          ( const UVec& v)       { return {-v.y, v.x } ;} // поворот на 90°
CE      UVec    operator - (              ) const { return {  -x,  -y } ;} // поворот на 180°

static const UVec î;
static const UVec ĵ;
};

CE const UVec UVec::î = { 1., 0. }; // единичный вектор вдоль оси 𝑋
CE const UVec UVec::ĵ = { 0., 1. }; // единичный вектор вдоль оси 𝑌

// Матрица поворота
class Rot
{
        UVec û1, û2;

CE      Rot( double cos, double sin ): û1( cos, -sin), û2( sin, cos) {};
public:
CE      Rot( double cos, Sing sing  ): Rot( cos   , sing * ce::sqrt( 1 - cos*cos) ) {};
CE      Rot( double θ               ): Rot( cos(θ), sin(θ)                        ) {};

CE       Vec operator * ( const  Vec& v̅ ) const { return { (û1, v̅), (û2, v̅) }; }
CE      UVec operator * ( const UVec& û ) const { return { (û1, û), (û2, û) }; }
};

struct Line
{
        double p; // растояние от начало координат до прямой
        UVec   n̂; // единичный (|n̂| = 1) перпедикуляр к прямой

        // Прямая, заданой нормальным уравнением (n̅, r̅) + 𝐶 = 0, |n̅| > 0
        // \param[in] n̅ - вектор, нормальный к прямой
        // \param[in] 𝐶 - скалярный параметр
CE      Line( const Vec& n̅, double 𝐶 )
        : p( -𝐶 ), n̂( n̅, &p )
        {}
        // Прямая, задана нормированным уравнением (n̂, r̅) = 𝑝, |n̂| = 1, 𝑝 ⩾ 0
        // \param[in] n̂ - единичный вектор, нормальный к прямой
        // \param[in] 𝑝 - растояние от начало координат до прямой
CE      Line( const UVec& n̂, double 𝑝 )
        : n̂( n̂ ), p( 𝑝 )
        {};
        // Прямая через точки a̅ и b̅
CE      Line( const Vec& a̅, const Vec& b̅ )
        : Line( L(b̅ - a̅), a̅ ^ b̅ )
        {};

        // растояние между прямой и точкой r̅
CE      double distance( const Vec& r̅) const   { return (n̂, r̅) - p; }
        // растояние между точкой r̅ и прямой l
CEfrnd  double ρ( const Vec& r̅, const Line& l) { return l.distance( r̅); }
        // растояние между прямой l и точкой r̅
//CEfrnd  double ρ( const Line& l, const Vec& r̅) { return l.distance( r̅); }

        // проекция точки r̅ на прямую l
CEfrnd  Vec operator >>( const Vec& r̅, const Line& l)
        {
                return r̅ - l.n̂ * ρ( r̅, l);
        }

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
CE      Vertical( double x0 ): Line( UVec::î, x0) {};
};
struct Horizontal: public Line
{
        // Горизонталь пересекающая ось 𝑌 в точке y0
CE      Horizontal( double y0 ): Line( UVec::ĵ, y0) {};
};

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
/*
        // центр проходящей через две точки окружность с радиусом R
CEstat  Vec center( const Vec& p1, const Vec& p2, double R, Sing case_)
        {
                double d² = ²(p1 - p2);
                //double h = sqrt( ²(R) - ²(d/2.));
                double h_div_d = case_ * ce::sqrt( (²(R)/d² - 1./4.) );

                return { (p2.x + p1.x)/2. + (p2.y - p1.y) * h_div_d
                       , (p2.y + p1.y)/2. - (p2.x - p1.x) * h_div_d
                       };                      
        };
        // проходящая через две точки окружность с радиусом R_
CE      Circle( const Vec& p1, const Vec& p2, double R_, Sing case_)
        : R( R_), o( center( p1, p2, R_, case_))
        {};
*/

CE      Circle( const Vec& center, double radius = 0. )
        : R( radius), o̅( center )
        {};

CE      Circle operator - () const
        {
                Circle a = *this;
                a.R = -a.R;
                return a;
        }

        void print( std::ostream &os, double α1, double α2 ) const
        {
                double Δα = α2 - α1; // угол поворота
                // кол. сегментов, на круг - 160 сегментов, примерно
                int segments = static_cast< int>( round( abs( Δα / 2π * 160)));

                Rot M = Δα / segments;          // матрица поворота на угол Δα/segments
                Vec r̅ = UVec( α1) * abs( R);    // радиус от центра окружности
                Vec s̅ = M*r̅ - r̅;                // сегментик, которым рисуем окружность
                r̅ += o̅;
                for( ; segments >= 0; --segments )
                {
                        os << r̅;
                        r̅ += s̅;
                        s̅ = M*s̅;
                }
        };
friend  std::ostream& operator<<( std::ostream &os, const Circle& _ )
        {
                _.print( os, 0, 2π );
                return os;
        };
};

struct Arc: public Circle
{
        Vec p̅1, p̅2;
        
CE      Arc( const Circle& circle, const Vec& start, const Vec& finish )
        : Circle( circle)
        , p̅1( start), p̅2( finish)
        {};

CE      Arc( const Vec& center, double R_, const Vec& start, const Vec& finish )
        : Circle( center, R_ )
        , p̅1( start), p̅2( finish)
        {};

        // сменить сектор откружности
CE      Arc operator - () const
        {
                Arc a = *this;
                a.R = -a.R;
                return a;
        }

friend  std::ostream& operator<<( std::ostream &os, const Arc& _ )
        {
                double α1 = angle( _.p̅1 - _.o̅ );
                double α2 = angle( _.p̅2 - _.o̅ );

                if( _.R < 0. )
                {
                        if( α1 >= α2) α2 += 2π;
                        else          α1 += 2π;
                }

                _.print( os, α1, α2 );
                return os;
        };
};

#pragma region // функции поиска пересечений
CE Vec cross( const Circle& c, const Line& l, Sing sing = plus )
{
        double h = ρ( c.o̅, l);                  // расстояние от центра окружности до прямой
        Vec    p̅ = c.o̅ - l.n̂ * h;               // точка проекции центра окружности на прямую
        UVec   û = L( l.n̂);                     // направляющий вектор прямой l
        Vec    v̅ = û * ce::sqrt(²(c.R) - ²(h)); // v̅ сонаправлен û но длиной √(𝘳² − 𝘩²)
        return p̅ - v̅ * sing;
}

CE Vec cross( const Circle& c1, const Circle& c2, Sing sing = plus )
{
        Line radical_line( (c2.o̅ - c1.o̅) * 2. 
                         , ²(c1.o̅) - ²(c2.o̅) + ²(c2.R) - ²(c1.R)
                         );

        return cross( c1, radical_line, -sing ); 
}

// сахарок синтаксический

// первая точка пересечения окружности и линии
CE Vec operator & ( const Circle& c,  const Line& l    ) { return cross( c, l, plus   ); }
// первая точка пересечения линии и окружности
CE Vec operator & ( const Line& l,    const Circle& c  ) { return cross( c, l, plus   ); }
// вторая точка пересечения окружности и линии
CE Vec operator ^ ( const Circle& c,  const Line& l    ) { return cross( c, l, minus  ); }
// вторая точка пересечения линии и окружности
CE Vec operator ^ ( const Line& l,    const Circle& c  ) { return cross( c, l, minus  ); }
// первая точка пересечения окружности и окружности
CE Vec operator & ( const Circle& c1, const Circle& c2 ) { return cross( c1, c2, plus ); }
// вторая точка пересечения окружности и окружности
CE Vec operator ^ ( const Circle& c1, const Circle& c2 ) { return cross( c1, c2, minus); }
#pragma endregion

#pragma region // функции поиска касательных

// касательная к окружностям c1 и c2
CE Line tangent_line( const Circle& c1, const Circle& c2, Sing sing = plus )
{
        Vec c̅ = c2.o̅ - c1.o̅;
        double cosθ = c2.R - c1.R;              // θ - угол между линией центров и касательной
        UVec n̂( c̅, &cosθ );                     // смаштабировать c̅ в n̂ и заодно cosθ
        UVec r̂ = Rot( cosθ, sing) * n̂;          // повернуть n̂ на θ в направлении sing
        return Line( r̂, (r̂, c1.o̅) - c1.R );     // TODO тут где-то баг... r̅ ∥ касат., а должен быть ⊥
}

// касательная к окружностям c1 и c2
CE Segment tangent_segment( const Circle& c1, const Circle& c2, Sing sing = plus )
{
        Line tl = tangent_line( c1, c2, -sing );
        return Segment( tl, c1.o̅ >> tl, c2.o̅ >> tl );
}

// поиск точки касания прямой исходящей из точки p и окружности c
CE Vec tangent_point( const Circle& c1, const Vec& p, Sing sing = plus )
{
        Line tl = tangent_line( c1, Circle( p, 0.), sing );
        return c1.o̅ >> tl;
};

// дуга радиусом R, исходящая из точки p и касательная к окружности c
CE Arc tangent_arc( const Circle& c, const Vec& p, double R, Sing sing = plus )
{
        Vec center = cross( Circle( p, R), Circle( c.o̅, R-c.R), sing );
        return Arc( center, R, p, cross( c, Line(c.o̅, center), -sing) );
}

// окруж. радиусом R, касательная к окружностям c1 и c2
CE Circle tangent_сircle( const Circle& c1, const Circle& c2, double R, Sing sing = plus )
{
        Vec center = cross( Circle( c1.o̅, R-c1.R), Circle( c2.o̅, R-c2.R), -sing );
        return Circle( center, R );
}
#pragma endregion

// очередная дуга из цепочки соприкасающихся окружностей
CE Arc chain_arc( const Vec& prev_arc_end, const Circle& current, const Circle& next, Sing sing = plus )
{
        return Arc( current, prev_arc_end, cross( next, Line( next.o̅, current.o̅), -sing) );
}

// зазор для соблюдения постулата Жуковского-Чаплыгина (Kutta condition)
CE const Vec TE1( 1.,  0.00001); // задняя кромка верх
CE const Vec TE2( 1., -0.00001); // задняя кромка низ

#ifndef NDEBUG
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
        CE double s½  = s / 2.         ; // относительная полутолщина стенки трубы
        CE double lr1 = lr1_abs / b_abs; // относительный радиус скругл. носка 1
        CE double lr2 = lr2_abs / b_abs; // относительный радиус скругл. носка 2

        CE Circle  c_lead({0, s½}, s½ );
        CE Circle  с_top = tangent_сircle( c_lead, TE1, R );
        CE Circle  c_bottom( с_top.o̅, R-s );
        CE Segment s_end = tangent_segment( c_bottom, TE2 );
        CE Circle  c_l1( Vertical( lr1) ^ Circle( с_top.o̅, c_bottom.R+lr1), lr1 );
        CE Circle  c_l2 = tangent_сircle( с_top, c_l1, lr2 );

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
        /*
        // круглая ПК немного сверху зализана и заостренная ЗК 
        // Cy/Cx = 25 при α ≈ 5° плоская полка α ≈ 0÷5°, срыв при α ≈ 7,5°
        CE double b_abs   =  50.0; // хорда профиля
        CE double D_abs   = 160.0; // диаметр трубы
        CE double s_abs   =   4.7; // толщина стенки трубы
        CE double lr1_abs =   1.0; // радиус скругл. носка 1
        CE double lr2_abs =  20.0; // радиус скругл. носка 2
        */
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
        CE Circle  с_top = tangent_сircle( c_lead, TE1, R );
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
        CE Circle  c_bottom = tangent_сircle( -c_l1, c_trail, R-s );
        CE Circle  с_top( c_bottom.o̅, R );
        CE Circle  c_l2  = tangent_сircle( c_l1, с_top, lr2, minus );
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
        CE Circle  c_bottom = tangent_сircle( -c_l1, c_trail, R-s );
        CE Circle  с_top( c_bottom.o̅, R );
        CE Circle  c_l2  = tangent_сircle( c_l1, с_top, lr2, minus );
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
}
int test5()
{
        // круглая ПК немного сверху зализана и заостренная ЗК
        // тупой пик, Cy/Cx > 30 α ≈ 4.5÷9.5°, Cy/Cx ≈ 35 α ≈ 6.6°, срыв α ≈ 15°
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
        CE Circle  c_bottom = tangent_сircle( -c_l1, c_trail, R-s );
        CE Circle  с_top( c_bottom.o̅, R );
        CE Circle  c_l2  = tangent_сircle( c_l1, с_top, lr2, minus );
        CE Segment s_end = tangent_segment( c_bottom, TE2 );

        CE Arc a_top = chain_arc( TE1     , с_top, c_l2     );
        CE Arc a_l2  = chain_arc( a_top.p̅2, c_l2 , c_l1     );
        CE Arc a_l1  = chain_arc( a_l2.p̅2 , c_l1 , c_bottom, minus );
        CE Arc a_bottom( c_bottom, a_l1.p̅2, s_end.p̅1 );

        std::cout << "PIPE " << D_abs << 'x' << s_abs << '-' << b_abs << " r" << lr1_abs << " f" << lr2_abs << '\n'
                << std::setprecision(5) << std::fixed
                << a_top << a_l2 << -a_l1 << a_bottom << s_end;

        return 0;
}
#endif

int main( unsigned argc, const char *argv[])
{
        static CE char *param_name[] =
        { "диаметр_трубы"
        , "толщина_стенки_трубы"
        , "хорда"
        , "радиус_передней_кромки"
        , "радиус_зализа_над_передней_кромкой"
        };

        double param[ std::size( param_name)];

        double &D_abs   = param[0]; // 160.0; // диаметр трубы
        double &s_abs   = param[1]; //   4.7; // толщина стенки трубы
        double &b_abs   = param[2]; //  50.0; // хорда профиля
        double &ler_abs = param[3]; //   1.5; // радиус передней кромки
        double &lef_abs = param[4]; //  18.0; // радиус скругл. передней кромки

#ifndef NDEBUG
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
#endif

        if( argc <= std::size( param_name) )
        {
                setlocale( LC_ALL, "" );
                std::cerr << "Параметры запуска:\npipefoil";
                for( unsigned i = 0; i < std::size( param_name); ++i)
                        std::cerr << ' ' << param_name[ i ];
                std::cerr << '\n';
                return 1;
        }

        for( unsigned i = 0; i < std::size( param); ++i)
        {
                char *last_char;
                param[i] = strtod( argv[i+1], &last_char);
                if( *last_char)
                {
                        setlocale( LC_ALL, "" );
                        std::cerr << "Ошибка параметра " << param_name[i] << "\n"
                                     "\"" << argv[i+1] << "\" распозналось как \"" << param[i] << "\"\n";
                        return 1;
                }
        }

        double R   = D_abs/2./ b_abs; // относительный внешний радиус трубы
        double s   = s_abs   / b_abs; // относительная толщина стенки трубы
        double ler = ler_abs / b_abs; // относительный радиус передней кромки
        double lef = lef_abs / b_abs; // относительный радиус скругл. передней кромки

        Circle  c_le   ( {ler, ler}, -ler );
        Circle  c_trail( TE1       , -s   );
        Circle  c_bottom = tangent_сircle( c_le, c_trail, R-s );
        Circle  с_top( c_bottom.o̅, R );
        Circle  c_lef = tangent_сircle( -c_le, с_top, lef, minus );
        Segment s_end = tangent_segment( c_bottom, TE2 );

        Arc a_top = chain_arc( TE1     , с_top, c_lef    );
        Arc a_lef = chain_arc( a_top.p̅2, c_lef, c_le     );
        Arc a_le  = chain_arc( a_lef.p̅2, c_le , c_bottom, minus );
        Arc a_bottom( c_bottom, a_le.p̅2, s_end.p̅1 );

        std::cout << "PIPE " << D_abs << 'x' << s_abs << '-' << b_abs << " r" << ler_abs << " f" << lef_abs << '\n'
                  << std::setprecision(5) << std::fixed
                  << a_top << a_lef << a_le << a_bottom << s_end;

        return 0;
}
