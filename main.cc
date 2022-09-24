#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <clocale>
#include <cassert>

#define CE      constexpr
#define CEstat  constexpr static
#define CExplct constexpr explicit
#define CEfrnd  constexpr friend

CE const double π = 3.14159265358979323846;

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
        return a < 1e-14; // std::numeric_limits< double>::epsilon() ?
}

struct Vec
{
        double x, y;
CE      Vec(                     ): x(NAN), y(NAN) {}
CE      Vec( double x_, double y_): x( x_), y( y_) {}

CE      double  operator , ( const Vec& v) const { return  x*v.x + y*v.y  ;     } // Скалярное произведение
CEfrnd  double  ²          ( const Vec& v)       { return  (v, v)         ;     } // Длина² (Скалярное произведение самого на себя)
CE      Vec     operator + ( const Vec& v) const { return { x+v.x, y+v.y };     }
CE      Vec     operator - ( const Vec& v) const { return { x-v.x, y-v.y };     }
CEfrnd  Vec     L          ( const Vec& v)       { return {  -v.y,   v.x };     } // поворот на 90°
CE      Vec     operator - (             ) const { return {-x    ,-y     };     } // поворот на 180°
CE      Vec     operator * ( double s    ) const { return { x*s  , y*s   };     } // Умножение на скаляр
CE      Vec     operator / ( double s    ) const { return { x/s  , y/s   };     } // Деление на скаляр

CE      Vec&    operator +=( const Vec& v) { x+=v.x, y+=v.y; return *this;  }
CE      Vec&    operator -=( const Vec& v) { x-=v.x, y-=v.y; return *this;  }
CE      Vec&    operator *=( double s    ) { x *= s, y *= s; return *this;  } // Умножение на скаляр
CE      Vec&    operator /=( double s    ) { x /= s, y /= s; return *this;  } // Деление на скаляр

CE      bool    operator ==( const Vec& v) const { return eq( x, v.x) && eq( y, v.y);   }
CE      bool    operator !=( const Vec& v) const { return !(*this == v);                }

friend  std::ostream& operator<<( std::ostream &os, const Vec& v )
        {
                static Vec last;
                if( last != v )
                {
                        os << std::setw(8) << v.x << std::setw(12) << v.y << '\n';
                        last = v;
                }
                return os;
        };

//CEfrnd  double cos( const Vec& l, const Vec& r) { return (l, r) / (~l * ~r); }

CEfrnd  double angle( const Vec& v )
        {
                if( v.x > +0. ) return atan( v.y / v.x );
                if( v.x < -0. ) return atan( v.y / v.x ) + π;

                if( v.y > +0. ) return π/2;
                if( v.y < -0. ) return π*3/2;

                return NAN;
        }
};

struct NVec: public Vec // Единичный вектор
{
private:
CExplct NVec( const Vec& v ): Vec( v ) {}
public:                                                  
CE      NVec( double x_, double y_): Vec( x_, y_) { /* TODO assert( (*this, *this) == 1.);*/ }
        // нормализующий всё подряд к-тор: нормализует вектор v и ещё что дадут (norm)
CE      NVec( const Vec& v, double* norm ): Vec( v)
        {
                double l = ce::sqrt((v, v));
                *this /= l;
                *norm /= l; 
        }

CEfrnd  double  ²          ( const NVec& v)       { return 1.           ;} // Скалярное произведение самого на себя = 1
                                                                           // при поворотах нормализованность сохраняется
CEfrnd  NVec    L          ( const NVec& v)       { return {-v.y, v.x } ;} // поворот на 90°
CE      NVec    operator - (              ) const { return {  -x,  -y } ;} // поворот на 180°

//CEfrnd  double cos( const NVec& l, const NVec& r) { return (l, r);         }
//CEfrnd  double cos( const  Vec& l, const NVec& r) { return (l, r) / ~l;    }
//CEfrnd  double cos( const NVec& l, const  Vec& r) { return (l, r) / ~r;    }
};

struct Matrix2x2
{
        NVec s1, s2;
CE      Matrix2x2( const NVec& v1, const NVec& v2 ): s1( v1), s2( v2) {};
CE       Vec operator * ( const  Vec& v) const { return { (s1, v), (s2, v) }; }
CE      NVec operator * ( const NVec& v) const { return { (s1, v), (s2, v) }; }
};

struct Line
{
        double p;
        NVec n;

        // Прямая, заданой нормальным уравнением (𝐧, 𝐫) + 𝐶 = 0, |𝐧| > 0
        // \param[in] 𝐧 - вектор, нормальный к прямой
        // \param[in] 𝐶 - скалярный параметр
CE      Line( const Vec& 𝐧, double 𝐶 )
        : p( -𝐶 ), n( 𝐧, &p )
        {}
        // Прямая, задана нормированным уравнением (𝐧, 𝐫) = 𝑝, |𝐧| = 1, 𝑝 ⩾ 0
        // \param[in] 𝐧 - единичный вектор, нормальный к прямой
        // \param[in] 𝑝 - растояние от начало координат до прямой
CE      Line( const NVec& 𝐧, double 𝑝 )
        : n( 𝐧 ), p( 𝑝 )
        {};
        // через точки p1 и p2
CE      Line( const Vec& p1, const Vec& p2 )
        : Line( L(p2 - p1), -(p1, L(p2)) )
        {};

friend  std::ostream& operator<<( std::ostream &os, const Line& obj )
        {
                os << Vec( 0.0, obj.p/obj.n.y)
                   << Vec( obj.p/obj.n.x, 0.0);
                return os;
        };
};
struct Vertical: public Line
{
        // Вертикаль пересекающая ось 𝑋 в точке x0
CE      Vertical( double x0 ): Line( NVec( 1., 0. ), x0) {};
};
struct Horizontal: public Line
{
        // Горизонталь пересекающая ось 𝑌 в точке y0
CE      Horizontal( double y0 ): Line( NVec( 0., 1. ), y0) {};
};

struct Segment: public Line
{
        Vec p1, p2;

CE      Segment( const Vec& p1_, const Vec& p2_ )
        : Line( p1_, p2_ )
        , p1( p1_), p2( p2_)
        {};

        // TODO стремный конст-ор, как бы его спрятать
CE      Segment( const Line& l, const Vec& p1_, const Vec& p2_ )
        : Line( l )
        , p1( p1_), p2( p2_)
        {};

        // обмен концов (рисоваться будет в другую сторону)
CE      Segment operator - () const { return Segment( this->p2, this->p1 ); }

friend  std::ostream& operator<<( std::ostream &os, const Segment& obj )
        {
                os << obj.p1 << obj.p2;
                return os;
        };
};

struct Circle
{
        Vec o;
        double R;

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

CE      Circle( const Vec& center, double R_ )
        : R( R_), o( center )
        {};
CE      Circle( const Vec& center            )
        : R( 0.), o( center )
        {};
        // проходящая через две точки окружность с радиусом R_
CE      Circle( const Vec& p1, const Vec& p2, double R_, Sing case_)
        : R( R_), o( center( p1, p2, R_, case_))
        {};

CE      Circle operator - () const
        {
                Circle a = *this;
                a.R = -a.R;
                return a;
        }

        void print( std::ostream &os, double a1, double a2 ) const
        {
                // вычисляем кол. сегментов
                int segments = static_cast< int>( round( abs( (a2 - a1) / (2.* π) * 160))); // на круг - 160 сегментов, примерно
                double Δa = (a2 - a1) / segments;
                double r = abs( R);
                for( ; segments >= 0; --segments )
                {
                        //os << std::setprecision(0) << std::fixed << a1*180/π << "*   " << std::setprecision(5) << std::fixed;
                        os << Vec( o.x + r * cos(a1), o.y + r * sin(a1));
                        a1 += Δa;
                }
        };
friend  std::ostream& operator<<( std::ostream &os, const Circle& c )
        {
                c.print( os, 0., 2.*π );
                return os;
        };
};

struct Arc: public Circle
{
        Vec p1, p2;
        
CE      Arc( const Circle& circle, const Vec& start, const Vec& finish )
        : Circle( circle)
        , p1( start), p2( finish)
        {};

CE      Arc( const Vec& center, double R_, const Vec& start, const Vec& finish )
        : Circle( center, R_ )
        , p1( start), p2( finish)
        {};

        // сменить сектор откружности
CE      Arc operator - () const
        {
                Arc a = *this;
                a.R = -a.R;
                return a;
        }

friend  std::ostream& operator<<( std::ostream &os, const Arc& arc )
        {
                double a1 = angle( arc.p1 - arc.o );
                double a2 = angle( arc.p2 - arc.o );

                if( arc.R < 0. )
                {
                        if( a1 >= a2) a2 += 2.* π;
                        else          a1 += 2.* π;
                }

                arc.print( os, a1, a2 );
                return os;
        };
};

// проекция точки (радиус-вектор 𝐫) на прямую l
CE Vec proj( const Vec& 𝐫, const Line& l)
{
        const NVec& 𝐧 = l.n;           // перпедикуляр к прямой
        return 𝐫 - 𝐧 * ((𝐧, 𝐫) - l.p);  // точка проекции на прямую
}

#pragma region // функции поиска пересечений
CE Vec cross( const Circle& c, const Line& l, Sing sing = plus )
{
        const Vec&   𝐨  = c.o;          // центр окружности
        const double 𝘳² = c.R*c.R;      // квадрат радиуса окружности

        Vec 𝐬 = proj( 𝐨, l);              // точка проекции центра окружности на прямую
        double 𝘩² = ²(𝐬 - 𝐨);           // квадрат расстояния от центра окружности до прямой
        NVec 𝐯 = L( l.n);               // вектор, параллельный прямой l
        Vec 𝐯₁ = 𝐯 * ce::sqrt( 𝘳² - 𝘩²); // 𝐯₁ паралелен 𝐯 но длиной √(𝘳² − 𝘩²)
        return 𝐬 - 𝐯₁ * sing;
}

CE Vec cross( const Circle& c1, const Circle& c2, Sing sing = plus )
{
        Line radical_line( (c2.o - c1.o) * 2. 
                         , ²(c1.o) - ²(c2.o) + ²(c2.R) - ²(c1.R)
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
        Vec c21 = c2.o - c1.o;
        double cosθ = (c2.R - c1.R); // /l; // θ - угол между линией центров и касательной
        NVec c21n( c21, &cosθ );
        double sinθ = ce::sqrt(1 - cosθ*cosθ);

        Matrix2x2 R( {      cosθ, -sing*sinθ }
                   , { sing*sinθ,       cosθ }
                   );

        NVec cr = R * c21n;

        return Line( cr, (cr, c1.o) - c1.R );
}

// касательная к окружностям c1 и c2
CE Segment tangent_segment( const Circle& c1, const Circle& c2, Sing sing = plus )
{
        Line tl = tangent_line( c1, c2, -sing );
        return Segment( tl, proj( c1.o, tl ), proj( c2.o, tl ) );
}

// поиск точки касания прямой исходящей из точки p и окружности c
CE Vec tangent_point( const Circle& c1, const Vec& p, Sing sing = plus )
{
        Line tl = tangent_line( c1, Circle( p, 0.), sing );
        return proj( c1.o, tl );
};

// дуга радиусом R, исходящая из точки p и касательная к окружности c
CE Arc tangent_arc( const Circle& c, const Vec& p, double R, Sing sing = plus )
{
        Vec center = cross( Circle( p, R), Circle( c.o, R-c.R), sing );
        return Arc( center, R, p, cross( c, Line(c.o, center), -sing) );
}

// окруж. радиусом R, касательная к окружностям c1 и c2
CE Circle tangent_сircle( const Circle& c1, const Circle& c2, double R, Sing sing = plus )
{
        Vec center = cross( Circle( c1.o, R-c1.R), Circle( c2.o, R-c2.R), -sing );
        return Circle( center, R );
}
#pragma endregion

// очередная дуга из цепочки соприкасающихся окружностей
CE Arc chain_arc( const Vec& prev_arc_end, const Circle& current, const Circle& next, Sing sing = plus )
{
        return Arc( current, prev_arc_end, cross( next, Line( next.o, current.o), -sing) );
}

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

        CE Vec TE1( 1.,  0.00001); // задняя кромка верх
        CE Vec TE2( 1., -0.00001); // задняя кромка низ, зазор для соблюдения  
                                     // постулата Жуковского-Чаплыгина (Kutta condition)
        CE Circle  c_lead({0, s½}, s½ );
        CE Circle  с_top = tangent_сircle( c_lead, TE1, R );
        CE Circle  c_bottom( с_top.o, R-s );
        CE Segment s_end = tangent_segment( c_bottom, TE2 );
        CE Circle  c_l1( Vertical( lr1) ^ Circle( с_top.o, c_bottom.R+lr1), lr1 );
        CE Circle  c_l2 = tangent_сircle( с_top, c_l1, lr2 );

        CE Arc a_top = chain_arc( TE1      , с_top , c_l2     );
        CE Arc a_l2  = chain_arc( a_top.p2 , c_l2  , c_l1     );
        CE Arc a_l1  = chain_arc( a_l2.p2  , c_l1  , c_bottom, minus );
        CE Arc a_bottom( c_bottom, a_l1.p2, s_end.p1 );

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

        CE Vec TE1( 1.,  0.00001); // задняя кромка верх
        CE Vec TE2( 1., -0.00001); // задняя кромка низ, зазор для соблюдения  
                                     // постулата Жуковского-Чаплыгина (Kutta condition)
        CE Circle  c_lead({lr1, lr1}, lr1 );
        CE Circle  с_top = tangent_сircle( c_lead, TE1, R );
        CE Circle  c_bottom( с_top.o, R-s );
        CE Segment s_start = tangent_segment( c_lead, -c_bottom, minus );
        CE Segment s_end = tangent_segment( c_bottom, TE2 );

        CE Arc a_top = chain_arc( TE1, с_top, c_lead      );
        CE Arc a_lead  ( c_lead  , a_top.p2  , s_start.p1 );
        CE Arc a_bottom( c_bottom, s_start.p2, s_end.p1   );

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

        CE Vec TE1( 1.,  0.00001); // задняя кромка верх
        CE Vec TE2( 1., -0.00001); // задняя кромка низ, зазор для соблюдения  
                                     // постулата Жуковского-Чаплыгина (Kutta condition)
        CE Circle  c_l1   ( {lr1, lr1}, lr1 );
        CE Circle  c_trail( TE1       , -s  );
        CE Circle  c_bottom = tangent_сircle( -c_l1, c_trail, R-s );
        CE Circle  с_top( c_bottom.o, R );
        CE Circle  c_l2  = tangent_сircle( c_l1, с_top, lr2, minus );
        CE Segment s_end = tangent_segment( c_l1, TE2, minus );

        CE Arc a_top = chain_arc( TE1     , с_top, c_l2 );
        CE Arc a_l2  = chain_arc( a_top.p2, c_l2 , c_l1 );
        CE Arc a_l1( c_l1, a_l2.p2, s_end.p1 );

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

        CE Vec TE1( 1.,  0.00001); // задняя кромка верх
        CE Vec TE2( 1., -0.00001); // задняя кромка низ, зазор для соблюдения  
                                     // постулата Жуковского-Чаплыгина (Kutta condition)
        CE Circle  c_l1   ( {lr1, lr1}, lr1 );
        CE Circle  c_trail( TE1       , -s  );
        CE Circle  c_bottom = tangent_сircle( -c_l1, c_trail, R-s );
        CE Circle  с_top( c_bottom.o, R );
        CE Circle  c_l2  = tangent_сircle( c_l1, с_top, lr2, minus );
        CE Segment s_end0 = tangent_segment( c_l1, TE2, minus );

        CE Arc a_top = chain_arc( TE1     , с_top, c_l2 );
        CE Arc a_l2  = chain_arc( a_top.p2, c_l2 , c_l1 );
        CE Arc a_l1( c_l1, a_l2.p2, s_end0.p1 );

        CE Segment s_start( s_end0.p1, s_end0 ^ c_bottom );
        CE Segment s_end( s_end0 & c_bottom, TE2 );
        CE Arc a_bottom( c_bottom, s_start.p2, s_end.p1 );

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

        CE Vec TE1( 1.,  0.00001); // задняя кромка верх
        CE Vec TE2( 1., -0.00001); // задняя кромка низ, зазор для соблюдения  
                                     // постулата Жуковского-Чаплыгина (Kutta condition)
        CE Circle  c_l1   ( {lr1, lr1}, lr1 );
        CE Circle  c_trail( TE1       , -s  );
        CE Circle  c_bottom = tangent_сircle( -c_l1, c_trail, R-s );
        CE Circle  с_top( c_bottom.o, R );
        CE Circle  c_l2  = tangent_сircle( c_l1, с_top, lr2, minus );
        CE Segment s_end = tangent_segment( c_bottom, TE2 );

        CE Arc a_top = chain_arc( TE1     , с_top, c_l2     );
        CE Arc a_l2  = chain_arc( a_top.p2, c_l2 , c_l1     );
        CE Arc a_l1  = chain_arc( a_l2.p2 , c_l1 , c_bottom, minus );
        CE Arc a_bottom( c_bottom, a_l1.p2, s_end.p1 );

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

        CE Vec TE1( 1.,  0.00001); // задняя кромка верх
        CE Vec TE2( 1., -0.00001); // задняя кромка низ, между ними зазор для соблюдения постулата Жуковского-Чаплыгина (Kutta condition)
        Circle  c_le   ( {ler, ler}, -ler );
        Circle  c_trail( TE1       , -s   );
        Circle  c_bottom = tangent_сircle( c_le, c_trail, R-s );
        Circle  с_top( c_bottom.o, R );
        Circle  c_lef = tangent_сircle( -c_le, с_top, lef, minus );
        Segment s_end = tangent_segment( c_bottom, TE2 );

        Arc a_top = chain_arc( TE1     , с_top, c_lef    );
        Arc a_lef = chain_arc( a_top.p2, c_lef, c_le     );
        Arc a_le  = chain_arc( a_lef.p2, c_le , c_bottom, minus );
        Arc a_bottom( c_bottom, a_le.p2, s_end.p1 );

        std::cout << "PIPE " << D_abs << 'x' << s_abs << '-' << b_abs << " r" << ler_abs << " f" << lef_abs << '\n'
                  << std::setprecision(5) << std::fixed
                  << a_top << a_lef << a_le << a_bottom << s_end;

        return 0;
}
