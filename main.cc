#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
//#include <limits>   

#define CE constexpr
CE double π = 3.14159265358979323846;
//enum Sing { minus = -1, plus = 1 };
CE int plus  =  1;
CE int minus = -1;
typedef int Sing;

namespace conexpr
{
        double CE sqrt_Newton_Raphson(double x, double curr, double prev)
        {
                return curr == prev
                        ? curr
                        : sqrt_Newton_Raphson(x, 0.5 * (curr + x / curr), curr);
        }

        /*
        * Constexpr version of the square root
        * Return value:
        *   - For a finite and non-negative value of "x", returns an approximation for the square root of "x"
        *   - Otherwise, returns NaN
        */
        double CE sqrt( double x)
        {
                return x >= 0 && x < std::numeric_limits<double>::infinity()
                        ? sqrt_Newton_Raphson(x, x, 0)
                        : std::numeric_limits<double>::quiet_NaN();
        }
}

CE double root_of_quadratic_equation( double A, double B, double C, Sing sing = plus )
{
        return (-B + sing * conexpr::sqrt( B*B - 4.*A*C ))/(2.*A);
}

CE double ²( double a)
{
        return a * a;
}

struct Point
{
        double x, y;
CE      Point(                     ): x(0.), y(0.) {};
CE      Point( double x_, double y_): x(x_), y(y_) {};
CE 
friend  double distance( const Point& p1, const Point& p2)
        {
                return conexpr::sqrt( ²( p2.x - p1.x) + ²( p2.y - p1.y));
        }
};

struct Line
{
        double a;
        double b;
        double c;

//CE      static double k ( const Point& p1, const Point& p2 ) { return (p2.y - p1.y)/(p2.x - p1.x); }
//CE      static double y0( const Point& p1, const Point& p2 ) { return p1.y - k( p1, p2 ) * p1.x;   }

        // ах + by + с = 0
CE      Line( double a_, double b_, double c_  ): a( a_), b( b_), c( c_) {};
        // y = kx + y0
CE      Line( double k, double y0              ): a(-k ), b( 1 ), c(-y0) {};
        // через точки p1 и p2
CE      Line( const Point& p1, const Point& p2 )
        : a( p1.y - p2.y           )
        , b( p2.x - p1.x           )
        , c( p1.x*p2.y - p2.x*p1.y )
        {};

//CE      double k()  const { return _k;   }
//CE      double y0() const { return _y0;  }

friend  std::ostream& operator<<( std::ostream &os, const Line& obj )
        {
                os << std::setprecision(5) << std::fixed
                   << std::setw(8) << 0.0 << std::setw(12) << (-obj.c/obj.b) << '\n'
                   << std::setw(8) << (-obj.c/obj.a) << std::setw(12) << 0.0 << '\n';

                return os;
        };
};
struct Vertical: public Line
{
        // Вертикальная линия пересекающая ось X в точке x0
CE      Vertical( double x0 ): Line( -1., 0., x0) {};
};
struct Horizontal: public Line
{
        // Вертикальная линия пересекающая ось Y в точке y0
CE      Horizontal( double y0 ): Line( 0., 1., -y0) {};
};

struct Segment: public Line
{
        Point p1;
        Point p2;

CE      Segment( const Point& p1_, const Point& p2_ )
        : Line( p1_, p2_ )
        , p1( p1_), p2( p2_)
        {};

        // обмен концов (рисоваться будет в другую сторону)
        void revers() { std::swap( p1, p2 ); }
        // обмен концов (рисоваться будет в другую сторону)
CE      Segment operator - () const { return Segment( this->p2, this->p1 ); }

friend  std::ostream& operator<<( std::ostream &os, const Segment& obj )
        {
                os << std::setprecision(5) << std::fixed
                   //<< std::setw(8) << obj.p1.x << std::setw(12) << obj.p1.y << '\n'
                   << std::setw(8) << obj.p2.x << std::setw(12) << obj.p2.y << '\n';

                return os;
        };
};

struct Circle: public Point
{
        double R;

CE static
        Point center( double R, const Point& p1, const Point& p2, Sing case_)
        {
                double d = case_ * distance( p1, p2);
                //double h = sqrt( ²(R) - ²(d/2.));
                double h_div_d = conexpr::sqrt( ²(R) - ²(d/2.)) / d;

                return Point( (p1.x + p2.x)/2. + (p2.y - p1.y) * h_div_d
                            , (p1.y + p2.y)/2. - (p2.x - p1.x) * h_div_d
                            );                      
        };

CE      Circle( double R_, double center_x, double center_y )
        : R( R_), Point( center_x, center_y )
        {};
CE      Circle( double R_, const Point& center )
        : R( R_), Point( center )
        {};
CE      Circle( double R_, const Point& p1, const Point& p2, Sing case_)
        : R( R_), Point( center( R_, p1, p2, case_))
        {};

CE      Circle& operator *= ( double scale)
        {
                R *= scale;
                return *this;
        };
};

struct Arc: public Circle
{
        Point p1;
        Point p2;
        bool sector;
        
CE      Arc( const Circle& circle, const Point& start, const Point& finish )
        : Circle( circle)
        , p1( start), p2( finish)
        , sector( false)
        {};

CE      Arc( double R_, const Point& center, const Point& start, const Point& finish )
        : Circle( R_, center )
        , p1( start), p2( finish)
        , sector( false)
        {};

CE      Arc( const Arc& obj, bool sector_ )
        : Circle( obj)
        , p1( obj.p1), p2( obj.p2)
        , sector( sector_ )
        {};

        // сменить сектор откружности
CE      Arc operator ~ () const { return Arc( *this, !(this->sector)  ); }

        // сменить напраление отрисовки
        void revers() { std::swap( p1, p2 ); }

        Arc& operator *= ( double scale)
        {
                //*((Circle *)this) *= scale;
                *static_cast< Circle *>(this) *= scale;
                p1.x = (p1.x - x) * scale + x;
                p1.y = (p1.y - y) * scale + y;
                p2.x = (p2.x - x) * scale + x;
                p2.y = (p2.y - y) * scale + y;
                return *this;
        };

static  double angle( double Δx, double Δy )
        {
                if( Δx > +0. ) return atan( Δy / Δx );
                if( Δx < -0. ) return atan( Δy / Δx ) + π;

                if( Δy > +0. ) return π/2;
                if( Δy < -0. ) return π*3/2;

                return NAN;
        }
        
friend  std::ostream& operator<<( std::ostream &os, const Arc& arc )
        {
                double a1 = angle( arc.p1.x - arc.x, arc.p1.y - arc.y );
                double a2 = angle( arc.p2.x - arc.x, arc.p2.y - arc.y );

                if( arc.sector )
                {
                        if( a1 > a2)
                                a2 += 2.* π;
                        else
                                a1 += 2.* π;
                }

                // вычисляем кол. сегментов
                //CE int n = 20;
                double n = abs( round( (a2 - a1) / (2.* π) * 160));
                double Δa = (a2 - a1) / n;
                os << std::setprecision(5) << std::fixed;
                //for( int i = 0; i <= n; i++)
                for( int i = 0; i < n; i++)
                {
                        double x = arc.x + arc.R * cos(a1);
                        double y = arc.y + arc.R * sin(a1);
                        //os << std::setprecision(0) << std::fixed << a1*180/π << "*   " << std::setprecision(5) << std::fixed;
                        os << std::setw(8) << x << std::setw(12) << y << '\n';
                        a1 += Δa;
                }

                return os;
        };
};

#pragma region // функции поиска пересечений
CE Point cross( const Circle& c, const Line& l, Sing sing = plus )
{
        /*
        ax + by + c = 0
        (x - x₀)² + (y - y₀)² = R²

        -(ax + c)/b = y
        (1 + (a/b)²)x² + 2(a(c/b + y₀)/b - x₀)x + ((c/b + y₀)² + x₀² - R²) = 0

        -(by + c)/a = x
        (1 + (b/a)²)y² + 2(b(c/a + x₀)/a - y₀)y + ((c/a + x₀)² + y₀² - R²) = 0
        */

        if( l.a != 0. )
        {
                double y = root_of_quadratic_equation
                        ( 1. + ²(l.b/l.a)
                        , 2.*(l.b*(l.c/l.a + c.x)/l.a - c.y)
                        , ²(l.c/l.a + c.x) + ²(c.y) - ²(c.R)
                        , sing
                        );
                return Point( -(l.b*y + l.c)/l.a, y );
        }
        
        //(b(c/a + X)/a - Y)² = (1 + (b/a)²)((c/a + X)² + Y² - R²)

        double x = root_of_quadratic_equation
                ( 1. + ²(l.a/l.b)
                , 2.*(l.a*(l.c/l.b + c.y)/l.b - c.x)
                , ²(l.c/l.b + c.y) + ²(c.x) - ²(c.R)
                , sing
                );
        return Point( x, -(l.a*x + l.c)/l.b );
}

CE Point cross( const Circle& c1, const Circle& c2, Sing sing = plus )
{

        Line radical_line( 2. * (c2.x - c1.x)
                         , 2. * (c2.y - c1.y)
                         , -(²(c1.R) - ²(c1.x) - ²(c1.y) - ²(c2.R) + ²(c2.x) + ²(c2.y))
                         );

        return cross( c1, radical_line, sing ); 
}

// сахарок синтаксический

// первая точка пересечения окружности и линии
CE Point operator & ( const Circle& c,  const Line& l    ) { return cross( c, l, plus   ); }
// первая точка пересечения линии и окружности
CE Point operator & ( const Line& l,    const Circle& c  ) { return cross( c, l, plus   ); }
// вторая точка пересечения окружности и линии
CE Point operator ^ ( const Circle& c,  const Line& l    ) { return cross( c, l, minus  ); }
// вторая точка пересечения линии и окружности
CE Point operator ^ ( const Line& l,    const Circle& c  ) { return cross( c, l, minus  ); }
// первая точка пересечения окружности и окружности
CE Point operator & ( const Circle& c1, const Circle& c2 ) { return cross( c1, c2, plus ); }
// вторая точка пересечения окружности и окружности
CE Point operator ^ ( const Circle& c1, const Circle& c2 ) { return cross( c1, c2, minus); }
#pragma endregion

// очередная дуга из цепочки соприкасающихся окружностей
CE Arc chain_arc( const Point& prev_arc_end, const Circle& current, const Circle& next )
{
        return Arc( current, prev_arc_end, cross( next, Line( next, current), plus) );
}

// поиск точки касания прямой исходящей из точки p и окружности c
CE Point tangent_point( const Circle& c, const Point& p, Sing sing = plus )
{
/*
https://www.cyberforum.ru/geometry/thread605358.html?ysclid=l77andw88d999612641 

x0=a,
y0=b - центр окр.
y и x- координаты ИСХОДНОЙ точки

(2ax-a²+R²-x²)k²+(2ab-2ay+2yx-2xb)k+(R²-b²-y²+2by) = 0
(2ax-a²+R²-x²)k²+2(ab-ay+yx-xb)k+(R²-b²-y²+2by) = 0

// решаем квадратное уравнение
A = (2ax-a²+R²-x²)
B = 2(ab-ay+yx-xb)
C = (R²-b²-y²+2by)
D = B² — 4AC

k = (-B±√D)/(2A)
нашли k

n=(y-kx)
X, Y - точка касания
X = (a+k(b-n))/(1+k²)
Y = kX-n
--------------------

/*
(x₀, y₀) исходная т.
(X, Y) центр окр.

R² = (aX + bY + c)²/(a² + b²)
ax₀ + by₀ + c = 0

a²(R² - (X - x₀)²) + b²(R² - (Y - y₀)²) = 2ab(Y - y₀)(X - x₀)
-(ax₀ + by₀) = c

c = 0
b = (-ax₀/y₀)
(R² - Y²)(x₀/y₀)² + 2XY(x₀/y₀) + R² - X² = 0

(-b/a) = x₀/y₀

(b/a)²(R² - Y²) - (b/a)2XY + R² - X² = 0





c = 1
b = (-(ax₀ + 1)/y₀)

a²(R²(1 + x₀²/y₀²) - (X - Yx₀/y₀)² ) + 
2a( x₀R²/y₀² - (1 - Y/y₀)(X - Yx₀/y₀)) + 
R²/y₀² - (1 - Y/y₀)² = 0
*/
        double Δx = c.x - p.x;
        double Δy = c.y - p.y;
        double R² = ²(c.R);

        // решаем квадратное уравнение
        double k = root_of_quadratic_equation( R²-²(Δx), 2.*Δx*Δy, R²-²(Δy), sing );

        // X, Y - точка касания
        double y0 = p.y - k*p.x;
        double X  = (c.x + k*(c.y - y0)) / (1.+²(k));
        double Y  = k*X + y0;

        return Point( X, Y );
};

// отрезок, исходящей из точки p и касающийся окружности c
CE Segment tangent_segment( const Circle& c, const Point& p, Sing sing = plus )
{
        return Segment( p, tangent_point( c, p, sing ) );
};

// дуга радиусом R, исходящая из точки p и касательная к окружности c
CE Arc tangent_arc( const Circle& c, const Point& p, double R, Sing sing = plus )
{
        Point center = cross( Circle( R, p), Circle( R - c.R, c), sing );
        return Arc( R, center, p, cross( c, Line(c, center), -sing) );
}

// окруж. радиусом R, через точку p и касательная к окружности c
CE Circle tangent_сircle( const Circle& c, const Point& p, double R, Sing sing = plus )
{
        Point center = cross( Circle( R, p), Circle( R - c.R, c), sing );
        return Circle( R, center );
}

// окруж. радиусом R, касательная к окружностям c1 и c2
CE Circle tangent_сircle( const Circle& c1, const Circle& c2, double R, Sing sing = plus )
{
        Point center = cross( Circle( c1.R-R, c1), Circle( c2.R-R, c2), sing );
        return Circle( R, center );
}

int main( int argc, const char *argv[])
{
        /*
        double b_abs   = atof( argv[1]); //  50.0; // хорда профиля
        double D_abs   = atof( argv[2]); // 160.0; // диаметр трубы
        double s_abs   = atof( argv[3]); //   4.7; // толщина стенки трубы
        double lr1_abs = atof( argv[4]); //   1.5; // радиус скругл. носка 1
        double lr2_abs = atof( argv[5]); //  18.0; // радиус скругл. носка 2
        */
        
        CE double b_abs   =  40.0; // хорда профиля
        CE double D_abs   = 110.0; // диаметр трубы
        CE double s_abs   =   2.7; // толщина стенки трубы
        CE double lr1_abs =   1.0; // радиус скругл. носка 1
        CE double lr2_abs =  18.0; // радиус скругл. носка 2
        
        CE double R   = D_abs   / b_abs / 2.; // относительный внешний радиус трубы
        CE double s   = s_abs   / b_abs     ; // относительная толщина стенки трубы
        CE double s2  = s / 2.              ; // относительная полутолщина стенки трубы
        CE double lr1 = lr1_abs / b_abs     ; // относительный радиус скругл. носка 1
        CE double lr2 = lr2_abs / b_abs     ; // относительный радиус скругл. носка 2

        CE Point TE1( 1.,  0.001); // задняя кромка верх
        CE Point TE2( 1., -0.001); // задняя кромка низ, разрывчик для соблюдения  
                                   // постулата Жуковского-Чаплыгина (Kutta condition)
        CE Circle  c_lead( s2, Point( 0, s2) );
        CE Circle  с_top = tangent_сircle( c_lead, TE1, R, minus );
        CE Circle  c_bottom( R-s, с_top );
        CE Segment s_end = -tangent_segment( c_bottom, TE2 );
        CE Circle  c_l1( lr1, Vertical( lr1) & Circle( c_bottom.R+lr1, с_top) );
        CE Circle  c_l2 = tangent_сircle( с_top, c_l1, lr2 );

        CE Arc a_top = chain_arc( TE1      , с_top , c_l2     );
        CE Arc a_l2  = chain_arc( a_top.p2 , c_l2  , c_l1     );
        CE Arc a_l1  = chain_arc( a_l2.p2  , c_l1  , c_bottom );
        CE Arc a_bottom( c_bottom, a_l1.p2, s_end.p1 );

        std::cout << "PIPE " << D_abs << 'x' << s_abs << '-' << b_abs << '\n'
                  << a_top << a_l2 << ~a_l1 << a_bottom << s_end;

        return 0;
}
