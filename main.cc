#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
//#include <limits>   

#define CE     constexpr
#define statCE constexpr static

CE double π = 3.14159265358979323846;

class Sing
{
        int value;
CE      Sing( int x): value(x) {};
public:
CE      Sing(): value(1) {};
CE      operator int   () const { return  value; };
CE      Sing operator- () const { return -value; };
statCE  Sing plus      ()       { return      1; };
statCE  Sing minus     ()       { return     -1; };
};

CE Sing plus  = Sing::plus();
CE Sing minus = Sing::minus();

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
        double a, b, c;

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

        // TODO стремный конст-ор, как бы его спрятать
CE      Segment( const Line& l, const Point& p1_, const Point& p2_ )
        : Line( l )
        , p1( p1_), p2( p2_)
        {};

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

statCE  Point center( double R, const Point& p1, const Point& p2, Sing case_)
        {
                double d = case_ * distance( p1, p2);
                //double h = sqrt( ²(R) - ²(d/2.));
                double h_div_d = conexpr::sqrt( ²(R) - ²(d/2.)) / d;

                return Point( (p1.x + p2.x)/2. + (p2.y - p1.y) * h_div_d
                            , (p1.y + p2.y)/2. - (p2.x - p1.x) * h_div_d
                            );                      
        };

CE      Circle( const Point& center, double R_ )
        : R( R_), Point( center )
        {};
CE      Circle( const Point& center            )
        : R( 0.), Point( center )
        {};
CE      Circle( const Point& p1, const Point& p2, double R_, Sing case_)
        : R( R_), Point( center( R_, p1, p2, case_))
        {};
};

struct Arc: public Circle
{
        Point p1;
        Point p2;
        
CE      Arc( const Circle& circle, const Point& start, const Point& finish )
        : Circle( circle)
        , p1( start), p2( finish)
        {};

CE      Arc( const Point& center, double R_, const Point& start, const Point& finish )
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

statCE  double angle( double Δx, double Δy )
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

                double R = arc.R;
                if( R < 0. )
                {
                        R = -R;
                        if( a1 > a2) a2 += 2.* π;
                        else         a1 += 2.* π;
                }

                // вычисляем кол. сегментов
                //CE int n = 20;
                double n = abs( round( (a2 - a1) / (2.* π) * 160));
                double Δa = (a2 - a1) / n;
                os << std::setprecision(5) << std::fixed;
                for( int i = 0; i < n; i++)
                {
                        double x = arc.x + R * cos(a1);
                        double y = arc.y + R * sin(a1);
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

#pragma region // функции поиска касательных и перпендикуляров

// касательная к окружностям c1 и c2
CE Line tangent_line( const Circle& c1, const Circle& c2, Sing sing = plus )
{
        double Δx = c2.x - c1.x;
        double Δy = c2.y - c1.y;
        double Δr = c2.R - c1.R;
        double d = conexpr::sqrt(Δx*Δx + Δy*Δy);
        double X = Δx/d;
        double Y = Δy/d;

        double cosθ = Δr/d; // θ - угол между линией центров и касательной
        double sinθ = conexpr::sqrt(1 - cosθ*cosθ);

        double a = X*cosθ -sing* Y*sinθ;
        double b = Y*cosθ +sing* X*sinθ;
        double c = c1.R - (a*c1.x + b*c1.y);

        return Line( a, b, c);
}

// основание перпендикуляра на прямую l из точки p
CE Point foot_of_perpendicular( const Line& l, const Point& p )
{
        double k = - (l.a*p.x + l.b*p.y + l.c) / (l.a*l.a + l.b*l.b);
        return { k*l.a + p.x, k*l.b + p.y };
}

// опустить перпендикуляр из точки p на прямую l
CE Segment perpendicular( const Line& l, const Point& p )
{
        return { p, foot_of_perpendicular( l, p ) };
}

// касательная к окружностям c1 и c2
CE Segment tangent_segment( const Circle& c1, const Circle& c2, Sing sing = plus )
{
        Line tl = tangent_line( c1, c2, sing );
        return Segment( tl, foot_of_perpendicular( tl, c1), foot_of_perpendicular( tl, c2) );
}

// поиск точки касания прямой исходящей из точки p и окружности c
CE Point tangent_point( const Circle& c1, const Point& p, Sing sing = plus )
{
        Line tl = tangent_line( c1, Circle( p, 0.), sing );
        return foot_of_perpendicular( tl, c1);
};

// дуга радиусом R, исходящая из точки p и касательная к окружности c
CE Arc tangent_arc( const Circle& c, const Point& p, double R, Sing sing = plus )
{
        Point center = cross( Circle( p, R), Circle( c, R-c.R), sing );
        return Arc( center, R, p, cross( c, Line(c, center), -sing) );
}

// окруж. радиусом R, касательная к окружностям c1 и c2
CE Circle tangent_сircle( const Circle& c1, const Circle& c2, double R, Sing sing = plus )
{
        Point center = cross( Circle( c1, R-c1.R), Circle( c2, R-c2.R), sing );
        return Circle( center, R );
}
#pragma endregion

// очередная дуга из цепочки соприкасающихся окружностей
CE Arc chain_arc( const Point& prev_arc_end, const Circle& current, const Circle& next )
{
        return Arc( current, prev_arc_end, cross( next, Line( next, current), plus) );
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
        
        CE double R   = D_abs/2./ b_abs; // относительный внешний радиус трубы
        CE double s   = s_abs   / b_abs; // относительная толщина стенки трубы
        CE double s½  = s / 2.         ; // относительная полутолщина стенки трубы
        CE double lr1 = lr1_abs / b_abs; // относительный радиус скругл. носка 1
        CE double lr2 = lr2_abs / b_abs; // относительный радиус скругл. носка 2

        CE Point TE1( 1.,  0.00001); // задняя кромка верх
        CE Point TE2( 1., -0.00001); // задняя кромка низ, зазор для соблюдения  
                                     // постулата Жуковского-Чаплыгина (Kutta condition)
        CE Circle  c_lead({0, s½}, s½ );
        CE Circle  с_top = tangent_сircle( c_lead, TE1, R, minus );
        CE Circle  c_bottom( с_top, R-s );
        CE Segment s_end = tangent_segment( c_bottom, TE2, minus );
        CE Circle  c_l1( Vertical( lr1) & Circle( с_top, c_bottom.R+lr1), lr1 );
        CE Circle  c_l2 = tangent_сircle( с_top, c_l1, lr2 );

        CE Arc a_top = chain_arc( TE1      , с_top , c_l2     );
        CE Arc a_l2  = chain_arc( a_top.p2 , c_l2  , c_l1     );
        CE Arc a_l1  = chain_arc( a_l2.p2  , c_l1  , c_bottom );
        CE Arc a_bottom( c_bottom, a_l1.p2, s_end.p1 );

        std::cout << "PIPE " << D_abs << 'x' << s_abs << '-' << b_abs << '\n'
                  << a_top << a_l2 << -a_l1 << a_bottom << s_end;

        return 0;
}
