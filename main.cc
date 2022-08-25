#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>   

namespace conexpr
{
        double constexpr sqrt_Newton_Raphson(double x, double curr, double prev)
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
        double constexpr sqrt( double x)
        {
                return x >= 0 && x < std::numeric_limits<double>::infinity()
                        ? sqrt_Newton_Raphson(x, x, 0)
                        : std::numeric_limits<double>::quiet_NaN();
        }
}

enum Sing { minus = -1, plus = 1 };

constexpr double ²( double a)
{
        return a * a;
}

struct Point
{
        double x, y;
constexpr Point(                     ): x(0.), y(0.) {};
constexpr Point( double x_, double y_): x(x_), y(y_) {};
constexpr 
friend  double distance( const Point& p1, const Point& p2)
        {
                return conexpr::sqrt( ²( p2.x - p1.x) + ²( p2.y - p1.y));
        }
};

struct Circle: public Point
{
        double R;

constexpr static
        Point circle_center( double R, const Point& p1, const Point& p2, Sing case_)
        {
                double d = case_ * distance( p1, p2);
                //double h = sqrt( ²(R) - ²(d/2.));
                double h_div_d = conexpr::sqrt( ²(R) - ²(d/2.)) / d;

                return Point( (p1.x + p2.x)/2. + (p2.y - p1.y) * h_div_d
                            , (p1.y + p2.y)/2. - (p2.x - p1.x) * h_div_d
                            );                      
        };

constexpr Circle( double R_, const Point& p1, const Point& p2, Sing case_)
        : Point( circle_center( R_, p1, p2, case_)), R( R_)
        {};

constexpr Circle& operator *= ( double scale)
        {
                R *= scale;
                return *this;
        };
};


// поиск точки касания прямой исходящей из точки src_p и окружности circle
constexpr Point tangent_point( const Circle& circle, const Point& src_p, Sing sing = plus )
{
/*
https://www.cyberforum.ru/geometry/thread605358.html?ysclid=l77andw88d999612641 

x0=a,
y0=b - центр окр.
y и x- координаты ИСХОДНОЙ точки

(2ax-a²+R²-x²)k²+(2ab-2ay+2yx-2xb)k+(R²-b²-y²+2by) = 0
(2ax-a²+R²-x²)k²+2(ab-ay+yx-xb)k+(R²-b²-y²+2by) = 0

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
*/

        // центр окр.
        double a = circle.x;
        double b = circle.y;
        double R²= ²(circle.R);
        //y и x - координаты ИСХОДНОЙ точки
        double x = src_p.x;
        double y = src_p.y;

        // (2ax-a²+R²-x²)k²+(2ab-2ay+2yx-2xb)k+(R²-b²-y²+2by) = 0
        // (2ax-a²+R²-x²)k²+2(ab-ay+yx-xb)k+(R²-b²-y²+2by) = 0

        // решаем квадратное уравнение
        double A = R² - ²(a-x);
        double B = 2.*(a*b-a*y+y*x-x*b);
        double C = R² - ²(b-y);
        double D = ²(B) - 4.*A*C;

        double k = (-B + sing * conexpr::sqrt(D))/(2.*A);

        // X, Y - точка касания
        double n = y - k*x;
        double X = (a+k*(b-n))/(1.+²(k));
        double Y = k*X + n;

        return Point( X, Y );
};

struct Line
{
        Point p1;
        Point p2;
        Sing direction;

constexpr Line( const Point& p1_, const Point& p2_, Sing direction_ = plus )
        : p1( p1_), p2( p2_), direction( direction_)
        {};
constexpr Line( const Circle& tangent_circle, const Point& p1_, Sing case_, Sing direction_ = plus )
        : p1( p1_), p2( tangent_point( tangent_circle, p1_, case_ )), direction( direction_)
        {};
};

double angle( double Δx, double Δy )
{
        constexpr double π = 3.14159265358979323846;

        if( Δx > +0. ) return atan( Δy / Δx );
        if( Δx < -0. ) return atan( Δy / Δx ) + π;

        if( Δy > +0. ) return π / 2.;
        if( Δy < -0. ) return 3./2. * π;

        return NAN;
}

struct Arc: public Circle
{
        Point p1;
        Point p2;
        Sing direction;

constexpr Arc( const Circle& circle_, const Point& p1_, const Point& p2_, Sing case_, Sing direction_ = plus )
                : Circle( circle_)
                , p1( p1_), p2( p2_), direction( direction_)
        {};

constexpr Arc( double r_, const Point& p1_, const Point& p2_, Sing case_, Sing direction_ = plus )
                : Circle( r_, p1_, p2_, case_)
                , p1( p1_), p2( p2_), direction( direction_)
        {};

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

friend  std::ostream& operator<<( std::ostream &os, const Arc& arc )
        {
                constexpr int n = 20;
 
                double a1 = angle( arc.p1.x - arc.x, arc.p1.y - arc.y );
                double a2 = angle( arc.p2.x - arc.x, arc.p2.y - arc.y );
                if( arc.direction == minus)
                        std::swap( a1, a2 );

                double Δa = (a2 - a1) / n;
                os << std::setprecision(5) << std::fixed;
                for( int i = 0; i <= n; i++)
                {
                        double x = arc.x + arc.R * cos(a1);
                        double y = arc.y + arc.R * sin(a1);
                        os << std::setw(8) << x << std::setw(12) << y << '\n';
                        a1 += Δa;
                }

                return os;
        };
};

int main( int argc, const char *argv[])
{
    //uint maxnum = argc < 2 ? 250000
    //                       : atoi( argv[1]);

    double b_abs =  70.0; // хорда профиля
    double D_abs = 160.0; // диаметр трубы
    double s_abs =   4.7; // толщина стенки трубы

    std::cout << "PIPE " << D_abs << 'x' << s_abs << '-' << b_abs << '\n';

    double D = D_abs / b_abs; // относительный диаметр трубы
    double s = s_abs / b_abs; // относительная толщина стенки трубы

    Arc a( D/2., Point( 1., 0.), Point( 0., 0.), minus );
    Arc b = a;
    b *= (1. - s_abs*2./D_abs);
    b.direction = minus;
    b.p1 = tangent_point( b, Point( 1., 0.), plus );
    std::cout << a << b;

    std::cout << " 1.00000     0.00000\n";

    return 0;
}
