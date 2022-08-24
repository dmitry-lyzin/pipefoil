#include <iostream>
#include <iomanip>
#include <cmath>

//#define DBL_MATb_EGO_MIN DBL_MIN
#define DBL_MATb_EGO_MIN 0.0000001
#define constex constexpr

double m( double x)
{
        return (fabs(x) < DBL_MATb_EGO_MIN) ? +0.0
                                            : x
                                            ;
}

enum Sing { minus = -1, plus = 1 };

double ²( double a)
{
        return a * a;
}

struct Point
{
        double x, y;
constex Point(                     ): x(0.), y(0.) {};
constex Point( double x_, double y_): x(x_), y(y_) {};
friend  double distance( const Point& p1, const Point& p2)
        {
                return sqrt( ²( p2.x - p1.x) + ²( p2.y - p1.y));
        }
};

struct Circle
{
        Point O;
        double R;

        Circle( double R_, const Point& p1, const Point& p2, Sing case_)
        : R( R_)
        {
                double d = case_ * distance( p1, p2);
                //double h = sqrt( ²(R) - ²(d/2.));
                double h_div_d = sqrt( ²(R) - ²(d/2.)) / d;

                O.x = (p1.x + p2.x)/2. + (p2.y - p1.y) * h_div_d;
                O.y = (p1.y + p2.y)/2. - (p2.x - p1.x) * h_div_d;
        };

        Circle& operator *= ( double scale)
        {
                R *= scale;
                return *this;
        };
};


// поиск точки касания прямой исходящей из точки src_p и окружности circle
Point tangent_point( const Circle& circle, const Point& src_p, Sing sing = plus )
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
        double a = circle.O.x;
        double b = circle.O.y;
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

        double k = (-B + sing*sqrt(D))/(2.*A);

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

constex Line( const Point& p1_, const Point& p2_, Sing direction_ = plus )
                : p1( p1_), p2( p2_), direction( direction_)
        {};
        Line( const Circle& tangent_circle, const Point& p1_, Sing case_, Sing direction_ = plus )
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

        Arc( const Circle& circle_, const Point& p1_, const Point& p2_, Sing case_, Sing direction_ = plus )
                : Circle( circle_)
                , p1( p1_), p2( p2_), direction( direction_)
        {};

        Arc( double r_, const Point& p1_, const Point& p2_, Sing case_, Sing direction_ = plus )
                : Circle( r_, p1_, p2_, case_)
                , p1( p1_), p2( p2_), direction( direction_)
        {};

        Arc& operator *= ( double scale)
        {
                //*((Circle *)this) *= scale;
                *static_cast< Circle *>(this) *= scale;
                p1.x = (p1.x - O.x) * scale + O.x;
                p1.y = (p1.y - O.y) * scale + O.y;
                p2.x = (p2.x - O.x) * scale + O.x;
                p2.y = (p2.y - O.y) * scale + O.y;
                return *this;
        };

friend  std::ostream& operator<<( std::ostream &os, const Arc& arc )
        {
                constexpr int n = 20;
 
                double a1 = angle( arc.p1.x - arc.O.x, arc.p1.y - arc.O.y );
                double a2 = angle( arc.p2.x - arc.O.x, arc.p2.y - arc.O.y );
                if( arc.direction == minus)
                        std::swap( a1, a2 );

                double Δa = (a2 - a1) / n;
                os << std::setprecision(5) << std::fixed;
                for( int i = 0; i <= n; i++)
                {
                        double x = arc.O.x + arc.R * cos(a1);
                        double y = arc.O.y + arc.R * sin(a1);
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
