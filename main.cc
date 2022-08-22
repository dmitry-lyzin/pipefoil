#include <iostream>
#include <iomanip>


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

//#define DBL_MATb_EGO_MIN DBL_MIN
#define DBL_MATb_EGO_MIN 0.0000001

double m( double x)
{
        return (fabs(x) < DBL_MATb_EGO_MIN) ? +0.0
                                            : x
                                            ;
}

enum Sing { minus = -1, plus = 1 };

double square( double a)
{
        return a * a;
}

struct Point
{
        double x, y;
        Point(                     ): x(0.), y(0.) {};
        Point( double x_, double y_): x(x_), y(y_) {};
friend  double distance( const Point& p1, const Point& p2)
        {
                return sqrt( square( p2.x - p1.x) + square( p2.y - p1.y));
        }
};

struct Circle
{
        Point center;
        double r;

        Circle( double r_, const Point& p1, const Point& p2, Sing case_)
        : r( r_)
        {
                double d = case_ * distance( p1, p2);
                //double h = sqrt( square( r) - square( d/2.));
                double h_div_d = sqrt( square( r) - square( d/2.)) / d;

                center.x = (p1.x + p2.x)/2. + (p2.y - p1.y) * h_div_d;
                center.y = (p1.y + p2.y)/2. - (p2.x - p1.x) * h_div_d;
        };

        Circle& operator *= ( double scale)
        {
                r *= scale;
                return *this;
        };
};

struct Arc: public Circle
{
        Point p1;
        Point p2;
        Sing direction;

        Arc( double r_, const Point& p1_, const Point& p2_, Sing case_, Sing direction_)
                : Circle( r_, p1_, p2_, case_)
                , p1( p1_), p2( p2_), direction( direction_)
        {};

        Arc& operator *= ( double scale)
        {
                *((Circle *)this) *= scale;
                p1.x = (p1.x - center.x) * scale + center.x;
                p1.y = (p1.y - center.y) * scale + center.y;
                p2.x = (p2.x - center.x) * scale + center.x;
                p2.y = (p2.y - center.y) * scale + center.y;
                return *this;
        };

friend  std::ostream& operator<<( std::ostream &os, const Arc& arc )
        {
                const int c = 40;
                double  x =  arc.p1.x;
                double dx = (arc.p2.x - x) / c;
                if( arc.direction == minus)
                {
                        x = arc.p2.x;
                        dx = -dx;
                }
                double r2 = square( arc.r);
                os << std::setprecision(5) << std::fixed;
                for( int i = 0; i <= c; i++)
                {
                        double y = arc.center.y + sqrt( r2 - square( x - arc.center.x) );
                        os << std::setw(8) << m(x) << std::setw(12) << m(y) << '\n';
                        x += dx;
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

    Arc a( D/2., Point( 1., 0.), Point( 0., 0.), minus, plus );
    Arc b = a;
    b *= (1. - s_abs*2./D_abs);
    b.direction = minus;
    std::cout << a << b;

    std::cout << " 1.00000     0.00000\n";

    return 0;
}

int main1( int argc, const char *argv[])
{
    //uint maxnum = argc < 2 ? 250000
    //                       : atoi( argv[1]);

    double b_abs =  70.0; // хорда профиля
    double D_abs = 160.0; // диаметр трубы
    double s_abs =   4.7; // толщина стенки трубы

    printf("PIPE %.0fx%.1f-%.0f\n", D_abs, s_abs, b_abs);

    double D = D_abs / b_abs; // относительный диаметр трубы
    double s = s_abs / b_abs; // относительная толщина стенки трубы

    // расчет координат центра
    double R2 = D*D/4.0; // квадрат радиуса трубы
    double x0 = 0.5;
    double y0 = - sqrt( R2 - x0*x0);

    for( double x = 1.0; x >= -DBL_MATb_EGO_MIN; x -= 0.02)
    {
        double y = y0 + sqrt( R2 - (x-x0)*(x-x0) );
        printf("% .5f    % .5f\n", x, y );
    }

    double k = 1. - s_abs*2/D_abs;
    //printf("%f\n", k);

    R2 = R2*k*k; // квадрат внутреннего радиуса трубы
    for( double x = 0.5*(1-k); x <= 0.5*(1+k); x += 0.02*k)
    {
        double y = y0 + sqrt( R2 - (x-x0)*(x-x0) );
        printf("% .5f    % .5f\n", x, y);
    }

    printf("% .5f    % .5f\n", 1.0, 0.0);

    return 0;
}
