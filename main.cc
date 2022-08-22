#include <iostream>
#include <iomanip>
#include <cmath>


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

        Arc( double r_, const Point& p1_, const Point& p2_, Sing case_, Sing direction_ = plus )
                : Circle( r_, p1_, p2_, case_)
                , p1( p1_), p2( p2_), direction( direction_)
        {};

        Arc& operator *= ( double scale)
        {
                //*((Circle *)this) *= scale;
                *static_cast< Circle *>(this) *= scale;
                p1.x = (p1.x - center.x) * scale + center.x;
                p1.y = (p1.y - center.y) * scale + center.y;
                p2.x = (p2.x - center.x) * scale + center.x;
                p2.y = (p2.y - center.y) * scale + center.y;
                return *this;
        };

friend  std::ostream& operator<<( std::ostream &os, const Arc& arc )
        {
                constexpr int n = 20;
 
                double a1 = angle( arc.p1.x - arc.center.x, arc.p1.y - arc.center.y );
                double a2 = angle( arc.p2.x - arc.center.x, arc.p2.y - arc.center.y );
                if( arc.direction == minus)
                        std::swap( a1, a2 );

                double Δa = (a2 - a1) / n;
                os << std::setprecision(5) << std::fixed;
                for( int i = 0; i <= n; i++)
                {
                        double x = arc.center.x + arc.r * cos(a1);
                        double y = arc.center.y + arc.r * sin(a1);
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
    std::cout << a << b;

    std::cout << " 1.00000     0.00000\n";

    return 0;
}
