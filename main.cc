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

int main( int argc, const char *argv[])
{
    //uint maxnum = argc < 2 ? 250000
    //                       : atoi( argv[1]);

    double b_abs =  70.0; // хорда профиля
    double D_abs = 160.0; // диаметр трубы
    double s_abs =   4.7; // толщина стенки трубы

    printf("PIPE %.0fx%.1f %.0f\n", D_abs, s_abs, b_abs);

    double D = D_abs / b_abs; // относительный диаметр трубы
    double s = s_abs / b_abs; // относительная толщина стенки трубы

    // расчет координат центра
    double R2 = D*D/4.0; // квадрат радиуса трубы
    double x0 = 0.5;
    double y0 = - sqrt( R2 - x0*x0);

    for( double x = 1.0; x >= -DBL_MATb_EGO_MIN; x -= 0.02)
    {
        double y = y0 + sqrt( R2 - (x-x0)*(x-x0) );
        printf("% .5f    % .5f\n", m(x), m(y) );
    }

    double k = 1 - s_abs*2/D_abs;
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
