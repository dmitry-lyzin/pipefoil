#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <clocale>
#include <cassert>
#include <type_traits>
#include <ciso646>
#include <cstdint>
#include <bit>
#include <limits>
#include <climits>
#include <cstddef>

#pragma warning( disable: 4455)

#define CE      constexpr
#define CEC     constexpr const
#define CEOP    constexpr operator
#define EXPL    constexpr explicit
#define FRIEND  constexpr friend const
#define STATIC  constexpr static
#define AUTO    constexpr auto
#define OP      operator
#define C       const
#define CSelf   const Self
#define self    (*this)
//#define ㄥ      *cis

using std::declval;
using std::ostream;
using std::false_type;
using std::true_type;
using std::enable_if;
using std::is_same;
using std::is_integral;
using std::is_floating_point;
using std::is_unsigned;
using std::is_signed;
using std::make_signed;
using std::make_unsigned;
using std::numeric_limits;
using std::cout;
using std::cerr;
#if __cplusplus >= 201411L
    using std::size;
#else
    template< class T, std::size_t N>
    CE std::size_t size( const T (&array)[N] ) { return N; }
#endif

// ищем bit_cast
#if __cplusplus >= 201806L
    using std::bit_cast;
#else
#   ifdef _MSC_VER
#       define bit_cast _Bit_cast
        using std::bit_cast;
#   else
#       ifdef __has_builtin
#           if __has_builtin( __builtin_bit_cast)

                template< typename Dest, typename Source>
                typename enable_if
                <          sizeof( Dest) == sizeof( Source)
                        && std::is_trivially_copyable< Dest  >::value
                        && std::is_trivially_copyable< Source>::value
                ,       Dest
                >::type CE bit_cast( const Source& source)
                {
                        return __builtin_bit_cast( Dest, source);
                }

#           endif
#       else
                template< class Dest, class Source>
                typename enable_if
                <          sizeof( Dest) == sizeof( Source)
                        && std::is_trivially_copyable< Source>::value
                        && std::is_trivially_copyable< Dest  >::value
                ,       Dest
                >::type CE bit_cast( const Source& source)
                {
                        static_assert( std::is_trivially_constructible_v< Dest>,
                                "This implementation additionally requires "
                                "destination type to be trivially constructible");

                        Dest dest;
                        std::memcpy( &dest, &source, sizeof( Dest)); // требуется constexpr
                        return dest;
                }
#       endif
#   endif
#endif

C class Ф{} ф; // МАМОЙ КЛЯНУСЬ (флаг для конструкторов из сырых данных)

#pragma region // template operators

#define FN( f) constexpr friend auto f( const Self& r) -> decltype( r.f()) { return r.f(); }

template< typename T>
AUTO lastbit = sizeof( T) * CHAR_BIT - 1;

template< typename X
        , typename = decltype( declval<X>().print( declval< ostream&>()) )
        >
ostream& operator <<( ostream& os, const X& x)
{
        x.print( os);
        return os;
}

#define TINT template< typename INT  , typename = typename enable_if< is_integral      < INT  >::value>::type> constexpr
#define TFLT template< typename FLOAT, typename = typename enable_if< is_floating_point< FLOAT>::value>::type> constexpr

namespace impl
{
        template<class T> struct Longer           {                        };
        template<       > struct Longer< char    >{ using type = short;    };
        template<       > struct Longer<  uint8_t>{ using type = uint16_t; };
        template<       > struct Longer< uint16_t>{ using type = uint32_t; };
        template<       > struct Longer< uint32_t>{ using type = uint64_t; };
        //template<       > struct Longer< uint64_t>{ using type = uint128_t;};
        template<       > struct Longer<   int8_t>{ using type =  int16_t; };
        template<       > struct Longer<  int16_t>{ using type =  int32_t; };
        template<       > struct Longer<  int32_t>{ using type =  int64_t; };
        //template<       > struct Longer<  int64_t>{ using type =  int128_t;};

        template<class T> struct Shorter           {                        };
        template<       > struct Shorter< uint16_t>{ using type =  uint8_t; };
        template<       > struct Shorter< uint32_t>{ using type = uint16_t; };
        template<       > struct Shorter< uint64_t>{ using type = uint32_t; };
        //template<       > struct Shorter<uint128_t>{ using type = uint64_t; };
        template<       > struct Shorter<  int16_t>{ using type =   int8_t; };
        template<       > struct Shorter<  int32_t>{ using type =  int16_t; };
        template<       > struct Shorter<  int64_t>{ using type =  int32_t; };
        //template<       > struct Shorter< int128_t>{ using type =  int64_t; };
}
template< typename T> using Longer  = typename impl::Longer < T>::type;
template< typename T> using Shorter = typename impl::Shorter< T>::type;
template< typename T> using Signed  = typename make_signed  < T>::type;
template< typename T> using Unsgned = typename make_unsigned< T>::type;

// пустой список
struct Ø {} ø;

namespace typelist_impl
{
        // предикат для построения списка (почти как в Prolog'е ;-)
        template< class H, class T    > struct ·{};

        template< class...            > struct List            {};
        template<                     > struct List<         > { using type = Ø;                                };
        template< class H, class... Ts> struct List< H, Ts...> { using type = ·< H, typename List< Ts...>::type>; };

        template< class H             > struct Head            { using type = H; };
        template< class H, class T    > struct Head< ·<H, T> > { using type = H; };

        template< class H             > struct Tail            { using type = Ø; };
        template< class H, class T    > struct Tail< ·<H, T> > { using type = T; };

        // немного тестов
        using Test_type_list = List< char, short, int, long>;
        static_assert( is_same< typename Test_type_list::type, ·<char, ·<short, ·<int, ·<long, Ø>>>> >::value	, "*");
}
template< class    T > using Head = typename typelist_impl::Head< T    >::type; // "голова" (первый элемент) списка типов
template< class    T > using Tail = typename typelist_impl::Tail< T    >::type; // "хвост" (все элементы кроме первого) списка типов
template< class... Ts> using List = typename typelist_impl::List< Ts...>::type; // преобразовать набор аргументов в список типов

template< typename T, typename Other = T>
struct Mul_div_assign
{
        FRIEND  T&    OP *=( T& t, C Other& other) { return t = t * other; }
        FRIEND  T&    OP /=( T& t, C Other& other) { return t = t / other; }
};

template< typename T, typename Other = T>
struct Arith_assign: Mul_div_assign< T, Other>
{
        FRIEND  T&   OP +=(  T& t, C Other& other) { return t = t +   other ; }
        FRIEND  T&   OP -=(  T& t, C Other& other) { return t = t + (-other); }
        FRIEND  auto OP - (C T& t, C Other& other) { return     t + (-other); }
};

// добавить операции сравнения с другим типом
template< typename T, typename Other = T>
struct Comparable
{
        FRIEND  bool OP !=( C T& t, C Other& other) { return !(t == other); }
        FRIEND	bool OP >=( C T& t, C Other& other) { return !(t <  other); }
        FRIEND	bool OP <=( C T& t, C Other& other) { return !(t >  other); }

        FRIEND  bool OP !=( C Other& other, C T& t) { return !(t == other); }
        FRIEND	bool OP >=( C Other& other, C T& t) { return !(t >  other); }
        FRIEND	bool OP <=( C Other& other, C T& t) { return !(t <  other); }

#if __cplusplus < 202000
        FRIEND  bool OP ==( C Other& other, C T& t) { return   t == other ; }
#endif
        FRIEND  bool OP < ( C Other& other, C T& t) { return   t >  other ; }
        FRIEND  bool OP > ( C Other& other, C T& t) { return   t <  other ; }
};

// добавить операции сравнения с самим собой
template< typename T>
struct Comparable< T, T>
{
        FRIEND	bool OP !=( C T& t1, C T& t2) { return !(t1 == t2); }
        FRIEND	bool OP >=( C T& t1, C T& t2) { return !(t1 <  t2); }

        FRIEND	bool OP <=( C T& t1, C T& t2) { return !(t2 <  t1); }
        FRIEND  bool OP > ( C T& t1, C T& t2) { return   t2 <  t1 ; }
};

template< class T, class Others>
struct Arith_binary_rec
     : Arith_binary_rec< T, Tail< Others>>

        , Arith_assign  < T, Head< Others>>
        , Comparable    < T, Head< Others>>
{
        using Other = Head< Others>;
        FRIEND  auto OP + ( C Other& other, C T& t) { return   t  + other; }
        FRIEND  auto OP - ( C Other& other, C T& t) { return (-t) + other; }
        FRIEND  auto OP * ( C Other& other, C T& t) { return   t  * other; }
        //FRIEND  auto OP / ( C Other& other, C T& t) { return T( other) / t;}
        FRIEND  auto OP / ( C Other& other, C T& t) { return (conj(t) * other) / abs²(t);}
};

template< class T>
struct Arith_binary_rec< T, Ø>
        : Arith_assign  < T>
        , Comparable    < T>
{};

template< class... Ts>
struct Arith_binary
        : Arith_binary_rec      < Head< List< Ts...>>
                                , Tail< List< Ts...>>
                                >
{};

template< class T, class Others>
struct Mul_div_binary_rec
     : Mul_div_binary_rec< T, Tail< Others>>

        , Mul_div_assign< T, Head< Others>>
{
        using Other = Head< Others>;
        FRIEND  auto OP * ( C Other& other, C T& t) { return t         * other; }
        FRIEND  auto OP / ( C Other& other, C T& t) { return t.recip() * other; }
        //FRIEND  auto OP / ( C Other& other, C T& t) { return (t.conj() * other) / t.abs²(); }
};

template< class T>
struct Mul_div_binary_rec< T, Ø>
        : Mul_div_assign  < T>
{};

template< class... Ts>
struct Mul_div_binary
        : Mul_div_binary_rec    < Head< List< Ts...>>
                                , Tail< List< Ts...>>
                                >
{};

template< typename T>
struct Arith_function
{
        FRIEND  auto    abs  ( C T& t) { return t.abs  (); }
        FRIEND  auto    abs² ( C T& t) { return t.abs² (); }
        FRIEND  auto    recip( C T& t) { return t.recip(); }
        FRIEND  auto    ⅟   ( C T& t) { return t.recip(); }
        FRIEND  auto    ²    ( C T& t) { return t.²    (); }
};

template< typename T>
class Near: Comparable< Near<T>, T>
{
        C T&	ref;	// ссылка на "родителя"
public:
        EXPL	Near		( C T& t): ref( t) {}
        CEC	bool	OP ==	( C T& t) C { return   ref.near		( t); }
        CEC	bool	OP <	( C T& t) C { return ! ref.near_greater	( t); }
        CEC	bool	OP >	( C T& t) C { return ! ref.near_less	( t); }
};

template< typename T>
struct Near_comparable
{
        FRIEND	Near< T> OP ~ ( C T& t) { return Near< T>(t); }
};

#pragma endregion

// constexpr'сный способ превращения T в signed T
template< typename T>
CE auto sign_cast( T x) -> typename
enable_if
<       is_integral< T>::value && is_unsigned< T>::value
,       Signed< T>
>::type
{
        using Signed_T = Signed< T>;
        CE Signed_T min = numeric_limits< Signed_T>::min();
        // этот if нужен для осчастливливания constexpr'а,
        // надеюсь, оптимизатор его выкинет
        if( x >= min )
                return Signed_T( x - min) + min;

        return x;
}

// абсолютное значение, но с минусом, чтоб избежать глюка с abs( самое_большое_отрицательное)
template< typename T>
CE auto nabs( T x) -> typename
enable_if
<       is_integral< T>::value && is_signed< T>::value
,       T
>::type
{
        using uT = Unsgned< T>;
        // x ⩾ 0 → mask = 0
        // x < 0 → mask = -1
        const uT mask = x >> lastbit< T>;
        return mask - (uT(x) ^ mask);
};

#pragma region // Double{}

struct Double: Near_comparable< Double>
{
        #define Self Double
private:
        double	val;

public:
CE	Self		(         ): val(  ) {}
CE	Self		( double x): val( x) {}
//CE      Self		( CSelf &x): val( x.val) {}
CEOP C	double		(         ) C { return  val;	}
CEC	Self	OP -	(         ) C { return -val;	}
CE	double* OP &	(         )   { return &val;	}
CEC	double* OP &	(         ) C { return &val;	}
//CEOP 	double&		(         )   { return  val;	}
//CEOP 	C double&	(         ) C { return  val;	}
CEC	bool    isfinite(         ) C { return raw_abs( val) < raw( infinity)	; }
CEC	double  ²       (         ) C { return val * val			; }
CEC	double	abs	(         ) C { return bit_cast< double>( raw_abs( val)); }
CEC	int64_t	signbit	(         ) C { return raw(val) & (uint64_t(1) << 63)	; }
                                        // val = (x >= 0. ? val : -val)
CEC	double	copysign( double x) C { return x < 0. ? -abs() : abs()		; }
//CEC	void	sgn	( double x)   { val = bit_cast< double>( raw(val) ^ raw(x) & (uint64_t(1) << 63)); }

CEC	double	root    (         ) C {
                                              if( val < 0        ) return quiet_NaN;
                                              if( not isfinite() ) return val;
                                                                   return root_rec( val, val, 0.);
                                      }
        FN( isfinite); FN( ²); FN( abs);

CE      Self&   OP +=	( double r)   { val += r; return self; }
CE      Self&   OP -=	( double r)   { val -= r; return self; }
CE      Self&   OP *=	( double r)   { val *= r; return self; }
CE      Self&   OP /=	( double r)   { val /= r; return self; }

static  CSelf   min             ; //   10000000000000
static  CSelf   max             ; // 7FEFFFFFFFFFFFFF
static  CSelf   lowest          ; // FFEFFFFFFFFFFFFF
static  CSelf   epsilon         ; // 3CB0000000000000 2.220446049250313081e-16
static  CSelf   round_error     ; // 3FE0000000000000 0.5
static  CSelf   denorm_min      ; //                1
static  CSelf   infinity        ; // 7FF0000000000000
static  CSelf   quiet_NaN       ; // 7FF8000000000000
static  CSelf   signaling_NaN   ; // 7FF8000000000001

CEC	bool	near_less	( CSelf &x) C { return sraw(val) < sraw(x) + ε	;} // почти равный и меньше
CEC	bool	near_greater	( CSelf &x) C { return sraw(val) > sraw(x) - ε	;} // почти равный и больше
CEC	bool	near		( CSelf &x) C // почти равный
        {
                return	   nabs( int64_t( uint64_t(raw(val)) - uint64_t(raw(x)))) >= -ε
                        || near_0(val) && near_0(x)
                        ;
        }
        
private:
STATIC	int64_t	ε = 3;  // епсилон (в точках на числовой прямой)
STATIC	int64_t	raw	( double x) { return bit_cast< int64_t>( x)			;}
STATIC	int64_t	sraw	( double x) { return raw(x) ^ (uint64_t( raw(x) >> 63) >> 1)	;} // заполнить все 64 бита копиями знакового бита, и обнуляет старший бит
STATIC	int64_t	raw_abs	( double x) { return raw(x) & ~(uint64_t(1) << 63)		;}
STATIC	bool	near_0	( double x) { return raw_abs(x) < raw( min)			;}
STATIC	double  root_rec( double r, double curr, double prev)
        {
                // sqrt Newton-Raphson
                // https://gist.github.com/alexshtf/eb5128b3e3e143187794
                if( curr == prev )
                        return curr;
                return root_rec( r, (curr + r/curr) / 2, curr);
        }
        #undef Self
};

CEC	double	root	( const Double& r) { return        r .root(); }
CEC	double	root	( double r       ) { return Double(r).root(); }
CEC	double	sgn	( double x, double y)
{
        return y < 0 ? -x : x;
        //return bit_cast< double>( bit_cast< uint64_t>(x) ^ (bit_cast< uint64_t>(y) & (uint64_t(1) << 63)));
}

CEC	Double	Double::min             = numeric_limits< double>::min          (); //   10000000000000
CEC	Double	Double::max             = numeric_limits< double>::max          (); // 7FEFFFFFFFFFFFFF
CEC	Double	Double::lowest          = numeric_limits< double>::lowest       (); // FFEFFFFFFFFFFFFF
CEC	Double	Double::epsilon         = numeric_limits< double>::epsilon      (); // 3CB0000000000000 2.220446049250313081e-16
CEC	Double	Double::round_error     = numeric_limits< double>::round_error  (); // 3FE0000000000000 0.5
CEC	Double	Double::denorm_min      = numeric_limits< double>::denorm_min   (); //                1
CEC	Double	Double::infinity        = numeric_limits< double>::infinity     (); // 7FF0000000000000
CEC	Double	Double::quiet_NaN       = numeric_limits< double>::quiet_NaN    (); // 7FF8000000000000
CEC	Double	Double::signaling_NaN   = numeric_limits< double>::signaling_NaN(); // 7FF8000000000001

CEC	Double	OP ""_d	(        long double x) { return x      ; }
CEC	Double	OP ""_d	( unsigned long long x) { return x      ; }
CEC	double	OP ""²	(        long double x) { return x * x  ; }
CEC	double	OP ""²	( unsigned long long x) { return x * x  ; }
CEC	double	OP ""⁰⁵	(        long double x) { return root(x); }
CEC	double	OP ""⁰⁵	( unsigned long long x) { return root(x); }

CE void Double_test()
{
        static_assert( sgn( 123.,  124.) ==  123., "*");
        static_assert( sgn( 123., -124.) == -123., "*");
        static_assert( sgn( 123.,  123.) ==  123., "*");
        static_assert( sgn(-123., -123.) ==  123., "*");
        static_assert( sgn( 123.,   -0.) ==  123., "*");
        static_assert( sgn( 123.,    0.) ==  123., "*");

        CE double a   = 0.00000000000000000001;
        CE Double a_m = a * 10;
        CE Double a_p = a+a+a+a+a+a+a+a+a+a;
        AUTO      vsv = a_m - a_p; // very small value
        CE int64_t a_m_ulp = bit_cast< int64_t>( a_m);
        CE int64_t a_p_ulp = bit_cast< int64_t>( a_p);
        CE int64_t vsv_ulp = bit_cast< int64_t>( vsv);
        CE int64_t eps_ulp = a_m_ulp - a_p_ulp;

        CE int64_t erwerwe = bit_cast< int64_t>( numeric_limits< double>::epsilon() / ( uint64_t(1) << 55));
        CE int64_t eeter = bit_cast< int64_t>( 4.440892098500626162e-16);
        AUTO       eter = bit_cast< double>( erwerwe );
        AUTO       ghjg = numeric_limits< double>::epsilon() / ( uint64_t(1) << 55);

        static_assert( -0.0_d==~ 0.0_d					,"*");
        static_assert( !(a_m ==  a_p)					,"*");
        static_assert(   a_m !=  a_p					,"*");
        static_assert(   a_m ==~ a_p					,"*");
        static_assert( !(a_m <=  a_p)					,"*");
        static_assert( !(a_p >=  a_m)					,"*");
        static_assert(   a_m <=~ a_p					,"*");
        static_assert(   a_m >=~ a_p					,"*");
        static_assert(   a_p <=~ a_m					,"*");
        static_assert(   a_p >=~ a_m					,"*");

        static_assert(     isfinite( Double::max		)	,"*");
        static_assert(     isfinite( Double::lowest		)	,"*");
        static_assert( not isfinite( Double::infinity		)	,"*");
        static_assert( not isfinite( Double::quiet_NaN		)	,"*");
        static_assert( not isfinite( Double::signaling_NaN	)	,"*");

        //static_assert( near_0( vsv)					,"*");
        //static_assert( vsv ==~ 0d					,"*");
        static_assert( abs( -.125_d) == .125				,"*");
        static_assert( root(  152.41383936_d ) ==~  12.3456_d		,"*");
        static_assert( root(²(5678.2023_d) ) == 5678.2023_d		,"*");
        static_assert( root(  5678.2023²   ) == 5678.2023		,"*");

        static_assert(   .122999999999999999 !=  .12300000000000001_d	,"*");
        static_assert(   .122999999999999999 ==~ .12300000000000001_d	,"*");
        static_assert(   .122999999999999999 <   .12300000000000001_d	,"*");
        static_assert( !(.122999999999999999 < ~ .12300000000000001_d)	,"*");
        static_assert(   .1229999999999999   < ~ .12300000000000001_d	,"*");
        static_assert(   .122999999999999999 <=~ .12300000000000001_d	,"*");

        static_assert(   .12300000000000001 !=  .122999999999999999_d	,"*");
        static_assert(   .12300000000000001 ==~ .122999999999999999_d	,"*");
        static_assert(   .12300000000000001 >   .122999999999999999_d	,"*");
        static_assert( !(.12300000000000001 > ~ .122999999999999999_d)	,"*");
        static_assert(   .12300000000000001 > ~ .1229999999999999_d	,"*");
        static_assert(   .12300000000000001 >=~ .122999999999999999_d	,"*");

        static_assert( ! 0.						,"*");
        static_assert(    double( vsv)					,"*");
        //static_assert( ! ~Double( vsv)				,"*");
        static_assert( ~1_d == .1 * 10					,"*");
        static_assert( 1. == Double(.1 * 10)				,"*");
}

#pragma endregion

#pragma region // Angle{}

struct Turn;

CEC	double	π	                        = 3.14159265358979323846;
CEC	double	OP ""π  (        long double a) { return π * a;}
CEC	double	OP ""π  ( unsigned long long a) { return π * a;}

// Угол
struct Angle    : Comparable	 < Angle>
                , Near_comparable< Angle>
{
        #define Self Angle
        #pragma warning( push)
        #pragma warning( disable: 4146)

        friend struct Turn; 

        using Val  = uint32_t;

private:
        Val     val;
STATIC  Val     semiturn = Val(1) << lastbit< Val>;
STATIC  Signed< Val> ε = 8;

TINT    Self         ( INT        v ): val( v) {}

public:
EXPL    Self         ( C Turn&      );
EXPL    Self         ( double radian): val( Val(-1) & Signed< Longer< Val>>( radian * (semiturn/π)) ) {}
EXPL OP C double     (              ) C { return val * (π/semiturn)		; }
CEC     bool    OP ==( CSelf&     r ) C { return val == r.val			; }
CEC     bool    OP < ( CSelf&     r ) C { return val <  r.val			; }
CEC	bool	near ( CSelf&     r ) C { return nabs( sign_cast(val-r.val))>-ε	; } // почти равный
CEC	bool	near_greater(CSelf&r) C { return       sign_cast(val-r.val) <-ε	; } // почти равный и больше
CEC	bool	near_less   (CSelf&r) C { return       sign_cast(val-r.val) > ε	; } // почти равный и меньше
CEC     Self    OP - (              ) C { return -val				; }
CEC     Val*    OP & (              ) C { return &val				; }

CEC     Self    OP + ( CSelf&     r ) C { return val + r.val                    ; }
CEC     Self    OP - ( CSelf&     r ) C { return val - r.val                    ; }
CEC     Self    OP * ( double     r ) C { return Val( val * r )                 ; }
CEC     Self    OP / ( double     r ) C { return Signed< Val>( val / r )        ; }
TINT C  Self    OP * ( INT        r ) C { return val * r                        ; }
TINT C  Self    OP / ( INT        r ) C { return sign_cast( val) / r            ; }
CEC     Self    OP / ( unsigned   r ) C { return            val  / r            ; }
CEC     Val     OP / ( CSelf&     r ) C { return val / r.val                    ; }
        void    print( ostream&   os) C { os << ( val * (180./semiturn)) << '°' ; }

FRIEND  Self    OP ""ᵒ( unsigned long long gradus);
FRIEND  Self    OP ""ᵒ(        long double gradus);
        #pragma warning( pop)
        #undef Self
};

CEC Angle OP ""ᵒ ( unsigned long long gradus)
{
        using Val = Angle::Val;
        return Val( Val(-1) & (Angle::semiturn * gradus / 180) );
};
CEC Angle OP ""ᵒ (        long double gradus)
{
        using Val = Angle::Val;
        return Val( Val(-1) & Longer< Val>( static_cast< long double>(Angle::semiturn) * gradus / 180) );
};

CE void test_Angle()
{
        #pragma warning( push)
        #pragma warning( disable: 4146)
        static_assert( Angle(-π   ) ==  180ᵒ, "***" );
        static_assert( Angle(-π/ 2) ==  270ᵒ, "***" );
        static_assert( Angle(5π   ) ==  180ᵒ, "***" );
        static_assert( Angle(5π/ 2) ==   90ᵒ, "***" );
        static_assert( Angle(3π/ 2) ==  270ᵒ, "***" );
        static_assert( Angle(3π/ 2) ==  -90ᵒ, "***" );
        static_assert( Angle( π   ) ==  180ᵒ, "***" );
        static_assert( Angle( π/ 2) ==   90ᵒ, "***" );
        static_assert( Angle( π/ 3) ==   60ᵒ, "***" );
        static_assert( Angle( π/ 4) ==   45ᵒ, "***" );
        static_assert( Angle( π/ 5) ==   36ᵒ, "***" );
        static_assert( Angle( π/ 6) ==   30ᵒ, "***" );
        static_assert( Angle( π/ 8) == 22.5ᵒ, "***" );
        static_assert( Angle( π/ 9) ==   20ᵒ, "***" );
        static_assert( Angle( π/10) ==   18ᵒ, "***" );
        static_assert( Angle( π/12) ==   15ᵒ, "***" );
        static_assert( Angle( π/15) ==   12ᵒ, "***" );
        static_assert( Angle( π/16) ==11.25ᵒ, "***" );
        static_assert( Angle( π/18) ==   10ᵒ, "***" );
        static_assert( Angle( π/20) ==    9ᵒ, "***" );
        static_assert( Angle( π/24) ==  7.5ᵒ, "***" );
        static_assert( Angle( π/25) ==  7.2ᵒ, "***" );
        static_assert( Angle(π/180) ==    1ᵒ, "***" );
        static_assert(        -180ᵒ ==  180ᵒ, "***" );
        static_assert(    1ᵒ - 359ᵒ ==    2ᵒ, "***" );
        static_assert(  1.5ᵒ + 2.5ᵒ ==    4ᵒ, "***" );
        static_assert(  1.5ᵒ *  10  ==   15ᵒ, "***" );
        static_assert(          15ᵒ == 1.5ᵒ * 10, "***" );
        static_assert(  -10ᵒ *  -3  ==   30ᵒ, "***" );
        static_assert(   30ᵒ /  -3  ==  -10ᵒ, "***" );
        static_assert(   30ᵒ /  -3. ==  -10ᵒ, "***" );
        static_assert(  -90ᵒ *  -3  ==  270ᵒ, "***" );
        static_assert(    1ᵒ *  10  ==~  10ᵒ, "***" );
        static_assert(   36ᵒ *  10  ==~   0ᵒ, "***" );
        static_assert(  359ᵒ *  10  ==~ -10ᵒ, "***" );
        static_assert(   90ᵒ /  10ᵒ ==    9 , "***" );
        static_assert(  -90ᵒ /  10ᵒ ==   27 , "***" );
        //static_assert( -90ᵒ  / -10ᵒ ==    9 , "***" );
        #pragma warning( pop)
}
#pragma endregion

#pragma region // Turn{}
// Оборот
struct Turn
{
        #define Self Turn
        friend struct Angle; 

        using Val  = Signed< Longer< Angle::Val>>;

private:
STATIC  Val     one_turn = Val(1) << (lastbit< Unsgned< Angle::Val>> + 1);
        Val     val;
CE      Self            ( Ф, Val v): val( v           ) {}
public:
CE      Self            ( Angle  x): val( x.val       ) {}
CE      Self            ( double x): val( x * one_turn) {}
TINT    Self            ( INT    x): val( x * one_turn) {}
CEOP  C double          (         ) C { return double(val) / one_turn; }
TINT OP C INT           (         ) C { return        val  / one_turn; }

CEC     Self    OP -    (         ) C { return { ф,-val         }; }
CEC     Self    OP +    ( CSelf& r) C { return { ф, val + r.val }; }
CEC     Self    OP -    ( CSelf& r) C { return { ф, val - r.val }; }
CEC     Self    OP *    ( double r) C { return { ф, Val( val*r) }; }
CEC     Self    OP /    ( double r) C { return Signed< Val>( double(val) / r ); }
TINT C  Self    OP *    ( INT    r) C { return { ф, val * r     }; }
TINT C  Self    OP /    ( INT    r) C { return { ф, val / r     }; }
CEC     Val     OP /    ( CSelf& r) C { return      val / r.val  ; }
        void    print   ( ostream& os) C
        {
                os << double( val) << " turn";
        }
        #undef Self
};

CEC Turn OP ""turn ( unsigned long long x) { return         x ; };
CEC Turn OP ""turn ( long double        x) { return double( x); };
AUTO     OP ""τ    ( unsigned long long x) { return OP ""turn ( x); };
AUTO     OP ""τ    ( long double        x) { return OP ""turn ( x); };

CE Angle::Angle( C Turn &x): val( x.val) {}

CE void test_Turn()
{
#pragma warning( push)
#pragma warning( disable: 4146)
        static_assert( Angle( Turn( 123ᵒ)) ==  123ᵒ, "***" );
        static_assert( Angle(         2τ ) ==    0ᵒ, "***" );
        static_assert( Angle(      1.25τ ) ==   90ᵒ, "***" );
        static_assert( Angle(     -1.25τ ) ==  -90ᵒ, "***" );
        static_assert( Angle(  0.1τ * 10 ) ==~   0ᵒ, "***" );
        static_assert( int( Turn(180ᵒ-120ᵒ)*160) == 26, "***" );
#pragma warning( pop)
}
#pragma endregion

#pragma region // Комплексные Числа Co{}, C͡o{}

struct  𝟏_t;	// комплексная единица
struct ˗𝟏_t;	// комплексная минус единица
struct  𝐢_t;	// мнимая единица
struct ˗𝐢_t;	// мнимая минус единица

template< typename Self>
struct Complex_function
{
        // унарный минус
        FRIEND  Self    OP -( CSelf& z) { return {ф, -z.re(), -z.im() }; }
        FRIEND  Self    conj( CSelf& z) { return {ф,  z.re(), -z.im() }; }
        FRIEND  Self    root( CSelf& z) // корень квадратный
        {
                const auto ǀzǀ = z.abs();
                return  { ф
                        ,      root( (ǀzǀ + z.re())/2 ) 
                        , sgn( root( (ǀzǀ - z.re())/2 ), z.im())
                        };
        }
};

// Комплексное Число
struct Co:
#define Self Co
          Arith_function	< Self	>
        , Complex_function	< Self	>
        , Near_comparable	< Self	>
        , Arith_binary		< Self, double, 𝐢_t, ˗𝐢_t, 𝟏_t, ˗𝟏_t>
{
protected:
        Double  r, i;

public:
CE      Self( double r, double i): r( r ), i( i ) {}
CE      Self( const Self& z     ): Self(z.r, z.i) {}
CE      Self(Ф,double r,double i): Self(  r,   i) {} // для унификации с C͡o
CE      Self( double s          ): Self(  s,   0) {}

CEC     double  abs²	() C { return  r*r + i*i     ; } // абсолютная величина в квадрате
CEC     double  abs	() C { return root( abs²())  ; } // абсолютная величина
CEC     Self    conj	() C { return {  r, -i }     ; } // Сопряжённое (conjugate) число
CEC     Self    recip	() C { return conj() / abs²(); } // 1/z - обратная величина (reciprocal, multiplicative inverse)
CEC     double  re	() C { return r              ; } // действительная часть
CEC     double  im	() C { return i              ; } // мнимая часть
CEC     double  ℜ	() C { return r              ; } // действительная часть
CEC     double  ℑ	() C { return i              ; } // мнимая часть
CEC     Angle   arg	() C
        {
                if( r > +0. ) return Angle( atan( i / r ));
                if( r < -0. ) return Angle( atan( i / r ) + π);

                if( i > +0. ) return Angle( π/2 );
                if( i < -0. ) return Angle( 3π/2);

                return Angle( NAN);
        }

        FN( re); FN( im); FN( ℜ); FN( ℑ); FN( arg);

CEC     bool    OP ==( CSelf& z) C { return r  ==  z.r  && i  ==  z.i     ; }
CEC	bool	near ( CSelf& z) C { return r.near(z.r) && i.near(z.i)    ; } // почти равный

CEC     Self    OP + ( CSelf& z) C { return {r     + z.r  , i     + z.i  }; }
CEC     Self    OP * ( CSelf& z) C { return {r*z.r - i*z.i, i*z.r + r*z.i}; }
CEC     Self    OP / ( CSelf& z) C { return (self * z.conj()) / z.abs²()  ; } // Деление КЧ на другое КЧ
/*
CEC     Self    OP * ( CSelf& z) C
        {
                double a = (r + i) * (z.r + z.i);
                double b = r * z.r;
                double c = i * z.i;

                return  { b - c
                        , a - b - c
                        };
        }
*/

CEC     Self    OP + ( double s) C { return {r + s, i    }; }
CEC     Self    OP * ( double s) C { return {r * s, i * s}; } // Умножение на скаляр
CEC     Self    OP / ( double s) C { return {r / s, i / s}; } // Деление на скаляр

CEC     Self    OP + ( C  𝐢_t& ) C { return {  r  , i + 1}; }
CEC     Self    OP * ( C  𝐢_t& ) C { return { -i  , r    }; } // Умножение на мнимую
CEC     Self    OP / ( C  𝐢_t& ) C { return {  i  ,-r    }; } // Деление на мнимую
CEC     Self    OP + ( C ˗𝐢_t& ) C { return {  r  , i - 1}; }
CEC     Self    OP * ( C ˗𝐢_t& ) C { return {  i  ,-r    }; } // Умножение на мнимую
CEC     Self    OP / ( C ˗𝐢_t& ) C { return { -i  , r    }; } // Деление на мнимую

CEC     Self    OP + ( C  𝟏_t& ) C { return  self + 1; }
CEC     Self    OP * ( C  𝟏_t& ) C { return  self    ; }
CEC     Self    OP / ( C  𝟏_t& ) C { return  self    ; }
CEC     Self    OP + ( C ˗𝟏_t& ) C { return  self - 1; }
CEC     Self    OP * ( C ˗𝟏_t& ) C { return -self    ; }
CEC     Self    OP / ( C ˗𝟏_t& ) C { return -self    ; }

// на время представим, что комплексные числа это векторы...
        // Скалярное произведение (Inner product)
CEC     double  OP , ( CSelf& v) C { return  r*v.r + i*v.i; }
CEC     double  OP ^ ( CSelf& v) C { return  r*v.i - i*v.r; } // Псевдоскалярное, Векторное, косое произведение (Outer, cross product)

        void    print( ostream& os) C
        {
                static Self last = { NAN, NAN };
                if( self !=~ last )
                {
                        os << std::setw(8) << r << std::setw(12) << i << '\n';
                        last = self;
                }
        };
#undef Self
};

CE void test_Co()
{
        CE Co z = {3, 4};
        static_assert( abs ( z) == 5		, "*");
        static_assert( abs²( z) == 25		, "*");
        static_assert( root( z) == Co( 2, 1)	, "*");
        static_assert( z + 1 == Co( 4	, 4   )	, "*");
        static_assert( 1 + z == Co( 4	, 4   )	, "*");
        static_assert( z * 2 == Co( 6	, 8   )	, "*");
        static_assert( 2 * z == Co( 6	, 8   )	, "*");
        static_assert( z / 2 == Co( 1.5	, 2   )	, "*");
        static_assert( 2 / z == Co( 0.24,-0.32)	, "*");
        static_assert( z + 1.== Co( 4	, 4   )	, "*");
        static_assert( 1.+ z == Co( 4	, 4   )	, "*");
        static_assert( z * 2.== Co( 6	, 8   )	, "*");
        static_assert( 2.* z == Co( 6	, 8   )	, "*");
        static_assert( z / 2.== Co( 1.5	, 2   )	, "*");
        static_assert( 2./ z == Co( 0.24,-0.32)	, "*");
        static_assert( z / z == 1		, "*");
        //static_assert( (Co( z) -= 1) == Co(2, 4), "*");
        //static_assert( (Co( z) /= 1) == z	, "*");
        //static_assert( (Co(-z) += z) == 0	, "*");
        static_assert(  -z != z			, "*");
        static_assert(   z != 1			, "*");
        static_assert(   1 != z			, "*");
        static_assert( !(1 == z)		, "*");
}

// Единичное Комплексное Число (Unit Complex Number). C͡o = {𝑧 ∈ ℂ: |𝑧| = 1}
struct C͡o: Co
#define Self C͡o
        , Arith_function	< Self	>
        , Complex_function	< Self	>
        , Mul_div_binary	< Self, Co, double, 𝟏_t, ˗𝟏_t, 𝐢_t, ˗𝐢_t>
{
protected:
CE      Self( double r, double i): Co( r, i) {}

public:
CE      Self( Ф, C Co& z ): Co( z                                    ) { /*assert( z.abs²()  == 1.);*/ } // мамой клянусь - |z| = 1
CE      Self( Ф, double r
               , double i): Co( r, i                                 ) { /*assert( r*r + i*i == 1.);*/ } // мамой клянусь - |r + i𝑖| = 1
EXPL    Self(    C Co& z ): Co( z / z.abs()                          ) {}
CE      Self( double sin ): Co( ::root( 1. - sin*sin), sin           ) {}
CE      Self( Angle    𝜑 ): Co( cos( double( 𝜑)), sin( double( 𝜑))  ) {} // фазор угла 𝜑 (единичное КЧ с аргументом 𝜑)

CEC     double  abs² () C { return     1.   ; } // абсолютная величина в квадрате = 1
CEC     double  abs  () C { return     1.   ; } // абсолютная величина = 1, этож единичное число ;-)
CEC     Self    conj () C { return { r, -i }; } // Сопряжённое (conjugate) число
CEC     Self    recip() C { return { r, -i }; } // 1/z - обратная величина (reciprocal, multiplicative inverse)
//CEC     double  ψarg () C { double x1 = r + 1.; return i>=0 ? -x1 : x1;} // псевдоугол (для сравнений)
CEC     double  ψarg () C { return sgn( r+1., -i);} // псевдоугол (для сравнений)
        FN( ψarg);

CEC     bool    OP < ( CSelf& z) C { return ψarg() < z.ψarg()   ; }
CEC     Self    OP * ( CSelf& z) C { return {ф, Co::OP *(z) }   ; } // ЕКЧ*ЕКЧ = ЕКЧ
CEC     Self    OP / ( CSelf& z) C { return self * z.conj()     ; } // ЕКЧ/ЕКЧ = ЕКЧ, при этом деление вырождено

CEC     Co      OP * ( C Co&  z) C { return Co::OP *(z); } // ЕКЧ*КЧ = КЧ
CEC     Co      OP / ( C Co&  z) C { return Co::OP /(z); } // ЕКЧ/КЧ = КЧ
CEC     Co      OP * ( double x) C { return Co::OP *(x); } // ЕКЧ* Ч = КЧ
CEC     Co      OP / ( double x) C { return Co::OP /(x); } // ЕКЧ/ Ч = КЧ

CEC     Self    OP * ( C  𝟏_t& ) C { return  self; }
CEC     Self    OP / ( C  𝟏_t& ) C { return  self; }
CEC     Self    OP * ( C ˗𝟏_t& ) C { return -self; }
CEC     Self    OP / ( C ˗𝟏_t& ) C { return -self; }

CEC     Self    OP * ( C  𝐢_t& ) C { return {-i, r}; }
CEC     Self    OP / ( C  𝐢_t& ) C { return { i,-r}; }
CEC     Self    OP * ( C ˗𝐢_t& ) C { return { i,-r}; }
CEC     Self    OP / ( C ˗𝐢_t& ) C { return {-i, r}; }
#undef Self
};

struct  𝟏_t: C͡o {	CE  𝟏_t	(): C͡o(  1.,  0.) {};  CEC ˗𝟏_t	OP - () C;   }; // комплексная единица
struct ˗𝟏_t: C͡o {	CE ˗𝟏_t	(): C͡o( -1.,  0.) {};  CEC  𝟏_t	OP - () C;   }; // комплексная минус единица
struct  𝐢_t: C͡o {	CE  𝐢_t	(): C͡o(  0.,  1.) {};  CEC ˗𝐢_t	OP - () C;   }; // мнимая единица
struct ˗𝐢_t: C͡o {	CE ˗𝐢_t	(): C͡o(  0., -1.) {};  CEC  𝐢_t	OP - () C;   }; // мнимая минус единица

CEC ˗𝟏_t  𝟏_t::OP - () C { return {}; };
CEC  𝟏_t ˗𝟏_t::OP - () C { return {}; };
CEC ˗𝐢_t   𝐢_t::OP - () C { return {}; };
CEC  𝐢_t  ˗𝐢_t::OP - () C { return {}; };

CEC 𝟏_t	root( C  𝟏_t& ) { return {}; }
CEC 𝐢_t	root( C ˗𝟏_t& ) { return {}; }

CEC 𝟏_t	𝟏;	// комплексная единица
CEC 𝐢_t	𝐢;	// мнимая единица
CEC Co	OP ""𝐢 ( unsigned long long a) { return { 0., double(a)}; }
CEC Co	OP ""𝐢 (        long double a) { return { 0., double(a)}; }

CE void test_C͡o()
{
        CE Co z( 3, 4);
        CE C͡o ẑ( ф, 0., 1.);

        AUTO minus_ẑ = -ẑ;
        static_assert( is_same< decltype(minus_ẑ), C C͡o>::value	, "*");
        static_assert( minus_ẑ == Co( 0, -1)			, "*");

        AUTO my_i = root( -𝟏 );
        static_assert( is_same< decltype( my_i	), C 𝐢_t>::value	, "*");
        static_assert( my_i == 𝐢					, "*");

        static_assert( root( -𝟏		) == 𝐢			, "*");
        static_assert( root( -1 + 0𝐢	) == 𝐢			, "*");
        static_assert( root( -𝐢		) == conj( root( 𝐢))	, "*");

        static_assert( is_same< decltype( ẑ * 𝐢	), C C͡o>::value	, "*");
        static_assert( is_same< decltype( ẑ / 𝐢	), C C͡o>::value	, "*");
        static_assert( is_same< decltype( 𝐢 * ẑ	), C C͡o>::value	, "*");
        static_assert( is_same< decltype( 𝐢 / ẑ	), C C͡o>::value	, "*");
        static_assert( is_same< decltype( ẑ * 𝟏	), C C͡o>::value	, "*");
        static_assert( is_same< decltype( ẑ / 𝟏	), C C͡o>::value	, "*");
        static_assert( is_same< decltype( 𝟏 * ẑ	), C C͡o>::value	, "*");
        static_assert( is_same< decltype( 𝟏 / ẑ	), C C͡o>::value	, "*");
        static_assert( is_same< decltype( ẑ * ẑ	), C C͡o>::value	, "*");
        static_assert( is_same< decltype( ẑ / ẑ	), C C͡o>::value	, "*");
        static_assert( is_same< decltype( z * ẑ	), C Co>::value	, "*");
        static_assert( is_same< decltype( z / ẑ	), C Co>::value	, "*");
        static_assert( is_same< decltype( ẑ * z	), C Co>::value	, "*");
        static_assert( is_same< decltype( ẑ / z	), C Co>::value	, "*");
        static_assert( ẑ * ẑ == -1				, "*");
        static_assert( ẑ / ẑ ==  1				, "*");
        static_assert( z * ẑ == Co(-4   , 3   )			, "*");
        static_assert( z / ẑ == Co( 4   ,-3   )			, "*");
        static_assert( ẑ * z == Co(-4   , 3   )			, "*");
        static_assert( ẑ / z == Co( 0.16, 0.12)			, "*");

        static_assert( ẑ + 1 == Co( 1	, 1   )			, "*");
        static_assert( 1 + ẑ == Co( 1	, 1   )			, "*");
        static_assert( ẑ * 2 == Co( 0	, 2   )			, "*");
        static_assert( 2 * ẑ == Co( 0	, 2   )			, "*");
        static_assert( ẑ / 2 == Co( 0	, 0.5 )			, "*");
        static_assert( 2 / ẑ == Co( 0	,-2   )			, "*");
        static_assert( ẑ + 1.== Co( 1	, 1   )			, "*");
        static_assert( 1.+ ẑ == Co( 1	, 1   )			, "*");
        static_assert( ẑ * 2.== Co( 0	, 2   )			, "*");
        static_assert( 2.* ẑ == Co( 0	, 2   )			, "*");
        static_assert( ẑ / 2.== Co( 0	, 0.5 )			, "*");
        static_assert( 2./ ẑ == Co( 0	,-2   )			, "*");
        /*
        AUTO ŵ = C͡o(ẑ) /= ẑ;
        AUTO û = C͡o(ẑ) *= ẑ;
        static_assert( is_same< decltype(ŵ), C C͡o>::value	, "*");
        static_assert( is_same< decltype(û), C C͡o>::value	, "*");
        static_assert( ŵ ==  1					, "*");
        static_assert( û == -1					, "*");

        static_assert( (C͡o( ẑ) -= 1) == Co(-1, 1)		, "*");
        static_assert( (C͡o( ẑ) /= 1) == ẑ			, "*");
        static_assert( (C͡o(-ẑ) += ẑ) == 0			, "*");
        */
        CE C͡o bad_ẑ( ф, 0.28, 100500);
        static_assert( root(bad_ẑ) == Co( 0.8, 0.6)		, "*");
        static_assert( abs (bad_ẑ) == 1				, "*");
        static_assert( abs²(bad_ẑ) == 1				, "*");
        static_assert( abs (    ẑ) == 1				, "*");
        static_assert( abs²(    ẑ) == 1				, "*");
        static_assert( 1 / bad_ẑ == conj(   bad_ẑ)		, "*");
        static_assert( ẑ / bad_ẑ == conj(-ẑ*bad_ẑ)		, "*");
}
CE void test_Co_1()
{
        CE Co z = {3, 4};

        static_assert( (5 + 𝐢)*(7 - 6𝐢)/(3+𝐢)		== (10 - 11𝐢)	, "*");
        static_assert( (4 + 𝐢)*(5 + 3𝐢)+(3 + 𝐢)*(3 - 2𝐢)	== (28 + 14𝐢)	, "*");

        static_assert( z + 𝐢 == Co( 3	, 5   )	, "*");
        static_assert( 𝐢 + z == Co( 3	, 5   )	, "*");
        static_assert( z - 𝐢 == Co( 3	, 3   )	, "*");
        static_assert( 𝐢 - z == Co(-3	,-3   )	, "*");
        static_assert( z * 𝐢 == Co(-4	, 3   )	, "*");
        static_assert( 𝐢 * z == Co(-4	, 3   )	, "*");
        static_assert( z / 𝐢 == Co( 4	,-3   )	, "*");
        static_assert( 𝐢 / z == Co( 0.16, 0.12)	, "*");
}

#pragma endregion

using ℂ  = Co;
using ℂ₁ = C͡o;

struct Line: Near_comparable< Line>
{
        ℂ₁    n̂; // единичный (|n̂| = 1) перпедикуляр к прямой
        Double p; // растояние от начало координат до прямой

        // Прямая, заданой нормальным уравнением (n̅, r̅) = 𝐶 т.е. 𝐴𝑥ᵣ + 𝐵𝑦ᵣ = 𝐶, n̅ = {𝐴, 𝐵}, |n̅| > 0
        // \param[in] n̅ - вектор, нормальный к прямой
        // \param[in] 𝐶 - скалярный параметр
CE      Line( C ℂ & n̅, double 𝐶 ): n̂( ф, n̅ ), p( 𝐶 )
        {
                double ǀn̅ǀ = abs( n̅ );
                *(ℂ*)&	n̂ /= ǀn̅ǀ;
                        p /= ǀn̅ǀ;
        }
        // Прямая, задана нормированным уравнением (n̂, r̅) = 𝑝, |n̂| = 1, 𝑝 ⩾ 0
        // \param[in] n̂ - единичный вектор, нормальный к прямой
        // \param[in] 𝑝 - растояние от начало координат до прямой
CE      Line( C ℂ₁&_n̂, double 𝑝 ): n̂(_n̂), p( 𝑝 )
        {
                if( p < 0 )
                {
                        p = -p;
                        n̂ = -n̂;
                }
        };
        // Прямая через точки a̅ и b̅
CE      Line( C ℂ& a̅, C ℂ& b̅   ): Line( (a̅ - b̅)*𝐢, (a̅ ^ b̅) ) {};

CE      bool    OP ==( C Line& l) C { return p ==  l.p && n̂ ==  l.n̂; }
CE	bool	near ( C Line& l) C { return p ==~ l.p && n̂ ==~ l.n̂; } // почти равный

CE      Line    OP - (          ) C { return Line( -n̂, -p)         ; }
CE      double  dist ( C   ℂ& r̅) C { return (n̂, r̅) - p            ; } // растояние между прямой и точкой r̅
CE      ℂ      proj ( C   ℂ& r̅) C { return r̅ - n̂ * dist( r̅)      ; } // проекция точки r̅ на прямую

        // растояние между точкой r̅ и прямой l
FRIEND  double dist( C ℂ& r̅, C Line& l) { return l.dist( r̅); }
        // растояние между прямой l и точкой r̅
FRIEND  double dist( C Line& l, C ℂ& r̅) { return l.dist( r̅); }
/*
friend  ostream& OP<<( ostream &os, C Line& obj )
        {
                os << ℂ( 0.0, obj.p/obj.n.i)
                   << ℂ( obj.p/obj.n.r, 0.0);
                return os;
        };
*/
};
struct Vertical: public Line
{
        // Вертикаль пересекающая ось 𝑋 в точке x0
CE      Vertical( double x0 ): Line( 𝟏, x0) {};
};
struct Horizontal: public Line
{
        // Горизонталь пересекающая ось 𝑌 в точке y0
CE      Horizontal( double y0 ): Line( 𝐢, y0) {};
};
CE void   Line_test()
{
        CE Vertical l_test( 2 );
        static_assert(       l_test.n̂		== 1,	"*");
        static_assert(       l_test.n̂ * 𝐢	== 𝐢,	"*");
        static_assert( dist( l_test, {3, 1})	== 1,	"*");
}

struct Segment: public Line
{
        ℂ p̅1, p̅2;

CE      Segment( C ℂ& p1_, C ℂ& p2_ )
        : Line( p1_, p2_ )
        , p̅1( p1_), p̅2( p2_)
        {};

        // TODO стремный конст-ор, как бы его спрятать
CE      Segment( C Line& l, C ℂ& p1_, C ℂ& p2_ )
        : Line( l )
        , p̅1( p1_), p̅2( p2_)
        {};

        // обмен концов (рисоваться будет в другую сторону)
CE      Segment OP - () C { return { p̅2, p̅1 }; }

        void print( ostream& os) C
        {
                os << p̅1 << p̅2;
        };
};

struct Circle: Near_comparable< Circle>
{
        ℂ     o̅; // центр окружности
        Double R; // радиус окружности

CE      Circle	     ( C ℂ& center, double radius = 0. )
        : R( radius), o̅( center )
        {};

CE      bool	OP ==( C Circle& c ) C { return R ==  c.R && o̅ ==  c.o̅; }
CE	bool	near ( C Circle& c ) C { return R ==~ c.R && o̅ ==~ c.o̅; } // почти равный
CE      Circle	OP - (             ) C { return { o̅, -R }             ; }

        // Вычислить единичный вектор сонаправленый со вторым катетом
        // прямоугольного треугольника по гипотенузе и длине первого катета
STATIC  ℂ₁ cathet2( C ℂ& hypotenuse, Double cathet1 )
        {
        /*
        // нормаль к касательной к окружностям self и c2
        CE      ℂ₁ tangent_norm( C Circle& c2 ) C
        {
                double	ΔR	= c2.R - R;
                ℂ	o͞o	= c2.o̅ - o̅;
                double ǀo͞oǀ²	= abs²( o͞o );
                //double ǀo͞oǀ	= abs( o͞o );
                //ℂ₁	o͡o	= ℂ₁( ф, o͞o /ǀo͞oǀ );	// направляющий вектор линии центров
                //double	sin𝜑	= ΔR /ǀo͞oǀ;		// 𝜑 - угол между линией центров и касательной
                //ℂ₁	o͡o_r	= o͡o * ℂ₁( sin𝜑) / 𝐢;	// повернуть o͡o на 𝜑-90°

                //		  o͡o     * ℂ₁( sin𝜑)                                     / 𝐢 =
                //		  o͡o     * { √(1            - sin²𝜑     ), sin𝜑   }       / 𝐢 =
                //		  o͞o/ǀo͞oǀ* { √(1            - (ΔR/ǀo͞oǀ)²), ΔR/ǀo͞oǀ}       / 𝐢 =
                //		  o͞o/ǀo͞oǀ* { √((ǀo͞oǀ/ǀo͞oǀ)² - (ΔR/ǀo͞oǀ)²), ΔR/ǀo͞oǀ}       / 𝐢 =
                //		  o͞o/ǀo͞oǀ* { √( ǀo͞oǀ²       -  ΔR²) /ǀo͞oǀ, ΔR/ǀo͞oǀ}       / 𝐢 =
                //		  o͞o/ǀo͞oǀ* { √( ǀo͞oǀ²       -  ΔR²)      , ΔR     } /ǀo͞oǀ / 𝐢 =
                //		  o͞o     * { √( ǀo͞oǀ²       -  ΔR²)      , ΔR     } /ǀo͞oǀ²/ 𝐢 =
                //		  o͞o     * { ΔR,                - √( ǀo͞oǀ²- ΔR²)  } /ǀo͞oǀ²
                ℂ	o͡o_r	= ℂ{ ΔR, -root(ǀo͞oǀ² - ²(ΔR)) } * o͞o /ǀo͞oǀ²;
                return ℂ₁( ф, o͡o_r );

                return cathet2( -o͞o, -ΔR); // почему минусы ???
        }
        */
                double ǀhypotenuseǀ² = abs²( hypotenuse );
                return ℂ₁( ф, ℂ{ cathet1, root( ǀhypotenuseǀ² - ²(cathet1)) } * hypotenuse /ǀhypotenuseǀ² );
        }

CE      ℂ intersect( C   Line& l ) C
        {
                Double cathet1 = dist( o̅, l);   // расстояние от центра окружности до прямой
                //ℂ    p̅ = (o̅ - l.n̂ * cathet1); // точка проекции центра окружности на прямую
                //ℂ₁   v̂ = l.n̂ * 𝐢;             // направляющий вектор прямой l
                //return p̅ + v̂ * root(²(R) - ²(cathet1));
                return o̅ - l.n̂ * ℂ{ cathet1, root( ²(R) - ²(cathet1)) };
        }
CE      ℂ intersect( C Circle& c2) C
        {
                C Circle& c1 = self;
                double cathet1 = (c2.o̅, c1.o̅) - (abs²(c2.o̅) + abs²(c1.o̅) - ²(c2.R) + ²(c1.R))/2;
                return c1.o̅ - c1.R * cathet2( (c2.o̅ - c1.o̅)*c1.R, cathet1);
        }

FRIEND  Circle tangent( C Circle& c1, C Circle& c2, double R )
        {
                ℂ center = Circle( c1.o̅, R+c1.R).intersect( Circle( c2.o̅, R+c2.R) );
                return Circle( center, R );
        }

        // касательная к окружностям self и c2
CE      Line tangent( C Circle& c2 ) C
        {
                ℂ₁ n̂ = cathet2( -(c2.o̅ - o̅), -(c2.R - R)); // TODO почему минусы ???
                return Line( n̂, (n̂, o̅) - R );
        }

        // точка касания касательной исходящей из точки p к окружности
CE      ℂ tangent_point( C ℂ& p̅ ) C
        {
                return o̅ - R * cathet2( o̅ - p̅, R);
        };

        void print( ostream &os) C
        {
                CE C int segs = 40;

                //CE C ℂ₁ m̂¹⁰ = cis( 2π/segs);
                CE C ℂ₁ m̂ = 2π/segs/10; // типа 𝛼 ≈ sin 𝛼
                CE C ℂ₁ m̂⁵  = m̂*m̂*m̂*m̂*m̂;
                CE C ℂ₁ m̂¹⁰ = m̂⁵*m̂⁵;

                ℂ₁ n̂ = 𝟏;
                for( int i = segs; i --> 0; )
                {
                        os << (o̅ + R*n̂);
                        n̂ *= m̂¹⁰;
                }
        };

        void print00( ostream &os, Angle α1, Angle α2 ) C
        {
                Angle Δα = α2 - α1; // угол поворота
                // кол. сегментов, на круг - 160 сегментов, примерно
                unsigned segments = Turn( Δα * 160);
                ℂ₁ m̂( Δα / segments);   // ед. век. повернутый на угол Δα/segments

                ℂ  r̅ = abs( R) * ℂ₁( α1);    // радиус от центра окружности
                ℂ  s̅ = m̂*r̅ - r̅;               // сегментик, которым рисуем окружность
                r̅ += o̅;
                for( ; segments --> 0;)
                {
                        os << r̅;
                        r̅ += ( s̅ *= m̂ );
                }
        };
};
CE void Circle_test()
{
        CE Circle c_test( {0, 0}, 3 );
        CE ℂ o1( 5, 0);
        CE ℂ o2( 5, 5);

        AUTO    b = c_test.tangent_point( ℂ( 5, 0)       );
        //AUTO    x = Line( ℂ( 0.6,  0.8), 5);

        static_assert( Circle( o1,  2).tangent( Circle( o2,  2)) == Line( ℂ( 1.0,   .0), 3), "***");
        static_assert( Circle( o1, -2).tangent( Circle( o2,  2)) == Line( ℂ( 0.6,  0.8), 5), "***");
        static_assert( Circle( o1,  2).tangent( Circle( o2, -2)) == Line( ℂ( 0.6, -0.8), 1), "***");
        static_assert( Circle( o1, -2).tangent( Circle( o2, -2)) == Line( ℂ( 1.0,   .0), 7), "***");

        static_assert( c_test.intersect( Line( {3, 0}, {0, 1.5})) == ℂ( 3,  0  ), "***");
        static_assert( c_test.intersect( Vertical( 3)           ) == ℂ( 3,  0  ), "***");
        static_assert( c_test.intersect( Circle( {4, 0}, 3)     ) == ℂ( 2, -5⁰⁵), "***");

        static_assert( c_test.tangent_point( ℂ( 5, 0)       ) ==~ ℂ(1.8, 2.4), "***");

        static_assert( c_test.tangent( Circle( {   2, 0}, 1) ) ==  Line( 𝟏                  ,  3), "***");
        static_assert( c_test.tangent( Circle( { 2⁰⁵, 0}, 2) ) ==~ Line( ℂ(  0.5⁰⁵,  0.5⁰⁵),  3), "***");
        static_assert( c_test.tangent( Circle( {   9, 0},-3) ) ==~ Line( ℂ(   2./3,  5⁰⁵/3),  3), "***");

        static_assert( tangent( c_test, Circle( {  6, 0}, 3),  3 ) == Circle( {  3, -3*3⁰⁵},  3), "***");
        static_assert( tangent( c_test, Circle( {  2, 0}, 1), 10 ) == Circle( { 13,      0}, 10), "***");
}

struct Arc: public Circle
{
        ℂ₁ n̂₁; // единич. вектор из центра дуги на первый конец
        ℂ₁ n̂₂; // единич. вектор из центра дуги на второй конец

CE      Arc( C Circle& c, C ℂ₁& _n̂₁, C ℂ₁& _n̂₂)
        : Circle( c), n̂₁(_n̂₁), n̂₂(_n̂₂)
        {}
CE      Arc( C Circle& c, C ℂ& start_point, C ℂ& next_point )
        : Circle( c)
        , n̂₁( start_point- c.o̅ )
        , n̂₂( next_point - c.o̅ )
        {}
CE      Arc( C Circle& c, C Arc& prev, C ℂ& next_point )
        : Arc( c, prev.o̅ + abs(prev.R)*prev.n̂₂, next_point )
        {
                // ставим знак вектору n̂₂ таким, чтоб конечная точка дуги была максимально близко к next_point
                //if( ((next_point - o̅), n̂₂) < 0)
                //        n̂₂ = -n̂₂;

                Double ǀRǀ = abs(R);
                R = sgn( ǀRǀ, prev.R * (prev.n̂₂, n̂₁)); //R = sgn( R, prev.R * R * (prev.n̂₂, n̂₁));
        }

CE      bool	OP ==( C Arc& a) C { return  Circle( self) ==  a && n̂₁ ==  a.n̂₁ && n̂₂ ==  a.n̂₂; }
CE	bool	near ( C Arc& a) C { return  Circle( self) ==~ a && n̂₁ ==~ a.n̂₁ && n̂₂ ==~ a.n̂₂; } // почти равный

        void	print( ostream &os) C
        {
                CE C int segs = 160;

                //CE C ℂ₁ m̂¹⁰ = cis( 2π/segs);
                CE C ℂ₁ m̂ = 2π/segs/10; // типа 𝛼 ≈ sin 𝛼
                CE C ℂ₁ m̂⁵ = m̂*m̂*m̂*m̂*m̂;

                ℂ₁ m̂¹⁰ = m̂⁵*m̂⁵;
                bool dir = (R >= 0); // направление рисования
                if( !dir )
                        m̂¹⁰ = conj( m̂¹⁰);

                ℂ r̅ = abs(R) * n̂₁;
                os << (o̅ + r̅);
                double a₂₁ = ψarg( n̂₂/n̂₁);

                for( ℂ₁ n̂ = m̂¹⁰; (ψarg(n̂) < a₂₁) == dir; n̂ *= m̂¹⁰ )
                        os << (o̅ + r̅*n̂);
        };
};

// зазор для соблюдения постулата Жуковского-Чаплыгина (Kutta condition)
CE C ℂ TE1( 1.,  0.00001); // задняя кромка верх
CE C ℂ TE2( 1., -0.00001); // задняя кромка низ

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
        CE Circle  c_top = tangent( c_lead, TE1, R );
        CE Circle  c_bottom( c_top.o̅, R-s );
        CE Segment s_end = tangent_segment( c_bottom, TE2 );
        CE Circle  c_l1( Vertical( lr1) ^ Circle( c_top.o̅, c_bottom.R+lr1), lr1 );
        CE Circle  c_l2 = tangent( c_top, c_l1, -lr2 );

        CE Arc a_top = chain_arc( TE1      , c_top , c_l2     );
        CE Arc a_l2  = chain_arc( a_top.p̅2 , c_l2  , c_l1     );
        CE Arc a_l1  = chain_arc( a_l2.p̅2  , c_l1  , c_bottom, minus );
        CE Arc a_bottom( c_bottom, a_l1.p̅2, s_end.p̅1 );

        cout    << "PIPE " << D_abs << 'r' << s_abs << '-' << b_abs << '\n'
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
        CE Circle  c_top = tangent( c_lead, TE1, R );
        CE Circle  c_bottom( c_top.o̅, R-s );
        CE Segment s_start = tangent_segment( c_lead, -c_bottom, minus );
        CE Segment s_end = tangent_segment( c_bottom, TE2 );

        CE Arc a_top = chain_arc( TE1, c_top, c_lead      );
        CE Arc a_lead  ( c_lead  , a_top.p̅2  , s_start.p̅1 );
        CE Arc a_bottom( c_bottom, s_start.p̅2, s_end.p̅1   );

        cout    << "PIPE " << D_abs << 'x' << s_abs << '-' << b_abs << '\n'
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
        CE Circle  c_top( c_bottom.o̅, R );
        CE Circle  c_l2  = tangent( c_l1, c_top, -lr2, minus );
        CE Segment s_end = tangent_segment( c_l1, TE2, minus );

        CE Arc a_top = chain_arc( TE1     , c_top, c_l2 );
        CE Arc a_l2  = chain_arc( a_top.p̅2, c_l2 , c_l1 );
        CE Arc a_l1( c_l1, a_l2.p̅2, s_end.p̅1 );

        cout    << "PIPE " << D_abs << 'x' << s_abs << '-' << b_abs << '\n'
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
        CE Circle  c_top( c_bottom.o̅, R );
        CE Circle  c_l2  = tangent( c_l1, c_top, -lr2, minus );
        CE Segment s_end0 = tangent_segment( c_l1, TE2, minus );

        CE Arc a_top = chain_arc( TE1     , c_top, c_l2 );
        CE Arc a_l2  = chain_arc( a_top.p̅2, c_l2 , c_l1 );
        CE Arc a_l1( c_l1, a_l2.p̅2, s_end0.p̅1 );

        CE Segment s_start( s_end0.p̅1, s_end0 ^ c_bottom );
        CE Segment s_end( s_end0 & c_bottom, TE2 );
        CE Arc a_bottom( c_bottom, s_start.p̅2, s_end.p̅1 );

        cout    << "PIPE " << D_abs << 'x' << s_abs << '-' << b_abs << '\n'
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
        CE Circle c_top    = { c_bottom.o̅, R };
        CE Circle c_lef    = tangent( c_top, c_le, -lef );
        CE ℂ     p_end    = c_bottom.tangent_point( TE2 );

        CE Arc a_top   ( c_top   , TE1  , c_lef.o̅   );
        CE Arc a_lef   ( c_lef   , a_top, c_le.o̅    );
        CE Arc a_le    ( c_le    , a_lef, c_bottom.o̅);
        CE Arc a_bottom( c_bottom, a_le , p_end     );
        
        // static_asserts
        {
                static_assert(    a_top == Arc{{{0.3656664547982899638, -1.468874254606819196 },  1.6  }, ℂ₁{ф, 0.3964584657510686894, 0.9180526591292618166}, {ф,-0.10593121205827985,   0.9943734601807630025}}, "*");
                static_assert(    a_lef == Arc{{{0.2385490003283541327, -0.2756261023899033713},  0.4  }, ℂ₁{ф,-0.1059312120582799194, 0.9943734601807632245}, {ф,-0.563645946833390421,  0.8260164929456861316}}, "*");
                static_assert(     a_le == Arc{{{0.03                 ,  0.03                 },  0.03 }, ℂ₁{ф,-0.5636459468333897549, 0.8260164929456866867}, {ф, 0.2185328481759700181,-0.9758295928429812083}}, "*");
                static_assert( a_bottom == Arc{{{0.3656664547982899638, -1.468874254606819196 }, -1.506}, ℂ₁{ф,-0.2185328481759699903, 0.9758295928429810973}, {ф, 0.0631673008458991242, 0.9980029519514678205}}, "*");
        }
        
        cout    << "PIPE " << D_abs << 'x' << s_abs << '-' << b_abs << " r" << ler_abs << " f" << lef_abs << '\n'
                << std::setprecision(5) << std::fixed
                << a_top << a_lef << a_le << a_bottom << p_end << TE2
                ;
        return 0;
}
#endif

int main( int argc, const char *argv[])
{
        ++argv;
        --argc;

#ifndef NDEBUG
        return test5();
#endif

        static CE const char *param_name[] =
        { "диаметр трубы"
        , "толщина стенки трубы"
        , "хорда профиля"
        , "радиус передней кромки"
        , "радиус зализа над передней кромкой"
        };

        static double param[ size( param_name)] = { 160.0, 4.7, 50.0, 1.5, 18.0};

        double &D_abs   = param[0];  // диаметр трубы
        double &s_abs   = param[1];  // толщина стенки трубы
        double &b_abs   = param[2];  // хорда профиля
        double &ler_abs = param[3];  // радиус передней кромки
        double &lef_abs = param[4];  // радиус скругл. передней кромки

#ifdef NONENONENONENONENONEDEF //#ifndef NDEBUG
        if( argc == 3 )
        {
                std::size_t n = 0;
                if( !strcmp( argv[1], "test") && ( n = atoi( argv[2])) )
                {
                        static CE int (*checks[])() = { &test1, &test2, &test3, &test4, &test5 };
                        if( (n > 0) && (n <= size( checks)) )
                                return checks[ n-1 ]();
                }
                setlocale( LC_ALL, "" );
                cerr << "Ошибка тестов!\n";
                return 1;
        }
#endif

        if( (unsigned)argc < size( param_name) )
        {
                setlocale( LC_ALL, "" );

                cerr << "Неверное количество параметров, их должно быть " << size( param_name) 
                     << " а введено " << argc << "\n"
                        "\n"
                        "Примерные параметры запуска:\n"
                        "  pipefoil";

                for( unsigned i = 0; i < size( param); ++i)
                        cerr << ' ' << param[ i ];

                cerr << " >1.dat\nгде:\n";
                for( unsigned i = 0; i < size( param); ++i)
                        cerr << "  " << param[ i ] << "\t- " << param_name[ i ] << '\n';

                cerr << "  1.dat\t- результат: файл с точками профиля, можно \"продуть\" его в xflr5\n";
                return 1;
        }

        for( unsigned i = 0; i < size( param); ++i)
        {
                char *last_char;
                param[i] = strtod( argv[i], &last_char);
                if( *last_char)
                {
                        setlocale( LC_ALL, "" );
                        cerr << "Ошибка " << (i+1) << "-ого параметра \"" << param_name[i] << "\"\n"
                                "\"" << argv[i] << "\" вроде как \"" << param[i] << "\"\n";
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
        Circle c_top    = { c_bottom.o̅, R   };
        Circle c_lef    = tangent( c_top, c_le, -lef );
        ℂ     p_end    = c_bottom.tangent_point( TE2 );

        Arc a_top   ( c_top   , TE1  , c_lef.o̅   );
        Arc a_lef   ( c_lef   , a_top, c_le.o̅    );
        Arc a_le    ( c_le    , a_lef, c_bottom.o̅);
        Arc a_bottom( c_bottom, a_le , p_end     );
        
        cout    << "PIPE " << D_abs << 'x' << s_abs << '-' << b_abs << " r" << ler_abs << " f" << lef_abs << '\n'
                << std::setprecision(5) << std::fixed
                << a_top << a_lef << a_le << a_bottom << p_end << TE2
                ;
        return 0;
}
