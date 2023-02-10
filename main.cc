#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <clocale>
#include <cassert>
#include <type_traits>
#include <ciso646>

#pragma warning( disable: 4455)

#define CE      constexpr
#define EXPL    constexpr explicit
#define FRIEND  constexpr friend
#define STATIC  constexpr static
#define AUTO    constexpr auto
#define OP      operator
#define C       const
#define Self    ​ // этот тип
#define CSelf   const Self
#define self    *this
#define ㄥ      *cis
#define bit_cast _Bit_cast

using std::bit_cast;
using std::declval;
using std::ostream;
using std::false_type;
using std::true_type;
using std::enable_if_t;
using std::is_same;
using std::is_integral;
using std::is_floating_point;
using std::make_signed;
using std::make_unsigned;
using std::numeric_limits;
using std::remove_reference_t;
using std::cout;
using std::cerr;
using std::size;

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

template< typename... >
using void_t = void;

// преобразовать строку в int
#define S2I( str ) ((str)[0] + 256*(str)[1])

// помогалка для проверки существования оператора oper для типов L и R
template< typename L, int oper, typename R, typename = void>
struct is_oper: false_type {};

#define MK_IS_OPER( oper ) \
template< typename L, typename R> \
struct is_oper< L, S2I(#oper), R, void_t< decltype( declval< L>() oper declval< R>())> >: true_type {};

MK_IS_OPER(==)
MK_IS_OPER(+=) MK_IS_OPER(-=) MK_IS_OPER(*=) MK_IS_OPER(/=)

template< typename L, int oper, typename R>
CE bool is_oper_v = is_oper< L, oper, R>::value;

// L o= R → L o R
#define OP_B( o, L )                                            \
template< typename R                                            \
        , typename = enable_if_t< is_oper_v< L, S2I( #o"="), R>>\
        >                                                       \
CE L operator o ( L l, const R& r)                              \
{ l o##= r; return l; };

// коммутативный оператор с каким-то типом
// L o= R → R o L
#define OP_C( o, L )                                            \
template< typename R                                            \
        , typename = enable_if_t< not is_oper_v< R, S2I( #o"="), L>> \
        >                                                       \
AUTO operator o ( const R& r, const L& l)                       \
        -> remove_reference_t< decltype( L(l) o##= r)>          \
{                                                               \
        return L(l) o##= r;                                     \
};

// некоммутативный оператор с каким-то типом
// L o= R → R o L
#define OP_Ȼ( o, L )                                            \
template< typename R                                            \
        , typename = enable_if_t< not is_oper_v< R, S2I( #o"="), L>> \
        >                                                       \
AUTO operator o ( const R& r, const L& l)                       \
        -> remove_reference_t< decltype( L(r) o##= l)>          \
{                                                               \
        return L(r) o##= l;                                     \
};

// оба варианта коммутативного оператора с каким-то типом
#define OPS_C( oper, l_type ) OP_B( oper, l_type ) OP_C( oper, l_type )

// оба варианта некоммутативного оператора с каким-то типом
#define OPS_Ȼ( oper, l_type ) OP_B( oper, l_type ) OP_Ȼ( oper, l_type )


// создать арифметические операторы (+-*/) для l_type
#define OPS_ARITH( l_type)\
OPS_C( +, l_type )        \
OPS_Ȼ( -, l_type )        \
OPS_C( *, l_type )        \
OPS_Ȼ( /, l_type )

#define TINT template< typename INT  , typename = enable_if_t< is_integral      < INT  >::value>> constexpr
#define TFLT template< typename FLOAT, typename = enable_if_t< is_floating_point< FLOAT>::value>> constexpr

namespace impl
{
        template<class T> struct longer           {                        };
        template<       > struct longer< char    >{ using type = short;    };
        template<       > struct longer<  uint8_t>{ using type = uint16_t; };
        template<       > struct longer< uint16_t>{ using type = uint32_t; };
        template<       > struct longer< uint32_t>{ using type = uint64_t; };
        //template<       > struct longer< uint64_t>{ using type = uint128_t;};
        template<       > struct longer<   int8_t>{ using type =  int16_t; };
        template<       > struct longer<  int16_t>{ using type =  int32_t; };
        template<       > struct longer<  int32_t>{ using type =  int64_t; };
        //template<       > struct longer<  int64_t>{ using type =  int128_t;};

        template<class T> struct shorter           {                        };
        template<       > struct shorter< uint16_t>{ using type =  uint8_t; };
        template<       > struct shorter< uint32_t>{ using type = uint16_t; };
        template<       > struct shorter< uint64_t>{ using type = uint32_t; };
        //template<       > struct shorter<uint128_t>{ using type = uint64_t; };
        template<       > struct shorter<  int16_t>{ using type =   int8_t; };
        template<       > struct shorter<  int32_t>{ using type =  int16_t; };
        template<       > struct shorter<  int64_t>{ using type =  int32_t; };
        //template<       > struct shorter< int128_t>{ using type =  int64_t; };
}
template< typename T> using longer  = typename impl::longer < T>::type;
template< typename T> using shorter = typename impl::shorter< T>::type;
template< typename T> using sgned   = typename make_signed  < T>::type;
template< typename T> using unsgned = typename make_unsigned< T>::type;

template< typename L, typename R
        , typename = enable_if_t< not is_same< L, R>::value
                               && not is_oper_v< L, S2I("=="), R>
                               &&     is_oper_v< R, S2I("=="), L>
                                >
        >
AUTO OP == ( C L& l, C R& r) -> decltype( r == l )
{ return r == l; }

template< typename L, typename R>
AUTO OP != ( C L& l, C R& r) -> decltype( !(l == r) )
{ return !(l == r); }

template< typename T>
class	Near
{
        C T&	ref;	// ссылка на "родителя"
public:
EXPL	Near		( C T& r): ref( r) {}
CE	bool	OP ==	( C T& r) C { return ref.near		( r); }
CE	bool	OP <=	( C T& r) C { return ref.less_near	( r); }
CE	bool	OP >=	( C T& r) C { return ref.greater_near	( r); }

FRIEND	bool	OP <=	( C T& l, C Near& r) { return r >= l; }
FRIEND	bool	OP >=	( C T& l, C Near& r) { return r <= l; }
};

template< typename X
        , typename = decltype( declval< X>().near( declval< X>()))
        >
CE	Near<X>	OP ~	( C X& x) { return Near<X>( x);}

#pragma endregion


// constexpr'сный способ превращения T в signed T
template< typename T>
enable_if_t< is_integral< T>::value && std::is_unsigned< T>::value,
sgned< T> > CE sign_cast( T x)
{
        using sT = sgned< T>;
        CE sT min = numeric_limits< sT>::min();
        // этот if нужен для осчастливливания constexpr'а,
        // надеюсь, оптимизатор его выкинет
        if( x >= min )
                return sT( x - min) + min;

        return x;
}

template< typename T>
enable_if_t< is_integral< T>::value && std::is_signed< T>::value,
T> CE nabs( T x)
{
        // x ⩾ 0 → mask = 0
        // x < 0 → mask = -1
        using uT = unsgned< T>;
        const uT mask = x >> lastbit< T>;
        return mask - (uT(x) ^ mask);
};

C class Ф{} ф; // МАМОЙ КЛЯНУСЬ (флаг для конструкторов из сырых данных)

struct Turn;

#pragma region // Double{}
struct Double
{
        using Self = Double;
private:
        double	val;

public:
CE	Self		( double x): val( x) {}
CE OP	double		(         ) C { return  val;	}
CE	Self	OP -	(         ) C { return -val;	}
CE	double* OP &	(         )   { return &val;	}
CE C	double* OP &	(         ) C { return &val;	}
//CE OP	double&		(         )   { return  val;	}
//CE OP	C double&	(         ) C { return  val;	}
CE	bool    isfinite(         ) C { return raw_abs( val) < raw( infinity)	; }
CE      double  ²       (         ) C { return val * val			; }
CE	double	abs	(         ) C { return bit_cast< double>( raw_abs( val)); }
CE      int64_t	signbit	(         ) C { return raw(val) & (uint64_t(1) << 63)	; }
                                        // val = (x >= 0. ? val : -val)
CE	void	revsign	( double x)   { val = bit_cast< double>( raw(val) ^ raw(x) & (uint64_t(1) << 63)); }
CE      double	root    (         ) C {
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

CE	bool	less_near	( CSelf  &x) C { return sraw(val) <= sraw(x) + ε	;}
CE	bool	greater_near	( CSelf  &x) C { return sraw(val) >= sraw(x) - ε	;}
CE	bool	near		( CSelf  &x) C // почти равный
        {
                return	   nabs( int64_t( uint64_t(raw(val)) - uint64_t(raw(x)))) >= -ε
                        || about_0(val) && about_0(x)
                        ;
        }

private:
STATIC	int64_t	ε = 3;  // епсилон (в точках на числовой прямой)
STATIC	int64_t	raw	( double x) { return bit_cast< int64_t>( x)			;}
STATIC	int64_t	sraw	( double x) { return raw(x) ^ (uint64_t( raw(x) >> 63) >> 1)	;} // заполнить все 64 бита копиями знакового бита, и обнуляет старший бит
STATIC	int64_t	raw_abs	( double x) { return raw(x) & ~(uint64_t(1) << 63)		;}
STATIC	bool	about_0	( double x) { return raw_abs(x) < raw( min)			;}
                                                          //raw( 1e-30)			;}
                                                          //raw( 1e-14)			;}
                                                          //raw( numeric_limits< double>::epsilon())		;}
                                                          //raw(epsilon/(uint64_t(1)<<55));}
                                                          //0x3CC0000000000000ul            ;}
STATIC	double  root_rec( double r, double curr, double prev)
        {
                // sqrt Newton-Raphson
                // https://gist.github.com/alexshtf/eb5128b3e3e143187794
                if( curr == prev )
                        return curr;
                return root_rec( r, (curr + r/curr) / 2, curr);
        }
};

CE double root( const Double& r) { return r.root(); }

CE const Double Double::min             = numeric_limits< double>::min          (); //   10000000000000
CE const Double Double::max             = numeric_limits< double>::max          (); // 7FEFFFFFFFFFFFFF
CE const Double Double::lowest          = numeric_limits< double>::lowest       (); // FFEFFFFFFFFFFFFF
CE const Double Double::epsilon         = numeric_limits< double>::epsilon      (); // 3CB0000000000000 2.220446049250313081e-16
CE const Double Double::round_error     = numeric_limits< double>::round_error  (); // 3FE0000000000000 0.5
CE const Double Double::denorm_min      = numeric_limits< double>::denorm_min   (); //                1
CE const Double Double::infinity        = numeric_limits< double>::infinity     (); // 7FF0000000000000
CE const Double Double::quiet_NaN       = numeric_limits< double>::quiet_NaN    (); // 7FF8000000000000
CE const Double Double::signaling_NaN   = numeric_limits< double>::signaling_NaN(); // 7FF8000000000001

CE Double operator ""d (        long double x) { return x; }
CE Double operator ""d ( unsigned long long x) { return x; }
CE double operator ""⁰⁵(        long double x) { return root( x); }
CE double operator ""⁰⁵( unsigned long long x) { return root( x); }

CE void Double_test()
{
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

        static_assert( -0.0d ==~ 0.0d, "***");
        static_assert( !(a_m ==  a_p), "***");
        static_assert(   a_m !=  a_p , "***");
        static_assert(   a_m ==~ a_p , "***");
        static_assert( !(a_m <=  a_p), "***");
        static_assert( !(a_p >=  a_m), "***");
        static_assert(   a_m <=~ a_p , "***");
        static_assert(   a_m >=~ a_p , "***");
        static_assert(   a_p <=~ a_m , "***");
        static_assert(   a_p >=~ a_m , "***");

        static_assert(     isfinite( Double::max           ), "***");
        static_assert(     isfinite( Double::lowest        ), "***");
        static_assert( not isfinite( Double::infinity      ), "***");
        static_assert( not isfinite( Double::quiet_NaN     ), "***");
        static_assert( not isfinite( Double::signaling_NaN ), "***");

        //static_assert( about_0( vsv)	, "***");
        //static_assert( vsv ==~ 0d	, "***");
        static_assert( abs( -.125d) == .125, "***");
        static_assert( root( 152.41383936d ) ==~  12.3456d, "***");
        static_assert( root( ²(5678.2023d) ) == 5678.2023d, "***");
        static_assert( .122999999999999999 !=  .12300000000000001d, "***");
        static_assert( .122999999999999999 ==~ .12300000000000001d, "***");
        static_assert( ! 0.		, "***");
        static_assert(    double( vsv)	, "***");
        //static_assert( ! ~Double( vsv)	, "***");
        static_assert( ~1d == .1 * 10	, "***");
        static_assert( 1. == Double(.1 * 10), "***");
}

#pragma endregion

#pragma region // Angle{}

CE C double π = 3.14159265358979323846;
CE C double operator ""π  (        long double a) { return π * a;                }
CE C double operator ""π  ( unsigned long long a) { return π * a;                }

// Угол
struct Angle
{
#pragma warning( push)
#pragma warning( disable: 4146)

        friend struct Turn; 

        using Self = Angle;
        using Val  = unsigned int;

private:
        Val     val;
STATIC  Val     semiturn = Val(1) << lastbit< Val>;
STATIC  sgned< Self::Val> ε = 8;

TINT    Self         ( INT        v ): val( v) {}

public:
EXPL    Self         ( C Turn&      );
EXPL    Self         ( double radian): val(   sgned< longer< Val>>(  radian * (semiturn/π))
                                            & Val( -1 )
                                          ) {}
EXPL OP double       (              ) C { return val * (π/semiturn)             ; }
CE      bool    OP ==( CSelf&     r ) C { return val == r.val                   ; }
CE	bool	near ( CSelf&     r ) C { return nabs( sign_cast( val - r.val)) > -ε; } // почти равный
CE      Self    OP - (              ) C { return -val                           ; }
CE      Val     OP / ( CSelf&     r ) C { return val / r.val                    ; }
CE      Self    OP / ( double     r ) C { return sgned< Val>( double( val) / r ); }
TINT    Self    OP / ( INT        r ) C { return           sign_cast( val) / r  ; }
CE      Self    OP / ( unsigned   r ) C { return                      val  / r  ; }

CE      Self&   OP +=( CSelf&     r )   { val += r.val; return self; }
CE      Self&   OP -=( CSelf&     r )   { val -= r.val; return self; }
TINT    Self&   OP *=( INT        r )   { val *= r    ; return self; }
TINT    Self&   OP /=( INT        r )   { val /= r    ; return self; }
CE      Self&   OP *=( double     r )   { val *= r    ; return self; }
        void    print( ostream&   os) C
        {
                os << ( val * (180./semiturn)) << '°';
        };

FRIEND  Self    OP ""ᵒ( unsigned long long);
FRIEND  Self    OP ""ᵒ( long double       );
#pragma warning( pop)
};

OPS_C( +, Angle)
OPS_Ȼ( -, Angle)
OPS_C( *, Angle)

CE Angle OP ""ᵒ ( unsigned long long x)
{
        return Angle::Val( x * Angle::semiturn / 180
                         & Angle::Val( -1)
                         );
};
CE Angle OP ""ᵒ (        long double x)
{
        return Angle::Val( longer< Angle::Val>( long double( Angle::semiturn) / 180 * x)
                         & Angle::Val( -1)
                         );
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
        friend struct Angle; 

        using Self = Turn;
        using Val  = sgned< longer< Angle::Val>>;

private:
STATIC  Val     one_turn = Val(1) << (lastbit< unsgned< Angle::Val>> + 1);
        Val     val;
CE      Self            ( Ф, Val v): val( v           ) {}
public:
CE      Self            ( Angle  x): val( x.val       ) {}
CE      Self            ( double x): val( x * one_turn) {}
TINT    Self            ( INT    x): val( x * one_turn) {}
CE   OP double          (         ) C { return double(val) / one_turn; }
TINT OP INT             (         ) C { return        val  / one_turn; }

CE      Self    OP -    (         ) C { return { ф, -val }; }
CE      Val     OP /    ( CSelf& r) C { return val / r.val; }
CE      Self    OP /    ( double r) C { return double( self) / r; }

CE      Self&   OP +=   ( CSelf& r)   { val += r.val; return self; }
CE      Self&   OP -=   ( CSelf& r)   { val -= r.val; return self; }
TINT    Self&   OP *=   ( INT    r)   { val *= r    ; return self; }
TINT    Self&   OP /=   ( INT    r)   { val /= r    ; return self; }
CE      Self&   OP *=   ( double r)   { val *= r    ; return self; }
        void    print   ( ostream& os) C
        {
                os << double( val) << " turn";
        };
};

OPS_C( +, Turn )
OPS_Ȼ( -, Turn )
OPS_C( *, Turn )

CE Turn OP ""turn ( unsigned long long x) { return         x ; };
CE Turn OP ""turn ( long double        x) { return double( x); };
AUTO    OP ""τ    ( unsigned long long x) { return OP ""turn ( x); };
AUTO    OP ""τ    ( long double        x) { return OP ""turn ( x); };

CE Angle::Angle( Turn C &x): val( x.val)
{}

void test_Turn()
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

#pragma region // Комплексные Числа ℂ{}, ℂ₁{}

C class 𝐈{} 𝐢; // самая мнимая единица на свете (настолько мнимая, что даже память не занимает!)

struct ℂ;
struct ℂ₁;

// Комплексное Число, но есть нюанс...
struct ℂ
{
        using Self = ℂ;

protected:
        Double  r, i;
CE      void    swap() { double r1 = r; r = i; i = r1; }

public:
CE      Self( 𝐈                 ): r( 0 ), i( 1 ) {}
EXPL    Self( double s          ): r( s ), i( 0 ) {}
CE      Self( double r, double i): r( r ), i( i ) {}

CE      double  abs² () C { return  r*r + i*i      ;} // абсолютная величина в квадрате
CE      double  abs  () C { return ::root( abs²()) ;} // абсолютная величина
CE      Self    OP - () C { return { -r, -i }      ;} // унарный минус
CE      Self    conj () C { return {  r, -i }      ;} // Сопряжённое (conjugate) число
CE      Self    recip() C { return conj() /= abs²();} // 1/z - обратная величина (reciprocal, multiplicative inverse)
CE      Self    ⅟   () C { return recip()         ;} // 1/z - обратная величина (reciprocal, multiplicative inverse)
CE      double  re   () C { return r               ;} // действительная часть
CE      double  im   () C { return i               ;} // мнимая часть
CE      double  ℜ   () C { return r               ;} // действительная часть
CE      double  ℑ    () C { return i               ;} // мнимая часть
CE      Angle   arg  () C
        {
                if( r > +0. ) return Angle( atan( i / r ));
                if( r < -0. ) return Angle( atan( i / r ) + π);

                if( i > +0. ) return Angle( π/2 );
                if( i < -0. ) return Angle( 3π/2);

                return Angle( NAN);
        }
CE      Self    root () C // корень квадратный
        {
                double abs1 = abs();
                Self res = { ::root( (abs1 + r)/2 ) 
                           , ::root( (abs1 - r)/2 )
                           };
                res.i.revsign( i);
                return res;
        }

        FN( abs²); FN( abs); FN( conj); FN( recip); FN( ⅟); FN( re); FN( im); FN( ℜ); FN( ℑ); FN( arg); FN( root);

CE      bool    OP ==( CSelf& z) C { return r ==  z.r && i ==  z.i;}
CE	bool	near ( CSelf& z) C { return r ==~ z.r && i ==~ z.i;} // почти равный

CE      Self&   OP +=( CSelf& z)   { r+=z.r; i+=z.i; return self; }
CE      Self&   OP -=( CSelf& z)   { r-=z.r; i-=z.i; return self; }
/*
CE      Self    OP * ( CSelf& z) C
        {
                double a = (r + i) * (z.r + z.i);
                double b = r * z.r;
                double c = i * z.i;

                return  { b - c
                        , a - b - c
                        };
        }
CE      Self&   OP *=( CSelf& z)   { return (self  = self * z );  }
*/
CE      Self&   OP *=( CSelf& z)
        {
                Self r90 = {-i, r}; // поворот на 90° *this;
                r90     *= z.i;
                self    *= z.r;
                self    += r90;
                return self;
        }
//CE      Self&   OP /=( CSelf& z)   { return (self *= z.recip()) ; } // Деление КЧ на другое КЧ
CE      Self    OP / ( CSelf& z) C // Деление КЧ на другое КЧ
        {
                Self res = { r*z.r + i*z.i
                           , i*z.r - r*z.i
                           };
                return res /= z.abs²();
        }
CE      Self&   OP /=( CSelf& z)   { return (self  = self / z );  }
CE      Self&   OP /=( C ℂ₁& n);                                   // Деление КЧ на единичное КЧ дает КЧ

CE      Self&   OP +=( double s)   { r += s;         return self; }
CE      Self&   OP -=( double s)   { r -= s;         return self; }
CE      Self&   OP *=( double s)   { r *= s; i *= s; return self; } // Умножение на скаляр
CE      Self&   OP /=( double s)   { r /= s; i /= s; return self; } // Деление на скаляр

CE      Self&   OP *=( C 𝐈&    )   { swap(); r = -r; return self; } // Умножение на мнимую
CE      Self&   OP /=( C 𝐈&    )   { swap(); i = -i; return self; } // Деление на мнимую

        // на время представим, что комплексные числа это векторы...
        // Скалярное произведение (Inner product)
CE      double  OP , ( CSelf& v) C { return  r*v.r + i*v.i; }
CE      double  OP ^ ( CSelf& v) C { return  r*v.i - i*v.r; } // Псевдоскалярное, Векторное, косое произведение (Outer, cross product)

        void    print( ostream& os) C
        {
                static Self last = { NAN, NAN };
                if( self !=~ last )
                {
                        os << std::setw(8) << r << std::setw(12) << i << '\n';
                        last = self;
                }
        };
};

OPS_ARITH( ℂ )

// Единичное Комплексное Число (Unit Complex Number). ℂ₁ = {𝑧 ∈ ℂ: |𝑧| = 1}
struct ℂ₁: public ℂ
{
        using Self  = ℂ₁;
        using Super = ℂ;

#ifndef NDEBUG
CE      Self( Ф, double r, double i): Super( r, i) {}
#endif
private:
CE      Self( double r, double i): Super( r, i) {}

public:
CE      Self( 𝐈          ): Super( 0., 1.                               ) {}
EXPL    Self( C Super& z ): Super( z / z.abs()                          ) {}
EXPL    Self(Ф,C Super&z ): Super( z                                    ) {} // мамой клянусь - |z| == 1
EXPL    Self( double sin ): Super( ::root( 1. - sin*sin), sin           ) {}        
CE      Self( Angle    𝜑 ): Super( cos( double( 𝜑)), sin( double( 𝜑))  ) {} // фазор угла 𝜑 (единичное КЧ с аргументом 𝜑)

CE      double  abs² () C { return     1.   ; } // абсолютная величина в квадрате = 1
CE      double  abs  () C { return     1.   ; } // абсолютная величина = 1, этож единичное число ;)
CE      Self    OP - () C { return {-r, -i }; } // унарный минус (ЕКЧ = -ЕКЧ)
CE      Self    conj () C { return { r, -i }; } // Сопряжённое (conjugate) число
CE      Self    recip() C { return { r, -i }; } // 1/z - обратная величина (reciprocal, multiplicative inverse)
CE      Self    ⅟   () C { return recip()  ; } // 1/z - обратная величина (reciprocal, multiplicative inverse)
CE      double  ψarg () C { double x1 = r + 1.; return i>=0 ? -x1 : x1;} // псевдоугол (для сравнений)
CE      Self    root () C
        {
                double i0 = ::root( (1. - r)/2 );
                if( i < 0 )
                        i0 = -i0;
                return { ::root( (1. + r)/2 )
                       , i0
                       };
        }
        FN( abs²); FN( abs); FN( conj); FN( recip); FN( ⅟); FN( ψarg); FN( root);

CE      bool    OP > ( CSelf& z) C { return ψarg() > z.ψarg();}
CE      bool    OP < ( CSelf& z) C { return ψarg() < z.ψarg();}

CE      Self&   OP *=( CSelf& z)   { Super::OP *=( z); return self; } // ЕКЧ*ЕКЧ = ЕКЧ
CE      Self&   OP /=( CSelf& z)   { return (self *= z.conj())    ; } // ЕКЧ/ЕКЧ = ЕКЧ, при этом само деление вырождено
CE      ℂ&     OP /=( C ℂ& z )   { return ℂ::OP /=( z)          ; } // КЧ/ЕКЧ = КЧ

CE      Self&   OP *=( C 𝐈&    )   { swap(); r = -r; return self; } // Умножение на мнимую
CE      Self&   OP /=( C 𝐈&    )   { swap(); i = -i; return self; } // Деление на мнимую

static  CSelf   𝟏;
static  CSelf   𝐢;
};

OPS_C( *, ℂ₁)
OPS_Ȼ( /, ℂ₁)

CE ℂ& ℂ::OP /=( C ℂ₁& n) 
{       return (self *= n.recip()); }

// комплексная единица (единичная, естественно)
CE C ℂ₁ ℂ₁::𝟏 = { 1., 0. };
// комплексная мнимая единица
CE C ℂ₁ ℂ₁::𝐢 = { 0., 1. };

// единица
CE C ℂ₁ 𝟏 = ℂ₁::𝟏;

CE ℂ OP ""𝐢 ( unsigned long long a) { return { 0., double(a)}; }
CE ℂ OP ""𝐢 (        long double a) { return { 0., double(a)}; }

AUTO OP -( 𝐈 ) { return      -ℂ₁::𝐢 ; }
AUTO root( 𝐈 ) { return root( ℂ₁::𝐢); }

void   ℂ_test()
{
        static_assert( root( -𝟏) == 𝐢                                   , "***");
        static_assert( root( -𝐢) == conj( root( 𝐢))                     , "***");
        static_assert( (5 + 𝐢)*(7 - 6𝐢) / (3 + 𝐢)          == (10 - 11𝐢), "***");
        static_assert( (4 + 𝐢)*(5 + 3𝐢) + (3 + 𝐢)*(3 - 2𝐢) == (28 + 14𝐢), "***");

        // это просто должно скомпилиться
        AUTO    a0 = ℂ₁( 0.1 ) * ℂ( 0, 1 );
        CE ℂ₁ u0 = 𝐢;
        CE ℂ₁ u1 = 𝟏 / 𝐢;
        CE ℂ  u2 = ℂ₁( 0.1 ) *= 𝐢;
        CE ℂ₁ u3 = ℂ₁( 0.1 );
        CE ℂ₁ u4 = u3 * ℂ₁( 0.1 ) / 𝐢;
        CE ℂ₁ u5 = 𝟏/ℂ₁( 0.1 );
        CE ℂ₁ u6 = 𝟏*ℂ₁( 0.1 );
        CE ℂ  u7 = ℂ( 0, 1 ) * ℂ( 0.995, 0.1 ) ;
        CE ℂ₁ u8 = 𝐢 * ℂ₁( 0.1 );
        CE ℂ  u9 = u7 / u8;
        CE ℂ  ua = u8 / u7;
        CE ℂ  ub = u7 + u9;
        CE ℂ  uc = u7 * 9;

        CE ℂ₁ a = ℂ₁( 0.1);
        static_assert( a / ℂ(3, -2) == ℂ( 0.2142278701015276898, 0.1761519134010184617), "***");
        static_assert(      a > 𝟏 && !(      a < 𝟏), "***");
        static_assert( conj(a)> 𝟏 && !( conj(a)< 𝟏), "***");
        static_assert( conj(a)> a && !( conj(a)< a), "***");
        static_assert( a*𝐢*𝐢*𝐢 > a && !( a*𝐢*𝐢*𝐢 < a), "***");

        static_assert( arg(0 + 𝐢) ==  90ᵒ, "***");
        static_assert( arg(0 - 𝐢) == 270ᵒ, "***");
}
#pragma endregion

struct Line
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
CE      Vertical( double x0 ): Line( ℂ₁::𝟏, x0) {};
};
struct Horizontal: public Line
{
        // Горизонталь пересекающая ось 𝑌 в точке y0
CE      Horizontal( double y0 ): Line( ℂ₁::𝐢, y0) {};
};
void   Line_test()
{
        CE Vertical l_test( 2 );
        static_assert(       l_test.n̂        == 𝟏, "***");
        static_assert(       l_test.n̂ * 𝐢    == 𝐢, "***");
        static_assert( dist( l_test, {3, 1}) == 1, "***");
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
CE      Segment OP - () C { return Segment( this->p̅2, this->p̅1 ); }

        void print( ostream& os) C
        {
                os << p̅1 << p̅2;
        };
};

struct Circle
{
        ℂ     o̅; // центр окружности
        Double R; // радиус окружности

CE      Circle	     ( C ℂ& center, double radius = 0. )
        : R( radius), o̅( center )
        {};

CE      bool	OP ==( C Circle& c ) C { return R ==  c.R && o̅ ==  c.o̅; }
CE	bool	near ( C Circle& c ) C { return R ==~ c.R && o̅ ==~ c.o̅; } // почти равный
CE      Circle	OP - (             ) C
        {
                Circle a = self;
                a.R = -a.R;
                return a;
        }

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
                ℂ₁ n̂ = cathet2( -(c2.o̅ - o̅), -(c2.R - R)); // почему минусы ???
                return Line( n̂, (n̂, o̅) - R );
        }

        // точка касания касательной к окружности исходящей из точки p
CE      ℂ tangent_point( C ℂ& p̅ ) C
        {
                return o̅ - R * cathet2( o̅ - p̅, R);
        };

        void print( ostream &os) C
        {
                CE C int segs = 40;

                //CE C ℂ₁ m̂¹⁰ = cis( 2π/segs);
                CE C ℂ₁ m̂( 2π/segs/10 ); // типа 𝛼 ≈ sin 𝛼
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
void   Circle_test()
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

        static_assert( c_test.tangent( Circle( {   2, 0}, 1) ) ==  Line( ℂ₁::𝟏             ,  3), "***");
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
                R.revsign( prev.R * R /*TODO R?*/ * (prev.n̂₂, n̂₁));
        }

CE      bool	OP ==( C Arc& a) C { return  Circle( self) ==  Circle( a) && n̂₁ ==  a.n̂₁ && n̂₂ ==  a.n̂₂; }
CE	bool	near ( C Arc& a) C { return  Circle( self) ==~ Circle( a) && n̂₁ ==~ a.n̂₁ && n̂₂ ==~ a.n̂₂; } // почти равный

        void	print( ostream &os) C
        {
                CE C int segs = 160;

                //CE C ℂ₁ m̂¹⁰ = cis( 2π/segs);
                CE C ℂ₁ m̂( 2π/segs/10 ); // типа 𝛼 ≈ sin 𝛼
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
        CE Circle  с_top = tangent( c_lead, TE1, R );
        CE Circle  c_bottom( с_top.o̅, R-s );
        CE Segment s_end = tangent_segment( c_bottom, TE2 );
        CE Circle  c_l1( Vertical( lr1) ^ Circle( с_top.o̅, c_bottom.R+lr1), lr1 );
        CE Circle  c_l2 = tangent( с_top, c_l1, -lr2 );

        CE Arc a_top = chain_arc( TE1      , с_top , c_l2     );
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
        CE Circle  с_top = tangent( c_lead, TE1, R );
        CE Circle  c_bottom( с_top.o̅, R-s );
        CE Segment s_start = tangent_segment( c_lead, -c_bottom, minus );
        CE Segment s_end = tangent_segment( c_bottom, TE2 );

        CE Arc a_top = chain_arc( TE1, с_top, c_lead      );
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
        CE Circle  с_top( c_bottom.o̅, R );
        CE Circle  c_l2  = tangent( c_l1, с_top, -lr2, minus );
        CE Segment s_end = tangent_segment( c_l1, TE2, minus );

        CE Arc a_top = chain_arc( TE1     , с_top, c_l2 );
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
        CE Circle  с_top( c_bottom.o̅, R );
        CE Circle  c_l2  = tangent( c_l1, с_top, -lr2, minus );
        CE Segment s_end0 = tangent_segment( c_l1, TE2, minus );

        CE Arc a_top = chain_arc( TE1     , с_top, c_l2 );
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
        CE Circle с_top    = { c_bottom.o̅, R };
        CE Circle c_lef    = tangent( с_top, c_le, -lef );
        CE ℂ     p_end    = c_bottom.tangent_point( TE2 );

        CE Arc a_top   ( с_top   , TE1  , c_lef.o̅   );
        CE Arc a_lef   ( c_lef   , a_top, c_le.o̅    );
        CE Arc a_le    ( c_le    , a_lef, c_bottom.o̅);
        CE Arc a_bottom( c_bottom, a_le , p_end     );
        
        // static_asserts
        {
                static_assert(    a_top == Arc{{{0.3656664547982899638, -1.468874254606819196 },  1.6  }, {ф, 0.3964584657510686894, 0.9180526591292618166}, {ф,-0.10593121205827985,    0.9943734601807630025}}, "***");
                static_assert(    a_lef == Arc{{{0.2385490003283541327, -0.2756261023899033713},  0.4  }, {ф,-0.1059312120582799194, 0.9943734601807632245}, {ф,-0.563645946833390421,   0.8260164929456861316}}, "***");
                static_assert(     a_le == Arc{{{0.03                 ,  0.03                 },  0.03 }, {ф,-0.5636459468333897549, 0.8260164929456866867}, {ф, 0.2185328481759700181, -0.9758295928429812083}}, "***");
                static_assert( a_bottom == Arc{{{0.3656664547982899638, -1.468874254606819196 }, -1.506}, {ф,-0.2185328481759699903, 0.9758295928429810973}, {ф, 0.0631673008458991242,  0.9980029519514678205}}, "***");
        }
        
        cout    << "PIPE " << D_abs << 'x' << s_abs << '-' << b_abs << " r" << ler_abs << " f" << lef_abs << '\n'
                << std::setprecision(5) << std::fixed
                << a_top << a_lef << a_le << a_bottom << p_end << TE2
                ;

        return 0;
}
#endif

int main( unsigned argc, C char *argv[])
{
        ++argv;
        --argc;

#ifndef NDEBUG
        return test5();
#endif

        static CE char *param_name[] =
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
/*
 #i fn def NDEBUG
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
 #e nd if
*/
        if( argc < size( param_name) )
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
        Circle с_top    = { c_bottom.o̅, R   };
        Circle c_lef    = tangent( с_top, c_le, -lef );
        ℂ     p_end    = c_bottom.tangent_point( TE2 );

        Arc a_top   ( с_top   , TE1  , c_lef.o̅   );
        Arc a_lef   ( c_lef   , a_top, c_le.o̅    );
        Arc a_le    ( c_le    , a_lef, c_bottom.o̅);
        Arc a_bottom( c_bottom, a_le , p_end     );
        
        cout    << "PIPE " << D_abs << 'x' << s_abs << '-' << b_abs << " r" << ler_abs << " f" << lef_abs << '\n'
                << std::setprecision(5) << std::fixed
                << a_top << a_lef << a_le << a_bottom << p_end << TE2;

        return 0;
}
