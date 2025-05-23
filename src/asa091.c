# include <math.h>
# include <stdlib.h>
# include <stdio.h>
# include <time.h>

# include "asa091.h"

/******************************************************************************/

double alnorm ( double x, int upper )

/******************************************************************************/
/*
  Purpose:

    ALNORM computes the cumulative density of the standard normal distribution.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 November 2010

  Author:

    Original FORTRAN77 version by David Hill.
    C version by John Burkardt.

  Reference:

    David Hill,
    Algorithm AS 66:
    The Normal Integral,
    Applied Statistics,
    Volume 22, Number 3, 1973, pages 424-427.

  Parameters:

    Input, double X, is one endpoint of the semi-infinite interval
    over which the integration takes place.

    Input, int UPPER, determines whether the upper or lower
    interval is to be integrated:
    1  => integrate from X to + Infinity;
    0 => integrate from - Infinity to X.

    Output, double ALNORM, the integral of the standard normal
    distribution over the desired interval.
*/
{
  double a1 = 5.75885480458;
  double a2 = 2.62433121679;
  double a3 = 5.92885724438;
  double b1 = -29.8213557807;
  double b2 = 48.6959930692;
  double c1 = -0.000000038052;
  double c2 = 0.000398064794;
  double c3 = -0.151679116635;
  double c4 = 4.8385912808;
  double c5 = 0.742380924027;
  double c6 = 3.99019417011;
  double con = 1.28;
  double d1 = 1.00000615302;
  double d2 = 1.98615381364;
  double d3 = 5.29330324926;
  double d4 = -15.1508972451;
  double d5 = 30.789933034;
  double ltone = 7.0;
  double p = 0.398942280444;
  double q = 0.39990348504;
  double r = 0.398942280385;
  int up;
  double utzero = 18.66;
  double value;
  double y;
  double z;

  up = upper;
  z = x;

  if ( z < 0.0 )
  {
    up = !up;
    z = - z;
  }

  if ( ltone < z && ( ( !up ) || utzero < z ) )
  {
    if ( up )
    {
      value = 0.0;
    }
    else
    {
      value = 1.0;
    }
    return value;
  }

  y = 0.5 * z * z;

  if ( z <= con )
  {
    value = 0.5 - z * ( p - q * y 
      / ( y + a1 + b1 
      / ( y + a2 + b2 
      / ( y + a3 ))));
  }
  else
  {
    value = r * exp ( - y ) 
      / ( z + c1 + d1 
      / ( z + c2 + d2 
      / ( z + c3 + d3 
      / ( z + c4 + d4 
      / ( z + c5 + d5 
      / ( z + c6 ))))));
  }

  if ( !up )
  {
    value = 1.0 - value;
  }

  return value;
}
/******************************************************************************/

void chi_square_cdf_values ( int *n_data, int *a, double *x, double *fx )

/******************************************************************************/
/*
  Purpose:

    CHI_SQUARE_CDF_VALUES returns some values of the Chi-Square CDF.

  Discussion:

    In Mathematica, the function can be evaluated by:

      Needs["Statistics`ContinuousDistributions`"]
      dist = ChiSquareDistribution [ df ]
      CDF [ dist, x ]

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 August 2004

  Author:

    John Burkardt

  Reference:

    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.

    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.

  Parameters:

    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.

    Output, int *A, the parameter of the function.

    Output, double *X, the argument of the function.

    Output, double *FX, the value of the function.
*/
{
# define N_MAX 21

  int a_vec[N_MAX] = { 
     1,  2,  1,  2, 
     1,  2,  3,  4, 
     1,  2,  3,  4, 
     5,  3,  3,  3, 
     3,  3, 10, 10, 
    10 };

  double fx_vec[N_MAX] = { 
     0.7965567455405796E-01,  
     0.4987520807317687E-02,   
     0.1124629160182849E+00,  
     0.9950166250831946E-02,  
     0.4729107431344619E+00,   
     0.1812692469220181E+00,   
     0.5975750516063926E-01,   
     0.1752309630642177E-01,   
     0.6826894921370859E+00,   
     0.3934693402873666E+00,   
     0.1987480430987992E+00,   
     0.9020401043104986E-01,   
     0.3743422675270363E-01,   
     0.4275932955291202E+00,   
     0.6083748237289110E+00,   
     0.7385358700508894E+00,   
     0.8282028557032669E+00,   
     0.8883897749052874E+00,   
     0.1721156299558408E-03,   
     0.3659846827343712E-02,   
     0.1857593622214067E-01 };

  double x_vec[N_MAX] = { 
     0.01E+00,   
     0.01E+00,    
     0.02E+00,   
     0.02E+00,   
     0.40E+00,   
     0.40E+00,   
     0.40E+00,   
     0.40E+00,   
     1.00E+00,   
     1.00E+00,   
     1.00E+00,   
     1.00E+00,   
     1.00E+00,   
     2.00E+00,   
     3.00E+00,   
     4.00E+00,   
     5.00E+00,   
     6.00E+00,   
     1.00E+00,   
     2.00E+00,   
     3.00E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *a = 0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *a = a_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
/******************************************************************************/

double gammad ( double x, double p, int *ifault )

/******************************************************************************/
/*
  Purpose:

    GAMMAD computes the Incomplete Gamma Integral

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 November 2010

  Author:

    Original FORTRAN77 version by B Shea.
    C version by John Burkardt.

  Reference:

    B Shea,
    Algorithm AS 239:
    Chi-squared and Incomplete Gamma Integral,
    Applied Statistics,
    Volume 37, Number 3, 1988, pages 466-473.

  Parameters:

    Input, double X, P, the parameters of the incomplete 
    gamma ratio.  0 <= X, and 0 < P.

    Output, int IFAULT, error flag.
    0, no error.
    1, X < 0 or P <= 0.

    Output, double GAMMAD, the value of the incomplete 
    Gamma integral.
*/
{
  double a;
  double an;
  double arg;
  double b;
  double c;
  double elimit = - 88.0;
  double oflo = 1.0E+37;
  double plimit = 1000.0;
  double pn1;
  double pn2;
  double pn3;
  double pn4;
  double pn5;
  double pn6;
  double rn;
  double tol = 1.0E-14;
  int upper;
  double value;
  double xbig = 1.0E+08;

  value = 0.0;
/*
  Check the input.
*/
  if ( x < 0.0 )
  {
    *ifault = 1;
    return value;
  }

  if ( p <= 0.0 )
  {
    *ifault = 1;
    return value;
  }

  *ifault = 0;

  if ( x == 0.0 )
  {
    value = 0.0;
    return value;
  }
/*
  If P is large, use a normal approximation.
*/
  if ( plimit < p )
  {
    pn1 = 3.0 * sqrt ( p ) * ( pow ( x / p, 1.0 / 3.0 ) 
      + 1.0 / ( 9.0 * p ) - 1.0 );

    upper = 0;
    value = alnorm ( pn1, upper );
    return value;
  }
/*
  If X is large set value = 1.
*/
  if ( xbig < x )
  {
    value = 1.0;
    return value;
  }
/*
  Use Pearson's series expansion.
  (Note that P is not large enough to force overflow in ALOGAM).
  No need to test IFAULT on exit since P > 0.
*/
  if ( x <= 1.0 || x < p )
  {
    arg = p * log ( x ) - x - lgamma ( p + 1.0 );
    c = 1.0;
    value = 1.0;
    a = p;

    for ( ; ; )
    {
      a = a + 1.0;
      c = c * x / a;
      value = value + c;

      if ( c <= tol )
      {
        break;
      }
    }

    arg = arg + log ( value );

    if ( elimit <= arg )
    {
      value = exp ( arg );
    }
    else
    {
      value = 0.0;
    }
  }
/*
  Use a continued fraction expansion.
*/
  else 
  {
    arg = p * log ( x ) - x - lgamma ( p );
    a = 1.0 - p;
    b = a + x + 1.0;
    c = 0.0;
    pn1 = 1.0;
    pn2 = x;
    pn3 = x + 1.0;
    pn4 = x * b;
    value = pn3 / pn4;

    for ( ; ; )
    {
      a = a + 1.0;
      b = b + 2.0;
      c = c + 1.0;
      an = a * c;
      pn5 = b * pn3 - an * pn1;
      pn6 = b * pn4 - an * pn2;

      if ( pn6 != 0.0 )
      {
        rn = pn5 / pn6;

        if ( fabs ( value - rn ) <= r8_min ( tol, tol * rn ) )
        {
          break;
        }
        value = rn;
      }

      pn1 = pn3;
      pn2 = pn4;
      pn3 = pn5;
      pn4 = pn6;
/*
  Re-scale terms in continued fraction if terms are large.
*/
      if ( oflo <= fabs ( pn5 ) )
      {
        pn1 = pn1 / oflo;
        pn2 = pn2 / oflo;
        pn3 = pn3 / oflo;
        pn4 = pn4 / oflo;
      }
    }

    arg = arg + log ( value );

    if ( elimit <= arg )
    {
      value = 1.0 - exp ( arg );
    }
    else
    {
      value = 1.0;
    }
  }

  return value;
}
/******************************************************************************/

double ppchi2 ( double p, double v, double g, int *ifault )

/******************************************************************************/
/*
  Purpose:

    PPCHI2 evaluates the percentage points of the Chi-squared PDF.

  Discussion

    Incorporates the suggested changes in AS R85 (vol.40(1),
    pages 233-5, 1991) which should eliminate the need for the limited
    range for P, though these limits have not been removed
    from the routine.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 June 2013

  Author:

    Original FORTRAN77 version by Donald Best, DE Roberts.
    C version by John Burkardt.

  Reference:

    Donald Best, DE Roberts,
    Algorithm AS 91:
    The Percentage Points of the Chi-Squared Distribution,
    Applied Statistics,
    Volume 24, Number 3, 1975, pages 385-390.

  Parameters:

    Input, double P,  value of the chi-squared cumulative
    probability density function.
    0.000002 <= P <= 0.999998.

    Input, double V, the parameter of the chi-squared probability
    density function.
    0 < V.

    Input, double G, the value of log ( Gamma ( V / 2 ) ).

    Output, int *IFAULT, is nonzero if an error occurred.
    0, no error.
    1, P is outside the legal range.
    2, V is not positive.
    3, an error occurred in GAMMAD.
    4, the result is probably as accurate as the machine will allow.

    Output, double PPCHI2, the value of the chi-squared random
    deviate with the property that the probability that a chi-squared random
    deviate with parameter V is less than or equal to PPCHI2 is P.
*/
{
  double a;
  double aa = 0.6931471806;
  double b;
  double c;
  double c1 = 0.01;
  double c2 = 0.222222;
  double c3 = 0.32;
  double c4 = 0.4;
  double c5 = 1.24;
  double c6 = 2.2;
  double c7 = 4.67;
  double c8 = 6.66;
  double c9 = 6.73;
  double c10 = 13.32;
  double c11 = 60.0;
  double c12 = 70.0;
  double c13 = 84.0;
  double c14 = 105.0;
  double c15 = 120.0;
  double c16 = 127.0;
  double c17 = 140.0;
  double c18 = 175.0;
  double c19 = 210.0;
  double c20 = 252.0;
  double c21 = 264.0;
  double c22 = 294.0;
  double c23 = 346.0;
  double c24 = 420.0;
  double c25 = 462.0;
  double c26 = 606.0;
  double c27 = 672.0;
  double c28 = 707.0;
  double c29 = 735.0;
  double c30 = 889.0;
  double c31 = 932.0;
  double c32 = 966.0;
  double c33 = 1141.0;
  double c34 = 1182.0;
  double c35 = 1278.0;
  double c36 = 1740.0;
  double c37 = 2520.0;
  double c38 = 5040.0;
  double ch;
  double e = 0.5E-06;
  int i;
  int if1;
  int maxit = 20;
  double pmax = 0.999998;
  double pmin = 0.000002;
  double p1;
  double p2;
  double q;
  double s1;
  double s2;
  double s3;
  double s4;
  double s5;
  double s6;
  double t;
  double value;
  double x;
  double xx;
/*
  Test arguments and initialize.
*/
  value = - 1.0;

  if ( p < pmin || pmax < p )
  {
    *ifault = 1;
    return value;
  }

  if ( v <= 0.0 )
  {
    *ifault = 2;
    return value;
  }

  *ifault = 0;
  xx = 0.5 * v;
  c = xx - 1.0;

/*
  Starting approximation for small chi-squared
*/
  if ( v < - c5 * log ( p ) )
  {
    ch = pow ( p * xx * exp ( g + xx * aa ), 1.0 / xx );

    if ( ch < e )
    {
      value = ch;
      return value;
    }
  }
/*
  Starting approximation for V less than or equal to 0.32
*/
  else if ( v <= c3 )
  {
    ch = c4;
    a = log ( 1.0 - p );

    for ( ; ; )
    {
      q = ch;
      p1 = 1.0 + ch * ( c7 + ch );
      p2 = ch * (c9 + ch * ( c8 + ch ) );

      t = - 0.5 + (c7 + 2.0 * ch ) / p1 - ( c9 + ch * ( c10 + 
        3.0 * ch ) ) / p2;

      ch = ch - ( 1.0 - exp ( a + g + 0.5 * ch + c * aa ) * p2 / p1) / t;

      if ( fabs ( q / ch - 1.0 ) <= c1 )
      {
        break;
      }
    }
  }
  else
  {
/*
  Call to algorithm AS 111 - note that P has been tested above.
  AS 241 could be used as an alternative.
*/
    x = ppnd ( p, ifault );
/*
  Starting approximation using Wilson and Hilferty estimate
*/
    p1 = c2 / v;
    ch = v * pow ( x * sqrt ( p1 ) + 1.0 - p1, 3 );
/*
  Starting approximation for P tending to 1.
*/
    if ( c6 * v + 6.0 < ch )
    {
      ch = - 2.0 * ( log ( 1.0 - p ) - c * log ( 0.5 * ch ) + g );
    }
  }

/*
  Call to algorithm AS 239 and calculation of seven term
  Taylor series
*/
  for ( i = 1; i <= maxit; i++ )
  {
    q = ch;
    p1 = 0.5 * ch;
    p2 = p - gammad ( p1, xx, &if1 );

    if ( if1 != 0 )
    {
      *ifault = 3;
      return value;
    }

    t = p2 * exp ( xx * aa + g + p1 - c * log ( ch ) );
    b = t / ch;
    a = 0.5 * t - b * c;
    s1 = ( c19 + a * ( c17 + a * ( c14 + a * ( c13 + a * ( c12 + 
      c11 * a ))))) / c24;
    s2 = ( c24 + a * ( c29 + a * ( c32 + a * ( c33 + c35 * a )))) / c37;
    s3 = ( c19 + a * ( c25 + a * ( c28 + c31 * a ))) / c37;
    s4 = ( c20 + a * ( c27 + c34 * a) + c * ( c22 + a * ( c30 + c36 * a ))) / c38;
    s5 = ( c13 + c21 * a + c * ( c18 + c26 * a )) / c37;
    s6 = ( c15 + c * ( c23 + c16 * c )) / c38;
    ch = ch + t * ( 1.0 + 0.5 * t * s1 - b * c * ( s1 - b * 
      ( s2 - b * ( s3 - b * ( s4 - b * ( s5 - b * s6 ))))));

    if ( e < fabs ( q / ch - 1.0 ) )
    {
       value = ch;
       return value;
    }
  }

 *ifault = 4;
 value = ch;

 return value;
}
/******************************************************************************/

double ppnd ( double p, int *ifault )

/******************************************************************************/
/*
  Purpose:

    PPND produces the normal deviate value corresponding to lower tail area = P.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 June 2013

  Author:

    Original FORTRAN77 version by J Beasley, S Springer.
    C version by John Burkardt.

  Reference:

    J Beasley, S Springer,
    Algorithm AS 111:
    The Percentage Points of the Normal Distribution,
    Applied Statistics,
    Volume 26, Number 1, 1977, pages 118-121.

  Parameters:

    Input, double P, the value of the cumulative probability
    densitity function.  0 < P < 1.

    Output, integer *IFAULT, error flag.
    0, no error.
    1, P <= 0 or P >= 1.  PPND is returned as 0.

    Output, double PPND, the normal deviate value with the property that
    the probability of a standard normal deviate being less than or
    equal to PPND is P.
*/
{
  double a0 = 2.50662823884;
  double a1 = -18.61500062529;
  double a2 = 41.39119773534;
  double a3 = -25.44106049637;
  double b1 = -8.47351093090;
  double b2 = 23.08336743743;
  double b3 = -21.06224101826;
  double b4 = 3.13082909833;
  double c0 = -2.78718931138;
  double c1 = -2.29796479134;
  double c2 = 4.85014127135;
  double c3 = 2.32121276858;
  double d1 = 3.54388924762;
  double d2 = 1.63706781897;
  double r;
  double split = 0.42;
  double value;

  *ifault = 0;
/*
  0.08 < P < 0.92
*/
  if ( fabs ( p - 0.5 ) <= split )
  {
    r = ( p - 0.5 ) * ( p - 0.5 );

    value = ( p - 0.5 ) * ( ( ( 
        a3   * r 
      + a2 ) * r 
      + a1 ) * r 
      + a0 ) / ( ( ( ( 
        b4   * r 
      + b3 ) * r 
      + b2 ) * r 
      + b1 ) * r 
      + 1.0 );
  }
/*
  P < 0.08 or P > 0.92,
  R = min ( P, 1-P )
*/
  else if ( 0.0 < p && p < 1.0 )
  {
    if ( 0.5 < p )
    {
      r = sqrt ( - log ( 1.0 - p ) );
    }
    else
    {
      r = sqrt ( - log ( p ) );
    }

    value = ( ( ( 
        c3   * r 
      + c2 ) * r 
      + c1 ) * r 
      + c0 ) / ( ( 
        d2   * r 
      + d1 ) * r 
      + 1.0 );

    if ( p < 0.5 )
    {
      value = - value;
    }
  }
/*
  P <= 0.0 or 1.0 <= P
*/
  else
  {
    *ifault = 1;
    value = 0.0;
  }

  return value;
}
/******************************************************************************/

double r8_min ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    R8_MIN returns the minimum of two R8's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 August 2004

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the quantities to compare.

    Output, double R8_MIN, the minimum of X and Y.
*/
{
  double value;

  if ( y < x )
  {
    value = y;
  } 
  else
  {
    value = x;
  }
  return value;
}
/******************************************************************************/

void timestamp ( )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    17 June 2014 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 June 2014

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}

