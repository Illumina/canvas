using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CanvasPartition
{
    internal static class R
    {
        // Not good! Depends on Math.NET Numerics
        public static double phyper(double x, double NR, double NB, double n,
           bool lowerTail, bool logP)
        {
            /* Sample of  n balls from  NR red  and	 NB black ones;	 x are red */
            double d, pd;
            if (Double.IsNaN(x) || Double.IsNaN(NR) || Double.IsNaN(NB) || Double.IsNaN(n))
            {
                return x + NR + NB + n;
            }
            x = Math.Floor(x + 1e-7);
            NR = R.R_D_forceint(NR);
            NB = R.R_D_forceint(NB);
            n = R.R_D_forceint(n);

            if (NR < 0 || NB < 0 || !R.R_finite(NR + NB) || n < 0 || n > NR + NB) { return Double.NaN; }

            if (x * (NR + NB) > n * NR)
            {
                /* Swap tails.	*/
                double oldNB = NB;
                NB = NR;
                NR = oldNB;
                x = n - x - 1;
                lowerTail = !lowerTail;
            }

            if (x < 0) { return R.R_DT_0(lowerTail, logP); }
            if (x >= NR || x >= n) { return R.R_DT_1(lowerTail, logP); }

            d = R.dhyper(x, NR, NB, n, logP);
            pd = R.pdhyper(x, NR, NB, n, logP);

            return logP ? R.R_DT_Log(d + pd, lowerTail) : R.R_D_Lval(d * pd, lowerTail);
        }

        public static double dhyper(double x, double r, double b, double n, bool giveLog)
        {
            if (Double.IsNaN(x) || Double.IsNaN(r) || Double.IsNaN(b) || Double.IsNaN(n))
            {
                return x + r + b + n;
            }
            if (R.R_D_negInonint(x))
            {
                return R.R_D__0(giveLog);
            }
            x = R.R_D_forceint(x);
            r = R.R_D_forceint(r);
            b = R.R_D_forceint(b);
            n = R.R_D_forceint(n);
            if (n < x || r < x || (n - x) > b) { return R.R_D__0(giveLog); }
            if (n == 0) { return (x == 0) ? R.R_D__1(giveLog) : R.R_D__0(giveLog); }
            double p = n / (r + b);
            double q = (r + b - n) / (r + b);

            double p1 = R.dbinom_raw(x, r, p, q, giveLog);
            double p2 = R.dbinom_raw(n - x, b, p, q, giveLog);
            double p3 = R.dbinom_raw(n, r + b, p, q, giveLog);

            return ((giveLog) ? p1 + p2 - p3 : p1 * p2 / p3);
        }

        public static double pdhyper(double x, double NR, double NB, double n, bool logP)
        {
            /*
             * Calculate
             *
             *	    phyper (x, NR, NB, n, TRUE, FALSE)
             *   [log]  ----------------------------------
             *	       dhyper (x, NR, NB, n, FALSE)
             *
             * without actually calling phyper.  This assumes that
             *
             *     x * (NR + NB) <= n * NR
             *
             */
            double sum = 0;
            double term = 1;

            while (x > 0 && term >= R.DBL_EPSILON * sum)
            {
                term *= x * (NB - n + x) / (n + 1 - x) / (NR + 1 - x);
                sum += term;
                x--;
            }

            double ss = (double)sum;
            return logP ? log1p(ss) : 1 + ss;
        }
        /*
                public static double dbinom_raw(double x, double n, double p, double q, bool log)
                {
                    int k = (int)(x);
                    var binom = new Binomial(p, (int)(n));
                    return log ? binom.ProbabilityLn(k) : binom.Probability(k);
                }
        */

        public static double dbinom_raw(double x, double n, double p, double q, bool giveLog)
        {
            double lf, lc;
            if (p == 0) { return (x == 0) ? R.R_D__1(giveLog) : R.R_D__0(giveLog); }
            if (q == 0) { return (x == n) ? R.R_D__1(giveLog) : R.R_D__0(giveLog); }
            if (x == 0)
            {
                if (n == 0) { return R.R_D__1(giveLog); }
                lc = (p < 0.1) ? (-R.bd0(n, n * q) - n * p) : (n * Math.Log(q));
                return R.R_D_exp(lc, giveLog);
            }
            if (x == n)
            {
                lc = (q < 0.1) ? -R.bd0(n, n * p) - n * q : n * Math.Log(p);
                return R.R_D_exp(lc, giveLog);
            }
            if (x < 0 || x > n) { return R.R_D__0(giveLog); }
            // n*p or n*q can underflow to zero if n and p or q are small.  This
            //   used to occur in dbeta, and gives NaN as from R 2.3.0.
            lc = R.stirlerr(n) - R.stirlerr(x) - R.stirlerr(n - x) - bd0(x, n * p) - bd0(n - x, n * q);

            // f = (M_2PI*x*(n-x))/n; could overflow or underflow
            // Upto R 2.7.1:
            // lf = log(M_2PI) + log(x) + log(n-x) - log(n);
            // -- following is much better for  x << n :
            lf = Math.Log(R.M_2PI) + Math.Log(x) + log1p(-x / n);

            return R.R_D_exp(lc - 0.5 * lf, giveLog);
        }

        public static double stirlerr(double n)
        {
            const double S0 = 0.083333333333333333333;  // 1/12
            const double S1 = 0.00277777777777777777778;    // 1/360
            const double S2 = 0.00079365079365079365079365;  // 1/1260 
            const double S3 = 0.000595238095238095238095238; // 1/1680 
            const double S4 = 0.0008417508417508417508417508;   // 1/1188 

            //  error for 0, 0.5, 1.0, 1.5, ..., 14.5, 15.0.

            double[] sferr_halves = {
	            0.0, // n=0 - wrong, place holder only 
	            0.1534264097200273452913848,  // 0.5 
	            0.0810614667953272582196702,  // 1.0 
	            0.0548141210519176538961390,  // 1.5 
	            0.0413406959554092940938221,  // 2.0 
	            0.03316287351993628748511048, // 2.5 
	            0.02767792568499833914878929, // 3.0 
	            0.02374616365629749597132920, // 3.5 
	            0.02079067210376509311152277, // 4.0 
	            0.01848845053267318523077934, // 4.5 
	            0.01664469118982119216319487, // 5.0 
	            0.01513497322191737887351255, // 5.5 
	            0.01387612882307074799874573, // 6.0 
	            0.01281046524292022692424986, // 6.5 
	            0.01189670994589177009505572, // 7.0 
	            0.01110455975820691732662991, // 7.5 
	            0.010411265261972096497478567, // 8.0 
	            0.009799416126158803298389475, // 8.5 
	            0.009255462182712732917728637, // 9.0 
	            0.008768700134139385462952823, // 9.5 
	            0.008330563433362871256469318, // 10.0 
	            0.007934114564314020547248100, // 10.5 
	            0.007573675487951840794972024, // 11.0 
	            0.007244554301320383179543912, // 11.5 
	            0.006942840107209529865664152, // 12.0 
	            0.006665247032707682442354394, // 12.5 
	            0.006408994188004207068439631, // 13.0 
	            0.006171712263039457647532867, // 13.5 
	            0.005951370112758847735624416, // 14.0 
	            0.005746216513010115682023589, // 14.5 
	            0.005554733551962801371038690  // 15.0 
            };
            double nn;

            if (n <= 15.0)
            {
                nn = n + n;
                if (nn == (int)nn) return (sferr_halves[(int)nn]);
                return (lgammafn(n + 1.0) - (n + 0.5) * Math.Log(n) + n - R.M_LN_SQRT_2PI);
            }

            nn = n * n;
            if (n > 500) return ((S0 - S1 / nn) / n);
            if (n > 80) return ((S0 - (S1 - S2 / nn) / nn) / n);
            if (n > 35) return ((S0 - (S1 - (S2 - S3 / nn) / nn) / nn) / n);
            // 15 < n <= 35 : 
            return ((S0 - (S1 - (S2 - (S3 - S4 / nn) / nn) / nn) / nn) / n);
        }

        public static double lgammafn_sign(double x, int[] sgn)
        {
            double ans, y, sinpiy;
            double xmax = 0.0;
            double dxrel = 0.0;

            if (xmax == 0)
            {/* initialize machine dependent constants _ONCE_ */
                xmax = d1mach(2) / Math.Log(d1mach(2));/* = 2.533 e305	 for IEEE double */
                dxrel = Math.Sqrt(d1mach(4));/* sqrt(Eps) ~ 1.49 e-8  for IEEE double */
            }

            if (sgn != null) sgn[0] = 1;
            if (Double.IsNaN(x)) return x;

            if (x < 0 && (Math.Floor(-x) % 2.0) == 0)
            {
                if (sgn != null) sgn[0] = -1;
            }

            if (x <= 0 && x == Math.Truncate(x))
            { // Negative integer argument
                ML_ERROR(ME_RANGE, "lgamma");
                return Double.PositiveInfinity; // +Inf, since lgamma(x) = log|gamma(x)|
            }

            y = Math.Abs(x);

            if (y < 1e-306) return -Math.Log(x); // denormalized range, R change
            if (y <= 10) return Math.Log(Math.Abs(gammafn(x)));
            //  ELSE  y = |x| > 10 ----------------------

            if (y > xmax)
            {
                ML_ERROR(ME_RANGE, "lgamma");
                return Double.PositiveInfinity;
            }

            if (x > 0)
            { // i.e. y = x > 10
                if (x > 1e17)
                    return (x * (Math.Log(x) - 1.0));
                else if (x > 4934720.0)
                    return (M_LN_SQRT_2PI + (x - 0.5) * Math.Log(x) - x);
                else
                    return M_LN_SQRT_2PI + (x - 0.5) * Math.Log(x) - x + lgammacor(x);
            }
            // else: x < -10; y = -x
            sinpiy = Math.Abs(Math.Sin(Math.PI * y));
            if (sinpiy == 0)
            { /* Negative integer argument ===
			  Now UNNECESSARY: caught above */
                MATHLIB_WARNING(" ** should NEVER happen! *** [lgamma.c: Neg.int, y={0}]\n", y.ToString());
                return ML_ERR_return_NAN();
            }

            ans = M_LN_SQRT_PId2 + (x - 0.5) * Math.Log(y) - x - Math.Log(sinpiy) - lgammacor(y);

            if (Math.Abs((x - Math.Truncate(x - 0.5)) * ans / x) < dxrel)
            {

                /* The answer is less than half precision because
                 * the argument is too near a negative integer. */
                ML_ERROR(ME_PRECISION, "lgamma");
            }

            return ans;
        }

        public static double gammafn(double x)
        {
            double[] gamcs = new double[42]{
	            +.8571195590989331421920062399942e-2,
	            +.4415381324841006757191315771652e-2,
	            +.5685043681599363378632664588789e-1,
	            -.4219835396418560501012500186624e-2,
	            +.1326808181212460220584006796352e-2,
	            -.1893024529798880432523947023886e-3,
	            +.3606925327441245256578082217225e-4,
	            -.6056761904460864218485548290365e-5,
	            +.1055829546302283344731823509093e-5,
	            -.1811967365542384048291855891166e-6,
	            +.3117724964715322277790254593169e-7,
	            -.5354219639019687140874081024347e-8,
	            +.9193275519859588946887786825940e-9,
	            -.1577941280288339761767423273953e-9,
	            +.2707980622934954543266540433089e-10,
	            -.4646818653825730144081661058933e-11,
	            +.7973350192007419656460767175359e-12,
	            -.1368078209830916025799499172309e-12,
	            +.2347319486563800657233471771688e-13,
	            -.4027432614949066932766570534699e-14,
	            +.6910051747372100912138336975257e-15,
	            -.1185584500221992907052387126192e-15,
	            +.2034148542496373955201026051932e-16,
	            -.3490054341717405849274012949108e-17,
	            +.5987993856485305567135051066026e-18,
	            -.1027378057872228074490069778431e-18,
	            +.1762702816060529824942759660748e-19,
	            -.3024320653735306260958772112042e-20,
	            +.5188914660218397839717833550506e-21,
	            -.8902770842456576692449251601066e-22,
	            +.1527474068493342602274596891306e-22,
	            -.2620731256187362900257328332799e-23,
	            +.4496464047830538670331046570666e-24,
	            -.7714712731336877911703901525333e-25,
	            +.1323635453126044036486572714666e-25,
	            -.2270999412942928816702313813333e-26,
	            +.3896418998003991449320816639999e-27,
	            -.6685198115125953327792127999999e-28,
	            +.1146998663140024384347613866666e-28,
	            -.1967938586345134677295103999999e-29,
	            +.3376448816585338090334890666666e-30,
	            -.5793070335782135784625493333333e-31
            };

            int i, n;
            double y;
            double sinpiy, value;

            int ngam = 0;
            double xmin = 0, xmax = 0.0, xsml = 0.0, dxrel = 0.0;

            // Initialize machine dependent constants, the first time gamma() is called.
            // FIXME for threads !
            if (ngam == 0)
            {
                ngam = chebyshev_init(gamcs, 42, DBL_EPSILON / 20); //was .1*d1mach(3)
                gammalims(ref xmin, ref xmax); //-> ./gammalims.c
                xsml = Math.Exp(fmax2(Math.Log(DBL_MIN), -Math.Log(DBL_MAX)) + 0.01);
                //   = exp(.01)*DBL_MIN = 2.247e-308 for IEEE
                dxrel = Math.Sqrt(DBL_EPSILON); //was sqrt(d1mach(4))
            }

            if (Double.IsNaN(x)) return x;

            // If the argument is exactly zero or a negative integer
            // then return NaN.
            if (x == 0 || (x < 0 && x == (long)x))
            {
                ML_ERROR(ME_DOMAIN, "gammafn");
                return Double.NaN;
            }

            y = Math.Abs(x);

            if (y <= 10)
            {
                // Compute gamma(x) for -10 <= x <= 10
                // Reduce the interval and find gamma(1 + y) for 0 <= y < 1
                // first of all. */

                n = (int)x;
                if (x < 0) --n;
                y = x - n;/* n = floor(x)  ==>	y in [ 0, 1 ) */
                --n;
                value = chebyshev_eval(y * 2 - 1, gamcs, ngam) + .9375;
                if (n == 0)
                    return value;/* x = 1.dddd = 1+y */

                if (n < 0)
                {
                    // compute gamma(x) for -10 <= x < 1 

                    // exact 0 or "-n" checked already above

                    // The answer is less than half precision
                    // because x too near a negative integer.
                    if (x < -0.5 && Math.Abs(x - (int)(x - 0.5) / x) < dxrel)
                    {
                        ML_ERROR(ME_PRECISION, "gammafn");
                    }

                    // The argument is so close to 0 that the result would overflow.
                    if (y < xsml)
                    {
                        ML_ERROR(ME_RANGE, "gammafn");
                        if (x > 0) return Double.PositiveInfinity;
                        else return Double.NegativeInfinity;
                    }

                    n = -n;

                    for (i = 0; i < n; i++)
                    {
                        value /= (x + i);
                    }
                    return value;
                }
                else
                {
                    // gamma(x) for 2 <= x <= 10

                    for (i = 1; i <= n; i++)
                    {
                        value *= (y + i);
                    }
                    return value;
                }
            }
            else
            {
                // gamma(x) for	 y = |x| > 10.

                if (x > xmax)
                {			// Overflow
                    ML_ERROR(ME_RANGE, "gammafn");
                    return Double.PositiveInfinity;
                }

                if (x < xmin)
                {			// Underflow
                    ML_ERROR(ME_UNDERFLOW, "gammafn");
                    return 0.0;
                }

                if (y <= 50 && y == (int)y)
                { // compute (n - 1)!
                    value = 1.0;
                    for (i = 2; i < y; i++) value *= i;
                }
                else
                { // normal case
                    value = Math.Exp((y - 0.5) * Math.Log(y) - y + M_LN_SQRT_2PI +
                        ((2 * y == (int)2 * y) ? stirlerr(y) : lgammacor(y)));
                }
                if (x > 0)
                    return value;

                if (Math.Abs((x - (int)(x - 0.5)) / x) < dxrel)
                {

                    /* The answer is less than half precision because */
                    /* the argument is too near a negative integer. */

                    ML_ERROR(ME_PRECISION, "gammafn");
                }

                sinpiy = Math.Sin(Math.PI * y);
                if (sinpiy == 0)
                {		/* Negative integer arg - overflow */
                    ML_ERROR(ME_RANGE, "gammafn");
                    return Double.PositiveInfinity;
                }

                return -Math.PI / (y * sinpiy * value);
            }
        }

        private static void gammalims(ref double xmin, ref double xmax)
        {
            // #ifdef IEEE_754
            xmin = -170.5674972726612;
            xmax = 171.61447887182298; // (3 Intel/Sparc architectures)
            // TODO: not sure if this will work
            // #else
            // ...
        }

        public static double fmax2(double x, double y)
        {
            if (Double.IsNaN(x) || Double.IsNaN(y))
                return x + y;
            return (x < y) ? y : x;
        }

        public static double lgammafn(double x)
        {
            return lgammafn_sign(x, null);
        }

        private static double lgammacor(double x)
        {
            double[] algmcs = new double[15]{
	            +.1666389480451863247205729650822e+0,
	            -.1384948176067563840732986059135e-4,
	            +.9810825646924729426157171547487e-8,
	            -.1809129475572494194263306266719e-10,
	            +.6221098041892605227126015543416e-13,
	            -.3399615005417721944303330599666e-15,
	            +.2683181998482698748957538846666e-17,
	            -.2868042435334643284144622399999e-19,
	            +.3962837061046434803679306666666e-21,
	            -.6831888753985766870111999999999e-23,
	            +.1429227355942498147573333333333e-24,
	            -.3547598158101070547199999999999e-26,
	            +.1025680058010470912000000000000e-27,
	            -.3401102254316748799999999999999e-29,
	            +.1276642195630062933333333333333e-30
            };

            double tmp;

            /* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
             *   xbig = 2 ^ 26.5
             *   xmax = DBL_MAX / 48 =  2^1020 / 3 */
            int nalgm = 5;
            double xbig = 94906265.62425156;
            double xmax = 3.745194030963158e306;

            if (x < 10)
                ML_ERR_return_NAN();
            else if (x >= xmax)
            {
                ML_ERROR(ME_UNDERFLOW, "lgammacor");
                /* allow to underflow below */
            }
            else if (x < xbig)
            {
                tmp = 10 / x;
                return chebyshev_eval(tmp * tmp * 2 - 1, algmcs, nalgm) / x;
            }
            return 1 / (x * 12);
        }

        public const double DBL_MIN = 2.22507e-308;
        public static double DBL_MAX = Double.MaxValue;
        public const double M_LOG10_2 = 0.301029995663981195213738894724; // log10(2)
        public const double M_LN_SQRT_PId2 = 0.225791352644727432363097614947; // log(sqrt(pi/2))

        private static double d1mach(int i)
        {
            switch (i)
            {
                case 1: return DBL_MIN;
                case 2: return DBL_MAX;

                case 3: /* = FLT_RADIX  ^ - DBL_MANT_DIG
	              for IEEE:  = 2^-53 = 1.110223e-16 = .5*DBL_EPSILON */
                    return 0.5 * DBL_EPSILON;

                case 4: /* = FLT_RADIX  ^ (1- DBL_MANT_DIG) =
	              for IEEE:  = 2^-52 = DBL_EPSILON */
                    return DBL_EPSILON;

                case 5: return M_LOG10_2;

                default: return 0.0;
            }
        }

        public static double bd0(double x, double np)
        {
            if (!R.R_finite(x) || !R.R_finite(np) || np == 0.0) { return Double.NaN; }
            if (Math.Abs(x - np) < 0.1 * (x + np))
            {
                double v = (x - np) / (x + np);
                double s = (x - np) * v;/* s using v -- change by MM */
                double ej = 2 * x * v;
                v = v * v;
                for (int j = 1; ; j++)
                { /* Taylor series */
                    ej *= v;
                    double s1 = s + ej / ((j << 1) + 1);
                    if (s1 == s) /* last term was effectively 0 */
                        return (s1);
                    s = s1;
                }
            }
            /* else:  | x - np |  is not too small */
            return (x * Math.Log(x / np) + np - x);
        }

        private const double DBL_EPSILON = 2.2204460492503131E-16; // dangerous...the smallest x such that 1 + x != 1
        private const double M_LN_SQRT_2PI = 0.918938533204672741780329736406;

        private const double M_2PI = Math.PI * 2;

        private const double M_LN2 = 0.693147180559945309417;

        public static double R_D_exp(double x, bool log)
        {
            return (log ? (x) : Math.Exp(x));
        }

        public static bool R_finite(double x)
        {
            return !(Double.IsNaN(x) || Double.IsInfinity(x));
        }

        public static bool R_D_nonint(double x)
        {
            return (Math.Abs(x - Math.Floor(x + 0.5)) > 1E-7);
        }

        public static bool R_D_negInonint(double x)
        {
            return (x < 0.0 || R.R_D_nonint(x));
        }

        public static double R_D_forceint(double x)
        {
            return Math.Floor(x + 0.5);
        }

        public static double R_D__0(bool log)
        {
            return log ? Double.NegativeInfinity : 0.0;
        }

        public static double R_D__1(bool log)
        {
            return log ? 0.0 : 1.0;
        }

        public static double R_DT_0(bool lowerTail, bool log)
        {
            return (lowerTail ? R.R_D__0(log) : R.R_D__1(log));
        }

        public static double R_DT_1(bool lowerTail, bool log)
        {
            return (lowerTail ? R.R_D__1(log) : R.R_D__0(log));
        }

        public static double R_DT_Log(double p, bool lowerTail)
        {
            return (lowerTail ? (p) : R.R_Log1_Exp(p));
        }

        /* log(1 - exp(x))  in more stable form than log1p(- R_D_qIv(x))) : */
        public static double R_Log1_Exp(double x)
        {
            return ((x) > -M_LN2 ? Math.Log(-expm1(x)) : log1p(-Math.Exp(x)));
        }

        // Computes exp(x) - 1
        public static double expm1(double x)
        {
            double y, a = Math.Abs(x);

            if (a < DBL_EPSILON) return x;
            if (a > 0.697) return Math.Exp(x) - 1;  /* negligible cancellation */

            if (a > 1E-8)
                y = Math.Exp(x) - 1;
            else /* Taylor expansion, more accurate in this range */
                y = (x / 2 + 1) * x;

            /* Newton step for solving   log(1 + y) = x   for y : */
            /* WARNING: does not work for y ~ -1: bug in 1.5.0 */
            y -= (1 + y) * (log1p(y) - x);
            return y;
        }

        public static double R_D_Lval(double p, bool lowerTail)
        {
            return (lowerTail ? (p) : (0.5 - (p) + 0.5));
        }

        // Computes log(1 + x)
        public static double log1p(double x)
        {
            /* series for log1p on the interval -.375 to .375
             *				     with weighted error   6.35e-32
             *				      log weighted error  31.20
             *			    significant figures required  30.93
             *				 decimal places required  32.01
             */
            double[] alnrcs = new double[43]{
	            +.10378693562743769800686267719098e+1,
	            -.13364301504908918098766041553133e+0,
	            +.19408249135520563357926199374750e-1,
	            -.30107551127535777690376537776592e-2,
	            +.48694614797154850090456366509137e-3,
	            -.81054881893175356066809943008622e-4,
	            +.13778847799559524782938251496059e-4,
	            -.23802210894358970251369992914935e-5,
	            +.41640416213865183476391859901989e-6,
	            -.73595828378075994984266837031998e-7,
	            +.13117611876241674949152294345011e-7,
	            -.23546709317742425136696092330175e-8,
	            +.42522773276034997775638052962567e-9,
	            -.77190894134840796826108107493300e-10,
	            +.14075746481359069909215356472191e-10,
	            -.25769072058024680627537078627584e-11,
	            +.47342406666294421849154395005938e-12,
	            -.87249012674742641745301263292675e-13,
	            +.16124614902740551465739833119115e-13,
	            -.29875652015665773006710792416815e-14,
	            +.55480701209082887983041321697279e-15,
	            -.10324619158271569595141333961932e-15,
	            +.19250239203049851177878503244868e-16,
	            -.35955073465265150011189707844266e-17,
	            +.67264542537876857892194574226773e-18,
	            -.12602624168735219252082425637546e-18,
	            +.23644884408606210044916158955519e-19,
	            -.44419377050807936898878389179733e-20,
	            +.83546594464034259016241293994666e-21,
	            -.15731559416479562574899253521066e-21,
	            +.29653128740247422686154369706666e-22,
	            -.55949583481815947292156013226666e-23,
	            +.10566354268835681048187284138666e-23,
	            -.19972483680670204548314999466666e-24,
	            +.37782977818839361421049855999999e-25,
	            -.71531586889081740345038165333333e-26,
	            +.13552488463674213646502024533333e-26,
	            -.25694673048487567430079829333333e-27,
	            +.48747756066216949076459519999999e-28,
	            -.92542112530849715321132373333333e-29,
	            +.17578597841760239233269760000000e-29,
	            -.33410026677731010351377066666666e-30,
	            +.63533936180236187354180266666666e-31
            };
            int nlnrel = 0;
            double xmin = 0.0;
            if (xmin == 0.0) xmin = -1 + Math.Sqrt(DBL_EPSILON);/*was sqrt(d1mach(4)); */
            if (nlnrel == 0) /* initialize chebychev coefficients */
                nlnrel = chebyshev_init(alnrcs, 43, DBL_EPSILON / 20);/*was .1*d1mach(3)*/

            if (x == 0.0) return 0.0;/* speed */
            if (x == -1) return Double.NegativeInfinity;
            if (x < -1) return ML_ERR_return_NAN();

            if (Math.Abs(x) <= .375)
            {
                /* Improve on speed (only);
               again give result accurate to IEEE double precision: */
                if (Math.Abs(x) < .5 * DBL_EPSILON)
                    return x;

                if ((0 < x && x < 1e-8) || (-1e-9 < x && x < 0))
                    return x * (1 - .5 * x);
                /* else */
                return x * (1 - x * chebyshev_eval(x / .375, alnrcs, nlnrel));
            }
            /* else */
            if (x < xmin)
            {
                /* answer less than half precision because x too near -1 */
                ML_ERROR(ME_PRECISION, "log1p");
            }
            return Math.Log(1 + x);
        }

        private static int chebyshev_init(double[] dos, int nos, double eta)
        {
            int i, ii;
            double err;

            if (nos < 1)
                return 0;

            err = 0.0;
            i = 0;			/* just to avoid compiler warnings */
            for (ii = 1; ii <= nos; ii++)
            {
                i = nos - ii;
                err += Math.Abs(dos[i]);
                if (err > eta)
                {
                    return i;
                }
            }
            return i;
        }

        private static double chebyshev_eval(double x, double[] a, int n)
        {
            double b0, b1, b2, twox;
            int i;

            if (n < 1 || n > 1000) return ML_ERR_return_NAN();

            if (x < -1.1 || x > 1.1) return ML_ERR_return_NAN();

            twox = x * 2;
            b2 = b1 = 0;
            b0 = 0;
            for (i = 1; i <= n; i++)
            {
                b2 = b1;
                b1 = b0;
                b0 = twox * b1 - b2 + a[n - i];
            }
            return (b0 - b2) * 0.5;
        }


        static double ML_ERR_return_NAN()
        {
            R.ML_ERROR(ME_DOMAIN, "");
            return Double.NaN;
        }

        /* For a long time prior to R 2.3.0 ML_ERROR did nothing.
           We don't report ME_DOMAIN errors as the callers collect ML_NANs into
           a single warning.
        */
        static void ML_ERROR(int x, string s)
        {
            if (x > ME_DOMAIN)
            {
                string msg = "";
                switch (x)
                {
                    case ME_DOMAIN:
                        msg = "argument out of domain in '{0}'\n";
                        break;
                    case ME_RANGE:
                        msg = "value out of range in '{0}'\n";
                        break;
                    case ME_NOCONV:
                        msg = "convergence failed in '{0}'\n";
                        break;
                    case ME_PRECISION:
                        msg = "full precision may not have been achieved in '{0}'\n";
                        break;
                    case ME_UNDERFLOW:
                        msg = "underflow occurred in '{0}'\n";
                        break;
                }
                MATHLIB_WARNING(msg, s);
            }
        }

        static void MATHLIB_WARNING(string fmt, string x)
        {
            Console.Error.Write(fmt, x);
        }

        private const int ME_NONE = 0;
        /*	no error */
        private const int ME_DOMAIN = 1;
        /*	argument out of domain */
        private const int ME_RANGE = 2;
        /*	value out of range */
        private const int ME_NOCONV = 4;
        /*	process did not converge */
        private const int ME_PRECISION = 8;
        /*	does not have "full" precision */
        private const int ME_UNDERFLOW = 16;
        /*	and underflow occured (important for IEEE)*/



    }
}
