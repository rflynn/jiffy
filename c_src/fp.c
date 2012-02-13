
#include <float.h>
#include <math.h>
#include <stdio.h>

#include "jiffy.h"

#define FLOAT_BIAS 1022
#define MIN_EXP -1074
#define BIG_POW 4503599627370496

#if defined(_M_IX86)
typedef unsigned __int64 US64;
#elif defined(__alpha) || defined(__x86_64)
typedef unsigned long US64;
#else
typedef unsigned long long US64
#endif

#define ADD_CHAR(p, b, c)   \
do {                        \
    (p)[b++] = (c);         \
    if(b >= 64) return 0;   \
} while(0)


static inline void
convert_double(double v, long* frac, long* expon)
{
    US64 *rep = (US64*) &v;
    *frac = (long) (*rep & 0x000FFFFFFFFFFFFFULL);
    *expon = (long) (((*rep) >> 52) & 0x7FF);
}


static inline int
int_pow(int x, int n)
{
    long r = 1;

    if(n == 0)
        return 1;

    while(n >= 2) {
        r = (n & 1) == 1 ? r * x : r;
        x = x * x;
        n = n >> 1;
    }

    return r * x;
}


static inline int
int_ceil(double d)
{
    double o = trunc(d);
    fprintf(stderr, "IC: %0.20lf %0.20lf\r\n", d, o);
    if(d - o > 0) return ((int) o) + 1;
    return (int) o;
}


int
fmt_double(double v, char* p)
{
    long frac;
    long expon;
    long bexp;
    long mplus;
    long mminus;
    long r;
    long s;
    long d;
    int scale;
    int k;
    int place;
    int est;
    int round;
    int highok;
    int lowok;
    int toolow;
    int tc1;
    int tc2;
    int count;
    size_t pos;

    fprintf(stderr, "FMT: %lf\r\n", v);

    convert_double(v, &frac, &expon);

    fprintf(stderr, "%ld %ld\r\n", frac, expon);

    // frexp_int
	if(expon == 0) {
		expon = MIN_EXP;
	} else {
        frac = frac + (((long) 1) << 52);
        expon = ((long) expon) - 53 - FLOAT_BIAS;
    }

    // digits1
    round = (frac & 1) == 0 ? 1 : 0;
    if(expon >= 0) {
        bexp = ((long) 1) << expon;
        if(frac != BIG_POW) {
            fprintf(stderr, "1\r\n");
            r = frac * bexp * 2;
            s = 2;
            mplus = bexp;
            mminus = bexp;
            highok = round;
            lowok = round;
        } else {
            fprintf(stderr, "2\r\n");
            r = frac * bexp * 4;
            s = 4;
            mplus = bexp * 2;
            mminus = bexp;
            highok = round;
            lowok = round;
        }
    } else {
        if(expon == MIN_EXP || frac != BIG_POW) {
            fprintf(stderr, "3 %ld\r\n", expon);
            r = frac * 2;
            s = ((long) 1) << (1 - expon);
            mplus = 1;
            mminus = 1;
            highok = round;
            lowok = round;
        } else {
            fprintf(stderr, "4\r\n");
            r = frac * 4;
            s = ((long) 1) << (2 - expon);
            mplus = 2;
            mminus = 1;
            highok = round;
            lowok = round;
        }
    }

    fprintf(stderr, "R S: %ld %ld\r\n", r, s);

    // scale
    est = int_ceil(log10(fabs(v)) - 1.0E-10);
    fprintf(stderr, "EST: %lf %d\r\n", log10(fabs(v)) - 1.0E-10, est);
    if(est > 0) {
        s *= int_pow(10, est);
        k = est;
    } else {
        scale = int_pow(10, -est);
        r *= scale;
        mplus *= scale;
        mminus *= scale;
        k = est;
    }

    // fixup
    if(highok) {
        toolow = (r + mplus) >= s ? 1 : 0;
    } else {
        toolow = (r + mplus) > s ? 1 : 0;
    }

    if(toolow) {
        place = k + 1;
    } else {
        place = k;
        r *= 10;
        mplus *= 10;
        mminus *= 10;
    }

    // generate
    tc1 = 0;
    tc2 = 0;
    pos = 0;
    count = 0;

    if(v < 0) {
        ADD_CHAR(p, pos, '-');
    }

    if(place > -6 && place < 0) {
        ADD_CHAR(p, pos, '0');
        ADD_CHAR(p, pos, '.');
        while(place < 0) {
            ADD_CHAR(p, pos, '0');
        }
    }

    while(tc1 == 0 && tc2 == 0) {
        d = r / s;
        fprintf(stderr, "%ld %ld %ld\r\n", r, s, d);
        r = r % s;

        if(lowok) {
            tc1 = (r <= mminus) ? 1 : 0;
        } else {
            tc1 = (r < mminus) ? 1 : 0;
        }
        if(highok) {
            tc2 = (r + mplus >= s) ? 1 : 0;
        } else {
            tc2 = (r + mplus > s) ? 1 : 0;
        }

        if(tc1 == 0 && tc2 == 0) {
            r *= 10;
            mplus *= 10;
            mminus *= 10;
            ADD_CHAR(p, pos, '0' + d);
            count++;
        } else if(tc1 == 0 && tc2 == 1) {
            ADD_CHAR(p, pos, '0' + (d+1));
            count++;
        } else if(tc1 == 1 && tc2 == 0) {
            ADD_CHAR(p, pos, '0' + d);
            count++;
        } else if(r * 2 < s) {
            ADD_CHAR(p, pos, '0' + d);
            count++;
        } else {
            ADD_CHAR(p, pos, '0' + (d+1));
            count++;
        }

        if((place <= -6 || place >= 6) && count == 1) {
            ADD_CHAR(p, pos, '.');
        } else if(count == place) {
            place = 0;
            ADD_CHAR(p, pos, '.');
        }
    }

    if(p[pos-1] == '.') {
        ADD_CHAR(p, pos, '0');
    }

    if(place > 0 && place < 6) {
        while(pos < place) {
            ADD_CHAR(p, pos, '0');
        }
        if(pos == place) {
            ADD_CHAR(p, pos, '.');
        }
    } else {
        ADD_CHAR(p, pos, 'e');
        if(place < 0) {
            ADD_CHAR(p, pos, '-');
            place = -place;
        } else {
            ADD_CHAR(p, pos, '+');
        }
        snprintf(p+pos, 31-pos, "%d", place);
    }

    return 1;
}

