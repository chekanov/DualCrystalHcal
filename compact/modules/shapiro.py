# Calculates the Shapiro-Wilk W test and its significance level
# Translated to Python by S.Chekanov (ANL)
import math


# The inverse of cdf.
def normalQuantile(p, mu, sigma):

    if (sigma < 0):
        return -1
    if (sigma == 0):
        return mu
    q = p - 0.5
    # print "p=",p
    if (0.075 <= p and p <= 0.925):
        r = 0.180625 - q * q
        val = q * (((((((r * 2509.0809287301226727 + 33430.575583588128105) * r + 67265.770927008700853) * r
                       + 45921.953931549871457) * r + 13731.693765509461125) * r + 1971.5909503065514427) * r + 133.14166789178437745) * r
                   + 3.387132872796366608) / (((((((r * 5226.495278852854561 + 28729.085735721942674) * r + 39307.89580009271061) * r
                                                  + 21213.794301586595867) * r + 5394.1960214247511077) * r + 687.1870074920579083) * r + 42.313330701600911252) * r + 1)

    else:  # closer than 0.075 from {0,1} boundary #
        # r = min(p, 1-p) < 0.075 #
        if (q > 0):
            r = 1 - p
        else:
            r = p  # R_DT_Iv(p) ^=  p #

        # r = sqrt(-log(r))  <==>  min(p, 1-p) = exp( - r^2 ) #
        r = math.sqrt(-math.log(r))

        if (r <= 5.):  # <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11  #
            r += -1.6
            val = (((((((r * 7.7454501427834140764e-4 + 0.0227238449892691845833) * r + .24178072517745061177) * r
                       + 1.27045825245236838258) * r + 3.64784832476320460504) * r + 5.7694972214606914055) * r
                    + 4.6303378461565452959) * r + 1.42343711074968357734) / (((((((r * 1.05075007164441684324e-9 + 5.475938084995344946e-4) * r
                                                                               + .0151986665636164571966) * r + 0.14810397642748007459) * r + 0.68976733498510000455) * r + 1.6763848301838038494) * r
                                                                               + 2.05319162663775882187) * r + 1)
        else:  # * very close to  0 or 1 #
            r += -5.
            val = (((((((r * 2.01033439929228813265e-7 + 2.71155556874348757815e-5) * r + 0.0012426609473880784386) * r
                       + 0.026532189526576123093) * r + .29656057182850489123) * r + 1.7848265399172913358) * r + 5.4637849111641143699) * r
                   + 6.6579046435011037772) / (((((((r * 2.04426310338993978564e-15 + 1.4215117583164458887e-7) * r
                                                 + 1.8463183175100546818e-5) * r + 7.868691311456132591e-4) * r + .0148753612908506148525) * r
                                                 + .13692988092273580531) * r + .59983220655588793769) * r + 1.)

        if (q < 0.0):
            val = -val
            # return (q >= 0.)? r : -r ; #
    return mu + sigma * val


def sign(x):
    if (x == 0):
        return 0
    return 1 if x > 0 else -1


def poly(cc, nord, x):
    # Algorithm AS 181.2    Appl. Statist.  (1982) Vol. 31, No. 2
    # Calculates the algebraic polynomial of order nord-1 with array of coefficients cc.
    # Zero order coefficient is cc(1) = cc[0] */
    ret_val = cc[0]
    if (nord > 1):
        p = x * cc[nord-1]
        for j in range(nord - 2, 0, -1):
            # for (j = nord - 2; j > 0; j--):
            p = (p + cc[j]) * x
        ret_val += p
    return ret_val


def ShapiroWilkW(x):
    # x = x.sort(function (a, b) { return a - b; });
    x.sort()
    n = len(x)
    if (n < 3):
        return undefined
    nn2 = int(math.floor(n / 2))
    # print nn2
    a = [None]*(int(math.floor(nn2)) + 1)
#	ALGORITHM AS R94 APPL. STATIST. (1995) vol.44, no.4, 547-551.
#	Calculates the Shapiro-Wilk W test and its significance level

    small = 1e-19
    #  polynomial coefficients
    g = [-2.273, 0.459]
    c1 = [0, 0.221157, -0.147981, -2.07119, 4.434685, -2.706056]
    c2 = [0, 0.042981, -0.293762, -1.752461, 5.682633, -3.582633]
    c3 = [0.544, -0.39978, 0.025054, -6.714e-4]
    c4 = [1.3822, -0.77857, 0.062767, -0.0020322]
    c5 = [-1.5861, -0.31082, -0.083751, 0.0038915]
    c6 = [-0.4803, -0.082676, 0.0030302]

    # Local variables */
    pw = 1
    an = n

    if (n == 3):
        a[1] = 0.70710678  # /* = sqrt(1/2) */
    else:
        an25 = an + 0.25
        summ2 = 0.0
        for i in range(1, nn2+1, 1):
            # for (i = 1; i <= nn2; i++) :
            #print (i - 0.375) / an25, 0, 1
            a[i] = normalQuantile((i - 0.375) / an25, 0, 1)  # p(X <= x),
            r__1 = a[i]
            summ2 += r__1 * r__1

        summ2 *= 2
        ssumm2 = math.sqrt(summ2)
        rsn = 1 / math.sqrt(an)
        a1 = poly(c1, 6, rsn) - a[1] / ssumm2

        # Normalize a[] */
        if (n > 5):
            i1 = 3
            a2 = -a[2] / ssumm2 + poly(c2, 6, rsn)
            fac = math.sqrt(
                (summ2 - 2 * (a[1] * a[1]) - 2 * (a[2] * a[2])) / (1 - 2 * (a1 * a1) - 2 * (a2 * a2)))
            a[2] = a2
        else:
            i1 = 2
            fac = math.sqrt((summ2 - 2 * (a[1] * a[1])) / (1 - 2 * (a1 * a1)))

        a[1] = a1
        # for (i = i1; i <= nn2; i++):
        for i in range(i1, nn2+1, 1):
            a[i] /= - fac

# Check for zero range

    hrange = x[n - 1] - x[0]
    if (hrange < small):
        print('range is too small!')
        return undefined


#  Check for correct sort order on range - scaled X */

    xx = x[0] / hrange
    sx = xx
    sa = -a[1]
    # for (i = 1, j = n - 1; i < n; j--):
    # print "Loop=",n
    for i in range(1, n, 1):
        j = n-i
        # print i, j
        xi = x[i] / hrange
        if (xx - xi > small):
            print("xx - xi is too big.", xx - xi)
            return undefined
        sx += xi
        i = i+1
        if (i != j):
            sa += sign(i - j) * a[min(i, j)]
        xx = xi

    if (n > 5000000):
        print("n is too big!")
        return undefined


#	Calculate W statistic as squared correlation
#	between data and coefficients */

    sa /= n
    sx /= n
    ssa = ssx = sax = 0.
    # print "New Loop=",n
    # for (i = 0, j = n - 1; i < n; i++, j--):
    for i in range(0, n, 1):
        j = n-i-1
        # print i,j
        if (i != j):
            asa = sign(i - j) * a[1 + min(i, j)] - sa
        else:
            asa = -sa
        xsx = x[i] / hrange - sx
        ssa += asa * asa
        ssx += xsx * xsx
        sax += asa * xsx

#	W1 equals (1-W) calculated to avoid excessive rounding error
#	for W very near 1 (a potential problem in very large samples)

    ssassx = math.sqrt(ssa * ssx)
    w1 = (ssassx - sax) * (ssassx + sax) / (ssa * ssx)
    w = 1 - w1

#	Calculate significance level for W

    if (n == 3):  # * exact P value :
        pi6 = 1.90985931710274  # = 6/pi */
        stqr = 1.04719755119660  # = asin(sqrt(3/4)) */
        pw = pi6 * (math.asin(math.sqrt(w)) - stqr)
        if (pw < 0.):
            pw = 0
            return w

    y = math.log(w1)
    xx = math.log(an)
    if (n <= 11):
        gamma = poly(g, 2, an)
        if (y >= gamma):
            pw = 1e-99  # an "obvious" value, was 'small' which was 1e-19f
            return w
        y = -math.log(gamma - y)
        m = poly(c3, 4, an)
        s = math.exp(poly(c4, 4, an))
    else:  # n >= 12
        m = poly(c5, 4, xx)
        s = math.exp(poly(c6, 3, xx))

    # Oops, we don't have pnorm
    # pw = pnorm(y, m, s, 0/* upper tail */, 0);

    return w


# print ShapiroWilkW([1, 1.3, 0.7, 1.2, 0.9, 1.1])
# print ShapiroWilkW([1, 20, 50, 100, 100, 200])
