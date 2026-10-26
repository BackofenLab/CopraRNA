"""Numerical conventions used by the original R statistical programs."""

from __future__ import annotations

import math
import numpy as np


def r_sum(values):
    """Sum with an extended-precision accumulator before returning a Python float.

    Retain R-compatible reduction rounding used by captured statistical fixtures.
    """
    return float(np.sum(np.asarray(values, dtype=np.longdouble)))


def r_mean(values):
    """Return the R-compatible corrected mean using an extended-precision residual sum."""
    a = np.asarray(values, dtype=np.longdouble)
    mean = np.sum(a) / len(a)
    return float(mean + np.sum(a - mean) / len(a))


def nelder_mead(
    function,
    initial,
    *,
    max_evaluations=500,
    relative_tolerance=math.sqrt(np.finfo(float).eps),
):
    """The simplex search used by legacy optim's default method.

    It contracts after any reflection that does not improve the best vertex;
    retaining a successful reflection if the subsequent contraction fails.
    Its convergence tolerance is relative to the initial function value.
    """
    initial = np.asarray(initial, dtype=float)
    dimension = len(initial)
    vertices = [initial.copy() for _ in range(dimension + 1)]
    step = max(abs(initial)) * 0.1 or 0.1
    for i in range(dimension):
        vertices[i + 1][i] += step
    values = [float(function(vertex)) for vertex in vertices]
    evaluations = len(values)
    tolerance = relative_tolerance * (abs(values[0]) + relative_tolerance)
    previous_size = math.inf
    while evaluations <= max_evaluations:
        low = min(range(dimension + 1), key=lambda i: values[i])
        high = max(range(dimension + 1), key=lambda i: values[i])
        if values[high] <= values[low] + tolerance:
            break
        centroid = np.array(
            [
                sum(vertices[j][i] for j in range(dimension + 1) if j != high)
                / dimension
                for i in range(dimension)
            ]
        )
        reflected = 2 * centroid - vertices[high]
        reflected_value = float(function(reflected))
        evaluations += 1
        if reflected_value < values[low]:
            expanded = 2 * reflected - centroid
            expanded_value = float(function(expanded))
            evaluations += 1
            if expanded_value < reflected_value:
                vertices[high], values[high] = expanded, expanded_value
            else:
                vertices[high], values[high] = reflected, reflected_value
        else:
            old_high = values[high]
            if reflected_value < values[high]:
                vertices[high], values[high] = reflected, reflected_value
            contracted = 0.5 * vertices[high] + 0.5 * centroid
            contracted_value = float(function(contracted))
            evaluations += 1
            if contracted_value < values[high]:
                vertices[high], values[high] = contracted, contracted_value
            elif reflected_value >= old_high:
                size = 0.0
                for i in range(dimension + 1):
                    if i != low:
                        vertices[i] = 0.5 * (vertices[i] + vertices[low])
                        size += sum(abs(vertices[i] - vertices[low]))
                        values[i] = float(function(vertices[i]))
                        evaluations += 1
                if size >= previous_size:
                    break
                previous_size = size
    low = min(range(dimension + 1), key=lambda i: values[i])
    return vertices[low]


def normal_upper_scalar(x):
    """Cody rational normal-tail approximation with split-exponential evaluation."""
    if math.isnan(x):
        return math.nan
    if math.isinf(x):
        return 0.0 if x > 0 else 1.0
    y = abs(x)
    if y <= 0.67448975:
        a = (
            2.2352520354606839287,
            161.02823106855587881,
            1067.6894854603709582,
            18154.981253343561249,
            0.065682337918207449113,
        )
        b = (
            47.20258190468824187,
            976.09855173777669322,
            10260.932208618978205,
            45507.789335026729956,
        )
        square = x * x
        numerator = a[4] * square
        denominator = square
        for j in range(3):
            numerator = (numerator + a[j]) * square
            denominator = (denominator + b[j]) * square
        return 0.5 - x * (numerator + a[3]) / (denominator + b[3])
    if y <= math.sqrt(32):
        c = (
            0.39894151208813466764,
            8.8831497943883759412,
            93.506656132177855979,
            597.27027639480026226,
            2494.5375852903726711,
            6848.1904505362823326,
            11602.651437647350124,
            9842.7148383839780218,
            1.0765576773720192317e-8,
        )
        d = (
            22.266688044328115691,
            235.38790178262499861,
            1519.377599407554805,
            6485.558298266760755,
            18615.571640885098091,
            34900.952721145977266,
            38912.003286093271411,
            19685.429676859990727,
        )
        numerator = c[8] * y
        denominator = y
        for j in range(7):
            numerator = (numerator + c[j]) * y
            denominator = (denominator + d[j]) * y
        ratio = (numerator + c[7]) / (denominator + d[7])
    else:
        p = (
            0.21589853405795699,
            0.1274011611602473639,
            0.022235277870649807,
            0.001421619193227893466,
            2.9112874951168792e-5,
            0.02307344176494017303,
        )
        q = (
            1.28426009614491121,
            0.468238212480865118,
            0.0659881378689285515,
            0.00378239633202758244,
            7.29751555083966205e-5,
        )
        square = 1 / (y * y)
        numerator = p[5] * square
        denominator = square
        for j in range(4):
            numerator = (numerator + p[j]) * square
            denominator = (denominator + q[j]) * square
        ratio = (
            1 / math.sqrt(2 * math.pi)
            - square * (numerator + p[4]) / (denominator + q[4])
        ) / y
    split = math.trunc(y * 16) / 16
    value = (
        math.exp(-split * split * 0.5)
        * math.exp(-((y - split) * (y + split)) * 0.5)
        * ratio
    )
    return value if x > 0 else 1 - value


def normal_upper(values):
    """Evaluate the retained standard-normal survival function elementwise.

    Preserve ndarray shape and return a float for scalar input. Using the same
    scalar approximation avoids changing the correlation fit through tail rounding.
    """
    values = np.asarray(values, dtype=float)
    result = np.array([normal_upper_scalar(float(x)) for x in values.ravel()]).reshape(
        values.shape
    )
    return float(result) if result.ndim == 0 else result


def normal_probit(values):
    """Convert upper-tail probabilities to standard-normal scores.

    Preserve input shape, propagate NaN, and map zero/one to positive/negative
    infinity. Retain the central rational approximation and tail conventions
    required by the original R p-value-combination fixtures.
    """
    from statistics import NormalDist

    normal = NormalDist()

    def quantile(p):
        if math.isnan(p):
            return math.nan
        if p == 0:
            return math.inf
        if p == 1:
            return -math.inf
        q = (0.5 - p + 0.5) - 0.5
        if abs(q) <= 0.425:
            r = 0.180625 - q * q
            numerator = (
                2509.0809287301226727,
                33430.575583588128105,
                67265.770927008700853,
                45921.953931549871457,
                13731.693765509461125,
                1971.5909503065514427,
                133.14166789178437745,
                3.3871328727963666080,
            )
            denominator = (
                5226.4952788528545610,
                28729.085735721942674,
                39307.895800092710610,
                21213.794301586595867,
                5394.1960214247511077,
                687.18700749205790830,
                42.313330701600911252,
                1.0,
            )
            n = numerator[0]
            d = denominator[0]
            for a, b in zip(numerator[1:], denominator[1:]):
                n = n * r + a
                d = d * r + b
            return q * n / d
        return -normal.inv_cdf(p)

    values = np.asarray(values, dtype=float)
    result = np.array([quantile(float(p)) for p in values.ravel()]).reshape(
        values.shape
    )
    return float(result) if result.ndim == 0 else result
