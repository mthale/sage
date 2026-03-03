r"""
Elliptic curves over a function field

As of now largely based upon ell_generic.py and ell_number_field.py

EXAMPLES:

We construct an elliptic curve over a function field::
    sage: p = 11
    sage: k = GF(p)
    sage: K.<t> = FunctionField(k)
    sage: E = EllipticCurve(K, [0,0,0,t,1]); E
    Elliptic Curve defined by y^2 = x^3 + t*x + 1 over Rational function field in t over Finite Field of size 11
"""
# Copyright etc

import math
import sage.rings.abc
from sage.categories.number_fields import NumberFields
from sage.categories.finite_fields import FiniteFields
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring import polygen
from sage.rings.rational_field import QQ
from sage.misc.misc_c import prod
from sage.schemes.elliptic_curves.ell_point import EllipticCurvePoint_field
from sage.schemes.curves.projective_curve import ProjectivePlaneCurve_field



from .ell_field import EllipticCurve_field
from .ell_point import EllipticCurvePoint_field
from .constructor import EllipticCurve
from .ell_curve_isogeny import EllipticCurveIsogeny, isogeny_codomain_from_kernel
from . import ell_generic

sqrt = math.sqrt
exp = math.exp

class EllipticCurve_rational_function_field_global(EllipticCurve_field):
    """
    Elliptic curve over rational function field with finite constant field
    """



    def __init__(self, K, ainvs):
        r"""
        EXAMPLES:

        We construct an elliptic curve over a function field::

            sage: p = 11
            sage: k = GF(p)
            sage: K.<t> = FunctionField(k)
            sage: E = EllipticCurve(K, [0,0,0,t,1]); E
            Elliptic Curve defined by y^2 = x^3 + t*x + 1 over Rational function field in t over Finite Field of size 11
        """
        # quick test
        from sage.rings.function_field.function_field_rational import RationalFunctionField_global
        
        if not isinstance(K, RationalFunctionField_global):
            raise TypeError("base field must be a rational function field with finite constant field")

        EllipticCurve_field.__init__(self, K, ainvs)

        def naive_height(self):
            """Calculate the naive height at point P"""
            raise NotImplementedError("not yet")
