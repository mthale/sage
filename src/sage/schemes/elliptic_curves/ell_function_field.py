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

from sage.rings.function_field.function_field_rational import RationalFunctionField
from sage.rings.function_field.place import FunctionFieldPlace
from sage.rings.function_field.function_field import FunctionField


from .ell_field import EllipticCurve_field
#from .ell_field import EllipticCurve_function_field
from .ell_point import EllipticCurvePoint_function_field
from .constructor import EllipticCurve
from .ell_curve_isogeny import EllipticCurveIsogeny, isogeny_codomain_from_kernel
from . import ell_generic

sqrt = math.sqrt
exp = math.exp

class EllipticCurve_function_field(EllipticCurve_field):
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
        """quick test
        from sage.rings.function_field.function_field import FunctionField
        K = self.base_ring()
        if not isinstance(K, FunctionField):
            raise TypeError("base field must be a rational function field with finite constant field")
        """
        EllipticCurve_field.__init__(self, K, ainvs)
        
        
    
    _point = EllipticCurvePoint_function_field
    
    
    def places(self):
        r"""
        Return the places of bad reduction on the elliptic curve    
        
        INPUT:

        - ``self`` -- elliptic curve defined over a rational function field
        
        OUTPUT: a list of
        
        - ``
        
        """
        #inf_place = self.base_field().places()[0]
        D = self.discriminant()
        _places = self.discriminant().poles()
        zero_places = self.discriminant().zeros()
        _places.extend(p for p in zero_places if p not in _places)
        return _places
    
    def is_local_integral(self, v) -> bool:
        r"""Determines if coeffecients have positive valuation at place v.
            In order to check minimality the model must be integral"""
        
        if all(a.valuation(v) >= 0 for a in self.a_invariants()):
            return True
        else:
            return False
            
    def local_integral_model(self, v):
        r"""Return a model integral at the place v.
        INPUT:
        - ``v`` -- a place of the base rational function field
        OUTPUT:
        An elliptic curve isomorphic to ``self`` whose Weierstrass coefficients
        all have non-negative valuation at ``v``.
        """

        vals = [a.valuation(v) for a in self.a_invariants()]

        # Compute scaling factor k = v(u):
        # v(a_i) - i*k >= 0 for all i  =>  k <= v(a_i)/i
        e = min(vals[i] / [1, 2, 3, 4, 6][i] for i in range(5)).floor()

        if e >= 0:
            return self
        # Uniformizer at v (element with valuation 1) u = 1/pi^{|k|}
        u = (v.local_uniformizer()) ** e

        return self.change_weierstrass_model([u, 0, 0, 0])    
        
    def is_local_minimal(self, v):
        r"""
        Determines if model is minimal at place v. 
        If no places are stated, a list of tuples is given with the boolean evaluation and the place
        """
        
        _E = self
        if _E.is_local_integral(v) == False:
            return False

        if _E.base_field().characteristic() != 2 or _E.base_field().characteristic() != 3:
            if _E.discriminant().valuation(v) < 12:
                return True
            elif _E.c4().valuation(v) < 4 or _E.c6().valuation(v) < 6:
                return True
            else:
                return False
        else:
            raise ValueError("Field characteristic must not be 2 or 3")
    
    def local_minimal_model(self, v):
        r"""Returns minimal model if self is not already minimal.
        
        INPUT:
        - ``v`` -- a place of the base rational function field
        OUTPUT:
        An elliptic curve isomorphic to ``self`` whose Weierstrass coefficients
        have minimal valuation at ``v``.
        """    
        _E = self.local_integral_model(v)
        
        if _E.is_local_minimal(v):
            return _E
        
        c4 = _E.c4()
        c6 = _E.c6()
        D = _E.discriminant()
        
        if _E.base_field().characteristic() == 2 or _E.base_field().characteristic() == 3:
            raise ValueError("Field characteristic must not be 2 or 3")
        
        if c4 != 0 and c6 != 0:
            q = max(math.ceil((c4.valuation(v))/4), math.ceil(c6.valuation(v)/6))
        elif c4 == 0:
            q = math.ceil(-c6.valuation(v)/6)
        elif c6 == 0:
            q = math.ceil(-c4.valuation(v)/4)
        else:
            q = 0
        
        return self.change_weierstrass_model([v.local_uniformizer()**-q, 0, 0, 0])
