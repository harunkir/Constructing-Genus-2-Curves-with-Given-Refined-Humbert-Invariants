# Checks if a given integral ternary quadratic form is geometric.
# Includes verification against the 199 examples from the author's PhD thesis. These examples were used for the main algorithm as well.
#IMPORTANT: WE DID NOT DISCUSS THE TERMINATION OF THIS CODE IN THE PAPER. 
#Indeed, the code is searching a perfect square n^2 represented by a given ternary form f (where gcd(n,disc(f))=1). We do not know how much we have to search for such a number.
#That's why this code was not discussed in the paper as an algorithm. However, one can still use this code 
#because if the code stops for the given (heuristic) bound, then it returns whether the given form is geometric or not.
#As a default, we search a suitable number for the vectors (x,y,z), where -10 \leq x,y,z \leq 10.
from sage.all import *
from itertools import product


# Global cache for search vectors (shared optimization)
_SEARCH_VECTORS = None
_CACHED_BOUND = 0

def get_search_vectors(B):
    global _SEARCH_VECTORS, _CACHED_BOUND
    if _SEARCH_VECTORS is None or _CACHED_BOUND != B:
        vecs = list(product(range(-B, B+1), repeat=3))
        # Sort by norm to check smaller vectors first
        vecs.sort(key=lambda v: v[0]**2 + v[1]**2 + v[2]**2)
        _SEARCH_VECTORS = vecs
        _CACHED_BOUND = B
    return _SEARCH_VECTORS

def is_geometric_form(f):
    """
    Decides if a given integral ternary form is a geometric form based on Theorem 4 
    of the paper as follows:
    
    A positive definite ternary form f is geometric if:
    1. It is Primitive (content=1) AND:
       a) f(x,y,z) = 0, 1 (mod 4) for all integers x,y,z.
       b) f represents a square n^2 coprime to disc(f).
       
    2. It is Imprimitive (content=4) AND:
       a) f/2 is improperly primitive.
       b) f represents a square (2n)^2 with gcd(n, disc(f)) = 1.
    """
    # 1. Positive Definite
    if not f.is_positive_definite():
        return False
    
    c = f.content()
    
    # Primitive Case (Content = 1)
    if c == 1:
        # Condition 3.3: f(x,y,z) = 0, 1 (mod 4)
        # A simple Fact from Lemma 5.3.6 of the author's PhD thesis: q=[a,b,c,r,s,t] = 0,1 mod 4 if and only if 
        # a,b,c = 0,1 mod 4 AND r=2bc mod 4, s=2ac mod 4, t=2ab mod 4.
        coeffs = f.coefficients() # [a, b, c, r, s, t]
        a, b, c_coeff, r, s, t = [Integer(x) for x in coeffs]
        
        if not (a % 4 in [0, 1] and b % 4 in [0, 1] and c_coeff % 4 in [0, 1]):
            return False
            
        if r % 4 != (2 * b * c_coeff) % 4:
            return False
        if s % 4 != (2 * a * c_coeff) % 4:
            return False
        if t % 4 != (2 * a * b) % 4:
            return False
                
        # Condition 3.4: Represents n^2 coprime to disc
        D = f.disc() 
        found_square = False
        vecs = get_search_vectors(10) # IMPORTANT: ONE MAY NEED TO INCREASE WHEN NO FOUND SQUARE CASE.
        for v in vecs:
            if v == (0,0,0): continue
            val = f(v)
            if val.is_square():
                n = ZZ(val.sqrt())
                if gcd(n, D) == 1:
                    found_square = True
                    break
        if not found_square:
            return False
            
        return True
        
    # Imprimitive Case (Content = 4)
    elif c == 4:
        # Condition 3.5: f/2 is improperly primitive
        f_div_4 = TernaryQF([x/4 for x in f.coefficients()])
        if f_div_4.content() != 1:
            return False
            
        # Condition 3.6: Represents (2n)^2 with gcd(n, disc(f)) = 1
        D = f.disc()
        found_square = False
        vecs = get_search_vectors(10)
        for v in vecs:
            if v == (0,0,0): continue
            val = f(v)
            if val.is_square():
                root = ZZ(val.sqrt())
                if root % 2 == 0:
                    n = root // 2
                    if gcd(n, D) == 1:
                        found_square = True
                        break
        if not found_square:
            return False
            
        return True
        
    else:
        # Geometric forms must have content 1 or 4
        return False

# ==============================================================================
# Verification for some examples
# ==============================================================================

# List of known geometric forms from the author's thesis
forms_to_check = [
    [4, 4, 5, 4, 4, 4], [4, 4, 9, 4, 4, 4], [4, 4, 5, 0, 0, -4], [4, 4, 9, 0, 0, -4],
    [4, 8, 61, 0, 0, -4], [4, 13, 32, 0, -4, 0], [4, 4, 89, 0, 0, -4], [4, 5, 48, 0, -4, 0],
    [4, 8, 29, 0, 0, 0], [4, 12, 17, 0, 0, -4], [4, 4, 37, 0, 0, 0], [4, 4, 49, 0, 0, -4],
    [4, 4, 41, 0, 0, -4], [4, 5, 24, 0, -4, 0], [4, 5, 28, 0, 0, -4], [4, 4, 25, 0, 0, 0],
    [4, 9, 12, 0, -4, 0], [4, 8, 13, 0, 0, -4], [4, 8, 13, -8, 0, 0], [4, 4, 25, 0, 0, -4],
    [4, 8, 9, 0, 0, 0], [4, 5, 17, 2, 4, 4], [4, 5, 12, 0, 0, 0], [4, 4, 13, 0, 0, 0],
    [4, 4, 17, 0, 0, -4], [4, 5, 12, 0, 0, -4], [4, 5, 8, 0, 0, 0], [4, 4, 9, 0, 0, 0],
    [4, 5, 8, 0, -4, 0], [4, 5, 9, 2, 4, 4], [4, 5, 5, -2, 0, 0], [4, 4, 5, 0, 0, 0],
    [4, 4, 5, 0, 0, -4], [4, 4, 60, -4, -4, 0], [4, 8, 32, -8, -4, 0], [4, 8, 16, 0, 0, -4],
    [4, 12, 12, 8, 4, 4], [4, 4, 24, -4, -4, 0], [4, 8, 12, 0, -4, 0], [4, 4, 20, -4, -4, 0],
    [4, 8, 12, -8, -4, 0], [4, 8, 8, -4, 0, 0], [4, 4, 20, 0, 0, -4], [4, 4, 16, 0, 0, -4],
    [4, 8, 8, 4, 4, 4], [4, 4, 12, -4, -4, 0], [4, 8, 8, 8, 4, 4], [4, 4, 8, -4, -4, 0],
    [4, 4, 8, 0, 0, -4], [5, 8, 17, -4, -2, 0], [5, 5, 28, 0, -4, -2], [5, 8, 20, -4, -4, -4],
    [5, 12, 12, -4, 0, -4], [5, 12, 13, 8, 2, 4], [8, 8, 13, 4, 8, 4], [5, 5, 12, 0, -4, -2],
    [5, 8, 8, -4, 0, -4], [5, 5, 8, 0, -4, -2], [5, 20, 20, 12, 4, 4], [12, 12, 13, -4, -4, -4],
    [5, 20, 20, -12, -4, -4], [8, 13, 16, 0, -4, 0], [5, 12, 21, -12, -2, 0], [5, 5, 45, -2, -2, -2],
    [8, 8, 21, 8, 8, 4], [9, 12, 13, -12, -2, 0], [12, 12, 13, 8, 12, 12], [5, 8, 24, -4, 0, 0],
    [5, 12, 16, -4, 0, 0], [8, 8, 17, -4, -4, -4], [9, 12, 12, 4, 8, 8], [5, 8, 24, 0, -4, 0],
    [5, 5, 40, 4, 4, 2], [8, 8, 17, -4, -8, 0], [8, 9, 13, -2, 0, 0], [8, 12, 13, -8, -8, 0],
    [5, 9, 20, -8, -4, -2], [8, 9, 13, 2, 8, 4], [5, 8, 17, -8, -2, 0], [5, 5, 25, 2, 2, 2],
    [8, 9, 12, -8, -8, 0], [8, 8, 13, 4, 4, 8], [5, 12, 13, -12, -2, 0], [5, 5, 25, -2, -2, -2],
    [8, 8, 13, 8, 8, 4], [5, 5, 21, -2, -2, -2], [5, 12, 12, 12, 4, 4], [8, 8, 9, 4, 4, 4],
    [5, 8, 12, -4, 0, 0], [8, 8, 9, -4, -4, -4], [5, 5, 17, 2, 2, 2], [8, 8, 9, 4, 4, 8],
    [5, 9, 12, -8, -4, -2], [5, 8, 13, 8, 2, 4], [5, 8, 9, 0, -2, 0], [5, 5, 16, 4, 4, 2],
    [5, 5, 13, -2, -2, -2], [5, 8, 8, 0, 0, -4], [5, 8, 8, 0, -4, -4], [5, 5, 9, 2, 2, 2],
    [5, 8, 8, 8, 4, 4], [5, 5, 9, -2, -2, -2], [5, 5, 8, 4, 4, 2], [8, 13, 48, 4, 8, 8],
    [8, 24, 25, 4, 8, 8], [5, 12, 32, -12, 0, 0], [8, 8, 29, 0, 0, -4], [12, 12, 29, 0, 0, -8],
    [12, 12, 33, -8, -8, -8], [5, 28, 48, -28, 0, 0], [12, 12, 41, 0, 0, -4], [5, 5, 68, 0, 0, -2],
    [8, 12, 17, 0, 0, 0], [5, 5, 64, -4, -4, -2], [8, 8, 33, 4, 8, 8], [5, 8, 36, -8, 0, 0],
    [8, 12, 17, 0, 0, -8], [5, 16, 16, -4, 0, 0], [9, 12, 12, -4, 0, 0], [5, 5, 52, 0, 0, -2],
    [8, 12, 13, 0, 0, 0], [9, 12, 12, -8, 0, 0], [12, 12, 13, -8, -8, -8], [5, 12, 56, -12, 0, 0],
    [8, 8, 53, 0, 0, -4], [5, 8, 28, 0, 0, 0], [8, 9, 20, 0, 0, -8], [5, 8, 76, 0, 0, 0],
    [8, 20, 21, 0, -8, 0], [5, 12, 16, 0, 0, 0], [9, 9, 12, 0, 0, -2], [5, 5, 40, -4, -4, -2],
    [8, 8, 21, 4, 8, 8], [5, 24, 24, -4, 0, 0], [13, 16, 16, -12, 0, 0], [5, 5, 120, -4, -4, -2],
    [8, 8, 61, 4, 8, 8], [5, 12, 16, -12, 0, 0], [8, 8, 13, 0, 0, -4], [5, 8, 20, -8, 0, 0],
    [8, 9, 12, 0, -8, 0], [5, 5, 28, 0, 0, -2], [8, 9, 12, 0, 0, -8], [5, 12, 12, -8, 0, 0],
    [9, 9, 12, -8, -8, -2], [5, 5, 24, -4, -4, -2], [8, 8, 13, 4, 8, 8], [5, 8, 12, 0, 0, 0],
    [5, 5, 20, 0, 0, -2], [5, 5, 16, -4, -4, -2], [8, 8, 9, 4, 8, 8], [5, 24, 24, -20, 0, 0],
    [12, 12, 17, 0, 0, -4], [5, 12, 40, -12, 0, 0], [8, 8, 37, 0, 0, -4], [8, 9, 40, 4, 8, 8],
    [8, 16, 21, 4, 8, 8], [5, 8, 52, 0, 0, 0], [8, 13, 20, 0, 0, 0], [8, 8, 20, 4, 4, 8],
    [8, 12, 12, -4, 0, -8], [8, 8, 8, 0, -4, -4], [8, 8, 8, 4, 4, 8], [8, 12, 20, -12, 0, 0],
    [8, 8, 36, 4, 8, 8], [8, 12, 16, -12, 0, 0], [8, 8, 28, 4, 8, 8], [8, 12, 12, -4, 0, 0],
    [8, 12, 16, 4, 8, 8], [8, 20, 24, -20, 0, 0], [8, 12, 40, 4, 8, 8], [8, 8, 16, 0, 0, -4],
    [12, 12, 12, 8, 12, 12], [8, 8, 16, 4, 8, 8], [8, 12, 12, 12, 8, 8], [8, 8, 8, 0, 0, -4],
    [8, 8, 12, 4, 8, 8], [8, 12, 28, 4, 8, 8], [8, 20, 20, 20, 8, 8], [8, 8, 9, 0, -8, 0],
    [8, 8, 57, 8, 8, 8], [8, 8, 25, 8, 8, 8], [8, 8, 17, 8, 8, 8], [12, 12, 13, 12, 12, 12],
    [8, 8, 9, 8, 8, 8], [9, 13, 13, -10, -6, -6], [8, 8, 41, -8, -8, 0], [8, 8, 25, 0, 0, 0],
    [8, 8, 17, -8, -8, 0], [8, 8, 9, 0, 0, 0], [8, 8, 9, -8, -8, 0], [8, 8, 89, 0, 0, -8],
    [12, 12, 49, 0, 0, -12], [8, 8, 49, 0, 0, -8], [8, 8, 41, 0, 0, -8], [12, 12, 25, 0, 0, -12],
    [8, 8, 25, 0, 0, -8], [8, 8, 17, 0, 0, -8], [12, 12, 16, 0, 0, -12]
]

def run_check():
    print(f"Checking {len(forms_to_check)} forms...")
    passed = 0
    failed = 0
    
    for i, coeffs in enumerate(forms_to_check):
        f = TernaryQF(coeffs)
        if is_geometric_form(f):
            # print(f"Form #{i+1} {coeffs}: OK")
            passed += 1
        else:
            print(f"Form #{i+1} {coeffs}: FAIL (Not geometric)")
            failed += 1
            
    print(f"\nResult: {passed} passed, {failed} failed.")

if __name__ == "__main__":
    run_check()
