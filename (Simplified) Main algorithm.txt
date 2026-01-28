# sage: Genus 2 Curve Construction - Algorithm 23 in the paper.

from sage.all import *
from itertools import product
import sys


def get_reciprocal_brandt(f):  
 # Instead of defining this, one can use an implemented code from SAGE (namely .reciprocal), 
 #but we did not want to do this to avoid any confusion of the aggrement of the reciprocal form.
 # Note that all these concepts with the same notations are discussed in Section 2 of the paper.
    """
    Computes F_f^B (Brandt's reciprocal form) and the content I1 (the basic genus invariant).
    Returns: (TernaryQF, I1)
    """
    M_adj = f.matrix().adjugate()
    
    coeffs = [M_adj[0,0], M_adj[1,1], M_adj[2,2], 
              2*M_adj[1,2], 2*M_adj[0,2], 2*M_adj[0,1]]
              
    I1 = gcd(coeffs)
    
    # Note: If f is positive definite, M_adj is positive definite.
    # The reciprocal form F should be positive definite.
    # We essentially return F = adj(f) / I1.
    return TernaryQF([c // I1 for c in coeffs]), I1

def step2_find_prime(f, D, B=100):
    """
    Step 2: Searches for a prime p represented by f.
    Uses shell search (expanding cube) to avoid generating all vectors at once.
    IMPORTANT: We just put a default bound B=100. For the larger calculations, one needs to increase this! The (implicit) bound is discussed in the paper.
    Conditions:
    1. p does not divide D.
    2. D is a quadratic residue mod p (kronecker != -1). (We want to check this to make sure we have a valid prime for the next step. 
    However, this is not necessary because we know this condition is guaranteed for the theoratical reasons, so this condition just helps if we did some mistakes before.
    """
    coeffs = f.coefficients()
    a, b, c, r, s, t = coeffs
    
    found_primes = set()
    
    # Iterate through shells of increasing radius `rad`
    for rad in range(1, B + 1):
        shell_vecs = []
        r_range = range(-rad, rad + 1)
        r_inner = range(-rad + 1, rad)
        
        # Build shell faces (vectors having at least one coord equal to +/- rad)
        # Face x = +/- rad
        shell_vecs.extend([(rad, y, z) for y in r_range for z in r_range])
        shell_vecs.extend([(-rad, y, z) for y in r_range for z in r_range])
        # Face y = +/- rad (excluding x = +/- rad)
        shell_vecs.extend([(x, rad, z) for x in r_inner for z in r_range])
        shell_vecs.extend([(x, -rad, z) for x in r_inner for z in r_range])
        # Face z = +/- rad (excluding x,y = +/- rad)
        shell_vecs.extend([(x, y, rad) for x in r_inner for y in r_inner])
        shell_vecs.extend([(x, y, -rad) for x in r_inner for y in r_inner])
        
        shell_vecs.sort(key=lambda v: v[0]**2 + v[1]**2 + v[2]**2)
        
        for x, y, z in shell_vecs:
            try:
                val = a*x**2 + b*y**2 + c*z**2 + r*y*z + s*x*z + t*x*y
                val = Integer(val)
                
                if val > 0 and val.is_prime() and val != 2: #of course, val is necessarily positive at nonzero vector.
                    if val in found_primes:
                        continue
                    found_primes.add(val)
                    
                    # Condition 1:
                    if gcd(val, D) == 1:
                         # Condition 2: Quadratic Residue (for the next step)
                         if kronecker(D, val) != -1:
                             return (val, x, y, z)
            except (TypeError, ValueError, ArithmeticError):
                continue
            
    return None

def step3_construct_binary(p, d):
    """
   Constructs a binary form [p, b, c] of discriminant d by following the proof of Lemma 19 of the paper.
    """
    # Check quadratic residue (we have already done this above).
    if kronecker(d, p) == -1:
        raise ValueError(f"{d} not quadratic residue mod {p}")
        
    # Solve b^2 = d (mod p)
    try:
        b = ZZ(GF(p)(d).sqrt())
    except:
        raise ValueError(f"Sqrt failed for {d} mod {p}")
    
    # Adjust parity: b must have same parity as d
    if b % 2 != d % 2:
        b = p - b
        
    # Compute c
    num = b**2 - d
    if num % (4*p) != 0:  # This is not possible by the theory!
         raise ValueError(f"b^2 - d not divisible by 4p (b={b}, d={d}, p={p})")

    c = num // (4*p)
    return BinaryQF([p, b, c])

def step_4_thm_38(f, vec):
    """
    Constructs primitive transformation T and binary form phi = f(Tx) by following the proof of Theorem 38 of Dickson's book; this is Algorithm 16 in our paper.
    """
    x0, y0, z0 = vec
    
    if z0 != 0:
        g = gcd(z0, x0 + y0)
        a1 = z0 // g
        a2 = z0 // g
        a3 = -(x0 + y0) // g
        
        d1, s1, t1 = xgcd(a1, a2)
        d2, s2, t2 = xgcd(d1, a3)
        
        sign = 1 if d2 == 1 else -1
        s2 *= sign; t2 *= sign
        
        l = s2 * s1; m = s2 * t1; n = t2
        
        b1 = a1 - m * z0 + n * y0
        b2 = a2 - n * x0 + l * z0
        b3 = a3 - l * y0 + m * x0
        
        T = matrix(ZZ, 3, 2, [a1, b1, a2, b2, a3, b3])
        
    else:
        # To deal with zero cases by permutation when z0=0.
        if y0 != 0:
            P = matrix(ZZ, 3, 3, [[1,0,0],[0,0,1],[0,1,0]])
            vec_perm = (x0, z0, y0)
            coeffs_perm = _permute_coeffs(f, P)
            T_perm, phi = step_4_thm_38(TernaryQF(coeffs_perm), vec_perm)
            T = P.transpose() * T_perm
        else:
            P = matrix(ZZ, 3, 3, [[0,0,1],[0,1,0],[1,0,0]])
            vec_perm = (0, 0, x0) # If z0=0 and y0=0, necessarily x0 is not 0.
            coeffs_perm = _permute_coeffs(f, P)
            T_perm, phi = step_4_thm_38(TernaryQF(coeffs_perm), vec_perm)
            T = P.transpose() * T_perm

    M_f = f.matrix()
    M_phi = T.transpose() * M_f * T
    phi = BinaryQF([M_phi[0,0]//2, M_phi[0,1], M_phi[1,1]//2])
    return T, phi

def _permute_coeffs(f, P):
    M = f.matrix()
    M_perm = P * M * P.transpose()
    return [M_perm[0,0]//2, M_perm[1,1]//2, M_perm[2,2]//2, 
            M_perm[1,2], M_perm[0,2], M_perm[0,1]]

def step_5_find_triple_s(phi, d):
    """
    Finds s = (n, m, k) in P(d) such that q_s ~ phi by following the proof of Theorem 17 of Kani's paper as in Algorithm 18.
    """
    N = None
    u, v = 0, 0
    
    # 1. Find N such that phi primitively represents N^2 and gcd(N,d)=1
    search_range = 100 #IMPORTANT: We put this bound as a default. One may need to increase this. By Proposition 17 in the paper, we have an absolute bound dependind on d.
    vecs = sorted(list(product(range(-search_range, search_range+1), repeat=2)), key=lambda p: p[0]**2 + p[1]**2)
    
    for x, y in vecs:
        if gcd(x, y) != 1: continue
        val = phi(x, y)
        if val > 0 and Integer(val).is_square():
            N_cand = ZZ(Integer(val).sqrt())
            if gcd(N_cand, d) == 1:
                N = N_cand; u, v = x, y
                break
                
    if N is None:
        raise ValueError("Could not find suitable N^2 represented by phi; try after increasing the search range.")
        
    # 2. Transform phi
    g, z, w = xgcd(u, -v)
    M = matrix(ZZ, 2, 2, [u, v, w, z])
    phi_bar = phi.matrix_action_left(M)
    
    A_bar, B_bar, C_bar = phi_bar[0], phi_bar[1], phi_bar[2]
    b, c = B_bar // 2, C_bar
    
    # 3. Calculate s=(n,m,k)
    modulus = N**2
    if N % 2 != 0:
        inv_2d = inverse_mod(2*d, modulus)
        k = (inv_2d * b) % modulus
        if k == 0: k = modulus
        m = (k**2 * d + 1) // N
        return (N, m, abs(k)) # note that k and -k lead to the equivalent form q_s.
    else:
        d_star = inverse_mod(d, modulus)
        c_bar = c % 4
        term = b//2 + (c_bar * N)//2
        k = (d_star * term)
        
        m_numerator = k**2 * d + 1
        if m_numerator % N != 0: #This is not possible.
             k_mod = k % modulus
             if (k_mod**2 * d + 1) % N == 0: k = k_mod
        
        m = (k**2 * d + 1) // N
        return (N, m, abs(k))

def verify_result(s, q_tilde, f_target): #Optional step. By using Proposition 6, we can make sure whether our output is correct or not.
    n, m, k = s
    a, b, c = q_tilde[0], q_tilde[1], q_tilde[2]
    mn = m * n
    coeffs = [n**2, (m**2)*(mn+3)*a, (b**2)*(k**2) + 4*c, 
              -2*b*m*(mn+1), -2*b*k*n, 2*k*a*(mn+2)]
    q_A = TernaryQF(coeffs)
    #To check the equivalence, we use the implemented code in SAGE, that is the Eisenstein reduced form.
    try:
        equiv = (q_A.reduced_form_eisenstein()[0] == f_target.reduced_form_eisenstein()[0]) 
    except:
        equiv = False
    return equiv, q_A

# ==============================================================================
# Main (Simplified) Algorithm for Algorithm 23 in the paper:
# ==============================================================================

def solve_primitive(f):
    """For geometric primitive forms."""
    # 1. D = disc(f) / 16
    sage_disc = f.disc() #our convention of the sign of the discriminant is not the same with the one in SAGE.
    if f.is_positive_definite() and sage_disc > 0: sage_disc = -sage_disc
    D = sage_disc // 16
    
    # 2. Find p using Brandt's reciprocal
    F_recip, I1 = get_reciprocal_brandt(f)
    
    # Calculate kappa (I1 = 16 * kappa for geometric primitive forms)
    # For geometric forms, I1 is divisible by 16, so kappa >= 1.
    kappa = I1 // 16
    
    D_prime = D // (kappa**2)

    # Since D = kappa^2 * D', checking D is a residue is equivalent to checking D' is a residue
    # provided p does not divide kappa (guaranteed by gcd(p, D) == 1).
    res = step2_find_prime(F_recip, D, B=100) #IMPORTANT: We put a search bound B=100. One may need larger than this as was mentioned above.
    if not res: return None, "Step 2 failed"
    p, x, y, z = res
    
    # 3. Construct binary form
    # Construct primitive form q' of discriminant D'
    try:
        q_tilde_prime = step3_construct_binary(p, D_prime)
    except ValueError: return None, "Step 3 failed"
    
    # The form on ExE' is found by kappa*q_tilde_prime 
    q_tilde = BinaryQF([kappa*q_tilde_prime[0], kappa*q_tilde_prime[1], kappa*q_tilde_prime[2]])
    
    # 4. Transform to phi_p
    T, phi_p = step_4_thm_38(f, (x, y, z))
    
    # 5. Find s (Target: phi_p, Type: kappa*p)
    try:
        s = step_5_find_triple_s(phi_p, kappa*p)
    except ValueError: return None, "Step 5 failed"
    
    # 6. Verify
    ver, _ = verify_result(s, q_tilde, f)
    return (q_tilde, p, s, ver), "Success"

def solve_imprimitive(f):
    """For geometric imprimitive forms."""
    f_div_4 = TernaryQF([c//4 for c in f.coefficients()])
    f_div_2 = TernaryQF([c//2 for c in f.coefficients()])
    
    # 1. Calc invariants
    # kappa = I1(f/4). Using get_reciprocal_brandt on f/4 to get correct I1 (kappa).
    F_recip_f4, kappa = get_reciprocal_brandt(f_div_4)
    
    sage_disc = f.disc()
    if f.is_positive_definite() and sage_disc > 0: sage_disc = -sage_disc
    D = sage_disc // 16
    D_prime = D // (kappa**2)
    
    # 2. Find p using Reciprocal of f/4
    # Condition: p doesn't divide D. D_prime must be quadratic residue mod p.
    # Note: D = kappa^2 * D_prime. Checking residue of D is equivalent to checking D_prime.
    res = step2_find_prime(F_recip_f4, D, B=100)
    if not res: return None, "Step 2 failed"
    p, x, y, z = res
    
    # 3. Construct q_tilde (using D_prime)
    try:
        q_tilde_prime = step3_construct_binary(p, D_prime)
    except ValueError: return None, "Step 3 failed"
    
    # 4. Transform f/2 to phi_p
    try:
        T, phi_p = step_4_thm_38(f_div_2, (x, y, z))
    except ValueError: return None, "Step 4 failed"
    
    # 5. Find s (Target: 2*phi_p, Type: kappa*p)
    target_phi = BinaryQF([2*phi_p[0], 2*phi_p[1], 2*phi_p[2]])
    type_d = kappa * p
    try:
        s = step_5_find_triple_s(target_phi, type_d)
    except ValueError: return None, "Step 5 failed"
    
    # 6. Verify
    # The form on ExE' is found by kappa*q_tilde_prime
    q_final_check = BinaryQF([kappa*q_tilde_prime[0], kappa*q_tilde_prime[1], kappa*q_tilde_prime[2]])
    ver, _ = verify_result(s, q_final_check, f)
    
    return (q_final_check, p, s, ver), "Success"

def algorithm_23(f_coeffs): #THE MAIN ALGORITHM
    f = TernaryQF(f_coeffs)
    cont = f.content()
    
    if cont == 1:
        type_str = "Primitive"
    elif cont == 4:
        type_str = "Imprimitive"
    else:
        type_str = f"Unknown (Content {cont})" #Note that the geometric form has content 1 or 4.
    
    print(f"\n>>> Given Ternary Form: {f_coeffs} [{type_str}]")
    
    result = None
    msg = ""
    
    try:
        if cont == 1:
            result, msg = solve_primitive(f)
        elif cont == 4:
            result, msg = solve_imprimitive(f)
        else:
            print(f"  Skipping: Content is {cont} (Expected 1 or 4).")
            return
            
        if result:
            q, p, s, ver = result
            # Calculate D manually to avoid the discriminant convention differences
            D = q[1]**2 - 4*q[0]*q[2]
            print(f"  D={D}, \\tilde{{q}}=[{q[0]}, {q[1]}, {q[2]}], s={s} and verified {ver}")
        else:
            print(f"  Failed: {msg}")
            
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"  Error: {e}")

def run_tests():
    examples = [
        # Example 6.1 from paper
        [4, 12, 28, 0, 4, 4],
        
        # Elliptic Subcovers Examples (from discussion on p.20)
        [4, 4, 5, 4, 4, 4],
        [4, 4, 9, 4, 4, 4],
        [4, 4, 5, 0, 0, -4],
        [4, 4, 9, 0, 0, -4]
    ]
    
    print(f"Running {len(examples)} Specific Test Cases from Paper...")
    for f in examples:
        algorithm_23(f)

if __name__ == "__main__":
    run_tests()
