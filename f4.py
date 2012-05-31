# -*- Mode: Python -*-
# # vi:si:et:sw=4:sts=4:ts=4
#

"""
Clean F4 Design

AUTHOR: Martin Albrecht <martinralbrecht@googlemail.com>
"""

from sage.misc.misc import exists
from sage.rings.ideal import is_Ideal

class Reduction:
    def __init__(self, F):
        F = Sequence([m*f for m,f in F])
        self.M = set(f.lm() for f in F)
        self.A, self.v = F.coefficient_matrix(sparse=False)
        self.called = False

    def __call__(self):
        self.called = True
        self.A.echelonize()

    def reduced_polynomials(self):
        if self.called is False:
            self()
        F = (self.A*self.v).list()
        return Sequence([f for f in F if (f and f.lm() not in self.M)])

def LM(F):
    if isinstance(F,(list,set,tuple)):
        return set([m*f.lm() for m,f in F])
    else:
        return F.lm()

class F4_orig:
    def __init__(self):
        pass

    def __call__(self, F, Sel=None):
        if Sel is None:
            Sel = self.normal_strategy

        if is_Ideal(F):
            F = F.gens()

        self.ring = F[0].parent()
        G = list(F)
        F0p = F
        i = 0
        P = set([self.pair(f,g) for f in G for g in G if f<g ] )

        while P != set():
            i += 1
            Pd, d = Sel(P)
            P = P.difference(Pd)
            Ld = set(self.left(Pd)).union(set(self.right(Pd)))
            Fdp = self.reduction(Ld,G)
            for h in Fdp:
                P = P.union(set([self.pair(h,g) for g in G ]))
                G.append(h)
        return Sequence(G)

    def reduction(self,L,G):
        F = self.symbolic_preprocessing(L, G)
        return self.row_echelon(F)

    def symbolic_preprocessing(self, L, G):
        F = L
        Done = LM(F)
        M = set([m*t for (m,f) in F for t in f.monomials()])
        R = self.ring
        while M != Done:
            m = M.difference(Done).pop()
            Done.add(m)
            t,g = self.ring.monomial_reduce(m,G)
            if t != 0:
                F.add( (t,g) )
            M = set([m*t for (m,f) in F for t in f.monomials()])
        return F

    def pair(self,f,g):
        lcm = self.ring.monomial_lcm(f.lm(), g.lm())
        return (lcm,f,g)

    def left(self,p):
        s = set()
        for f in p:
            s.add((self.ring.monomial_quotient(f[0],f[1].lm()),f[1]))
        return s

    def right(self,p):
        s = set()
        for f in p:
            s.add((self.ring.monomial_quotient(f[0],f[2].lm()),f[2]))
        return s

    def row_echelon(self, F):
        R = Reduction(F)
        R()
        return R.reduced_polynomials()

    def normal_strategy(self,P):
        d = min(set([ lcm.total_degree() for (lcm,fi,fj) in P ]))
        return set([ (lcm,fi,fj) for (lcm,fi,fj) in P if lcm.total_degree()==d]), d

    def update_pairs(self,G,B,h):
        R = self.ring

        # if G is a set then C only contains unique elements
        C = [self.pair(h,g) for g in G]
        D = list() # only adding elements of C, thus unique

        # Criterion M

        while C!=list():
            (lcmhg1,h,g1) = C.pop()

            lcm_divides = lambda lcmhg2: R.monomial_divides(  lcmhg2[0], lcmhg1 )

            # if LM(h) and LM(g_1) are disjoint
            if R.monomial_pairwise_prime(h.lm(),g.lm()) or \
               (\
                   not exists(C, lcm_divides )[0] \
                   and \
                   not exists(D, lcm_divides )[0]\
                ):
                D.append((lcmhg1,h,g1))

        E = list() #only adding elements of D, thus unique

        # Criterion F

        while D != list():
            (lcmhg,h,g) = D.pop()
            # if LM(h) and LM(g) are not disjoint
            if not R.monomial_pairwise_prime(h.lm(),g.lm()):
                E.append((lcmhg,h,g))

        B_new = set()

        # Criterion B_k

        while B != set():
            lcmg1g2,g1,g2 = B.pop()
            if not self.ring.monomial_divides( h.lm(), lcmg1g2 ) or \
                   self.ring.monomial_lcm(g1.lm(), h.lm()) == lcmg1g2 or \
                   self.ring.monomial_lcm( h.lm(),g2.lm()) == lcmg1g2 :
                B_new.add((lcmg1g2,g1,g2))

        B_new = B_new.union(E)

        G_new = list()

        while G != list():
            g = G.pop()
            if not R.monomial_divides( h.lm(), g.lm() ):
                G_new.append(g)

        G_new.append(h)

        return G_new,B_new

    def update_simple(self,G, P, h):
        return G+[h],P.union([self.pair(g,h) for g in G])


class F4(F4_orig):
    def __call__(self, F, Sel=None, Update=None):
        if is_Ideal(F):
            F = F.gens()

        # pretty looking code
        Left = self.left
        Right = self.right
        Reduction = self.reduction
        first = self.first
        if Sel is None:
            Sel = self.normal_strategy
        if Update is None:
            Update = self.update_pairs

        self.ring = F[0].parent()
        self.term_order = self.ring.term_order()

        # We maintain a list of dictionaries which contain f.lm() => f
        # maps for the sets $F_j^~$ to allow O(1) lookups for this code:
        #"$F_j^~$ is the row echelon form of F_j w.r.t. < there exists a
        # (unique) $p \in F_j^~ such that LM(p) = LM(u*f)"
        self.Ftd = [[]]

        F = list(F)
        Fd = dict()

        G = list()
        P = set()

        while F:
            f = first(F)
            F.remove(f)
            G,P = Update(G,P,f)

        while P:
            Pd, d = Sel(P)
            P = P.difference(Pd)
            Ld = Left(Pd).union( Right(Pd) )
            Fdp, Fd[i] = Reduction(Ld, G, Fd)
            for h in Fdp:
                G,P = Update(G,P,h)
        return Sequence(G)

    def reduction(self, L, G, Fset):
        F = self.symbolic_preprocessing(L,G,Fset)

        Ftp = self.row_echelon(F)

        return Ftp,F

    def symbolic_preprocessing(self,L,G,Fset):
        Simplify = self.simplify
        R = self.ring

        F = set([Simplify(m,f,Fset) for (m,f) in L])

        Done = LM(F)

        M = set([m*t for (m,f) in F for t in f.monomials()])

        MdivDone = M.difference(Done)
        G = tuple(G)

        while MdivDone:
            m = MdivDone.pop()
            Done.add(m)
            t,g = self.ring.monomial_reduce(m,G)
            if t != 0:
                t,g = Simplify(t,g,Fset)
                F.add((t,g))
                for tgm in (t*g).monomials():
                    M.add(tgm)
                    if tgm not in Done:
                        MdivDone.add(tgm)
        return F

    def simplify(self,t,f,F):
        for u in sorted(self.ring.monomial_all_divisors(t), reverse=True):
            uf = u*f
            for j in F:
                if uf in F[j]:
                    # F~_j is the row echelon form of F_j w.r.t. <
                    # there exists a (unique) p \in F~_j such that LM(p) = LM(u*f)
                    p = self.Ftd[j][uf.lm()]
                    if u != t:
                        return self.simplify(self.ring.monomial_quotient(t,u),p,F) #t/u
                    else:
                        return (self.ring(1),p)
        return (t,f)

    def first(self,G):
        mg = G[0]
        mm = mg.lm()
        for g in G:
            if g.lm() > mm:
                mm = g.lm()
                mg = g
        return mg

f4 = F4()
