// Quadratic modelings of binary syndrome decoding: given H in F_2^{(n-k) x n},
// s in F_2^{n-k} and t, find x in F_2^n with H x^T = s^T and wt(x) = t.
//
// Each modeling takes an instance as read by ReadInstance (instance.m) and
// returns a system of polynomials over GF(2) in grevlex order whose common
// zeros give exactly the solutions:
//   CCMMP2, RCI, CSA  a weight counter of bit registers, checked against t;
//   RST               the support of x as t ordered positions (also returns X);
//   LMTC, LMTQFC      the C-Hamming and QFC-Hamming modelings of La Scala,
//                     Marchesin and Tiwari.
// For RST the error vector is read from the X it returns; for the others it
// is the first n variables of the ring.


// Equations H * x^T = s^T, one per row of H. X is a sequence (or vector)
// of polynomials in the system's ring.
function ParityCheckEquations(X, H, s)
  X := Eltseq(X);
  Fqx := Universe(X);
  H := ChangeRing(H, Fqx);
  s := ChangeRing(s, Fqx);

  return Eltseq(Vector(Fqx, X) * Transpose(H) - s);
end function;

function FieldEquations(X)
  return [x^2 - x : x in Eltseq(X)];
end function;

// The "exact weight comparison" equations are {y_{n,j} + t_j : 1 <= j <= ell_t}
// where t_j is bit j-1 of t (1 <= j <= ell_t), and y_{n,j} = 0 for j > ell_t.
function ExactWeightEquations(Y, t)
    Yn := Y[#Y];
    // this is ell or ell_t depending on Y
    ell := #Yn;
    error if t ge 2^ell, "t does not fit in the weight counter";

    tbits := Intseq(t, 2, ell);
    return [Yn[j] + tbits[j] : j in [1..ell]];
end function;

////////// CCMMP2 modeling //////////

// Write Y_i = (y_{i,1}, ..., y_{i, \ell_t}) for 1 <= i <= n.
// The CCMMP2 Hamming weight computation equations are
// {x_1 + y_{1, 1}, y_{1, 2}, ..., y_{1, ell_t}}
// cup {x_i + y_{i,1} + y_{i-1, 1} : 2 <= i <= n}
// cup {y_{i,j-1}y_{i-1,j-1} + y_{i,j} + y_{i-1,j-1} + y_{i-1,j} : 2 <= i <= n, 2 <= j <= ell_t}
// cup {y_{i,ell_t}y_{i-1, ell_t} + y_{i-1,ell_t} : 2 <= i <= n}.
//
// X = [x_1, ..., x_n], Y = [Y_1, ..., Y_n] with Y[i][j] = y_{i,j}.
function CCMMP2_WeightComputationEquations(X, Y)
    n := #X;
    ell := #Y[1];

    eqs := [X[1] + Y[1][1]] cat [Y[1][j] : j in [2..ell]];
    for i in [2..n] do
        Append(~eqs, X[i] + Y[i][1] + Y[i-1][1]);
        eqs cat:= [Y[i][j-1]*Y[i-1][j-1] + Y[i][j] + Y[i-1][j-1] + Y[i-1][j]
                   : j in [2..ell]];
        Append(~eqs, Y[i][ell]*Y[i-1][ell] + Y[i-1][ell]);
    end for;
    return eqs;
end function;

function CCMMP2(sd_inst)
    n := sd_inst`n;
    H := sd_inst`H;
    s := sd_inst`s;
    t := sd_inst`t;
    ell_t := Floor(Log(2, t)) + 1;
    nvars := n*(ell_t + 1);

    // Variables: x_1, ..., x_n, then y_{1,1}, ..., y_{1,ell_t}, y_{2,1}, ...
    R := PolynomialRing(GF(2), nvars, "grevlex");
    AssignNames(~R, [Sprintf("x%o", i) : i in [1..n]] cat
        &cat[[Sprintf("y%o_%o", i, j) : j in [1..ell_t]] : i in [1..n]]);
    X := [R.i : i in [1..n]];
    Y := [[R.(n + (i - 1)*ell_t + j) : j in [1..ell_t]] : i in [1..n]];

    eqs := ParityCheckEquations(X, H, s)
       cat FieldEquations([R.i : i in [1..nvars]])
       cat ExactWeightEquations(Y, t)
       cat CCMMP2_WeightComputationEquations(X, Y);

    return eqs;
end function;

////////// RST modeling //////////

// Write N = n - t + 1 and Z = [Z_1, ..., Z_t], Z_i = (z_{i,1}, ..., z_{i,N}),
// with Z[i][j] = z_{i,j}. z_{i,j} = 1 iff the i-th support position is
// p_i = i + j - 1 (i.e. its shifted position is a_i = j).

// The error vector as linear forms in Z:
// x_p = sum_{1 <= i <= t, 1 <= p-i+1 <= N} z_{i,p-i+1}, 1 <= p <= n.
function RST_ErrorVector(Z, n)
    t := #Z;
    N := #Z[1];
    return [&+[Universe(Z[1]) | Z[i][p - i + 1] : i in [1..t] | 1 le p - i + 1 and p - i + 1 le N]
            : p in [1..n]];
end function;

// The ordering equations: {z_{i,j} z_{l,m} : 1 <= i <= l <= t, 1 <= m < j <= N}.
function RST_OrderingEquations(Z)
    t := #Z;
    N := #Z[1];

    eqs := [];
    for i in [1..t] do
        for l in [i..t] do
            for j in [2..N] do
                eqs cat:= [Z[i][j]*Z[l][m] : m in [1..j-1]];
            end for;
        end for;
    end for;
    return eqs;
end function;

// The weight equations: {sum_{j=1}^{N} z_{i,j} - 1 : 1 <= i <= t}
function RST_WeightEquations(Z)
    return [&+Z[i] - 1 : i in [1..#Z]];
end function;

// Returns the equations and the error vector X = [x_1, ..., x_n] as linear
// forms in Z (there are no x variables), for
//   IdealSolutions(inst, Ideal(eqs) : X := X).
function RST(sd_inst)
    n := sd_inst`n;
    H := sd_inst`H;
    s := sd_inst`s;
    t := sd_inst`t;
    error if t lt 1, "RST needs t >= 1";
    N := n - t + 1;
    nvars := t*N;

    // Variables: z_{1,1}, ..., z_{1,N}, z_{2,1}, ..., z_{t,N}
    R := PolynomialRing(GF(2), nvars, "grevlex");
    AssignNames(~R, &cat[[Sprintf("z%o_%o", i, j) : j in [1..N]] : i in [1..t]]);
    Z := [[R.((i - 1)*N + j) : j in [1..N]] : i in [1..t]];
    X := RST_ErrorVector(Z, n);

    // No field equations: z^2 - z is implied by the ordering and weight equations.
    eqs := ParityCheckEquations(X, H, s)
       cat RST_OrderingEquations(Z)
       cat RST_WeightEquations(Z);

    return eqs, X;
end function;

////////// LMT modelings //////////

// The C-Hamming and QFC-Hamming modelings of La Scala, Marchesin and Tiwari
// compute elementary symmetric polynomials of x along the balanced binary
// decomposition tree T_n of [1, n]: every interval I = [a, b] with a < b has
// the children I' = [a, m] and I'' = [m + 1, b], m = floor((a + b)/2), and
// the leaves are the intervals [a, a]. Both use
//   L = ceil(log2(max(t, n - t) + 1)) - 1 = floor(log2(max(t, n - t))),
// since the vectors of weight t are the common zeros of e_{2^h}(x) + bit_h(t),
// 0 <= h <= L, and the field equations. As in the variable counts of La Scala
// et al., the variable of a leaf [a, a] is x_a itself, and the root keeps only
// the variables y_{2^h}, 0 <= h <= L, used by the weight equations.

// The nodes of the subtree of T_n rooted at [a, b], children before parents.
function LMT_Subtree(a, b)
    if a eq b then
        return [[a, b]];
    end if;
    m := (a + b) div 2;
    return LMT_Subtree(a, m) cat LMT_Subtree(m + 1, b) cat [[a, b]];
end function;

function LMT_Children(I)
    m := (I[1] + I[2]) div 2;
    return [I[1], m], [m + 1, I[2]];
end function;

function LMT_L(n, t)
    return Ilog2(Max([t, n - t, 1]));
end function;

// C-Hamming. Node I has the variables y^I_d, meant to be e_d(x_I), for
// 1 <= d <= min(|I|, 2^L), except that the root has only d = 2^h; and
// y^I_0 = 1, y^I_d = 0 for d > |I|. The node and root polynomials are
// {y^I_d + sum_{j=0}^{d} y^{I'}_j y^{I''}_{d-j} : I not a leaf}
// cup {y^{[1,n]}_{2^h} + bit_h(t) : 0 <= h <= L}.
function LMTC(sd_inst)
    n := sd_inst`n;
    H := sd_inst`H;
    s := sd_inst`s;
    t := sd_inst`t;
    L := LMT_L(n, t);
    root := [1, n];
    internal := [I : I in LMT_Subtree(1, n) | I[1] lt I[2]];
    degrees := func<I | I eq root select [2^h : h in [0..L]]
                                     else [1..Min(I[2] - I[1] + 1, 2^L)]>;

    // Variables: x_1, ..., x_n, then the y^I_d of each internal node I,
    // children before parents
    names := [Sprintf("x%o", i) : i in [1..n]];
    for I in internal do
        names cat:= [Sprintf("y%o_%o_%o", I[1], I[2], d) : d in degrees(I)];
    end for;
    R := PolynomialRing(GF(2), #names, "grevlex");
    AssignNames(~R, names);
    X := [R.i : i in [1..n]];

    // Y[I][d] = y^I_d
    Y := AssociativeArray();
    for a in [1..n] do
        Y[[a, a]] := AssociativeArray();
        Y[[a, a]][1] := X[a];
    end for;
    v := n;
    for I in internal do
        Y[I] := AssociativeArray();
        for d in degrees(I) do
            v +:= 1;
            Y[I][d] := R.v;
        end for;
    end for;
    y := func<I, d | d eq 0 select R!1 else (d gt I[2] - I[1] + 1 select R!0 else Y[I][d])>;

    tbits := Intseq(t, 2, L + 1);
    eqs := ParityCheckEquations(X, H, s)
       cat FieldEquations([R.i : i in [1..Rank(R)]])
       cat [y(root, 2^h) + tbits[h + 1] : h in [0..L]];
    for I in internal do
        I1, I2 := LMT_Children(I);
        eqs cat:= [y(I, d) + &+[y(I1, j)*y(I2, d - j) : j in [0..d]] : d in degrees(I)];
    end for;

    return eqs;
end function;

// QFC-Hamming. Node I has the variables y^I_{2^g}, meant to be bit g of
// wt(x_I), for 0 <= g < rho(I) := min(floor(log2 |I|), L) + 1. For j >= 0
// write c^I_j for the product of the y^I_{2^g} over the bits g of j, so
// c^I_0 = 1, and c^I_j = 0 if j >= 2^rho(I). The node polynomials (the
// factorized convolutions) are
// {y^I_{2^h} + sum_{j=0}^{2^h} c^{I'}_j c^{I''}_{2^h-j} : I not a leaf, 0 <= h < rho(I)}.
// Every c^J_j with at least two factors that occurs in them, i.e. with
// j <= 2^(rho(I)-1) for the parent I of J, is a new variable, with the product
// polynomial c^J_j + c^J_{j - 2^g} y^J_{2^g}, g the highest bit of j. The root
// polynomials are as for C-Hamming.
function LMTQFC(sd_inst)
    n := sd_inst`n;
    H := sd_inst`H;
    s := sd_inst`s;
    t := sd_inst`t;
    L := LMT_L(n, t);
    root := [1, n];
    nodes := LMT_Subtree(1, n);
    internal := [I : I in nodes | I[1] lt I[2]];
    rho := func<I | Min(Ilog2(I[2] - I[1] + 1), L) + 1>;
    weight := func<j | &+Intseq(j, 2)>;

    // products[J]: the j of the product variables c^J_j
    products := AssociativeArray();
    products[root] := [Integers() | ];
    for I in internal do
        I1, I2 := LMT_Children(I);
        for J in [I1, I2] do
            products[J] := [j : j in [3..Min(2^(rho(I) - 1), 2^rho(J) - 1)] | weight(j) ge 2];
        end for;
    end for;

    // Variables: x_1, ..., x_n, then for each internal node, children before
    // parents, its y^I_{2^g} and then its c^I_j
    names := [Sprintf("x%o", i) : i in [1..n]];
    for I in internal do
        names cat:= [Sprintf("y%o_%o_%o", I[1], I[2], 2^g) : g in [0..rho(I) - 1]];
        names cat:= [Sprintf("c%o_%o_%o", I[1], I[2], j) : j in products[I]];
    end for;
    R := PolynomialRing(GF(2), #names, "grevlex");
    AssignNames(~R, names);
    X := [R.i : i in [1..n]];

    // C[I][j] = c^I_j for j a power of two or a product variable.
    C := AssociativeArray();
    for a in [1..n] do
        C[[a, a]] := AssociativeArray();
        C[[a, a]][1] := X[a];
    end for;
    v := n;
    for I in internal do
        C[I] := AssociativeArray();
        for j in [2^g : g in [0..rho(I) - 1]] cat products[I] do
            v +:= 1;
            C[I][j] := R.v;
        end for;
    end for;
    c := func<I, j | j eq 0 select R!1 else (j ge 2^rho(I) select R!0 else C[I][j])>;

    tbits := Intseq(t, 2, L + 1);
    eqs := ParityCheckEquations(X, H, s)
       cat FieldEquations([R.i : i in [1..Rank(R)]])
       cat [c(root, 2^h) + tbits[h + 1] : h in [0..L]];
    for I in internal do
        I1, I2 := LMT_Children(I);
        eqs cat:= [c(I, 2^h) + &+[c(I1, j)*c(I2, 2^h - j) : j in [0..2^h]]
                   : h in [0..rho(I) - 1]];
    end for;
    for J in nodes do
        eqs cat:= [c(J, j) + c(J, j - 2^g)*c(J, 2^g) where g := Ilog2(j)
                   : j in products[J]];
    end for;

    return eqs;
end function;

////////// RCI modeling //////////

// Write Y_i = (y_{i,1}, ..., y_{i,ell_t}) for 1 <= i <= n, and let c_{i,j} be
// the carry entering bit j of the register while adding x_i, with
// c_{i,1} := x_i. The RCI (ripple-carry incrementer) Hamming weight
// computation equations are
// {x_1 + y_{1,1}, y_{1,2}, ..., y_{1,ell_t}}
// cup {y_{i,j} + y_{i-1,j} + c_{i,j} : 2 <= i <= n, 1 <= j <= ell_t}
// cup {c_{i,j+1} + y_{i-1,j} c_{i,j} : 2 <= i <= n, 1 <= j <= ell_t - 1}
// cup {y_{i-1,ell_t} c_{i,ell_t} : 2 <= i <= n},
// where the last set is the no-overflow condition.
//
// X = [x_1, ..., x_n], Y = [Y_1, ..., Y_n] with Y[i][j] = y_{i,j}, and
// C = [C_1, ..., C_n] with C[i][j] = c_{i,j} for i >= 2 (so C[i][1] = x_i).
// C[1] is not used.
function RCI_WeightComputationEquations(X, Y, C)
    n := #X;
    ell := #Y[1];

    eqs := [X[1] + Y[1][1]] cat [Y[1][j] : j in [2..ell]];
    for i in [2..n] do
        eqs cat:= [Y[i][j] + Y[i-1][j] + C[i][j] : j in [1..ell]];
        eqs cat:= [C[i][j+1] + Y[i-1][j]*C[i][j] : j in [1..ell-1]];
        Append(~eqs, Y[i-1][ell]*C[i][ell]);
    end for;
    return eqs;
end function;

function RCI(sd_inst)
    n := sd_inst`n;
    H := sd_inst`H;
    s := sd_inst`s;
    t := sd_inst`t;
    ell_t := Floor(Log(2, t)) + 1;
    nxy := n*(ell_t + 1);
    nvars := nxy + (n - 1)*(ell_t - 1);

    // Variables: x_1, ..., x_n, then y_{1,1}, ..., y_{1,ell_t}, y_{2,1}, ...,
    // then c_{2,2}, ..., c_{2,ell_t}, c_{3,2}, ..., c_{n,ell_t}
    R := PolynomialRing(GF(2), nvars, "grevlex");
    names := [Sprintf("x%o", i) : i in [1..n]] cat
        &cat[[Sprintf("y%o_%o", i, j) : j in [1..ell_t]] : i in [1..n]];
    for i in [2..n] do
        names cat:= [Sprintf("c%o_%o", i, j) : j in [2..ell_t]];
    end for;
    AssignNames(~R, names);
    X := [R.i : i in [1..n]];
    Y := [[R.(n + (i - 1)*ell_t + j) : j in [1..ell_t]] : i in [1..n]];
    C := [[X[1]]] cat
         [[X[i]] cat [R.(nxy + (i - 2)*(ell_t - 1) + j - 1) : j in [2..ell_t]]
          : i in [2..n]];

    eqs := ParityCheckEquations(X, H, s)
       cat FieldEquations([R.i : i in [1..nvars]])
       cat ExactWeightEquations(Y, t)
       cat RCI_WeightComputationEquations(X, Y, C);

    return eqs;
end function;

////////// CSA modeling //////////

// The carry-save adder circuit computing wt(x) for x in F_2^n, with
// ell = floor(log2(n)) + 1 columns, where column j holds bits of weight
// 2^(j-1). Wires are numbered: 1, ..., n are the inputs x_1, ..., x_n, and
// gate g outputs the sum wire n + 2g - 1 and the carry wire n + 2g.
// For j = 1, ..., ell in order, column j is compressed first-in-first-out:
//   1. while column j has at least three wires, remove the three oldest
//      a, b, c and apply a full adder: its sum is appended to column j and
//      its carry to column j + 1;
//   2. if exactly two wires a, b remain, apply a half-adder.
// The single wire w_j left in column j is bit j-1 of wt(x).
//
// Returns the gates as sequences [a, b, c] of input wires (c = 0 for a
// half-adder) and the output wires [w_1, ..., w_ell].
function CSA_Circuit(n)
    ell := Floor(Log(2, n)) + 1;

    cols := [[1..n]] cat [[Integers() | ] : j in [1..ell]];
    gates := [];
    w := [];
    nwires := n;
    for j in [1..ell] do
        col := cols[j];
        while #col ge 2 do
            if #col ge 3 then
                abc := col[1..3];
                col := col[4..#col];
            else
                abc := col cat [0];
                col := [];
            end if;
            Append(~gates, abc);
            Append(~col, nwires + 1);
            Append(~cols[j + 1], nwires + 2);
            nwires +:= 2;
        end while;
        Append(~w, col[1]);
    end for;
    // Column ell receives exactly one wire, so no carry leaves it.
    assert IsEmpty(cols[ell + 1]);

    return gates, w;
end function;

// The CSA Hamming weight computation equations: for a gate with inputs
// a, b, c, sum sigma and carry kappa,
// {sigma + a + b + c, kappa + maj(a, b, c)},  maj(a, b, c) = ab + ac + bc,
// where c = 0 for a half-adder (maj(a, b, 0) = ab).
//
// V = [v_1, ..., v_N] are the variables indexed by wire number, and gates is
// as returned by CSA_Circuit.
function CSA_WeightComputationEquations(V, gates)
    n := #V - 2*#gates;
    zero := Universe(V)!0;

    eqs := [];
    for g in [1..#gates] do
        a := V[gates[g][1]];
        b := V[gates[g][2]];
        c := gates[g][3] eq 0 select zero else V[gates[g][3]];
        sigma := V[n + 2*g - 1];
        kappa := V[n + 2*g];
        eqs cat:= [sigma + a + b + c, kappa + a*b + a*c + b*c];
    end for;
    return eqs;
end function;

function CSA(sd_inst)
    n := sd_inst`n;
    H := sd_inst`H;
    s := sd_inst`s;
    t := sd_inst`t;
    gates, w := CSA_Circuit(n);
    nvars := n + 2*#gates;

    // Variables: x_1, ..., x_n, then s_1, k_1, s_2, k_2, ..., the sum and
    // carry of each gate in the order the gates are applied
    R := PolynomialRing(GF(2), nvars, "grevlex");
    names := [Sprintf("x%o", i) : i in [1..n]];
    for g in [1..#gates] do
        names cat:= [Sprintf("s%o", g), Sprintf("k%o", g)];
    end for;
    AssignNames(~R, names);
    V := [R.i : i in [1..nvars]];
    X := V[1..n];
    // The output wires form the final register: W[j] = w_j is bit j-1 of wt(x).
    W := [V[i] : i in w];

    eqs := ParityCheckEquations(X, H, s)
       cat FieldEquations(V)
       cat ExactWeightEquations([W], t)
       cat CSA_WeightComputationEquations(V, gates);

    return eqs;
end function;
