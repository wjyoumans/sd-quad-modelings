// Planted syndrome decoding instances: a parity-check matrix H in
// F_2^{(n-k) x n}, a syndrome s, a weight t, and the planted solution, a
// vector of weight t with H x^T = s^T.
RF_SDInstance := recformat<n, k, t, seed, H, s, solution>;

// The instance in a file written by generate.py, with the sections "# n",
// "# k", "# seed", "# w" (= t), "# H (each line is a row of H)",
// "# s^transpose" and "# solution".
function ReadInstance(file)
    // Section name -> its non-header lines, e.g. "# H (each ...)" is section "H".
    sections := AssociativeArray();
    key := "";
    for l in Split(Read(file), "\n") do
        if #l eq 0 then continue; end if;
        if l[1] eq "#" then
            key := Split(l, " ")[2];
            sections[key] := [];
        else
            error if key eq "", "data before the first section in " cat file;
            sections[key] cat:= [l];
        end if;
    end for;

    n := StringToInteger(sections["n"][1]);
    k := StringToInteger(sections["k"][1]);
    m := n - k;

    F := GF(2);
    bits := func< l | [F | StringToInteger(l[i]) : i in [1..#l]] >;
    rows := sections["H"];
    error if #rows ne m or exists{r : r in rows | #r ne n}, "unexpected H in " cat file;

    inst := rec<RF_SDInstance | >;
    inst`n := n;
    inst`k := k;
    inst`t := StringToInteger(sections["w"][1]);
    inst`seed := StringToInteger(sections["seed"][1]);
    inst`H := Matrix(F, m, n, &cat[bits(r) : r in rows]);
    inst`s := Vector(F, bits(sections["s^transpose"][1]));
    inst`solution := Vector(F, bits(sections["solution"][1]));
    return inst;
end function;

// Whether x solves the instance: H * x^T = s^T and wt(x) = t.
function IsSDSolution(sd_inst, x)
    return x * Transpose(sd_inst`H) eq sd_inst`s and Weight(x) eq sd_inst`t;
end function;

// The solutions of the ideal I of a modeling of sd_inst as error vectors, and
// whether every one of them solves sd_inst. The error vector is read from X, a
// sequence of n polynomials giving x_1, ..., x_n; by default X is the first n
// variables of the ring (pass the X returned by RST). Computes the Groebner
// basis of I first unless it is already known.
function IdealSolutions(sd_inst, I : X := [])
    n := sd_inst`n;
    F := BaseRing(sd_inst`H);

    R := Generic(I);
    if #X eq 0 then
        X := [R.i : i in [1..n]];
    end if;

    // Variety reuses the basis stored on I. [1] (no solutions) gives [].
    sols := [Vector(F, [Evaluate(x, [p[i] : i in [1..Rank(R)]]) : x in X])
             : p in Variety(I)];
    ok := &and[Booleans() | IsSDSolution(sd_inst, x) : x in sols];

    return sols, ok;
end function;
