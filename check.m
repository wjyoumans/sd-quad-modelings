// Checks the six modelings against brute force. Run from the repository root:
//   magma -b check.m
// For n = 16, 18, 20, k = n/2 and n/4, and seeds 0 and 1, checks that
//   1. the instance is as generate.py describes it: H has full rank and the
//      planted solution has weight t_GV(n, k);
//   2. every modeling's solutions are exactly the solutions found by brute
//      force, each once.

SetQuitOnError(true);
SetColumns(0);
load "instance.m";
load "modelings.m";

function TGV(n, k)
    return Max([t : t in [0..n div 2] | Binomial(n, t) le 2^(n - k)]);
end function;

function BruteForceSolutions(inst)
    n := inst`n;
    sols := {};
    for S in Subsets({1..n}, inst`t) do
        x := Vector(GF(2), [i in S select 1 else 0 : i in [1..n]]);
        if IsSDSolution(inst, x) then
            Include(~sols, x);
        end if;
    end for;
    return sols;
end function;

ncases := 0;
for n in [16, 18, 20], k in [n div 2, n div 4], seed in [0, 1] do
    file := Sprintf("instances/SD_%o_%o_%o", n, k, seed);
    inst := ReadInstance(file);
    error if inst`k ne k or Rank(inst`H) ne n - k or inst`t ne TGV(n, k)
        or not IsSDSolution(inst, inst`solution), "unexpected instance", file;
    expected := BruteForceSolutions(inst);
    for model in ["CCMMP2", "RST", "RCI", "CSA", "CHAM", "QFCHAM"] do
        X := [];
        case model:
            when "CCMMP2": eqs := CCMMP2(inst);
            when "RST":    eqs, X := RST(inst);
            when "RCI":    eqs := RCI(inst);
            when "CSA":    eqs := CSA(inst);
            when "CHAM":   eqs := LMTC(inst);
            when "QFCHAM": eqs := LMTQFC(inst);
        end case;
        sols, ok := IdealSolutions(inst, Ideal(eqs) : X := X);
        error if not ok or #sols ne #Set(sols) or Set(sols) ne expected,
            "wrong solutions", model, file;
        ncases +:= 1;
    end for;
end for;
printf "ok (%o systems)\n", ncases;
quit;
