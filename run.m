// Times the Groebner basis computation of one modeling on one instance.
// Run from the repository root:
//   magma -b inst:=instances/SD_32_16_0 model:=CSA mem:=40 run.m
// where model is one of CCMMP2, RST, RCI, CSA, CHAM (C-Hamming) and QFCHAM
// (QFC-Hamming), and mem is the memory limit in GiB (optional).
// It prints, in order,
//   JOB ...     the instance, the modeling and the largest degree of its
//               polynomials,
//   the verbose output of Magma's F4, with the degree of each step and the
//   size of each matrix,
//   RESULT ...  the CPU and real time of GroebnerBasis in seconds, the number
//               of solutions, and whether all of them solve the instance and
//               the planted solution is among them.

SetQuitOnError(true);
SetColumns(0);
// The F4 matrices can depend on Magma's random seed, so fix it to make runs
// repeatable.
SetSeed(1);
load "instance.m";
load "modelings.m";

if assigned mem then
    SetMemoryLimit(StringToInteger(mem) * 2^30);
end if;

sd := ReadInstance(inst);
X := [];
case model:
    when "CCMMP2": eqs := CCMMP2(sd);
    when "RST":    eqs, X := RST(sd);
    when "RCI":    eqs := RCI(sd);
    when "CSA":    eqs := CSA(sd);
    when "CHAM":   eqs := LMTC(sd);
    when "QFCHAM": eqs := LMTQFC(sd);
    else error "unknown model", model;
end case;
printf "JOB n=%o k=%o t=%o seed=%o model=%o maxdeg=%o\n",
    sd`n, sd`k, sd`t, sd`seed, model, Max([TotalDegree(f) : f in eqs]);

I := Ideal(eqs);
SetVerbose("Groebner", 1);
t0 := Cputime();
r0 := Realtime();
G := GroebnerBasis(I);
gb_cpu := Cputime(t0);
gb_real := Realtime(r0);
SetVerbose("Groebner", 0);

sols, ok := IdealSolutions(sd, I : X := X);
ok and:= sd`solution in sols;
printf "RESULT gb_cpu=%o gb_real=%o nsols=%o ok=%o\n", gb_cpu, gb_real, #sols, ok;
quit;
