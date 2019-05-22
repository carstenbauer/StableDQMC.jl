# Takes two .bench files as input
length(ARGS) >= 2 || error("Please provide two .bench files as input args.")

using Pkg
Pkg.activate(@__DIR__)

using PkgBenchmark

f1 = ARGS[1]
f2 = ARGS[2]

r1 = readresults(f1)
r2 = readresults(f2)

j = judge(r1, r2)

name1 = replace(f1, ".bench" => "")
name2 = replace(f2, ".bench" => "")

export_markdown(string("judge_",name1,"_",name2,".md"), j)