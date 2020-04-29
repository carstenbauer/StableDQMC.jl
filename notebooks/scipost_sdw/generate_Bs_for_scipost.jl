# load DQMC framework
include(joinpath(ENV["JULIA_DQMC"], "src/dqmc_framework.jl"))

using HDF5
c = h5read("conf_near_qcp_L_10_B_40.h5", "conf")

p = Params()
p.output_file = "live.out.h5.running"
xml2parameters!(p, "dqmc.in.xml")
mc = DQMC(p)

h5open("scipost_sdw.h5", isfile("scipost_sdw.h5") ? "r+" : "w") do file
    for β in (40)#range(5,40,step=5)
        mc.p.beta = β
        mc.p.slices = Int(mc.p.beta/mc.p.delta_tau)
        init!(mc)
        p.hsfield .= c
        for l in 1:mc.p.slices
            file["B_$β/$l"] = slice_matrix(mc,l)
        end
    end
end
