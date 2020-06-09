"Basic Monte Carlo function, parameter distribution sampling in Main script (so far)"
n = 8 # how many runs
using NamedArrays
MC_Results = DataFrame()

function mc()
    for i = 1:n
        include("CGE-J-EL_Main.jl")
        MC_Results[Symbol("Unempl$Unempli")] =  [JuMP.value(w),JuMP.value(r),JuMP.value(HHI),JuMP.value(Ld_f[1]),JuMP.value(Ld_f[2]),
        JuMP.value(Ld_f[3]),JuMP.value(Ld_f[4]),JuMP.value(Ld_f[5]),JuMP.value(Ld_f[6]),JuMP.value(Kd_f[1]),JuMP.value(Kd_f[2]),JuMP.value(Kd_f[3]),
        JuMP.value(Kd_f[4]),JuMP.value(Kd_f[5]),JuMP.value(Kd_f[6]) ,JuMP.value(Pr_Commods[1]) ,JuMP.value(Pr_Commods[2]),JuMP.value(Commodsd_HH[1]),
        JuMP.value(Commodsd_HH[2]),JuMP.value(YOut[1]),JuMP.value(YOut[2]), HHU,JuMP.value(TotSav),Walras,GovBudg]
        println(i)
    end
end

mc()

MC_Results = NamedArray(convert(Matrix,MC_Results),(["w","r","HHI","Ld_f[1]","Ld_f[2]","Ld_f[3]","Ld_f[4]","Ld_f[5]","Ld_f[6]","Kd_f[1]","Kd_f[2]",
"Kd_f[3]","Kd_f[4]","Kd_f[5]","Kd_f[6]","Pr_Commods[1]","Pr_Commods[2]","Commodsd_HH[1]","Commodsd_HH[2]","YOut[1]","YOut[2]","HHU","TotSav","Walras",
"GovBudg",],names(MC_Results)))