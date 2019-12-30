using CSV, JuMP, Ipopt #, DataFrames
"An open source Computerised General Equlibrium model"
# ```jldoctest
# a = 1
# b = 2
# c = a + b
# # output
# 3
# ```

CGE_EL = Model(with_optimizer(Ipopt.Optimizer))
# Input Output Table
IOdata = CSV.read("IOdata.csv") #Careful Eli :: columns & rows transposed

sectors = [1,2]
commods = [1,2]

# Initial Values
rKi = 1                     #initial return to Kapital
wLi = 1                     #initial return to Labor (initial wage)
Kdi = IOdata[:Kdi][sectors]          #initial Kapital demand
Ldi = IOdata[:Ldi][sectors]          #initial Labour demand
Cdi = convert(Matrix,IOdata[[3],sectors.+1])       #initial Consumption demand
# This is a complicated way to get the Cdi as an array from IOdata, just so I can keep the data in the
## same structure as EcoMod
Kei = sum(Kdi)              #initial Kapital endowment
Lei = sum(Ldi)              #initial Labour endowment
IOi = [IOdata[sectors,:Sec1] IOdata[sectors,:Sec2]] # Input Output table (as array)
#?? A bit tricky here, unclear what dimensions this shld have; I THINK the Totals are right
YiOut =  sum(IOi, dims=2) + Kdi + Ldi  #initial gross Total income (by sector)
YiIn = rKi * Kei + wLi * Lei    #initial Income level
TechCf = IOi./YiOut

#this is probably not efficient...
CPi = []
for i in 1:length(sectors)
    push!(CPi, 1)
end                    #initial Commodity Price Level (1 for each sector commod)

HHUexp = CPi[sectors] .* Cdi / YiIn
HHUi = prod(Cdi.^HHUexp)        #initial HouseHold Utility level
Le = Lei * 1.1     #Labour (endowment) grows at 1.1
Ke = Kei * 1.2    #Kapital (endowment) grows at 1.2
w = wLi * 5

## Parameters per sector
Kexp = rKi * Kdi[sectors] ./ (rKi * Kdi[sectors] + wLi * Ldi[sectors])  #exponent on K in Cobb-Douglas Production Function
Lexp = 1 .- Kexp                     ##exponent on L in Cobb-Douglas Production Function
PFeFs = YiOut ./ (Kdi.^Kexp .* Ldi.^Lexp)  #Production Function efficiency of Firms ???shld this be dot product or matrix product???

#Initial (equilibrium) levels for the endogenous variables and lower bounds to preven numerical problems in opt
@variable(CGE_EL, r, lower_bound = 0.001 * rKi, start = rKi)#, lower_bound = 0)   
@NLconstraint(CGE_EL, r >=0) # This does NOT have the same behaviour as incorporating in variable ...
@variable(CGE_EL, w, lower_bound = 0.001 * wLi, start = wLi)#, lower_bound = 0)  
@NLconstraint(CGE_EL, w >=0) # This does NOT have the same behaviour as incorporating in variable ...
@variable(CGE_EL, Kd_f[i = sectors], start = Kdi[i], lower_bound=0.001*Kdi[i], lower_bound=0)
#@NLconstraint(CGE_EL, Kd_f >=0) # Is that K? # This doesn't run
@variable(CGE_EL, Ld_f[i = sectors], start = Ldi[i], lower_bound=0.001*Ldi[i], lower_bound=0)
#@NLconstraint(CGE_EL, Ld_f >=0) # This doesn't run
@variable(CGE_EL, HHi, lower_bound = 0.001 * YiIn, start = YiIn, lower_bound=0)
#@NLconstraint(CGE_EL, HHi >=0)
@variable(CGE_EL, Pr_Commods[i = sectors], start = CPi[i], lower_bound=0.001*CPi[i], lower_bound=0) 
#@NLconstraint(CGE_EL, Pr_Commods[i] >=0) #This doesn't run, of course
@variable(CGE_EL, YOut[i = sectors], start = YiOut[i], lower_bound=0.001*YiOut[i], lower_bound=0) 
#@NLconstraint(CGE_EL, YOut >=0) # This doesn't run anyway
@variable(CGE_EL, Commodsd_HH[i = sectors], start = Cdi[i], lower_bound=0.001*Cdi[i], lower_bound=0) 
#@NLconstraint(CGE_EL, Commodsd_HH >=0) # This doesn't run anyway

#@NLconstraint(CGE_EL, Ke >= 0)
#@NLconstraint(CGE_EL, Le >= 0)

@NLconstraint(CGE_EL, ECd[i = sectors], Pr_Commods[i] * Commodsd_HH[i] == HHUexp[i] * HHi) #Consumer commodity demand function (income =spending)
@NLconstraint(CGE_EL, EKd_f[i = sectors], Kd_f[i] == (YOut[i] / PFeFs[i]) * (Kexp[i] * w / (Lexp[i] * r))^Lexp[i]) #Kapital demand function (Kap = ?...)
@NLconstraint(CGE_EL, ELd_f[i = sectors], Ld_f[i] == (YOut[i] / PFeFs[i]) * (Lexp[i] * r / (Kexp[i] * w))^Kexp[i]) #Labour demand function (Kap = ?...)
@NLconstraint(CGE_EL, ENoProf[i = sectors], Pr_Commods[i] * YOut[i] == r * Kd_f[i] + w * Ld_f[i] + 
sum(TechCf[i,j]*YOut[j]*Pr_Commods[i] for j in sectors)) #Competitive Equlibrium, no profit (value of output=production)
# Might need  -- sum(TechCf,dims=2) -- I'm not sure how the JuMP indexing will play out.
######## I'm up to here (2019-12-29)
@NLconstraint(CGE_EL, ELs_f, sum(Ld_f[i] for i in sectors) == Le) # Competitve Eq: L demand = L supply
@NLconstraint(CGE_EL, EKs_f, sum(Kd_f[i] for i in sectors) == Ke) # Competitve Eq: K demand = K supply
@NLconstraint(CGE_EL, EC[i = sectors], Commodsd_HH[i] == YOut[i]) # 
@NLconstraint(CGE_EL, EHHi, HHi == r * Ke + w * Le) #Income Definition
#@NLconstraint(CGE_EL, ETrick, Trick == 1) # This just makes everything stay at initial values

Trick = 1
@NLobjective(CGE_EL, Max, Trick) #Numeraire from EcoMod...not clear exactly how to translate
@time optimize!(CGE_EL)

#Look at some results
CGE_EL
# This should calclule Utility, but doesn't run :: print("HHU=           ", prod(JuMP.value.(Commodsd_HH)^HHUexp))
print("w=             ",JuMP.value(w),"\n")
print("r=             ",JuMP.value(r),"\n")
print("HHi=           ",JuMP.value(HHi),"\n")
print("Ld_f[1]=       ",JuMP.value(Ld_f[1]),"\n")
print("Ld_f[2]=       ",JuMP.value(Ld_f[2]),"\n")
print("Kd_f[1]=       ",JuMP.value(Kd_f[1]),"\n")
print("Kd_f[2]=       ",JuMP.value(Kd_f[2]),"\n")
print("Pr_Commods[1]= ",JuMP.value(Pr_Commods[1]),"\n")
print("Pr_Commods[2]= ",JuMP.value(Pr_Commods[2]),"\n")
print("Commodsd_HH[1]=",JuMP.value(Commodsd_HH[1]),"\n")
print("Commodsd_HH[2]=",JuMP.value(Commodsd_HH[2]),"\n")
print("YOut[1]=       ",JuMP.value(YOut[1]),"\n")
print("YOut[2]=       ",JuMP.value(YOut[2]),"\n")
