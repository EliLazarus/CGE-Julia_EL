using CSV, JuMP, Ipopt
"An open source Computerised General Equlibrium model"
# ```jldoctest
# a = 1
# b = 2
# c = a + b
# # output
# 3
# ```

US_CGE = Model(with_optimizer(Ipopt.Optimizer))
# Input Output Table
IOdata = CSV.read("IOdata.csv") #Careful Eli :: columns & rows transposed

sectors = [1,2]
commods = [1,2]

# Initial Values
rKi = 1                     #initial return to Kapital
wLi = 1                     #initial return to Labor (initial wage)
Kdi = IOdata[:Kdi]          #initial Kapital demand
Ldi = IOdata[:Ldi]          #initial Labour demand
Cdi = IOdata[:Cdi]          #initial Consumption demand
Kei = sum(Kdi)              #initial Kapital endowment
Lei = sum(Ldi)              #initial Labour endowment
YiOut =  rKi * Kdi + wLi * Ldi  #initial gross Total income (by sector)
YiIn = rKi * Kei + wLi * Lei    #initial Income level

CPi = []
for i in 1:length(IOdata[1])
    push!(CPi, 1)
end                    #initial Commodity Price Level (1 for each sector commod)

HHUexp = CPi .* Cdi / YiIn
HHUi = prod(Cdi.^HHUexp)        #initial HouseHold Utility level
Le = Lei * 1.1     #Labour (endowment) grows at 1.1
Ke = Kei * 1.2    #Kapital (endowment) grows at 1.2
w = wLi * 5

## Parameters per sector
Kexp = rKi * Kdi ./ (rKi * Kdi + wLi * Ldi)  #exponent on K in Cobb-Douglas Production Function
Lexp = 1 .- Kexp                     ##exponent on L in Cobb-Douglas Production Function
PFeFs = YiOut ./ (Kdi.^Kexp .* Ldi.^Lexp)  #Production Function efficiency of Firms ???shld this be dot product or matrix product???

#Initial (equilibrium) levels for the endogenous variables and lower bounds to preven numerical problems in opt
@variable(US_CGE, r, lower_bound = 0.001 * rKi, start = rKi)#, lower_bound = 0)   
@NLconstraint(US_CGE, r >=0) # This does NOT have the same behaviour...
@variable(US_CGE, w, lower_bound = 0.001 * wLi, start = wLi)#, lower_bound = 0)  
@NLconstraint(US_CGE, w >=0) # This does NOT have the same behaviour...
@variable(US_CGE, Kd_f[i = sectors], start = Kdi[i], lower_bound=0.001*Kdi[i], lower_bound=0)
#@NLconstraint(US_CGE, Kd_f >=0) # Is that K? # This doesn't run
@variable(US_CGE, Ld_f[i = sectors], start = Ldi[i], lower_bound=0.001*Ldi[i], lower_bound=0)
#@NLconstraint(US_CGE, Ld_f >=0) # This doesn't run
@variable(US_CGE, HHi, lower_bound = 0.001 * YiIn, start = YiIn, lower_bound=0)
#@NLconstraint(US_CGE, HHi >=0)
@variable(US_CGE, Pr_Commods[i = sectors], start = CPi[i], lower_bound=0.001*CPi[i], lower_bound=0) 
#@NLconstraint(US_CGE, Pr_Commods[i] >=0) #This doesn't run, of course
@variable(US_CGE, YOut[i = sectors], start = YiOut[i], lower_bound=0.001*YiOut[i], lower_bound=0) 
#@NLconstraint(US_CGE, YOut >=0) # This doesn't run anyway
@variable(US_CGE, Commodsd_HH[i = sectors], start = Cdi[i], lower_bound=0.001*Cdi[i], lower_bound=0) 
#@NLconstraint(US_CGE, Commodsd_HH >=0) # This doesn't run anyway

#@NLconstraint(US_CGE, Ke >= 0)
#@NLconstraint(US_CGE, Le >= 0)

@NLconstraint(US_CGE, ECd[i = sectors], Pr_Commods[i] * Commodsd_HH[i] == HHUexp[i] * HHi) #Consumer commodity demand function (income =spending)
@NLconstraint(US_CGE, EKd_f[i = sectors], Kd_f[i] == (YOut[i] / PFeFs[i]) * (Kexp[i] * w / (Lexp[i] * r))^Lexp[i]) #Kapital demand function (Kap = ?...)
@NLconstraint(US_CGE, ELd_f[i = sectors], Ld_f[i] == (YOut[i] / PFeFs[i]) * (Lexp[i] * r / (Kexp[i] * w))^Kexp[i]) #Labour demand function (Kap = ?...)
@NLconstraint(US_CGE, ENoProf[i = sectors], Pr_Commods[i] * YOut[i] == r * Kd_f[i] + w * Ld_f[i]) #Competitive Equlibrium, no profit (value of output=production)
@NLconstraint(US_CGE, ELs_f, sum(Ld_f[i] for i in sectors) == Le) # Competitve Eq: L demand = L supply
@NLconstraint(US_CGE, EKs_f, sum(Kd_f[i] for i in sectors) == Ke) # Competitve Eq: K demand = K supply
@NLconstraint(US_CGE, EC[i = sectors], Commodsd_HH[i] == YOut[i]) # 
@NLconstraint(US_CGE, EHHi, HHi == r * Ke + w * Le) #Income Definition
#@NLconstraint(US_CGE, ETrick, Trick == 1) # This just makes everything stay at initial values

Trick = 1
@NLobjective(US_CGE, Max, Trick) #Numeraire from EcoMod...not clear exactly how to translate
@time optimize!(US_CGE)

#Look at some results
US_CGE
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
