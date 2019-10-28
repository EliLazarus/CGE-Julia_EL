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
rKi = 1                 #initial return to Kapital
wLi = 1                 #initial return to Labor (initial wage)
Kdi = IOdata[:Kdi]          #initial Kapital demand
Ldi = IOdata[:Ldi]          #initial Labour demand
Cdi = IOdata[:Cdi]          #initial Consumption demand
Kei = sum(Kdi)          #initial Kapital endowment
Lei = sum(Ldi)          #initial Labour endowment
YiOut =  rKi * Kdi + wLi * Ldi  #initial gross Total income (by sector)
YiIn = rKi * Kei + wLi * Lei    #initial Income level

CPi = []
for i in 1:length(IOdata[1])
    push!(CPi, 1)
end                    #initial Commodity Price Level (1 for each sector commod)

HHUexp = CPi .* Cdi / YiIn
HHUi = Cdi.^HHUexp        #initial HouseHold Utility level
Le = Lei * 1.1     #Labour grows at 1.1
Ke = Kei * 1.2     #Kapital grows at 1.2
w = wLi * 5

## Parameters per sector
Kexp = rKi * Kdi / (rKi * Kdi + wLi * Ldi)  #exponent on K in Cobb-Douglas Production Function
Lexp = 1 .- Kexp                     ##exponent on L in Cobb-Douglas Production Function
PFeFs = YiOut ./ (Kdi.^Kexp .* Ldi.^Lexp)  #Production Function efficiency of Firms ???shld this be dot product or matrix product???

#Initial (equilibrium) levels for the endogenous variables and lower bounds to preven numerical problems in opt
@variable(US_CGE, r, lower_bound = 0.001 * rKi, start = rKi) 
#Kd_f=Kdi
@variable(US_CGE, Kd_f[i = sectors], start = Kdi[i]) #??? lower_bound not working...lower_bound=0.001*Kdi, start not working , start=Kdi
@variable(US_CGE, w, lower_bound = 0.001 * wLi, start = wLi) 
#Ld_f = wdi
@variable(US_CGE, Ld_f[i = sectors], start = Ldi[i]) #??? lower_bound not working...lower_bound=0.001*Ldi, start=Ldi
@variable(US_CGE, HHi, lower_bound = 0.001 * YiIn, start = YiIn)
#Pr_Commods = CPi
@variable(US_CGE, Pr_Commods[i = sectors], start = CPi[i]) #??? lower_bound not working..., lower_bound=0.001*CPi, start=CPi
#YOut = YiOut
@variable(US_CGE, YOut[i = sectors], start = YiOut[i]) #??? lower_bound not working..., lower_bound=0.001*YiOut, start=YiOut
#Commodsd_HH = Cdi
@variable(US_CGE, Commodsd_HH[i = sectors], start = Cdi[i])  #??? lower_bound not working...lower_bound=0.001*Cdi, start=Cdi

#?????? Ut = Commodsd_HH.^HHUexp
@NLconstraint(US_CGE, EKd_f[i = sectors], Kd_f[i] == (YOut[i] / PFeFs[i]) * (Kexp[i] * w / Lexp[i] * r)^Lexp[i]) #Kapital demand function (Kap = ?...)
@NLconstraint(US_CGE, ELd_f[i = sectors], Ld_f[i] == (YOut[i] / PFeFs[i]) * (Lexp[i] * r / Kexp[i] * w)^Kexp[i]) #Labour demand function (Kap = ?...)
@NLconstraint(US_CGE, ELs_f, sum(Ld_f[i] for i in sectors) == Le) # Competitve Eq: L demand = L supply
@NLconstraint(US_CGE, EKs_f, sum(Kd_f[i] for i in sectors) == Ke) # Competitve Eq: K demand = K supply
@NLconstraint(US_CGE, EHHi, HHi == r * Ke + w * Le) #Income Definition
@NLconstraint(US_CGE, ENoProf[i = sectors], Pr_Commods[i] * YOut[i] == r * Kd_f[i] + w * Ld_f[i]) #Competitive Equlibrium, no profit (value of output=production)
@NLconstraint(US_CGE, EC[i = sectors], Commodsd_HH[i] == YOut[i]) # 
@NLconstraint(US_CGE, ECd[i = sectors], Pr_Commods[i] * Commodsd_HH[i] == HHUexp[i] * HHi) #Consumer commodity demand function (income =spending)

Trick = 1
@NLobjective(US_CGE, Max, Trick) #Numeraire from EcoMod...not clear exactly how to translate
@time optimize!(US_CGE)

#Look at some results
US_CGE
JuMP.value(w)
JuMP.value(r)
JuMP.value(HHi)
JuMP.value(Ld_f[1])
JuMP.value(Ld_f[2])
JuMP.value(Kd_f[1])
JuMP.value(Kd_f[2])
JuMP.value(Pr_Commods[1])
JuMP.value(Pr_Commods[2])
JuMP.value(Commodsd_HH[1])
JuMP.value(Commodsd_HH[2])
JuMP.value(YOut[1])
JuMP.value(YOut[2])
