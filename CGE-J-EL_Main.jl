using CSV, NamedArrays, JuMP, Ipopt 
"An open source Computerised General Equlibrium model"

CGE_EL = Model(with_optimizer(Ipopt.Optimizer))
# Input Output Table
IOdata = CSV.read("IOdata.csv") #SAM
IOdata = NamedArray(convert(Matrix,IOdata[1:size(IOdata,1),2:size(IOdata,2)]),
    (IOdata[:,1],names(IOdata[2:size(IOdata,2)]))) #Named Array to use row names
#set # sectors from data
numsectors = count(x -> occursin("Sec", string(x)), names(IOdata[:1,:])[1])
sectors = Array{Int64}(undef, 1, numsectors); for i in 1:numsectors; sectors[i] = i; end
# # commodities set manually so far
numcommonds = 2
commods = Array{Int64}(undef, 1, numcommonds); for i in 1:numcommonds; commods[i] = i; end

# Initial Values
rKi = 1                     #initial return to Kapital
wLi = 1                     #initial return to Labor (initial wage)
PrIndi = 1                  #Baseline for Price change/inflation
Kdi = IOdata["Kdi",sectors][][:]  #initial Kapital demand
Ldi = IOdata["Ldi",sectors][][:]  #initial Labour demand
#create base nominal initial Commodity Price = 1 for each sector ##if is just for running single lines
CPi = []; if length(CPi)<length(sectors); for i in 1:length(sectors); push!(CPi, 1); end; end #initial Commodity Price Level (1 for each sector commod)
Cdi = IOdata[sectors,:Cdi][][:]  #initial Investment Demand for Commodities (from data)
Invi = IOdata[sectors,:Invi][][:]#Initial Investment (from data)
Kei = sum(Kdi)              #initial Kapital endowment
Unempli = 10.               #initial level of unemployment
Lei = sum(Ldi) + Unempli    #initial Labour endowment
# Loop to build IO square array from csv data with n sectors
IOi = Array{Int64}(undef,length(sectors),length(sectors))
for i in 1:length(sectors); 
    for j in 1:length(sectors)
        IOi[i,j] = IOdata[i,j]; end; end

YiOut =  sum(IOi,dims=1)[:] + Kdi + Ldi                 #initial gross Total income (by sector)
YiIn = rKi * Kei + wLi * (Lei - Unempli)              #initial Income level
ConsBudgi = sum(CPi.*Cdi)                             #Initial Consumption Budget (all initial income)
HHSavi = YiIn - ConsBudgi                             #Initial Savings set by the difference -> in the data
TotSavi = HHSavi                                      #Savings starts

#Factors
mps = HHSavi/YiIn #marginal propensity to save:Fixed: (initial savings as a proportion of initial income)
IOtechCf = IOi./transpose(YiOut) #Input Output coefficiencts of transformation note: this is now NOT a transpose of EcoMod...
KLsubselasi = IOdata[sectors,:KLsubselasi][][:] #eg. [.8,1.2] # initial Kapital/Labor substitution elasticities
YinelasCommodsi = IOdata[sectors,:YinelasCommodsi][]#eg. [.9,1.1] # inititial income elasticity of commodities demand #(How responsive demand for each commod to changes in income)
frisch = -1.1 #  expenditure elasticity of the marginal utility of expenditure #how response the changes in utility of expenditure
Phili = -.1 # (rate of) change in wages to (rate of) change in unemployment
#Uh oh...
HHUlesexpi = YinelasCommodsi[:] .* CPi.* Cdi / ConsBudgi #Initial marginal budget shares note:gams updates assignment but I added a (i) variable
#HHUlesexpTot = sum(HHUlesexp) #gratuitous? Yes, semms like it
HHUlesexp = HHUlesexpi/sum(HHUlesexpi) #LES exponents of commods for HH utility (from initial marginal budget shares)
HHCsubsist = Cdi + HHUlesexp*ConsBudgi./(CPi*frisch) #(Stone-Geary!) subsistence quantity of consumption of each good 
HHUi = prod((Cdi-HHCsubsist).^HHUlesexp)        #initial HouseHold Utility level (mod 5 = 108.14)
BankUexp = Invi .* CPi / TotSavi  #(Cobb-Douglas) exponenent for Bank's Utility function 
CESdist = 1 ./(1 .+(wLi/rKi)*(Kdi ./Ldi) .^(-1 ./KLsubselasi))  #maybe not dot (Kdi /Ldi) #Constant Elasticity of substitution
# parameters in the production function (ie. how much labor to Kapital)
PFeFs = YiOut ./(CESdist .*Kdi .^((KLsubselasi .- 1) ./KLsubselasi) +
    (1 .-CESdist) .*Ldi .^((KLsubselasi .- 1) ./KLsubselasi)) .^
    (KLsubselasi./(KLsubselasi .-1))  #Production Effiency Factor (CES production)

#Variables at Initial (equilibrium) levels for the endogenous variables and lower bounds to prevent numerical problems in opt
#Labor and Kapital
@variable(CGE_EL, w, start = wLi, lower_bound = 0.001 * wLi, lower_bound = 0) #wages (return to labor)
@variable(CGE_EL, r, start = rKi, lower_bound = 0.001 * rKi, lower_bound = 0) #interest/rent (return to Kapital)
@variable(CGE_EL, Le, start = Lei, lower_bound = 0, lower_bound = 0.001 * Lei) #Labor supply
@variable(CGE_EL, Ke, lower_bound = 0.001 * Kei, start = Kei, lower_bound = 0) #Kapital supply
@variable(CGE_EL, Unempl, start = Unempli, lower_bound = 0.001* Unempli, lower_bound = 0) #Level of Unemployment

#Commodities
@variable(CGE_EL, Pr_Commods[i = sectors], start = CPi[i], lower_bound=0.001*CPi[i], lower_bound=0) #Commodity Prices
@variable(CGE_EL, Commodsd_HH[i = sectors], start = Cdi[i], lower_bound=0.001*Cdi[i], lower_bound=0) #Commodity Demand
@variable(CGE_EL, YOut[i = sectors], start = YiOut[i], lower_bound=0.001*YiOut[i], lower_bound=0) #Output (per sector)

#Household savings and income
@variable(CGE_EL, HHSav, start = HHSavi, lower_bound = 0.001* HHSavi, lower_bound = 0) #HH Savings
@variable(CGE_EL, TotSav, start = TotSavi, lower_bound = 0.001* TotSavi, lower_bound = 0) #What is Total Savings?
@variable(CGE_EL, Inv[i = sectors], lower_bound = 0.001* Invi[i], start = Invi[i], lower_bound = 0) #Investment Demand for Commodities
@variable(CGE_EL, HHI, lower_bound = 0.001 * YiIn, start = YiIn, lower_bound=0) #HH Income
@NLconstraint(CGE_EL, EHHI, HHI == r * Ke + w * (Le - Unempl)) #Total HH Income Definition
@variable(CGE_EL, ConsBudg, start = ConsBudgi, lower_bound = 0.001* ConsBudgi, lower_bound = 0) #Consumer Budget
@NLconstraint(CGE_EL, EConsBudg, ConsBudg == HHI - TotSav) # Consumption Budget, as remainder (Income - Savings)

@NLconstraint(CGE_EL, ECd[i = sectors], Pr_Commods[i] * Commodsd_HH[i] == Pr_Commods[i] * HHCsubsist[i] +
 HHUlesexp[i] * (ConsBudg - sum(HHCsubsist[j] * Pr_Commods[j] for j in commods))) #Consumer commodity demand function (income =spending)

#K and L demand
@variable(CGE_EL, Ld_f[i = sectors], start = Ldi[i], lower_bound=0.001*Ldi[i], lower_bound=0) #Labor demand (per sector)
@variable(CGE_EL, Kd_f[i = sectors], start = Kdi[i], lower_bound=0.001*Kdi[i], lower_bound=0) #Kapital demand (per sector)
@NLconstraint(CGE_EL, EKd_f[i = sectors], Kd_f[i] ==
    (YOut[i] / PFeFs[i]) * (CESdist[i] /r)^KLsubselasi[i] *
    (CESdist[i]^KLsubselasi[i] * r^(1 -KLsubselasi[i]) +
    (1 - CESdist[i])^KLsubselasi[i] * w^(1 -KLsubselasi[i]))^
    (KLsubselasi[i]/(1-KLsubselasi[i]))) #Kapital demand function (Kap = ?...)
@NLconstraint(CGE_EL, ELd_f[i = sectors], Ld_f[i] ==
    (YOut[i] / PFeFs[i]) * ((1 -CESdist[i]) /w)^KLsubselasi[i] *
    (CESdist[i]^KLsubselasi[i] * r^(1 -KLsubselasi[i]) + 
    (1 - CESdist[i])^KLsubselasi[i] * w^(1 -KLsubselasi[i]))^
    (KLsubselasi[i]/(1-KLsubselasi[i]))) #Labor demand function (L = ?...)

@NLconstraint(CGE_EL, EHHSav, HHSav == mps * HHI) # Household Savings, proportion of Income
@NLconstraint(CGE_EL, ESav, TotSav == HHSav) # Total Savings .....
#@NLconstraint(CGE_EL, ETotSav, TotSav == mps*HHI) # Total Savings (Fixed fraction of Income)

#Markets Clearing
@NLconstraint(CGE_EL, ENoProf[i = sectors], Pr_Commods[i] * YOut[i] == r * Kd_f[i] + w * Ld_f[i] +
 sum(IOtechCf[j,i] * YOut[i] * Pr_Commods[j] for j in sectors )) #Competitive Equlibrium, no profit (value of output [for each sector])=production [of each sector])
@NLconstraint(CGE_EL, ELs_f, sum(Ld_f[i] for i in sectors) == Le - Unempl) # Competitve Eq: Total L demand = Total L supply
@NLconstraint(CGE_EL, EKs_f, sum(Kd_f[i] for i in sectors) == Ke) # Competitve Eq: Total K demand = Total K supply
@NLconstraint(CGE_EL, EC[i = sectors], Commodsd_HH[i] + Inv[i] + sum(IOtechCf[i,j] * YOut[j] for j in sectors) ==
 YOut[i]) #Market clearing for commodities (Sum of consumption, investment [and production uses] equals Total Output)
@NLconstraint(CGE_EL, EInv[i = sectors], Pr_Commods[i] * Inv[i] == BankUexp[i]*TotSav) # Investment Demand = Investment Supply

#Price Index
@variable(CGE_EL, PrInd, start = PrIndi, lower_bound = 0.001* PrIndi, lower_bound = 0) #Price/Inflation Index
@NLconstraint(CGE_EL, EPrInd[i = sectors], PrInd == sum(Pr_Commods[i] * Cdi[i] for i in sectors)/
sum(CPi[i]*Cdi[i] for i in sectors)) #Inflation Index

@NLconstraint(CGE_EL, EPhil, ((w/PrInd)/(wLi/PrIndi)-1) == Phili * ((Unempl/Le)/(Unempli/Lei)-1)) # Wage(/unemployment) Curve

@variable(CGE_EL, Trick, start = 1)

#Utility
@variable(CGE_EL, HHU, start = HHUi)

fix(Trick, 1, force = true)
fix(Ke, Kei, force = true)
fix(Le, Lei*1.1, force = true)
fix(w, wLi, force = true)

@NLobjective(CGE_EL, Max, Trick) #Numeraire from EcoMod...not clear exactly how to translate
@time optimize!(CGE_EL)

HHU = prod((JuMP.value(Commodsd_HH[i])-HHCsubsist[i])^HHUlesexp[i] for i in sectors) #Final HH Utility calculated AFTER model solved

#Look at some results
CGE_EL
print("w=             ",JuMP.value(w),"\n")
print("r=             ",JuMP.value(r),"\n")
print("HHI=           ",JuMP.value(HHI),"\n")
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
print("HHU=           ",HHU,"\n")
print("TotSav=        ",JuMP.value(TotSav),"\n")