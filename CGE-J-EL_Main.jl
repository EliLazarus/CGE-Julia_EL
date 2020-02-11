using CSV, JuMP, Ipopt #, DataFrames
"An open source Computerised General Equlibrium model"

CGE_EL = Model(with_optimizer(Ipopt.Optimizer))
# Input Output Table
IOdata = CSV.read("IOdata.csv") #Careful Eli :: columns & rows transposed

sectors = Array{Int64}(undef,length(IOdata[1]))
for i in 1:length(IOdata[1])
    sectors[i] = i
end
commods = [1,2]

# Initial Values
rKi = 1                     #initial return to Kapital
wLi = 1                     #initial return to Labor (initial wage)
PrIndi = 1                  #Baseline for Price change/inflation
Kdi = IOdata[sectors,:Kdi]          #initial Kapital demand
Ldi = IOdata[sectors,:Ldi]          #initial Labour demand
Cdi = IOdata[sectors,:Cdi]
Invi = IOdata[sectors,:Invi]
Kei = sum(Kdi)              #initial Kapital endowment
Unempli = 10. # initial level of unemployment
Lei = sum(Ldi) + Unempli             #initial Labour endowment
# Loop to build IO square array from csv data with n sectors
IOi = Array{Int64}(undef,length(sectors),length(sectors))  
    for i in 1:length(sectors)
        for j in 2:length(sectors)+1
            IOi[i,j-1] = IOdata[i,j]
        println(i)
        end
    end
#[IOdata[sectors,:Sec1] IOdata[sectors,:Sec2]] # Input Output table (as array)

CPi = []
for i in 1:length(sectors)
    push!(CPi, 1)
end                    #initial Commodity Price Level (1 for each sector commod)

YiOut =  sum(IOi, dims=2) + Kdi + Ldi  #initial gross Total income (by sector)
YiIn = rKi * Kei + wLi * (Lei - Unempli)    #initial Income level
ConsBudgi = sum(CPi.*Cdi)
HHSavi = YiIn - ConsBudgi
TotSavi = HHSavi
mps = HHSavi/YiIn

TechCf = IOi./YiOut #this is transpose of EcoMod...

KLsubselasi = IOdata[sectors,:KLsubselasi] #[.8,1.2] # initial Kapital/Labor substitution elasticities
YinelasCommodsi = IOdata[sectors,:YinelasCommodsi]# [.9,1.1] # inititial income elasticity of commodities demand
frisch = -1.1 #  expenditure elasticity of the marginal utility of expenditure
Phili = -.1

HHUlesexp = YinelasCommodsi .* CPi.* Cdi / ConsBudgi
#HHUlesexpTot = sum(HHUlesexp) #gratuitous?
HHUlesexp = HHUlesexp/sum(HHUlesexp) #rescale marginal budgetshares by sum
HHCsubsist = Cdi + HHUlesexp*ConsBudgi./(CPi*frisch)
HHUi = prod((Cdi-HHCsubsist).^HHUlesexp)        #initial HouseHold Utility level (mod 2 = 125.373)
BankUexp = Invi .* CPi / TotSavi

CESdist = 1 ./(1 .+(wLi/rKi)*(Kdi ./Ldi) .^(-1 ./KLsubselasi))  #maybe not dot (Kdi /Ldi)
PFeFs = YiOut ./(CESdist .*Kdi .^((KLsubselasi .- 1) ./KLsubselasi) +
    (1 .-CESdist) .*Ldi .^((KLsubselasi .- 1) ./KLsubselasi)) .^
    (KLsubselasi./(KLsubselasi .-1))

    #Initial (equilibrium) levels for the endogenous variables and lower bounds to preven numerical problems in opt
@variable(CGE_EL, r, start = rKi, lower_bound = 0.001 * rKi, lower_bound = 0)   
@variable(CGE_EL, w, start = wLi, lower_bound = 0.001 * wLi, lower_bound = 0)  
@variable(CGE_EL, HHSav, start = HHSavi, lower_bound = 0.001* HHSavi, lower_bound = 0)
@variable(CGE_EL, TotSav, start = TotSavi, lower_bound = 0.001* TotSavi, lower_bound = 0)
@variable(CGE_EL, Kd_f[i = sectors], start = Kdi[i], lower_bound=0.001*Kdi[i], lower_bound=0)
@variable(CGE_EL, Ld_f[i = sectors], start = Ldi[i], lower_bound=0.001*Ldi[i], lower_bound=0)
@variable(CGE_EL, Pr_Commods[i = sectors], start = CPi[i], lower_bound=0.001*CPi[i], lower_bound=0) 
@variable(CGE_EL, YOut[i = sectors], start = YiOut[i], lower_bound=0.001*YiOut[i], lower_bound=0) 
@variable(CGE_EL, Commodsd_HH[i = sectors], start = Cdi[i], lower_bound=0.001*Cdi[i], lower_bound=0) 
@variable(CGE_EL, Inv[i = sectors], lower_bound = 0.001* Invi[i], start = Invi[i], lower_bound = 0)

@variable(CGE_EL, ConsBudg, start = ConsBudgi, lower_bound = 0.001* ConsBudgi, lower_bound = 0)
@variable(CGE_EL, PrInd, start = PrIndi, lower_bound = 0.001* PrIndi, lower_bound = 0)
@variable(CGE_EL, Unempl, start = Unempli, lower_bound = 0.001* Unempli, lower_bound = 0)
@variable(CGE_EL, Trick, start = 1)

@variable(CGE_EL, Le, lower_bound = 0.001 * Lei, start = Lei, lower_bound = 0)   
#@NLconstraint(CGE_EL, Le == Lei)
@variable(CGE_EL, Ke, lower_bound = 0.001 * Kei, start = Kei, lower_bound = 0)   
#@NLconstraint(CGE_EL, Ke == Kei) # Kapital Endowment (supply)
@variable(CGE_EL, HHi, lower_bound = 0.001 * YiIn, start = YiIn, lower_bound=0)
@NLconstraint(CGE_EL, EHHi, HHi == r * Ke + w * (Le - Unempl)) #Income Definition
@NLconstraint(CGE_EL, ETotSav, TotSav == mps*HHi) # Total Savings (fraction of Income)
@NLconstraint(CGE_EL, EConsBudg, ConsBudg == HHi - TotSav) # Consumption (Income - Savings)
@NLconstraint(CGE_EL, EPrInd[i = sectors], PrInd == sum(Pr_Commods[i] * Cdi[i] for i in sectors)/
sum(CPi[i]*Cdi[i] for i in sectors)) #Inflation Index

@NLconstraint(CGE_EL, ECd[i = sectors], Pr_Commods[i] * Commodsd_HH[i] == Pr_Commods[i] * HHCsubsist[i] +
 HHUlesexp[i] * (ConsBudg - sum(HHCsubsist[j] * Pr_Commods[j] for j in commods))) #Consumer commodity demand function (income =spending)
@NLconstraint(CGE_EL, EPhil, ((w/PrInd)/(wLi/PrIndi)-1) == Phili * ((Unempl/Le)/(Unempli/Lei))) # Wage Curve

#K and L demand
@NLconstraint(CGE_EL, EKd_f[i = sectors], Kd_f[i] == (YOut[i] / PFeFs[i]) * (CESdist[i]/r)^KLsubselasi[i] *
    (CESdist[i]^KLsubselasi[i] * r ^(1 -KLsubselasi[i]) + (1 - CESdist[i])^KLsubselasi[i] * w^(1 -KLsubselasi[i]))^
    (KLsubselasi[i]/(1-KLsubselasi[i]))) #Kapital demand function (Kap = ?...)
 @NLconstraint(CGE_EL, ELd_f[i = sectors], Ld_f[i] == (YOut[i] / PFeFs[i]) * ((1 -CESdist[i]) /w) ^KLsubselasi[i] *
    (CESdist[i] ^KLsubselasi[i] * r ^(1 -KLsubselasi[i]) + (1 - CESdist[i]) ^KLsubselasi[i] * w ^(1 -KLsubselasi[i])) ^
    (KLsubselasi[i]/(1-KLsubselasi[i]))) #Labour demand function (Kap = ?...)

@NLconstraint(CGE_EL, ESav, TotSav == HHSav) # Savings .....
@NLconstraint(CGE_EL, EInv[i = sectors], Pr_Commods[i] * Inv[i] == BankUexp[i]*TotSav) # Investment Demand?

#Market Clearing
@NLconstraint(CGE_EL, ENoProf[i = sectors], Pr_Commods[i] * YOut[i] == r * Kd_f[i] + w * Ld_f[i] +
 sum(TechCf[i,j] * YOut[i] * Pr_Commods[i] for j in sectors )) #Competitive Equlibrium, no profit (value of output=production)
@NLconstraint(CGE_EL, ELs_f, sum(Ld_f[i] for i in sectors) == Le - Unempl) # Competitve Eq: L demand = L supply
@NLconstraint(CGE_EL, EKs_f, sum(Kd_f[i] for i in sectors) == Ke) # Competitve Eq: K demand = K supply
@NLconstraint(CGE_EL, EC[i = sectors], Commodsd_HH[i] + Inv[i] + sum(TechCf[j,i] * YOut[j] for j in sectors) == YOut[i]) # 

#Utility
@variable(CGE_EL, HHU, start = HHUi)
@NLconstraint(CGE_EL, EHHU, HHU == prod((Cdi[i]-HHCsubsist[i])^HHUlesexp[i] for i in sectors))


#@variable(CGE_EL, Trick == 1)
fix(Trick, 1, force = true)
fix(Ke, Kei, force = true)
fix(Le, Lei, force = true )
fix(w, wLi, force = true)

@NLobjective(CGE_EL, Max, Trick) #Numeraire from EcoMod...not clear exactly how to translate
@time optimize!(CGE_EL)

#Look at some results
CGE_EL
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
print("HHU=           ",JuMP.value(HHU),"\n")
