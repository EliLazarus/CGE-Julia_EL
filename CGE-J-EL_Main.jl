using CSV, NamedArrays, JuMP, Ipopt, DataFrames, XLSX
"An open source Computerised General Equlibrium model"

"*COLUMN* accounts record *SPENDING*" # *ROW* acounts record *INCOME*

CGE_EL = Model(with_optimizer(Ipopt.Optimizer))#with_optimizer()
"A. Wächter and L. T. Biegler, On the Implementation of a Primal-Dual Interior Point Filter Line Search Algorithm for Large-Scale Nonlinear Programming, Mathematical Programming 106(1), pp. 25-57, 2006 (preprint)"

# SAM Table
SAMdataX = XLSX.readxlsx("SAMdataPlus.xlsx")
SAMdataY = XLSX.gettable(SAMdataX["SAMdataPlus"], header=false)
SAMdata = DataFrames.DataFrame(SAMdataY[1])
SAMdata = NamedArray(convert(Matrix,SAMdata[2:size(SAMdata,1),2:size(SAMdata,2)]),(SAMdata[2:end,1],values.(SAMdata[1,2:end]))) #Named Array to use row names

#set # sectors from data
numsectors = count(x -> occursin("Sec", string(x)), names(SAMdata[:1,:])[1])
sectors = Array{Int64}(undef, 1, numsectors); for i in 1:numsectors; sectors[i] = i; end
numcommonds = convert(Int64,SAMdata[1,"numcommonds"]) #So far as factor in data sheet (until commods get own SAM row/columns)
commods = Array{Int64}(undef, 1, numcommonds); for i in 1:numcommonds; commods[i] = i; end

# Initial Values
rKi = 1                     #initial return to Kapital
wLi = 1                     #initial return to Labor (initial wage)
PrIndi = 1                  #Baseline for Price change/inflation
Pr_W_Im_foreignCurri = SAMdata[sectors,"Pr_W_Im_foreignCurri"][][:] # Initial prices of Imports at World price in foreign currency
Pr_W_Exp_foreignCurri = SAMdata[sectors,"Pr_W_Exp_foreignCurri"][][:] # Initial prices of Exports at World price in foreign currency
Pr_Commods_toHomei = SAMdata[sectors,"Pr_Commods_toHomei"][][:] # Initial prices of Domestic commodity/output of firm to domestic market
Pr_CombCommods_toHomei = SAMdata[sectors,"Pr_CombCommods_toHomei"][][:] # Initial Prices of Combined Domestic and Foreign commodities to home market
Kdi = SAMdata["Kdi",sectors][][:]  #initial Kapital demand
Ldi = SAMdata["Ldi",sectors][][:]  #initial Labour demand
Cdi = SAMdata["Cdi",sectors][][:]  #initial Consumer Demand for Commodities (from data)
GKdi = SAMdata["Kdi","Gov"] #initial Govt Kapital demand (from data)
GLdi = SAMdata["Ldi","Gov"] #initial Govt Labour demand (from data)
GSavi = SAMdata["SavInv","Gov"] # Government Savings (0: balanced budget, >0: surplus, <0: deficit )
TaxInRevi = SAMdata["TaxIn","HH"] #Inititial Income tax
UnemplBenRate = SAMdata[1,"UnemplBenRate"] #What proportion wage HH gets if unemployed
Transfi = SAMdata["HH","Gov"] #Government direct tranfers to Households
TransfOthi = 15. # Other Transfers... I don't understand where this fits in the SAM...
GovCdi = SAMdata[sectors,"Gov"][][:] # Government purchased from firms 
TaxCRevi = SAMdata["TaxC", sectors][][:]    #initial Consumption tax
TaxKRevi = SAMdata["TaxK",sectors][][:]    #initial Kapital use tax
TaxLRevi = SAMdata["TaxL",sectors][][:]   #initial Payroll tax? (labour use)

ForeignSavi = SAMdata["SavInv","RoW"]
#create base nominal initial Commodity Price = 1 for each sector ##the if is just for running single lines during debugging etc
CPi = []; if length(CPi)<length(sectors); for i in 1:length(sectors); push!(CPi, 1); end; end #initial Commodity Price Level (1 for each sector commod)
Invi = SAMdata[sectors,"InvSav"][][:]        #Initial Investment (from data)
Expi = SAMdata["Expi",sectors][][:]         #Initial Export demand (from data)
Impi = SAMdata["Impi",sectors][][:]         #Initial Import demand (from data)
TaxImpRevi = SAMdata["TaxImp",sectors][][:] #initial Import tax (tariff) revenue
Kei = sum(Kdi) + GKdi                       #initial Kapital endowment (assume supply = demand?)
Unempli = SAMdata[1,"Unempli"]               #initial level of unemployment
XRatei = SAMdata[1,"Xchangei"]               #initial exchange rate (1 obviously)
Lei = sum(Ldi) + GLdi + Unempli    #initial Labour endowment

Pr_Exp_DomCurri = Pr_W_Exp_foreignCurri * XRatei # Export Prices

# Loop to build IO square sector by sector array from csv data with n sectors
# IMPORTANT: Only works so far if sectors are the first columns and rows
IOi = Array{Int64}(undef,length(sectors),length(sectors))
for i in 1:length(sectors); 
    for j in 1:length(sectors)
        IOi[i,j] = SAMdata[i,j]; end; end

YOuti =  sum(IOi,dims=1)[:] + Kdi + TaxKRevi + Ldi + TaxLRevi  #initial gross Total income (by sector)
Yout_toHomei = YOuti - Expi                         # Domestic consumption from domestic production (everything - exports)
YCombOut_toHomei = Yout_toHomei + Impi + TaxImpRevi # Total Domestic consumption (dom cons from dom prod + imports + tarrif revenue)
YIni = rKi * Kei + wLi * (Lei - Unempli) + Transfi  # Initial Income level
ConsBudgi = sum(CPi.*Cdi) + sum(TaxCRevi)           # Initial Consumption Budget (all initial income)
HHSavi = YIni - ConsBudgi - TaxInRevi               # Initial Savings set by the difference -> in the data
TotSavi = HHSavi + GSavi * PrIndi + ForeignSavi * XRatei # Savings starts
TaxTotRevi = sum(TaxCRevi+TaxKRevi+TaxLRevi+TaxImpRevi) + TaxInRevi # Initial total government revenue
TaxCRate = TaxCRevi./Cdi.*Pr_CombCommods_toHomei # Consumption tax rate
TaxKRate = TaxKRevi./Kdi.*rKi                    # Kapital use tax rate from data
TaxLRate = TaxLRevi./Ldi.*wLi                    # Labour use tax rate from data
TaxCRatei = TaxCRate                             # initial Consumption tax rate (for Price Index)
TaxInRate = TaxInRevi/YIni                       #Income tax rate from data (Revenue/Total Income)
TaxImpRate = TaxImpRevi./(Impi.*Pr_W_Im_foreignCurri * XRatei) #Imports tax rate (tarrifs) from data (Revenue/Total Imports)
Pr_Im_DomCurri = (1 .+ TaxImpRate) .* Pr_W_Im_foreignCurri * XRatei # Import prices at domestic currency

#Factors
frisch = SAMdata[1,"Frisch"]      # expenditure elasticity of the marginal utility of expenditure #how response the changes in utility of expenditure
Phili = SAMdata[1,"Philli"]       # (rate of) change in wages to (rate of) change in unemployment
mps = HHSavi/(YIni - TaxInRevi)  # marginal propensity to save:Fixed: (initial savings as a proportion of initial income)
IOtechCf = IOi./transpose(YOuti) # Input Output coefficiencts of transformation note: this is now NOT a transpose of EcoMod...
ArmSubElasi = SAMdata[sectors,"ArmSubElasi"][][:]       # Initial sustitution elasticities of the Armington function (between foreign and domestic)
TransformElasi = SAMdata[sectors,"TransformElasi"][][:] # Initital elasticities of transformation in CET funtion
KLsubselasi = SAMdata[sectors,"KLsubselasi"][][:]#eg. [.8,1.2 - values MUST not be =1.] # initial Kapital/Labor substitution elasticities
YinelasCommodsi = SAMdata[sectors,"YinelasCommodsi"][][:]#eg. [.9,1.1] # inititial income elasticity of commodities demand #(How responsive demand for each commod to changes in income)
HHUlesexpi = YinelasCommodsi[:] .*(1 .+ TaxCRate).* Pr_CombCommods_toHomei.* Cdi / ConsBudgi #Initial marginal budget shares note:gams updates assignment but I added a (i) variable
HHUlesexp = HHUlesexpi./sum(HHUlesexpi) # nested ELES exponents for HH utility (from initial marginal budget shares)
HHCsubsist = Cdi + HHUlesexp*ConsBudgi./(Pr_CombCommods_toHomei*frisch.*(1 .+ TaxCRate)) #(Stone-Geary!) subsistence quantity of consumption of each good 
HHUi = prod((Cdi-HHCsubsist).^HHUlesexp) # initial HouseHold Utility level
BankUexp = Invi .* Pr_CombCommods_toHomei/ TotSavi # (Cobb-Douglas) exponenent for Bank's Utility function 
CESdist = 1 ./(1 .+((1 .+TaxLRate)*wLi)./((1 .+TaxKRate)*rKi).*(Kdi ./Ldi).^(-1 ./KLsubselasi)) # Constant Elasticity of substitution
# parameters in the production function (ie. how much labor compared to Kapital)
PFeFs = YOuti ./(CESdist .*Kdi .^((KLsubselasi .- 1) ./KLsubselasi) +
    (1 .-CESdist) .*Ldi .^((KLsubselasi .- 1) ./KLsubselasi)) .^
    (KLsubselasi./(KLsubselasi .-1))  #Production Effiency Factor (CES production)
Armexp = 1 ./(1 .+(Pr_Commods_toHomei./Pr_Im_DomCurri).*(Impi./Yout_toHomei) .^(-1 ./ArmSubElasi)) # Exponent for Armington function (how much imports to domestic)
ArmEF = YCombOut_toHomei./(Armexp.*Impi.^ ((ArmSubElasi .-1)./ArmSubElasi) .+
    (1 .- Armexp).*Yout_toHomei.^ ((ArmSubElasi .-1)./ArmSubElasi)).^ (ArmSubElasi./(ArmSubElasi .-1)) #Efficiency (coefficient) for Armington
TransElasexp_Yout  = 1 ./(1 .+(Pr_Commods_toHomei./Pr_Exp_DomCurri).*(Expi ./Yout_toHomei) .^(-1 ./TransformElasi)) 
ShiftParTransElas_Yout = YOuti./(TransElasexp_Yout.*Expi.^((TransformElasi .-1)./TransformElasi) .+
 (1 .- TransElasexp_Yout).*Yout_toHomei .^ ((TransformElasi .-1)./TransformElasi)).^(TransformElasi ./(TransformElasi .-1)) #...

GovUExpC = Pr_CombCommods_toHomei.*GovCdi/(TaxTotRevi-Transfi-PrIndi*GSavi)
GovUExpL = wLi * GLdi/(TaxTotRevi-Transfi-PrIndi*GSavi)
GovUExpK = rKi * GKdi/(TaxTotRevi-Transfi-PrIndi*GSavi)

#Variables at Initial (equilibrium) levels for the endogenous variables and lower bounds to prevent numerical problems in opt
#Labor and Kapital
@variable(CGE_EL, w, start = wLi, lower_bound = 0.001 * wLi, lower_bound = 0) #wages (return to labor)
@variable(CGE_EL, r, start = rKi, lower_bound = 0.001 * rKi, lower_bound = 0) #interest/rent (return to Kapital)
@variable(CGE_EL, Le, start = Lei, lower_bound = 0.001 * Lei, lower_bound = 0) #Labor supply
@variable(CGE_EL, Ke, start = Kei, lower_bound = 0.001 * Kei, lower_bound = 0) #Kapital supply
@variable(CGE_EL, Unempl, start = Unempli, lower_bound = 0.001* Unempli, lower_bound = 0) #Level of Unemployment

#Commodities
@variable(CGE_EL, Pr_CombCommods_toHome[i = sectors], start = Pr_CombCommods_toHomei[i], lower_bound = 0.001 * Pr_CombCommods_toHomei[i], lower_bound = 0) #
@variable(CGE_EL, Pr_Commods[i = sectors], start = CPi[i], lower_bound=0.001*CPi[i], lower_bound=0) #Commodity Prices for domestically produced
@variable(CGE_EL, Pr_Commods_toHome[i = sectors], start = Pr_Commods_toHomei[i], lower_bound = 0.001 * Pr_Commods_toHomei[i], lower_bound = 0) #Prices of domestic commods for domestic consumption
@variable(CGE_EL, Pr_Exp_DomCurr[i = sectors], start = Pr_Exp_DomCurri[i], lower_bound = 0.001 * Pr_Exp_DomCurri[i], lower_bound = 0) #Prices of exports at port (f.o.b.?)
@variable(CGE_EL, Pr_Im_DomCurr[i = sectors], start = Pr_Im_DomCurri[i], lower_bound = 0.001 * Pr_Im_DomCurri[i], lower_bound = 0) #Princes of Imports (in domestic currency [tarrifs and at exchange])
@variable(CGE_EL, Commodsd_HH[i = sectors], start = Cdi[i], lower_bound=0.001*Cdi[i], lower_bound=0) #Commodity Demand
@variable(CGE_EL, YCombOut_toHome[i = sectors], start = YCombOut_toHomei[i], lower_bound=0.001*YCombOut_toHomei[i], lower_bound=0) #Combined Output -dom and imports - for domestic consumption (per sector)
@variable(CGE_EL, Yout_toHome[i = sectors], start = Yout_toHomei[i], lower_bound=0.001*Yout_toHomei[i], lower_bound=0) #Domestic Output - for dom consumption (per sector)
@variable(CGE_EL, Exp[i = sectors], start = Expi[i], lower_bound=0.001*Expi[i], lower_bound=0) #Export supply (per sector)
@variable(CGE_EL, Imp[i = sectors], start = Impi[i], lower_bound=0.001*Impi[i], lower_bound=0) #Import supply (per sector)
@variable(CGE_EL, YOut[i = sectors], start = YOuti[i], lower_bound=0.001*YOuti[i], lower_bound=0) #Total Output (per sector)

# Government
@variable(CGE_EL, GovCd[i = sectors], start = GovCdi[i], lower_bound=0.001*GovCdi[i], lower_bound=0) #Governmend commodity demand (per sector)
@variable(CGE_EL, GKd, start = GKdi, lower_bound=0.001*GKdi, lower_bound=0) #Government Kapital demand
@variable(CGE_EL, GLd, start = GLdi, lower_bound=0.001*GLdi, lower_bound=0) #Government Labour demand
@variable(CGE_EL, TaxTotRev, start = TaxTotRevi, lower_bound=0.001*TaxTotRevi, lower_bound=0) # Tax Revenue
@variable(CGE_EL, Transf, start = Transfi, lower_bound=0) # Total transfers (unemployment + other)     #lower_bound=0.001*Transfi,
@variable(CGE_EL, TransfOth, start = TransfOthi, lower_bound=0.001*TransfOthi, lower_bound=0) # Other (non-unemployment) transfers
@variable(CGE_EL, GSav, start = GSavi, lower_bound=0.001*GSavi, lower_bound=0) # Government Savings

#Household savings and income
@variable(CGE_EL, HHSav, start = HHSavi, lower_bound = 0.001* HHSavi, lower_bound = 0) #HH Savings
@variable(CGE_EL, TotSav, start = TotSavi, lower_bound = 0.001* TotSavi, lower_bound = 0) # Total Savings (HH+Gov+Foreign)
@variable(CGE_EL, Inv[i = sectors], lower_bound = 0.001* Invi[i], start = Invi[i], lower_bound = 0) #Investment Demand for Commodities
@variable(CGE_EL, HHI, start = YIni, lower_bound = 0.001 * YIni, lower_bound=0) #HH Income
@NLconstraint(CGE_EL, EHHI, HHI == r * Ke + w * (Le - Unempl) + Transf) #Total HH Income Definition ('interest', wages, and gov transfers)
@variable(CGE_EL, ConsBudg, start = ConsBudgi, lower_bound = 0.001* ConsBudgi, lower_bound = 0) #Consumer Budget
@NLconstraint(CGE_EL, EConsBudg, ConsBudg == (1-TaxInRate)*HHI - HHSav) # Consumption Budget = after tax income - Savings (balanced)

@NLconstraint(CGE_EL, ECd[i = sectors], (1 + TaxCRate[i]) * Pr_CombCommods_toHome[i] * Commodsd_HH[i] ==
 (1 + TaxCRate[i]) * Pr_CombCommods_toHome[i] * HHCsubsist[i] + HHUlesexp[i] * 
 (ConsBudg - sum(HHCsubsist[j] * (1 + TaxCRate[j]) * Pr_CombCommods_toHome[j] for j in commods))) #Consumer commodity demand function (income =spending)

#K and L demand
@variable(CGE_EL, Ld_f[i = sectors], start = Ldi[i], lower_bound=0.001*Ldi[i], lower_bound=0) #Labor demand (per sector)
@variable(CGE_EL, Kd_f[i = sectors], start = Kdi[i], lower_bound=0.001*Kdi[i], lower_bound=0) #Kapital demand (per sector)
@NLconstraint(CGE_EL, EKd_f[i = sectors], Kd_f[i] == (YOut[i] / PFeFs[i]) *
    (CESdist[i] /((1+TaxKRate[i])*r))^KLsubselasi[i] *
    (CESdist[i]^KLsubselasi[i] *((1+TaxKRate[i])*r)^(1 -KLsubselasi[i]) +
    (1 - CESdist[i])^KLsubselasi[i] *((1+TaxLRate[i])*w)^(1 -KLsubselasi[i]))^
    (KLsubselasi[i]/(1-KLsubselasi[i]))) #Kapital demand function for firms (Kap = ?...)
@NLconstraint(CGE_EL, ELd_f[i = sectors], Ld_f[i] == (YOut[i] / PFeFs[i]) *
    ((1 -CESdist[i]) /((1+TaxLRate[i])*w))^KLsubselasi[i] *
    (CESdist[i]^KLsubselasi[i] * ((1+TaxKRate[i])*r)^(1 -KLsubselasi[i]) + 
    (1 - CESdist[i])^KLsubselasi[i] * ((1+TaxLRate[i])*w)^(1 -KLsubselasi[i]))^
    (KLsubselasi[i]/(1-KLsubselasi[i]))) #Labor demand function for firms (L = ?...)

#Price Index
@variable(CGE_EL, PrInd, start = PrIndi, lower_bound = 0.001* PrIndi, lower_bound = 0) #Price/Inflation Index
@NLconstraint(CGE_EL, EPrInd, PrInd == sum((1+TaxCRate[i]) * Pr_CombCommods_toHome[i] * Cdi[i] for i in sectors)/
    sum((1+TaxCRatei[i])* Pr_CombCommods_toHomei[i] * Cdi[i] for i in sectors)) #Inflation Index
@variable(CGE_EL, XRate, start = XRatei, lower_bound = 0.001 * XRatei, lower_bound = 0) # Exchange rate btw domestic and RoW (based on trade balance?)

#Savings
@NLconstraint(CGE_EL, EHHSav, HHSav == mps * (HHI - TaxInRate * HHI)) # Household Savings, proportion of Income
@variable(CGE_EL, ForeignSav, start = ForeignSavi) # Foreign Savings      # , lower_bound = ForeignSavi*0.001, lower_bound = 0
@NLconstraint(CGE_EL, ESav, TotSav == HHSav + GSav * PrInd + ForeignSav * XRate) # Total Savings .....

#Government
@NLconstraint(CGE_EL, EGCd[i = sectors], Pr_CombCommods_toHome[i] * GovCd[i] == GovUExpC[i]*(TaxTotRev - Transf - GSav * PrInd)) #Gov consumer purchases = Gov budget for commods
@NLconstraint(CGE_EL, EGKd, r * GKd == GovUExpK*(TaxTotRev - Transf - GSav * PrInd)) # Gov spending on Kapital = Govt budget for Kapital
@NLconstraint(CGE_EL, EGLd, w * GLd == GovUExpL*(TaxTotRev - Transf - GSav * PrInd)) # Gov spending on Labour = Govt budget for Labour
@NLconstraint(CGE_EL, ETaxRev[i=sectors], TaxTotRev == TaxInRate*HHI + sum(Pr_CombCommods_toHome[i]*TaxCRate[i]*Commodsd_HH[i] +
    TaxKRate[i]*Kd_f[i]*r + TaxLRate[i]*Ld_f[i]*w + TaxImpRate[i]*Imp[i]*Pr_W_Im_foreignCurri[i]*XRate for i in sectors)) #Tax Revenue = tax rates*amount of taxed quantites 
@NLconstraint(CGE_EL, Etransf, Transf == UnemplBenRate*w*Unempl + TransfOth*PrInd) # Govt transfers for unemployment and 'other'

#Trade
@NLconstraint(CGE_EL, EExp[i = sectors], Exp[i] == (YOut[i]/ShiftParTransElas_Yout[i])*
(TransElasexp_Yout[i]/Pr_Exp_DomCurr[i])^TransformElasi[i]*
((TransElasexp_Yout[i]^TransformElasi[i])*
(Pr_Exp_DomCurr[i]^(1-TransformElasi[i])) +
((1-TransElasexp_Yout[i])^TransformElasi[i])*
(Pr_Commods_toHome[i]^(1-TransformElasi[i])))^(TransformElasi[i]/(1-TransformElasi[i]))) # Exports =....

@NLconstraint(CGE_EL, EYOut_toHome[i = sectors], Yout_toHome[i] == (YOut[i]/ShiftParTransElas_Yout[i])*
((1-TransElasexp_Yout[i])/Pr_Commods_toHome[i])^TransformElasi[i]*
((TransElasexp_Yout[i]^TransformElasi[i])*
(Pr_Exp_DomCurr[i]^(1-TransformElasi[i])) +
((1-TransElasexp_Yout[i])^TransformElasi[i])*
(Pr_Commods_toHome[i]^(1-TransformElasi[i])))^(TransformElasi[i]/(1-TransformElasi[i]))) # Domestic output for domestic consumption = ...

@NLconstraint(CGE_EL, EPr_Imp[i = sectors], Pr_Im_DomCurr[i] == (1+TaxImpRate[i])*XRate*Pr_W_Im_foreignCurri[i]) # Price of imports is world price (at exchange rate to domestic currency) and tarrifs
@NLconstraint(CGE_EL, EPr_Exp[i = sectors], Pr_Exp_DomCurr[i] == Pr_W_Exp_foreignCurri[i]*XRate) # Price of exports at port (at world prices with exchange rates to make domestic currency)
@NLconstraint(CGE_EL, ETransform[i = sectors], Pr_Commods[i] * YOut[i] == Pr_Exp_DomCurr[i]*Exp[i]+Pr_Commods_toHome[i]*Yout_toHome[i]) #Total commodities produced = exports and domestic consumption (at domestic currency) 
@NLconstraint(CGE_EL, EArmImpDemand[i = sectors], Imp[i] == (YCombOut_toHome[i]/ArmEF[i])*
(Armexp[i]/Pr_Im_DomCurr[i])^ArmSubElasi[i]*
((Armexp[i]^ArmSubElasi[i])*
(Pr_Im_DomCurr[i]^(1-ArmSubElasi[i])) +
((1-Armexp[i])^ArmSubElasi[i])*
(Pr_Commods_toHome[i]^(1-ArmSubElasi[i]))
)^(ArmSubElasi[i]/(1-ArmSubElasi[i]))) # Imports = ... (some function that determines the proportion of imports to domestic production for demand...)

@NLconstraint(CGE_EL, EArmDomDemand[i = sectors],Yout_toHome[i] == (YCombOut_toHome[i]/ArmEF[i])*
((1-Armexp[i])/Pr_Commods_toHome[i])^ArmSubElasi[i]*
((Armexp[i]^ArmSubElasi[i])*
(Pr_Im_DomCurr[i]^(1-ArmSubElasi[i])) +
((1-Armexp[i])^ArmSubElasi[i])*
(Pr_Commods_toHome[i]^(1-ArmSubElasi[i]))
)^(ArmSubElasi[i]/(1-ArmSubElasi[i]))) # Domestic consumption = ... (some function that determines the proportion of imports to domestic production for demand...)

@NLconstraint(CGE_EL, EArmington[i = sectors], Pr_CombCommods_toHome[i] * YCombOut_toHome[i] == Pr_Im_DomCurr[i] * Imp[i] + Pr_Commods_toHome[i] * Yout_toHome[i]) #Total combined commods = combined imports and domestically produced at domestic curreny

#Markets Clearing
@NLconstraint(CGE_EL, ENoProf[i = sectors], Pr_Commods[i] * YOut[i] == (1 +TaxKRate[i])*r * Kd_f[i] +  (1 +TaxLRate[i])*w * Ld_f[i] +
 sum(IOtechCf[j,i] * YOut[i] * Pr_CombCommods_toHome[j] for j in commods )) #Competitive Equlibrium, no profit (value of output [for each sector])=production [of each sector])
@NLconstraint(CGE_EL, ELs_f, sum(Ld_f[i] for i in sectors) + GLd == Le - Unempl) # Competitve Eq: Total L demand = Total L supply
@NLconstraint(CGE_EL, EKs_f, sum(Kd_f[i] for i in sectors) + GKd == Ke) # Competitve Eq: Total K demand = Total K supply
@NLconstraint(CGE_EL, EC[i = sectors], Commodsd_HH[i] + Inv[i] + sum(IOtechCf[i,j] * YOut[j] for j in commods) + GovCd[i]  ==
YCombOut_toHome[i]) #Market clearing for commodities (Sum of consumption, investment [and production uses] equals Total Output)
@NLconstraint(CGE_EL, EInv[i = sectors], Pr_CombCommods_toHome[i] * Inv[i] == BankUexp[i]*TotSav) # Investment Demand = Investment Supply
@NLconstraint(CGE_EL, EPhil, ((w/PrInd)/(wLi/PrIndi)-1) == Phili * ((Unempl/Le)/(Unempli/Lei)-1)) # Wage:unemployment Curve
@NLconstraint(CGE_EL, ETradeBal, sum(Imp[i] * Pr_W_Im_foreignCurri[i] for i in sectors) == 
sum(Pr_W_Exp_foreignCurri[i]*Exp[i] for i in sectors) + ForeignSav) # Trade Balance: Imports = Exports + Foreign Savings (in foreign currency [just cos it's easier not to use exchange rates I guess])

@variable(CGE_EL, Trick, start = 1)

#Utility
@variable(CGE_EL, HHU, start = HHUi)

fix(Ke, Kei, force = true)
fix(Le, Lei, force = true)
fix(TransfOth, TransfOthi, force = true)
fix(GSav, GSavi, force = true)
fix(ForeignSav, ForeignSavi, force= true)
fix(w, wLi, force = true)
fix(Trick, 1, force = true)

@NLobjective(CGE_EL, Min, Trick) #(Not) Numeraire from EcoMod...not clear exactly how to translate
@time optimize!(CGE_EL)

HHU = prod((JuMP.value(Commodsd_HH[i])-HHCsubsist[i])^HHUlesexp[i] for i in sectors) #Final HH Utility calculated AFTER model solved
Walras = JuMP.value(sum(Ld_f[i] for i in sectors) + GLd + Unempl - Le)
GovBudg = JuMP.value(TaxTotRev - Transf - sum(Pr_Commods[i]*GovCd[i] for i in sectors) - w*GLd - r*GKd)

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
print("Walras=           ",Walras,"\n")
print("GovBudg=           ",GovBudg,"\n")