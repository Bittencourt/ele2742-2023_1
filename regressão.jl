using JuMP, GLPK, Ipopt, CSV, DataFrames, Plots

caminho_arq = raw"C:\Users\pedro\Documents\Scanned Documents\Documents\PUC\progLin\ele2742-2023_1\dados_Geração.csv" 

dados = CSV.read(caminho_arq, DataFrame)

dados_Cumaru = dados.Cumaru

scatter(dados_Cumaru[1:end-1], dados_Cumaru[2:end])

#######
# PNL #
#######

regressao_NL = Model(Ipopt.Optimizer)

@variable(regressao_NL, β_0)
@variable(regressao_NL, β_1)

@objective(regressao_NL, Min, sum((dados_Cumaru[t] - β_0 - β_1*dados_Cumaru[t-1])^2 for t in 2:length(dados_Cumaru)))

optimize!(regressao_NL)

β_0_NL = value(β_0)
β_1_NL = value(β_1)

reta_NL = [β_0_NL + β_1_NL*x for x in 1:maximum(dados_Cumaru)]

plot!(reta_NL)

##############
# Quantílica #
##############

regressao_Quantilica = Model(GLPK.Optimizer)

@variable(regressao_Quantilica, β_0)
@variable(regressao_Quantilica, β_1)
@variable(regressao_Quantilica, δ[t=2:744])

@constraint(regressao_Quantilica, restricao_1[t=2:744], δ[t] >= dados_Cumaru[t] - β_0 - β_1*dados_Cumaru[t-1])
@constraint(regressao_Quantilica, restricao_2[t=2:744], δ[t] >= -(dados_Cumaru[t] - β_0 - β_1*dados_Cumaru[t-1]))

@objective(regressao_Quantilica, Min, sum(δ[t] for t in 2:length(dados_Cumaru)))

optimize!(regressao_Quantilica)

β_0_Quantilica = value(β_0)
β_1_Quantilica = value(β_1)

reta_Quantilica = [β_0_Quantilica + β_1_Quantilica*x for x in 1:maximum(dados_Cumaru)]

plot!(reta_Quantilica)