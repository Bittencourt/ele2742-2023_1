using Pkg;

Pkg.add("JuMP")
Pkg.add("GLPK")
Pkg.add("Plots")

using JuMP, GLPK, Plots
# declaração de um modelo

modeloProducao = Model(GLPK.Optimizer)

#variaveis

@variable(modeloProducao, x_1 >= 0)
@variable(modeloProducao, x_2 >= 0)

#restrições

@constraint(modeloProducao, 2*x_1 + x_2 <= 4)
@constraint(modeloProducao, x_1 + 2*x_2 <= 4)

print(modeloProducao)

optimize!(modeloProducao)

##################################
# Despacho economico bem simples #
##################################

despachoEconomico_simples = Model(GLPK.Optimizer)

# variaveis

@variable(despachoEconomico_simples, g_1 >= 0)
@variable(despachoEconomico_simples, g_2 >= 0)

# restricoes

@constraint(despachoEconomico_simples, g_1 + g_2 == 10)

#objetivo

@objective(despachoEconomico_simples, Min, 1*g_1 + 10*g_2)


optimize!(despachoEconomico_simples)
g_1_otimo = value(g_1)
g_2_otimo = value(g_2)

funcao_objetivo = objective_value(despachoEconomico_simples)

##################################
# Despacho economico com restricao da geracao #
##################################

despachoEconomico_rest = Model(GLPK.Optimizer)

# variaveis

@variable(despachoEconomico_rest, g_1 >= 0)
@variable(despachoEconomico_rest, g_2 >= 0)

# restricoes

@constraint(despachoEconomico_rest, g_1 <=5)
@constraint(despachoEconomico_rest, g_1 + g_2 == 10)

#objetivo

@objective(despachoEconomico_rest, Min, 1*g_1 + 10*g_2)


optimize!(despachoEconomico_rest)
g_1_otimo = value(g_1)
g_2_otimo = value(g_2)

funcao_objetivo = objective_value(despachoEconomico_rest)

##################################
# Despacho economico com tres barras #
##################################

despachoEconomico_tres_barras = Model(GLPK.Optimizer)

# variaveis

@variable(despachoEconomico_tres_barras, 0 <= g_1 <= 10)
@variable(despachoEconomico_tres_barras, 0 <= g_2 <= 10)
@variable(despachoEconomico_tres_barras, -5 <= f_1 <= 5)
@variable(despachoEconomico_tres_barras, -10 <= f_2 <= 10)
@variable(despachoEconomico_tres_barras, -3 <= f_3 <= 3)

# restricoes

@constraint(despachoEconomico_tres_barras, g_1 + g_2 == 10)
@constraint(despachoEconomico_tres_barras, g_1 - f_1 - f_3 == 0)
@constraint(despachoEconomico_tres_barras, f_1 + f_2 == 10)
@constraint(despachoEconomico_tres_barras, g_2 + f_3 - f_2 == 0)

#objetivo

@objective(despachoEconomico_tres_barras, Min, 1*g_1 + 10*g_2)


optimize!(despachoEconomico_tres_barras)
g_1_otimo = value(g_1)
g_2_otimo = value(g_2)
f_1_otimo = value(f_1)
f_2_otimo = value(f_2)
f_3_otimo = value(f_3)

funcao_objetivo = objective_value(despachoEconomico_tres_barras)

###################################################
# Despacho economico notacao generica (conjuntos) #
###################################################

# Inputs

B = [1; 2; 3] #conjunto de barras
G = [1; 2] #conjunto de geradores
n_gereadores = length(G)
c = [1; 10] #custos de geração
g_max = [10; 10] #gerações máximas
g_min = [0; 0] #gerações mínimas
G_b = Array[[], [1], [2]] #geradores por barra 
F = [1; 2; 3] #conjunto dos Fluxos
n_fluxos = length(F)
f_max = [5; 10; 3]
f_min = [-5; -10; -3]
F_b_c = Array[[1;2], [],[3]]
F_b_s = Array[[], [1;3],[2]]
d = [10;0;0]

#modelo

despachoEconomico_conjuntos = Model(GLPK.Optimizer)

# variaveis

@variable(despachoEconomico_conjuntos, g_min[i] <= g[i=1:n_gereadores] <= g_max[i])
@variable(despachoEconomico_conjuntos, f_min[i] <= f[i=1:n_fluxos] <= f_max[i])

# restricoes

@constraint(despachoEconomico_conjuntos, balanco[b in B], sum(g[i] for i in G_b[b])+sum(f[i] for i in F_b_c[b])-sum(f[i] for i in F_b_s[b])== d[b])

#objetivo

@objective(despachoEconomico_conjuntos, Min, c'*g)


optimize!(despachoEconomico_conjuntos)

termination_status(despachoEconomico_conjuntos)
g_otimo = value.(g)
f_otimo = value.(f)

funcao_objetivo = objective_value(despachoEconomico_conjuntos)
