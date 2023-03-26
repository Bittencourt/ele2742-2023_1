using JuMP, GLPK, Plots

##################################
# Despacho economico bem simples #
##################################

despachoEconomico_Simples = Model(GLPK.Optimizer)

# variáveis

@variable(despachoEconomico_Simples, g_1 >= 0)
@variable(despachoEconomico_Simples, g_2 >= 0)

# restrições

@constraint(despachoEconomico_Simples, g_1 + g_2 == 10)

# objetivo

@objective(despachoEconomico_Simples, Min, 1*g_1 + 10*g_2)

# minimiza

optimize!(despachoEconomico_Simples)

g_1_otimo = value(g_1)
g_2_otimo = value(g_2)

funcao_objetivo = objective_value(despachoEconomico_Simples)

###############################################
# Despacho economico com restrição na geração #
###############################################

despachoEconomico_RestricaoGeracao = Model(GLPK.Optimizer)

# variáveis

@variable(despachoEconomico_RestricaoGeracao, 0 <= g_1 <= 5)
@variable(despachoEconomico_RestricaoGeracao, 0 <= g_2 <= 10)

# restrições

@constraint(despachoEconomico_RestricaoGeracao, g_1 + g_2 == 10)

# objetivo

@objective(despachoEconomico_RestricaoGeracao, Min, 1*g_1 + 10*g_2)

# minimiza

optimize!(despachoEconomico_RestricaoGeracao)

g_1_otimo = value(g_1)
g_2_otimo = value(g_2)

funcao_objetivo = objective_value(despachoEconomico_RestricaoGeracao)

###################################
# Despacho economico com 3 barras #
###################################

despachoEconomico_3barras = Model(GLPK.Optimizer)

# variáveis

@variable(despachoEconomico_3barras, 0 <= g_1 <= 10)
@variable(despachoEconomico_3barras, 0 <= g_2 <= 10)
@variable(despachoEconomico_3barras, -5 <= f_1 <= 5)
@variable(despachoEconomico_3barras, -10 <= f_2 <= 10)
@variable(despachoEconomico_3barras, -3 <= f_3 <= 3)

# restrições

@constraint(despachoEconomico_3barras, f_1 + f_2 == 10)
@constraint(despachoEconomico_3barras, g_1 - f_1 - f_3 == 0)
@constraint(despachoEconomico_3barras, g_2 - f_2 + f_3 == 0)

# objetivo

@objective(despachoEconomico_3barras, Min, 1*g_1 + 10*g_2)

# minimiza

optimize!(despachoEconomico_3barras)

g_1_otimo = value(g_1)
g_2_otimo = value(g_2)

f_1_otimo = value(f_1)
f_2_otimo = value(f_2)
f_3_otimo = value(f_3)

funcao_objetivo = objective_value(despachoEconomico_3barras)

####################################
# Despacho economico via conjuntos #
####################################

# inputs

B = [1; 2; 3]

d = [10; 15; 0]

G = [1; 2; 3]
n_geradores = length(G)
c = [1; 10; 5]
g_max = [10; 10; 10]
g_min = [0; 0; 0]
G_b = Array[[3], [1], [2]]

F = [1; 2; 3]
n_fluxos = length(F)
f_max = [5; 10; 3]
f_min = [-5; -10; -3]
F_b_c = Array[[1; 2], [], [3]]
F_b_s = Array[[], [1; 3], [2]]

# Modelo

despachoEconomico_Conjuntos = Model(GLPK.Optimizer)

# Variáveis

@variable(despachoEconomico_Conjuntos, g_min[i] <= g[i=1:n_geradores] <= g_max[i])
@variable(despachoEconomico_Conjuntos, f_min[i] <= f[i=1:n_fluxos] <= f_max[i])

# restrições

@constraint(despachoEconomico_Conjuntos, balanco[b in B], sum(g[i] for i in G_b[b]) + sum(f[i] for i in F_b_c[b]) - sum(f[i] for i in F_b_s[b]) == d[b])

# objetivo

@objective(despachoEconomico_Conjuntos, Min, c'*g)

# otimiza

optimize!(despachoEconomico_Conjuntos)

# resultados

g_otimo = value.(g)

f_otimo = value.(f)

funcao_objetivo = objective_value(despachoEconomico_Conjuntos)


#########################################################
# Despacho economico via conjuntos com demanda variando #
#########################################################

# inputs

B = [1; 2; 3]

d = Array[rand(24).*10, [4*cos(2*pi*t/24) + 10 for t in 1:24], zeros(24)]

plot([5*sin(2*pi*t/24) + 10 for t in 1:24])

G = [1; 2; 3]
n_geradores = length(G)
c = [1; 10; 5]
g_max = [10; 10; 10]
g_min = [0; 0; 0]
G_b = Array[[3], [1], [2]]

F = [1; 2; 3]
n_fluxos = length(F)
f_max = [5; 10; 3]
f_min = [-5; -10; -3]
F_b_c = Array[[1; 2], [], [3]]
F_b_s = Array[[], [1; 3], [2]]

# Modelo

despachoEconomico_t = Model(GLPK.Optimizer)

# Variáveis

@variable(despachoEconomico_t, g_min[i] <= g[i=1:n_geradores, t=1:24] <= g_max[i])
@variable(despachoEconomico_t, f_min[i] <= f[i=1:n_fluxos, t=1:24] <= f_max[i])

# restrições

@constraint(despachoEconomico_t, balanco[b in B, t in 1:24], sum(g[i,t] for i in G_b[b]) + sum(f[i,t] for i in F_b_c[b]) - sum(f[i,t] for i in F_b_s[b]) == d[b][t])

# objetivo

@objective(despachoEconomico_t, Min, sum(c'*g[:,t] for t in 1:24))

# otimiza

optimize!(despachoEconomico_t)

# resultados

termination_status(despachoEconomico_t)

g_otimo = value.(g)

f_otimo = value.(f)

funcao_objetivo = objective_value(despachoEconomico_t)

plot(d[1][:], label = "Demanda barra 1", ylabel = "MWh", xlabel = "Hora do dia")
plot!(d[2][:], label = "Demanda barra 2", ylabel = "MWh", xlabel = "Hora do dia")

plot!(g_otimo[1,:], label = "Gerador 1")
plot!(g_otimo[2,:], label = "Gerador 2")
plot!(g_otimo[3,:], label = "Gerador 3")