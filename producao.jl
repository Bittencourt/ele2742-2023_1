using Pkg

Pkg.add("JuMP")
Pkg.add("GLPK")

using JuMP, GLPK

# modelo 

modeloProducao = Model(GLPK.Optimizer)

# variáveis 

@variable(modeloProducao, x_1 >= 0)
@variable(modeloProducao, x_2 >= 0)

# restrições

@constraint(modeloProducao, 2*x_1 + x_2 <= 4)
@constraint(modeloProducao, x_1 + 2*x_2 <= 4)

print(modeloProducao)

# função objetivo

@objective(modeloProducao, Max, 4*x_1 + 3*x_2)

# otimizar

optimize!(modeloProducao)

# resultados

termination_status(modeloProducao)

x_1_otimo = value(x_1)
x_2_otimo = value(x_2)

funcao_objetivo = objective_value(modeloProducao)

# adicionando restrições

for i in 1:4
    @constraint(modeloProducao, x_1 + x_2 <= i)
    optimize!(modeloProducao)
    println(modeloProducao)
end