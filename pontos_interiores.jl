# Lista 3 - Pedro Bittencourt
# Pontos Interiores
#

using Pkg
Pkg.add("MathOptInterface")
Pkg.add("NLPModelsJuMP")

import MathOptInterface as MOI
using NLPModelsJuMP
using JuMP, GLPK, Cbc, SparseArrays
using LinearAlgebra, Plots

function PontosInteriores(A, b, c, ρ, α, ϵ, maxIter, show_iter = true)
    # A: matriz de coeficientes
    # b: vetor de termos independentes
    # c: vetor de custos
    # ρ: parâmetro de redução de gap
    # ϵ: tolerância
    # maxIter: número máximo de iterações
    
    # Dimensões
    m, n = size(A)

    A_ = sparse(A)
    
    # Inicialização
    dim_x = 1:n
    dim_y = n+1:n+m
    dim_s = n+m+1:2*n+m

    x = ones(n)  
    y = ones(m)    
    s = ones(n)  

    iter = 1
    gap = x' * s
    μ = gap / n
    
    # Vetores de armazenamento
    zk = [[x; y; s]]
    Fk = []
    
    converge = 0

    # Iteração
    while converge == 0
        F = [A_*x - b ; A_'*y + s - c ; diagm(x)*diagm(s)*ones(n) - μ*ones(n)]
        push!(Fk, F)
        # Teste de convergência
        if maximum(abs.(F[1:n+m])) < ϵ && gap < ϵ
            converge = 1
            if( show_iter == true)
                println("Iteração ", iter)
                println("----------------")
                println("x = ", x)
                println("y = ", y)
                println("s = ", s)
                println()
                println("f(x) =", c'*x)
            end
            println("Encontrada solução ótima")
            break
        end

        # Teste de número máximo de iterações
        if iter >= maxIter
            converge = 2
            println("Número máximo de iterações atingido")
            break
        end
        
        # Passo de Newton
        J = spzeros(2*n+m, 2*n+m)
        J[1:m, 1:n] = A_
        J[m+1:m+n, n+1:n+m] = A_'
        J[m+1:m+n, n+m+1:m+2*n] = sparse(I, n, n)
        J[n+m+1:m+2*n, 1:n] = sparse(diagm(s))
        J[n+m+1:m+2*n, n+m+1:m+2*n] = sparse(diagm(x))
        
        # Resolução do sistema linear
        #d = try
          d = J\-F
        # catch e
        #   converge = 3
        #   println("Problema ilimtado")
        #   break;
        # end
        #println("d = ", d)

        # Atualização
        dx = d[1:n]
        dy = d[n+1:n+m]
        ds = d[n+m+1:end]

        B_p = minimum([1 ; α.*(-x[dx .< 0] ./ dx[dx .< 0])])
        B_d = minimum([1 ; α.*(-s[ds .< 0] ./ ds[ds .< 0])])
        #println("B_p = ", B_p, " B_d = ", B_d)

        x += B_p*dx
        y += B_d*dy
        s += B_d*ds

         if(maximum(s) > 1e15 )
             println("Problema inviável")
             converge = 3
             break;
         end
        if(any(x->x>1e9,x))
            println("Problema ilimitado")
            converge = 3
            break;
        end

        push!(zk, [x; y; s])
        

        gap = x' * s
        if( show_iter == true)
            println("Iteração ", iter)
            println("----------------")
            println("x = ", x)
            println("y = ", y)
            println("s = ", s)
            println()
            println("f(x) =", c'*x)
            println()
            println("gap = ", gap)
            println("----------------")
            println()
        end
        μ = ρ*gap / n
        iter += 1
    end
    return x, y, s, iter, gap, converge, zk, Fk
end

########################################
# problema viável (problema da produção)
########################################

A = [-2 -1 -1 0 ; -1 -2 0 -1]
b = [-4 ; -4]
c = [-4 ; -3 ; 0; 0]
α = 0.95
ρ = 0.1
ϵ = 1e-6

x, y, s, iter, gap, converge, zk, Fk = PontosInteriores(A, b, c, ρ, α, ϵ, 1000)

m, n = size(A)
if converge == 1
    y_vet1 = [zk[i][n+1] for i in 1:iter]
    y_vet2 = [zk[i][n+2] for i in 1:iter]
    y_vet3 = [zk[i][n+3] for i in 1:iter]
    x_vet1 = [zk[i][1] for i in 1:iter]
    x_vet2 = [zk[i][2] for i in 1:iter]
    x_vet3 = [zk[i][3] for i in 1:iter]
    f_vet = [c'*zk[i][1:n] for i in 1:iter]
    # Plot da trajetória de x e y e do gap primal dual
    gap_k = [zk[i][1:n]'*zk[i][n+m+1:end] for i in 1:iter]
    plot(-f_vet, markersize = 5, markershape = :hexagon, label = "Obj.")
    plot(y_vet1, y_vet2, y_vet3, markersize = 5, markershape = :hexagon, label = "y")
    annotate!(y_vet1, y_vet2, y_vet3, text.(1:iter, :bottom, :left))
    plot(x_vet1, x_vet2,x_vet3, markersize = 5, markershape = :hexagon, label = "x")
    annotate!(x_vet1, x_vet2,x_vet3, text.(1:iter, :bottom, :left))
    plot(log10.(gap_k), markersize = 5, markershape = :hexagon, label = "log gap")
end

####################
# problema ilimitado
####################
A = -[-1 1 -1 1 0 0; 2 -1 0 0 1 0; 0 -1 2 0 0 1]
b = [-5 ; -3; -5]
c = -[2 ; 0 ; 1 ; 0; 0 ;0]
α = 0.95
ρ = 0.1
ϵ = 1e-6 

x, y, s, iter, gap, converge, zk, Fk = PontosInteriores(A, b, c, ρ, α, ϵ, 1000)

###################
# problema inviável
###################
A = -[2 -1 -2 1 0 0;2 -3 -1 0 1 0;-1 1 1 0 0 1]
b = -[4; -5; -1]
c = -[1; -1; 1; 0; 0; 0]
α = 0.95
ρ = 0.1
ϵ = 1e-6

x, y, s, iter, gap, converge, zk, Fk = PontosInteriores(A, b, c, ρ, α, ϵ, 1000)


###########################
#problema da produção maior
###########################
N = [2;5;10;50;100;200] # número de produtos

results_time_a = zeros(6,20)
results_time_b = zeros(6,20)
results_obj_a = zeros(6,20)
results_obj_b = zeros(6,20)
results_iter_a = zeros(6,20)
results_iter_b = zeros(6,20)
index = 0
for n in N
    index += 1
    for k in 1:20
        A = -[rand((1:10),2,n) I] # quanto de cada insumo para cada produto rand((1:4),2,n) !!não pode ser zero!! + variavies Slack (I)
        b = -[400;400]
        c = [-rand((0:10),n);zeros(2)]
        α = 0.95
        ρ = 0.1
        ϵ = 1e-4
        #results[index,k] = @elapsed PontosInteriores(A, b, c, ρ, α, ϵ, 1000, false)
        x, y, s, iter, gap, converge, zk, Fk = PontosInteriores(A, b, c, ρ, α, ϵ, 1000,false)
        results_obj_a[index,k] = c'*x
        results_iter_a[index,k] = iter
        model = Model(GLPK.Optimizer)
        @variable(model, x[1:n+2] >= 0)
        @objective(model, Min, c'*x)
        @constraint(model, A*x .== b)
        optimize!(model)
        #println("n = ", n, " iter = ", iter, " converge = ", converge, " f(x) = ", c'*x)
    end
end


####################
#problema da bateria
####################

model = Model(GLPK.Optimizer)

T = 2

master_g = 2.5*[-0.026;-0.026;-0.026;-0.026;-0.026;-0.026;1.293;22.447;56.507;86.463;87.749;115.92;98.812;115.92;86.249;44.803;9.683;-0.187;-0.029;-0.029;-0.027;-0.028;-0.027;-0.027]
g = master_g[1:T]
# demand from a real scenario in kWh
master_d = [188.475;179.527;203.068;198.753;166.385;190.016;225.074;229.42;212.373;233.244;243.932;253.594;237.235;234.595;236.405;274.686;290.507;304.974;287.483;290.268; 242.027; 252.202; 227.184; 236.652]
d = master_d[1:T]
# dynamic princing in R$/kWh (arbitrary values)
master_c = [1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1.8; 1.8; 1.8; 1.8; 1; 1; 1; 1]
c = master_c[1:T]
#g = zeros(24)
Max_Discharge = -100
Max_Charge = 100
Max_Capacity = 550

# Define variables
@variable(model, F[t=1:T])
@variable(model, Max_Discharge <= b[t=1:T] <= Max_Charge)

# Define constraints
@constraint(model, balaco[t in 1:T], F[t] == d[t] - g[t] + b[t])
@constraint(model, baterycharge[t in 1:T], 0 <= sum(b[1:t]) <= Max_Capacity)

# Define objective function
@objective(model, Min, sum(c[t]*F[t] for t in 1:T))

# Solve model
optimize!(model)

# Print results
println("Objective function value: ", objective_value(model))


#plot B

plot(d, label="demand")
plot!(value.(g), label="generation")
plot!(value.(F), label="flux")
plot!(value.(b), label="battery")
plot!(title="Battery optimization", xlabel="hour", ylabel="MWh")


termination_status(model)

function extract_matrices(model::Model)
    num_vars = length(all_variables(model))
    constraints = all_constraints(model; include_variable_in_set_constraints = false)
    quad_constraints = all_constraints(model, JuMP.GenericQuadExpr{Float64, JuMP.VariableRef}, MOI.LessThan{Float64})
    num_cons = length(constraints)
    num_quad_cons = length(quad_constraints)

    A = spzeros(num_cons, num_vars)
    B = spzeros(num_cons)
    C = spzeros(num_vars)
    Q = spzeros(num_vars, num_vars)
    Q_constraints = Vector{SparseMatrixCSC{Float64, Int}}()
    b_constraints = spzeros(num_quad_cons)

    # Extract the constraint matrix (A) and the right-hand side vector (B)
    for (i, con) in enumerate(constraints)
        constraint_function = JuMP.constraint_object(con).func
        for var in keys(constraint_function.terms)
            A[i, JuMP.index(var).value] = constraint_function.terms[var]
        end
        try
            B[i] = JuMP.constraint_object(con).set.upper
        catch
            B[i] = JuMP.constraint_object(con).set.value
        end
    end

    # Extract quadratic constraints
    for (i, con) in enumerate(quad_constraints)
        quad_constraint_function = JuMP.constraint_object(con).func
        quad_constraint = spzeros(num_vars, num_vars)
        b_constraints[i] = JuMP.constraint_object(con).set.upper - quad_constraint_function.aff.constant
        for (quad_var, coeff) in quad_constraint_function.terms
            quad_constraint[index(getfield(quad_var, 1)).value, index(getfield(quad_var, 2)).value] = coeff
        end
        push!(Q_constraints, quad_constraint)
    end

    # Extract the cost vector (C) and the quadratic matrix (Q)
    objective_function = JuMP.objective_function(model)
    if isa(objective_function, JuMP.QuadExpr)
        for var in keys(objective_function.aff.terms)
            C[index(var).value] = objective_function.aff.terms[var]
        end
        for (quad_var, coeff) in objective_function.terms
            Q[index(getfield(quad_var, 1)).value, index(getfield(quad_var, 2)).value] = coeff
            Q[index(getfield(quad_var, 2)).value, index(getfield(quad_var, 1)).value] = coeff
        end
    else
        for var in keys(objective_function.terms)
            C[JuMP.index(var).value] = objective_function.terms[var]
        end
    end

    return A, B, C, Q, Q_constraints, b_constraints
end

# Extract standard form
nlp = MathOptNLPModel(model) # the input "< model >" is the name of the model you created by JuMP before with variables and constraints (and optionally the objective function) attached to it.
x = zeros(nlp.meta.nvar)
c = NLPModelsJuMP.grad(nlp, x)
A = Matrix(NLPModelsJuMP.jac(nlp, x))

lb = nlp.meta.lcon
ub = nlp.meta.ucon
m = size(A)[1]
b = deepcopy(lb)
for i = 1:m
    v_aux = zeros(m)
    l = lb[i]
    u = ub[i]
    if l==u && !isinf(l)
        continue
    end
    if isinf(l)
        b[i] = ub[i]
        v_aux[i] = 1
        A = hcat(A, v_aux)
    end
    if isinf(u)
        v_aux[i] = -1
        A = hcat(A, v_aux)
    end
    display(A)
end
#c = vcat(c, zeros(size(A)[2]-n))

A,B,C = extract_matrices(model)
x, y, s, iter, gap, converge, zk, Fk = PontosInteriores(A, collect(B), collect(C), ρ, α, ϵ, 1000)
