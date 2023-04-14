using JuMP, GLPK, Plots,ForwardDiff

domainMax = [2,2,1]
domainMin = [-1/2,-1/2,-2]
problemDimension = length(domainMax)
ϵi = 10^-3

function testFunction(x)
    sum(x[i]^2 for i=1:lastindex(x))
end
function testFunction2(x)
    max.(log.(abs.(x)),x,exp.(x))
end
function diffTestFunction(x)
    ForwardDiff.gradient(testFunction,x)
end
function benders(f,df,ϵ)
    #cria o modelo JuMP
    problem = Model(GLPK.Optimizer)

    #incializa os valores inciais do algoritmo
    maxIter = 1000                #maximo de iterações
    local k = 1                   #contador de iterações
    xc = Array[domainMax]         #chute inicial do valor de entrada
    LB = -Inf                     #lower bound inicial -∞
    UB = f(xc[k])                 #upper bound inicial f(x0) valor da função no x inicial
    local optimal = xc[k]         #variável de saída com o valor de X ótimo, inicializada com X inicial
    
    #define a variavel de otimização X dentro dos limites dados (domainMin e domainMax)
    @variable(problem,domainMin[i] <= x[i=1:problemDimension] <= domainMax[i])

    #define variável δ
    @variable(problem,δ)

    #primeiro corte definido pela restrição δ ≥ f(x0) + ∇f(x0)ᵀ(x-x0)
    @constraint(problem,δ>=testFunction(xc[k])+df(xc[k])'*(x-xc[k]))

    #definição do problema como min δ
    @objective(problem,Min,δ)

    #executa a otimização
    optimize!(problem)

    #atualiza os valores de lower e upper bounds
    LB = value(δ)                       #lower bound recebe o δ resultante da otimização com primeiro corte
    UB = min(UB,f(JuMP.value.(x)))      #upper bound é o menor valor entre o UB anterior e f(xₖ₊₁) 

    #enquanto a diferença entre UB e LB for maior que o ϵ dado itero o método de Benders
    while(UB-LB>ϵ || k<maxIter)
        k = k + 1                       #incrementa o contator de iterações
        insert!(xc,k,JuMP.value.(x))    #inclui o próximo X candidato no vetor de candidatos

        #incluo o novo corte como restrição de δ atualizado para o novo X candidato 
        @constraint(problem,δ>=testFunction(xc[k])+df(xc[k])'*(x-xc[k]))

        #executa a otimização
        optimize!(problem)

        #atualiza os valores de LB e UB
        LB = value(δ)                    #LB recebe o δ resultante depois de inseridos K cortes
        UB = min(UB,f(JuMP.value.(x)))   #UB é atualizado se for melhor candidato
    end
    optimal = JuMP.value.(x)
    println(optimal)
    println("error is ",UB-LB)
    return optimal
end

benders(testFunction,diffTestFunction,ϵi)
