using JuMP, GLPK, Plots


function genCost(G,d)
    g_min =[0,0,0]
    c = [10,50,100]

    m = Model(GLPK.Optimizer)

    @variable(m, g_min[i] <= g[i=1:length(g_min)] <= G[i])
    @constraint(m, sum(g) == d)

    @objective(m, Min, sum(c'*g))

    optimize!(m)
    return objective_value(m)
end

function minInvest(d)
    G_max = [1000,1000,1000]
    G_min = [0,0,0]
    I = [200,200,50]
    m_i = Model(GLPK.Optimizer)

    @variable(m_i,G_min[i] <= G[i=1:length(G_max)] <= G_max[i])
    @objective(m_i,Min, I'*G + genCost(G,d))

    optimize!(m_i)
end