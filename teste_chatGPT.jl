using JuMP
using GLPK

function cutting_planes(f, subgradient, X, x0, ε; max_iter=1000)
    k = 0
    x_k = x0
    n = length(x0)

    # Define the linear programming model
    model = Model(GLPK.Optimizer)

    @variable(model, x[1:n])
    @variable(model, z)
    @constraint(model, x .∈ X) # Add constraints to ensure x ∈ X

    while k < max_iter
        # Compute the subgradient of f at x_k
        g_k = subgradient(x_k)

        # Define the cutting plane at x_k
        h_k(y) = f(x_k) + g_k' * (y - x_k)

        # Add the cutting plane constraint to the model
        @constraint(model, h_k(x) <= z)

        # Solve the linear programming problem
        @objective(model, Min, z)
        optimize!(model)

        # Get the new point x_(k+1) from the LP solution
        x_kp1 = value.(x)

        # Check the stopping criterion
        if abs(f(x_kp1) - f(x_k)) <= ε
            break
        end

        x_k = x_kp1
        k += 1
    end

    return x_k, f(x_k)
end

function testF(x)
    x^2
end

function gradTestF(x)
    2*x 
end

cutting_planes(testF,gradTestF,[-10:0.1:10],2.1,0.01)