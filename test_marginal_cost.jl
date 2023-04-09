using Plots

function C(p, c, Gmax)
    return sum(c ./ Gmax) * min.(Gmax, p / (sum(c) * Gmax))
end

# Example problem
p = 1000
c = [10, 15, 20]
Gmax = [300, 400, 500]

# Generate points for plotting
x = 1:1:p
y = [C(i, c, Gmax) for i in x]

# Plot the functions
plot(x, y, label="C(p)")
xlabel!("Total demand (p)")
ylabel!("Generation cost")
title!("Generation cost and its derivative")
