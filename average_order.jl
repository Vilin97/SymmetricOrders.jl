using Combinatorics

function num_elem_of_cycle_struct(n, cycles)
    cycles_sparse = to_sparse_cycles(cycles)
    denom = prod([factorial(big(cycles_sparse[d])) * d^cycles_sparse[d] for d in 1:n])
    round(BigInt, factorial(big(n))/denom)
end

function to_sparse_cycles(cycles)
    result = zeros(Int, sum(cycles))
    for k in cycles
        result[k] += 1
    end
    return result
end

function remove_trailing_zeros(list)
    non_zero_length = 0
    for (i, elem) in enumerate(list)
        if elem != 0
            non_zero_length = i
        end
    end
    list[1:non_zero_length]
end

function order_distribution(n)
    result = zeros(BigInt, round(BigInt, MathConstants.e^(n/MathConstants.e))+1)
    for partition in partitions(n)
        result[lcm(partition)] += num_elem_of_cycle_struct(n, partition)
    end
    remove_trailing_zeros(result)
end

function avg_order(order_dist)
    sum([order*count for (order, count) in enumerate(order_dist)])/sum(order_dist)
end

dist = order_distribution(52)

function sparse_dist(dist)
    sparse_dist = Dict()
    for (order, count) in enumerate(dist)
        if count != zero(typeof(count))      
            sparse_dist[order] = count
        end
    end
    sparse_dist
end

function cumulative_sums(list)
    cum_sum = copy(list)
    for (k, count) in collect(enumerate(list))[2:end]
        cum_sum[k] = cum_sum[k-1] + count
    end
    cum_sum
end

function cumulative_fractions(list)
    cum_sums = cumulative_sums(list)
    cum_sums ./ cum_sums[end]
end

function get_x_y_from_pairs(pairs)
    return [x for (x, y) in pairs], [y for (x, y) in pairs]
end

using Plots
using PyPlot
pyplot()
p1 = plot(dist[1:1000], xaxis = "order", yaxis = "number of permutations", title = "Permutations in S_52 for order 1 to 1000", label = false); # captures 73% of all permutations in S₅₂
p2 = plot(dist[1000:10000], xaxis = "order", yaxis = "number of permutations", title = "Permutations in S_52 for order 1000 to 10,000", label = false);
p3 = plot(cumulative_fractions(dist), xaxis = "order, log axis", yaxis = "fraction of permutations", title = "Cumulative fraction of permutations", xscale = :log10, xticks = [10^i for i in 1:5], label = false);
descending_counts = reverse(sort(collect(enumerate(dist)), by = x -> x[2])[end-1216:end])
orders, counts = get_x_y_from_pairs(descending_counts)
cutoff = 20
p4 = bar(string.(orders[1:cutoff]), counts[1:cutoff], xaxis = "order", yaxis = "number of permutations", title = "Top 20 orders", xticks = (1:cutoff, string.(orders[1:cutoff])), xrotation = -45, label = false);

plot(p1, p2, p3, p4, layout = 4, size=(1200,800))
savefig("collage_plot")