using AlgebraicOptimization
using Catlab
using Test
using Plots

export threeTrials, test_num_nodes, test_num_clusters, test_lots, test_small

dd_runtime = []
dd_memory_usage = []
hdd_runtime = []
hdd_memory_usage = []

function threeTrials(cluster_amount, connectivities, ss, iters)
    dd_time = 0
    dd_mem = 0
    hdd_time = 0
    hdd_mem = 0

    g, gs, d = SBM_flowgraph(cluster_amount, connectivities)
    dd_optimizer = dual_decomposition(g, ss)
    os = [dual_decomposition(g, ss) for g in gs]
    hdd_optimizer = oapply(OpenDiscreteOpt(), d, os)
    λ0 = zeros(nvertices(data(g)))
    res1 = simulate(dd_optimizer, λ0, iters)
    res2 = simulate(hdd_optimizer, λ0,iters)

    for _ in 1:3
        #@test res1 ≈ res2
        λ0 = zeros(nvertices(data(g)))
        t1 = @timed simulate(dd_optimizer, λ0, iters)
        t2 = @timed simulate(hdd_optimizer, λ0, iters)
        dd_time = dd_time + t1.time
        dd_mem = dd_mem + t1.bytes
        hdd_time = hdd_time + t2.time
        hdd_mem = hdd_mem + t2.bytes
    end
    push!(dd_runtime, dd_time / 3)
    push!(dd_memory_usage, dd_mem / 3)
    push!(hdd_runtime, hdd_time / 3)
    push!(hdd_memory_usage, hdd_mem / 3)
end

#fixed number of clusters and connectivity and increase number of nodes per cluster
function test_num_nodes(cluster_amount, connectivities, ss, iters)
    n = 10
    nodes = []
    for i in 1:n
        threeTrials(cluster_amount, connectivities, ss, iters)
        cluster_amount = [cluster_amount[1] + 5, cluster_amount[2] + 5, cluster_amount[3] + 5]
        push!(nodes, cluster_amount[1])
    end
    plt1 = plot_runtime(nodes, dd_runtime, nodes, hdd_runtime, "Nodes per cluster")
    plt2 = plot_memory_usage(nodes, dd_memory_usage, nodes, hdd_memory_usage, "Nodes per cluster")
    return plt1, plt2
end

#fixed number of nodes and connectivity, increase number of clusters
function test_num_clusters(cluster_amount, connectivities, ss, iters)
    n = 20
    clusters = []
    for i in 1:n
        threeTrials(cluster_amount, connectivities, ss, iters)
        push!(cluster_amount, cluster_amount[i])
        connectivities = repeat([connectivities[i]], length(cluster_amount),length(cluster_amount))
        push!(clusters, length(cluster_amount)-1)
    end
    plt1 = plot_runtime(clusters, dd_runtime, clusters, hdd_runtime, "Number of clusters")
    plt2 = plot_memory_usage(clusters, dd_memory_usage, clusters, hdd_memory_usage, "Number of clusters")
    return plt1, plt2
end

#vary connectivities and lots of cluster and nodes
function test_lots()
    cluster_amount = []
    num_clusters = 10
    connectivities = [0.1]
    connectivity = []
    ss = 0.1
    iters = 10
    n = 15

    for i in 1:num_clusters
        push!(cluster_amount, 15)
    end
    connectivities = repeat([connectivities[1]], length(cluster_amount),length(cluster_amount))
    for j in 1:n
        threeTrials(cluster_amount, connectivities, ss, iters)
        push!(connectivity, connectivities[1])
        connectivities = repeat([connectivities[1,1] + 0.05], length(cluster_amount),length(cluster_amount))
    end
    plt1 = plot_runtime(connectivity, dd_runtime, connectivity, hdd_runtime, "Connectivity percentage")
    plt2 = plot_memory_usage(connectivity, dd_memory_usage, connectivity, hdd_memory_usage, "Connectivity percentage")
    return plt1, plt2
end

#vary connectivities and small amount of clusters and nodes
function test_small()
    cluster_amount = []
    num_clusters = 3
    connectivities = [0.1]
    connectivity = []
    ss = 0.1
    iters = 10
    n = 15

    for i in 1:num_clusters
        push!(cluster_amount, 10)
    end
    connectivities = repeat([connectivities[1]], length(cluster_amount),length(cluster_amount))
    for j in 1:n
        threeTrials(cluster_amount, connectivities, ss, iters)
        push!(connectivity, connectivities[1])
        connectivities = repeat([connectivities[1,1] + 0.05], length(cluster_amount),length(cluster_amount))
    end
    plt1 = plot_runtime(connectivity, dd_runtime, connectivity, hdd_runtime, "Connectivity percentage")
    plt2 = plot_memory_usage(connectivity, dd_memory_usage, connectivity, hdd_memory_usage, "Connectivity percentage")
    return plt1, plt2
end

function plot_runtime(dd_x, dd_y, hdd_x, hdd_y, x_label)
    f = "Computer Modern"
    p1 = plot(dd_x, dd_y, label="Standard DD",
    size = (1200,800),
    marker=:circle,
    title=x_label * " vs. Runtime",
    titlefont = (14,f),
    linewidth = 2,
    xlabel=x_label,
    ylabel ="Execution time (s)",
    thickness_scaling = 2,
    tickfont = (10,f),
    legend = :topleft,
    legend_font_family = f,
    #smooth = true,
    legendfontsize=10,
    seriescolor=palette(:default)[1],
    ms=4,
    guidefont=(f,12)
)
    plot!(hdd_x, hdd_y, label="Hierarchical DD", marker=:square, linewidth=2, seriescolor=palette(:default)[2])
end

function plot_memory_usage(dd_x, dd_y, hdd_x, hdd_y, x_label)
    f = "Computer Modern"
    p1 = plot(dd_x, dd_y, label="Standard DD",
    size = (1200,800),
    marker=:circle,
    title=x_label * " vs. Memory Usage",
    titlefont = (14,f),
    linewidth = 2,
    xlabel=x_label,
    ylabel ="Memory Usage (bytes)",
    thickness_scaling = 2,
    tickfont = (10,f),
    legend = :topleft,
    legend_font_family = f,
    #smooth = true,
    legendfontsize=10,
    seriescolor=palette(:default)[1],
    ms=4,
    guidefont=(f,12)
)
    plot!(hdd_x, hdd_y, label="Hierarchical DD", marker=:square, linewidth=2, seriescolor=palette(:default)[2])
end

#p1, p2 = test_num_nodes([10,10,10], [0.3 0.2 0.2; 0.2 0.3 0.2; 0.2 0.2 0.3], 0.1, 10)
#p1, p2 = test_num_clusters([10,10,10], [0.2 0.2 0.2; 0.2 0.2 0.2; 0.2 0.2 0.2], 0.1, 10)
#p11, p12 = test_lots()
p21, p22 = test_small()