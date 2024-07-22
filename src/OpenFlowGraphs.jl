module OpenFlowGraphs

export FlowGraph, underlying_graph, FG, OpenFG, to_problem,
    node_incidence_matrix, dual_decomposition, nvertices, nedges, random_open_flowgraph, SBM_flowgraph, test,
    add_wire_with_junction!, helper, bipartite_random_graph, SBM_UWD, SBM_flowgraph

using ..FinSetAlgebras
import ..FinSetAlgebras: hom_map, laxator
using ..Objectives
using ..Optimizers
using Catlab
import Catlab: oapply, dom, src, tgt, nvertices
using Test
using Optim
using ForwardDiff
using StatsBase
using Random
using Distributions

struct FlowGraph
    nodes::FinSet
    edges::FinSet
    src::FinFunction
    tgt::FinFunction
    edge_costs::Vector{Function}
    flows::Vector{Float64}
end

function FlowGraph(g::Graph, costs, flows)
    return FlowGraph(FinSet(nv(g)), FinSet(ne(g)), 
        FinFunction(src(g), nv(g)), 
        FinFunction(tgt(g), nv(g)), 
        costs, flows)
end

src(g::FlowGraph, e::Int) = g.src(e)
tgt(g::FlowGraph, e::Int) = g.tgt(e)

dom(g::FlowGraph) = FinSet(g.nodes)
function underlying_graph(g::FlowGraph)
    res = Graph(length(g.nodes))
    add_edges!(res, g.src.func, g.tgt.func)
    return res
end

nvertices(g::FlowGraph) = length(g.nodes)
nedges(g::FlowGraph) = length(g.edges)

struct FG <: FinSetAlgebra{FlowGraph} end

hom_map(::FG, ϕ::FinFunction, g::FlowGraph) =
    FlowGraph(codom(ϕ), g.edges, 
        g.src⋅ϕ, g.tgt⋅ϕ, 
        g.edge_costs, pushforward_matrix(ϕ)*g.flows)

function laxator(::FG, gs::Vector{FlowGraph})
    laxed_src = reduce(⊕, [g.src for g in gs])
    laxed_tgt = reduce(⊕, [g.tgt for g in gs])
    laxed_flows = vcat([g.flows for g in gs]...)
    laxed_costs = vcat([g.edge_costs for g in gs]...)
    return FlowGraph(codom(laxed_src), dom(laxed_src),
        laxed_src, laxed_tgt, laxed_costs, laxed_flows)
end

struct OpenFG <: CospanAlgebra{Open{FlowGraph}} end

function oapply(d::AbstractUWD, gs::Vector{Open{FlowGraph}})
    return oapply(OpenFG(), FG(), d, gs)
end

# Flow graphs to min cost net flow objective functions
function node_incidence_matrix(g::FlowGraph)
    V = nvertices(g)
    E = nedges(g)
    A = zeros(V,E)
    for (v,e) in Iterators.product(1:V, 1:E)
        if src(g, e) == tgt(g, e) && tgt(g, e) == v
            continue
        elseif src(g,e) == v
            A[v,e] = -1
        elseif tgt(g,e) == v
            A[v,e] = 1
        end
    end
    return A
end

function to_problem(og::Open{FlowGraph})
    g = data(og)
    S = og.S
    m = og.m
    A = node_incidence_matrix(g)
    function obj(x,λ)
        return sum([g.edge_costs[i](x[i]) for i in 1:nedges(g)]) + λ'*(A*x-g.flows)
    end

    return Open{SaddleObjective}(S, SaddleObjective(FinSet(nedges(g)), S, obj), m)
end

function dual_decomposition(og::Open{FlowGraph}, γ)
    g = data(og)
    A = node_incidence_matrix(g)
    N = nedges(g)
    
    L(i) = (x,λ) -> g.edge_costs[i](x[1]) + λ'*(A[:,i]*x[1] - 1/N*g.flows)
    L(x,λ) = sum([L(i)(x[i], λ) for i in 1:N])

    function dual_decomp_dynamics(λ)
        x = zeros(N)
        #=Threads.@threads=# for i in 1:N
            x[i] = (optimize(x -> L(i)(x,λ), [0.0], LBFGS(), autodiff=:forward).minimizer)[1]
        end
        return λ + γ*ForwardDiff.gradient(λ->L(x,λ), λ)
    end
    return Open{Optimizer}(og.S, dual_decomp_dynamics, og.m)
end

function random_connected_graph(nv, p)
    g = erdos_renyi(Graph, nv, p)
    while(length(connected_components(g))>1)
        g = erdos_renyi(Graph, nv, p)
    end
    return g
end

g = random_connected_graph(10, .2)
fg = FlowGraph(g, [], [])
@test underlying_graph(fg) == g

function random_quadratic()
    a = rand()
    b = rand()*rand([-1,1])
    c = rand()*rand([-1,1])
    return x -> a*x^2 + b*x + c
    #return x -> x^2
end

function random_flow(n::Int, n_nonzeros::Int)
    u = 2*rand(n_nonzeros) .- 1
    x = [(u[1]-u[n_nonzeros])/2; diff(u) ./ 2]
    res = vcat(x, zeros(n - n_nonzeros))
    return shuffle(res)
end

function random_flow_graph(N::Int, connectivity)
    g = random_connected_graph(N, connectivity)
    E = ne(g)
    flow_costs = [random_quadratic() for i in 1:E]
    #flows = zeros(N)
    flows = random_flow(N, 4)
    return FlowGraph(g, flow_costs, flows)
end

function random_injection(dom::Int, codom::Int)
    f = sample(1:codom, dom, replace=false)
    return FinFunction(sort(f), codom)
end

function random_open_flowgraph(n_vertices, p, n_boundary)
    return Open{FlowGraph}(FinSet(n_vertices), 
        random_flow_graph(n_vertices, p), 
        random_injection(n_boundary, n_vertices))
end

function nports(d, b)
    return length(incident(d, b, :box))
end

function add_wire_with_junction!(d::AbstractUWD, port1::Int, port2::Int)
    j = add_junction!(d)
    set_junction!(d, port1, j)
    set_junction!(d, port2, j)
    
end

function SBM_flowgraph(num_clusters, inter_connectivity, intra_connectivity, nodes_per_cluster)
    d = UndirectedWiringDiagram(0)
    for i = 1:num_clusters
        add_box!(d, nodes_per_cluster)
    end

    for i in 1:nboxes(d)
        for j in i+1:nboxes(d)
            for k in 1:nports(d, i)
                for l in 1:nports(d, j)
                    coin = Bernoulli(inter_connectivity)
                    flip = rand(coin)
                    if flip == true
                        add_wire_with_junction!(d, (i, k) , (j, l))
                    end
                end
            end
        end
    end
    return d
end

function test(num_clusters, inter_connectivity, intra_connectivity, nodes_per_cluster)
    d = UndirectedWiringDiagram(0)

    for i in 1:2
        add_box!(d, nodes_per_cluster)
    end

    for k in ports(d,1)
        for l in ports(d,2)
            coin = Bernoulli(inter_connectivity)
            flip = rand(coin)
            if flip == true #&& !has_wire(d, 1, k) && !has_wire(d, 2, l)
                add_wire_with_junction!(d, k, l)
            end
        end
    end
    return d
end


#=function test(num_clusters, inter_connectivity, intra_connectivity, nodes_per_cluster)
    d = UndirectedWiringDiagram(0)

    for i in 1:2
        add_box!(d, nodes_per_cluster)
    end

    for k in 1:nports(d, 1)
        for l in k+1:nports(d, 2)
            coin = Bernoulli(inter_connectivity)
            flip = rand(coin)
            if flip == true && !has_wire(d, 1, k) && !has_wire(d, 2, l)
                add_wire_with_junction!(d, (1, k) , (2, l))
            end
        end
    end
    return d
end
=#

function helper!(d, box_One, box_Two)
    n = nports(d, box_One)
    m = nports(d, box_Two)
    n_plus_m = add_box!(d, n+m)

    for k in 1:n
       add_wire!(d, (box_One, k), (n_plus_m, k))
    end
    for l in 1:m
       add_wire!(d, (box_Two, l), (n_plus_m, l+n))
    end
end
function bipartite_random_graph(n, m, p)
    g = Graph(n+m)
    for i in 1:n
        for j in n+1:m+n
            c = Bernoulli(p)
            f = rand(c)
            if f == true
                coin = Bernoulli(p)
                flip = rand(coin)
                if flip == true && i != j
                    add_edge!(g, i, j)
                    
                elseif flip == false && i != j
                    add_edge!(g, j, i)
                end
            end     
        end
    end

    E = ne(g)
    
    flow_costs = [random_quadratic() for i in 1:E]
    flows = zeros(nv(g))
    return FlowGraph(g, flow_costs, flows)
end

function random_bipartite_flowgraph(n, m, p)
    return Open{FlowGraph}(FinSet(n+m), 
        bipartite_random_graph(n, m, p), 
        id(FinSet(n+m)))
end

function SBM_UWD(cluster_amount)
    d = UndirectedWiringDiagram(0)
    for i in 1:length(cluster_amount)
        num_ports = cluster_amount[i]
        add_box!(d, num_ports)
    end
    n = nboxes(d)
    for i in 1:n
        for j in i+1:n
            helper!(d, i, j)
        end
    end
    return d
end

function SBM_flowgraph(cluster_amount, connectivities)
    d = SBM_UWD(cluster_amount)
    num_cluster = length(cluster_amount)
    flowgraphs = Open{FlowGraph}[]
    for i in 1:length(cluster_amount)
        g = random_open_flowgraph(cluster_amount[i], connectivities[i, i], cluster_amount[i])
        push!(flowgraphs, g)
    end
    n = num_cluster
    for i in 1:n
        for j in i+1:n
            g = random_bipartite_flowgraph(nports(d, i), nports(d, j), connectivities[i, j])
            push!(flowgraphs, g)
        end
    end
    g = oapply(d, flowgraphs)
    return g, flowgraphs, d
end

#=function helper(num_clusters, inter_connectivity, intra_connectivity, nodes_per_cluster, num_ports_N, num_ports_M)
    d = UndirectedWiringDiagram(0)

    num_ports_N_plus_M = num_ports_N + num_ports_M
    n = add_box!(d, num_ports_N)
    n_plus_m = add_box!(d, num_ports_N_plus_M)
    m = add_box!(d, num_ports_M)
    #only add a wire 3 times
    for k in 1:num_ports_N
        
        add_wire!(d, (n, k), (n_plus_m, k))
    end
    for l in 1:num_ports_M
        add_wire!(d, (m, l), (n_plus_m, l+num_ports_N))
    end
    return d
end=#






end