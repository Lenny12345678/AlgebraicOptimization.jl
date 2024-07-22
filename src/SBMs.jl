using Catlab
using Distributions
function random_graph(n, p)
    g = Graph(n)
    for i in 1:n
        for j  in i+1:n
            coin = Bernoulli(p)
            flip = rand(coin)
            if flip == true #&& i != j
                add_edge!(g, i, j)
            end
            
            #another way of doing the same thing.
            #if has_edge(g, i, j) == true && has_edge(g, j, i) == true
                #rem_edge!(g, i, j)
                #rem_edge!(g, j, i)
            #end
        end
    end
    return g
end

 g = random_graph(10, 0.4)

to_graphviz(g)
