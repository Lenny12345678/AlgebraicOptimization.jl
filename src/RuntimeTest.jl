using Plots

runtimeList = []
listSize = []

function threeTrials(x)
    s = 0
    for i = 1:3
        list = rand(x)
        sort!(list)
        t = @elapsed sort!(list)
        s = s + t
    end
    push!(runtimeList, s / 3)
end

function test()
    n = 1000
    while n <= 1000000     
        threeTrials(n)
        push!(listSize, n)
        n += 1000
    end
end
test()
scatter(listSize, runtimeList)