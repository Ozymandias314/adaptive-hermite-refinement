using LinearAlgebra, JLD2, Plots

function plotData(filename)
    data = load_object(filename)
    print(data[1:5])
    field_type = split(filename,"_")[1]
    timestamp = split(filename,['_','.'])[2]
    colormap = heatmap(data,c = :thermal)
    plot!(colormap,title = field_type*" at time "*timestamp)

end

filename = "apar_0.jld2"
plotData(filename)