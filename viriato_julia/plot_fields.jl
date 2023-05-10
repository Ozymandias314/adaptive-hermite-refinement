using LinearAlgebra, JLD2, Plots, LaTeXStrings

function plotData(filename)
    data = load_object(filename)
    print(data[1:5])
    field_type = split(filename,"_")[1]
    timestamp = split(filename,['_','.'])[2]
    colormap = heatmap(data,c = :thermal)
    plot!(colormap,title = L"A_∥"* " at time "*timestamp,xaxis=L"x",yaxis=L"y",right_margin=12Plots.mm,aspect_ratio=:equal,grid=false,lims=(0,256))

end

filename = "apar_200.jld2"
plotData(filename)