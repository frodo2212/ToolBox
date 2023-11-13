module PlottingExt 

    using ToolBox
    using HDF5
    using KM3io
    using Dates
    using DataFrames 
    using FileIO
    using CairoMakie 
    using Statistics

    function ylimits!(ax, bound)
        if bound == (nothing,nothing)
            return 
        elseif bound[2] === nothing
            ylims!(ax, low=bound[1])
        elseif bound[1] === nothing
            ylims!(ax, high=bound[2])
        else
            ylims!(ax, bound[1], bound[2])
        end
    end

    include("plot_DomData.jl")


end