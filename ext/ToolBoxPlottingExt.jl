module ToolBoxPlottingExt 

    using ToolBox
    using HDF5
    using KM3io
    using Dates
    using DataFrames 
    using FileIO
    using CairoMakie 

    #export plot_DomDataV3_Rings, plot_DomDataV3_PMT, plot_DataV3_Events, plot_DataV3_Event, plot_DomDataV3_Floors

    include("plot_DomDataV3.jl")


end