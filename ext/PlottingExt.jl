module PlottingExt 

    using ToolBox
    using HDF5
    using KM3io
    using Dates
    using DataFrames 
    using FileIO
    using CairoMakie 

    #export plot_DomDataV3_Rings, plot_DomDataV3_PMT, plot_DataV3_Events, plot_DataV3_Event, plot_DomDataV3_Floors

    #funktionen ohne parameter gehen nicht, da er dann die aus den emptyplotfunctions nimmt
    #bzw gehen schon, wenn man in emptyplotfunctions halt andere parameter angibt
    #TODO: Tamas fragen wie das richtig gehen sollte, das kanns so ja eigentlichnicht sein
    function ToolBox.test()
        println("Juhuu es funktioniert")
    end

    include("plot_DomData_alt.jl")


end