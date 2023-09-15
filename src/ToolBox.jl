module ToolBox

#die ganzen Module, die verwendet werden
using HDF5
using KM3io
using Statistics
using Dates
using CairoMakie
using DataFrames #brauch ich das überhaupt?
using FileIO

#export Data, Data_folder, DomData, DomData_folder, load_Data
#export plot_File, plotDomData1, plotDomDataRings


#ein Artefakt aus alter Zeit, möglichst bald ausbauen
module config 
    Color = [:orange, :green1, :purple4, :navy, :red, :navajowhite3, :blue, :darkred,:darkgreen, :cyan, :maroon, :dimgrey, :purple, :yellow]
    markers = [:circle, :utriangle, :rect, :diamond, :vline, :star5, :hexagon, :cross, :dtriangle, :star4, :xcross, :ltriangle, :pentagon, :star6, :star8, :hline, :rtriangle, 'a', 'B', '↑', '✈']
    Detector_Strings = Int32[5,9,10,11,12,13,14,15,16,19,20,21,22,23,24,25,26,27,28,30,32]
    Detector_Floors = Int32[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]  
    Detector_PMT_Ring = Int32[6,5,5,5,6,6,5,6,6,6,5,5,4,3,2,4,4,3,2,2,3,3,1,4,2,2,2,4,3,3,4]
    Ring_Color = [:orange, :green1, :purple4, :navajowhite3, :navy, :red]
    Detector_PMT_Ringe = [Int32[23], Int32[15,20,26,25,27,19], Int32[14,22,30,29,21,18], Int32[13,16,24,31,28,17], Int32[11,7,4,3,2,12] ,Int32[10,9,5,1,6,8]]
    #dz als abhängigkeit des Rings: dz = [-1,-0.850,-0.555,-0.295,0.295,0.555]
end 

include("smallfunctions.jl")
# include("Data.jl")
include("DomData.jl")
# include("DomDataV1.jl")
# include("DomDataV2.jl")
include("DomDataV3.jl")
include("plot_DomDataV3.jl")

#die ganzen export Functions - kommt das vo oder nach die includes?

end # module ToolBox
