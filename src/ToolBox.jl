module ToolBox

#all the used Modules
using HDF5
using KM3io
using Statistics
using Dates
using DataFrames #brauch ich das überhaupt?
using FileIO
using CurveFit

export DomData, DomDataV3, DomDataV3_Floors
export Search_DomDataV3, linfit_DomDataV3, linfit_DomDataV3_intervalls, intensiveSearch_DomDataV3
export optical_DomIds, maskTime, PMT_Direction, pos_Strings, Floors
export DateTime_autocorrect, T_intervall, T_intervall2

module config 
    #Arrays with random colors/markers, to better apply a specific color/marker to plots
    Color = [:orange, :green1, :purple4, :navy, :red, :navajowhite3, :blue, :darkred,:darkgreen, :cyan, :maroon, :dimgrey, :purple, :yellow]
    markers = [:circle, :utriangle, :rect, :diamond, :vline, :star5, :hexagon, :cross, :dtriangle, :star4, :xcross, :ltriangle, :pentagon, :star6, :star8, :hline, :rtriangle, 'a', 'B', '↑', '✈']
    Ring_Color = [:orange, :green1, :purple4, :navajowhite3, :navy, :red]
    #The Ring-Positions from the Document of the geometry of the Detector
        #Better implementation via DetectorData possible?
    Detector_PMT_Ring = Int32[6,5,5,5,6,6,5,6,6,6,5,5,4,3,2,4,4,3,2,2,3,3,1,4,2,2,2,4,3,3,4]
    Detector_PMT_Ringe = [Int32[23], Int32[15,20,26,25,27,19], Int32[14,22,30,29,21,18], Int32[13,16,24,31,28,17], Int32[11,7,4,3,2,12] ,Int32[10,9,5,1,6,8]]
end 
PMT_count = 31

 
include("smallfunctions.jl")
include("DomData.jl")
include("DomDataV3.jl")


end 
