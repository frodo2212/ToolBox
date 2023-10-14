module ToolBox

using HDF5
using KM3io
using Statistics
using Dates
using DataFrames
using FileIO
using LsqFit
using Glob

export Data, DomData, DomData_Floors, Search_DomData, linfit_DomData, linfit_DomData_intervalls, intensiveSearch_DomData
export optical_DomIds, PMT_Direction, pos_Strings, Floors, DateTime_autocorrect, T_intervall, T_intervall2
export Event, linfitData

module config 
    Color = [:orange, :green1, :purple4, :navy, :red, :navajowhite3, :blue, :darkred,:darkgreen, :cyan, :maroon, :dimgrey, :purple, :yellow]
    markers = [:circle, :utriangle, :rect, :diamond, :vline, :star5, :hexagon, :cross, :dtriangle, :star4, :xcross, :ltriangle, :pentagon, :star6, :star8, :hline, :rtriangle, 'a', 'B', '↑', '✈']
    Ring_Color = [:orange, :green1, :purple4, :navajowhite3, :navy, :red]
    Detector_PMT_Ring = Int32[6,5,5,5,6,6,5,6,6,6,5,5,4,3,2,4,4,3,2,2,3,3,1,4,2,2,2,4,3,3,4]
    Detector_PMT_Ringe = [Int32[23], Int32[15,20,26,25,27,19], Int32[14,22,30,29,21,18], Int32[13,16,24,31,28,17], Int32[11,7,4,3,2,12] ,Int32[10,9,5,1,6,8]]
end 
PMT_count = 31
 
struct Event
    Typ::String
    Dom_Id::UInt32
    pmt::UInt32
    array_start::UInt32
    array_end::UInt32
    array_length::UInt32
    time_start::UInt32
    time_end::UInt32
    time_length::UInt32 
    missing_timestamps::UInt32
    slice_length::Int32
end

struct linfitData
    Typ::String
    Dom_Id::UInt32
    pmt::UInt32
    time::UInt32
    time_intervall::UInt32 
    params::Tuple{Float64,Float64}
    rel_values::Float64
    slice_length::Int32
end

include("smallfunctions.jl")
include("Data.jl")
include("DomData_search.jl")
include("emptyplotfunctions.jl")


end 
