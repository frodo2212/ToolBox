module ToolBox

using HDF5
using KM3io
using Statistics
using Dates
using DataFrames
using FileIO
using LsqFit
using Glob

export Data, DomData, DomData_Floors
export Search_DomData, intSearch_DomData, linfit_DomData, linfit_DomData_intervalls, Search_linFitData, Search_linFitData_intervalls, DomData_Dom_mean, DomData_PMT_mean
export optical_DomIds, DateTime_autocorrect, T_intervall, T_intervall2, save_Events, load_Events, save_linfitData, load_linfitData
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

plot_DomData_Floors(test::String) = error("Extension not loaded")
plot_DomData_PMT() = error("Extension not loaded")
plot_DomData_Rings() = error("Extension not loaded")
plot_DomData_linFit_array() = error("Extension not loaded")
plot_DomData_linFit() = error("Extension not loaded")
plot_DomData_Event_array() = error("Extension not loaded")
plot_DomData_Event() = error("Extension not loaded")
plot_Strings() = error("Extension not loaded")
plot_Doms() = error("Extension not loaded")
plot_DomData() = error("Extension not loaded")
plot_lf_Intervalls() = error("Extension not loaded")
plot_surrounding_Doms() = error("Extension not loaded")
plot_multiple_Doms() = error("Extension not loaded")
plot_DomData_Drift() = error("Extension not loaded")
plot_Det() = error("Extension not loaded")


plot_Doms_3D() = error("Extension not loaded")
plot_allDoms_3D() = error("Extension not loaded")
plot_DomData_3D() = error("Extension not loaded")
plot_closeDoms_3D() = error("Extension not loaded")
plot_DomData_PMT_interactive() = error("Extension not loaded")
plot_DomData_Rings_interactive() = error("Extension not loaded")
plot_closeDoms_3D_interactive() = error("Extension not loaded")

include("smallfunctions.jl")
include("Data.jl")
include("DomData_search.jl")


end 
