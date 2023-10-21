module PlottingExt3D

    using ToolBox
    using HDF5
    using KM3io
    using Dates
    using DataFrames 
    using FileIO
    using GLMakie
    include("interactivePlots.jl")
    # brauche ich das??
    # GLMakie.activate!()
    """
    plots the 3D konfiguration of a Detector
    """
    function ToolBox.plot_allDoms_3D(detector::Detector; stringnumbers::Bool=false)
        ToolBox.plot_Doms_3D(ToolBox.optical_DomIds(detector), detector, stringnumbers=stringnumbers)
    end
    """
    plots specific Doms of a Detector 
    if plot_allDoms is activated, it shows the whole detector and highlights given Dom Ids
    """
    function ToolBox.plot_Doms_3D(Dom_Ids::Vector{Int32}, detector::Detector; stringnumbers::Bool=false, plot_allDoms::Bool=false, main_Dom::Int32=Int32(0))
        aspect=(1, 1, 1)
        perspectiveness=0.5
        fig = Figure(; resolution=(1200, 400))
        ax1 = Axis3(fig[1, 1]; aspect, perspectiveness)
        Dom_positions = ToolBox.pos_Doms(detector)
        if plot_allDoms
            all_positions = [v for (k,v) in Dom_positions]
            scatter!(ax1, all_positions; markersize=10)
            main_Dom!=0 && (scatter!(Dom_positions[main_Dom]; markersize=18, color=:darkred))
        end
        positions = [Dom_positions[Dom] for Dom in Dom_Ids]
        scatter!(ax1, positions; markersize=15, color=:darkgreen) #meshscatter oder scatter ändern die Texture
        if stringnumbers
            str_pos = ToolBox.pos_Strings(detector)
            Strings = collect(keys(str_pos))
            string_positions = [(str_pos[i][1],str_pos[i][2],0) for i in Strings]
            scatter!(ax1, string_positions)
            for i in (1:length(Strings))
                text!(ax1, string_positions[i], text = string(Strings[i]), offset = (4, 0)) #, align = (:left, :top))
            end
        end
        fig
    end
    """
    takes a DomId, range and detector the plots the whole detector, and marks all doms in the range of given reference Dom
    """
    function ToolBox.plot_closeDoms_3D(Dom_Id::Int32, range::Real, detector::Detector; stringnumbers::Bool=false, plot_allDoms::Bool=true)
        close_Doms = ToolBox.close_Doms(Dom_Id, range, detector)
        return ToolBox.plot_Doms_3D(close_Doms, detector, stringnumbers=stringnumbers, plot_allDoms=plot_allDoms, main_Dom=Dom_Id)
    end
    """
    takes either a Dictionary of Doms => values, or two Vectors of Doms and values 
    it then plots The Doms in a 3D Grid and can incorporate the values either via colors, or sizes (default: colored)
    """
    function ToolBox.plot_DomData_3D(Dom_Ids::Vector{Int32},values::Vector{Float64}, detector::Detector; valuetyp::String="", add_markersizes::Bool=false, colored::Bool=true, size_bounds::Tuple{Real,Real}=(5,25), stringnumbers::Bool=false)
        Dom_positions = ToolBox.pos_Doms(Dom_Ids, detector)
        positions = [v for (k,v) in Dom_positions]
        aspect=(1, 1, 1)
        perspectiveness=0.5
        sizes = values .* ((size_bounds[2]-size_bounds[1])/(maximum(values)-minimum(values)))
        sizes = sizes .+ (size_bounds[1]-minimum(sizes))
        fig = Figure(; resolution=(1200, 400))
        ax1 = Axis3(fig[1, 1]; aspect, perspectiveness, title="position of the Doms inside the Detector")
        if add_markersizes && colored
            scatter!(ax1, positions; markersize=sizes, color = values, colormap = :thermal, colorrange = (minimum(values), maximum(values)))
        elseif colored
            scatter!(ax1, positions; markersize=16, color = values, colormap = :thermal, colorrange = (minimum(values), maximum(values)))
        elseif add_markersizes
            scatter!(ax1, positions; markersize=sizes)
        else
            scatter!(ax1, positions; markersize=16)
        end
        colored && Colorbar(fig[1, 2], limits = (minimum(values), maximum(values)), colormap = :thermal, label=valuetyp)
        if stringnumbers
            str_pos = ToolBox.pos_Strings(detector)
            Strings = collect(keys(str_pos))
            string_positions = [(str_pos[i][1],str_pos[i][2],0) for i in Strings]
            scatter!(ax1, string_positions)
            for i in (1:length(Strings))
                text!(ax1, string_positions[i], text = string(Strings[i]), offset = (4, 0)) #, align = (:left, :top))
            end
        end
        fig
    end
    function ToolBox.plot_DomData_3D(Dom_data::Dict{Int32,Float64}, detector::Detector; valuetyp::String="", add_markersizes::Bool=false, colored::Bool=true, size_bounds::Tuple{Real,Real}=(5,25), stringnumbers::Bool=false)
        Dom_Ids = collect(keys(Dom_data))
        values = [Dom_data[Dom] for Dom in Dom_Ids]
        ToolBox.plot_DomData_3D(Dom_Ids,values, detector, valuetyp=valuetyp, add_markersizes=add_markersizes, colored=colored, size_bounds=size_bounds, stringnumbers=stringnumbers)
    end


end #End of module write nothing after that