module PlottingExt3D

    using ToolBox
    using HDF5
    using KM3io
    using Dates
    using DataFrames 
    using FileIO
    using GLMakie

    # brauche ich das??
    # GLMakie.activate!()

    function ToolBox.plot_Doms_3D(detector::Detector)
        ToolBox.plot_Doms_3D(ToolBox.optical_DomIds(detector), detector)
    end
    function ToolBox.plot_Doms_3D(Dom_Ids::Vector{Int32}, detector::Detector)
        Dom_positions = ToolBox.pos_Doms(Dom_Ids, detector)
        positions = [v for (k,v) in Dom_positions]
        aspect=(1, 1, 1)
        perspectiveness=0.5
        fig = Figure(; resolution=(1200, 400))
        ax1 = Axis3(fig[1, 1]; aspect, perspectiveness)
        scatter!(ax1, positions; markersize=15) #meshscatter oder scatter ändern die Texture
        fig
    end
    function ToolBox.plot_Doms_3D(Dom_Ids::Vector{Int32},values::Vector{Float64}, detector::Detector, valuetyp::String="")
        Dom_positions = ToolBox.pos_Doms(Dom_Ids, detector)
        positions = [v for (k,v) in Dom_positions]
        aspect=(1, 1, 1)
        perspectiveness=0.5
        fig = Figure(; resolution=(1200, 400))
        ax1 = Axis3(fig[1, 1]; aspect, perspectiveness, title="position of the Doms inside the Detector")
        scatter!(ax1, positions; markersize=15, color = values, colormap = :thermal, colorrange = (minimum(values), maximum(values)))
        Colorbar(fig[1, 2], limits = (minimum(values), maximum(values)), colormap = :thermal, label=valuetyp)
        fig
    end
end

function ToolBox.plot_DomData_Rings(DomID::Integer; loadpath::String="../Data/DomData_Doms", alpha::Float64=1.0, slice_length::Integer=6000)
    fig = Figure() # resolution = (3840, 2160))
    ax1 = Axis(fig[1,1],title ="DomData Rings")
    ringmenu = Menu(fig[1,2], options = [1, 2, 3, 4, 5, 6], fontsize = 30)
    button = Button(fig, label = "refresh")
    fig[1, 2] = vgrid!(
        Label(fig, "Ring:", fontsize = 30, width = 150), ringmenu,
        button, tellheight = false, width = 200)
    file = h5open(string(loadpath, "/Dom_", DomID,"_",Int32(slice_length/600),".h5"), "r")
    Times = read(file["Time"])
    pmt_values = Vector{Matrix{Float64}}(undef, 6)
    Zeiten, Typ = ToolBox.autoscale_time(minimum(Times), maximum(Int64, Times), intervalls=3)
    ax1.xticks[] = (datetime2unix.(Zeiten.dates) , Dates.format.(Zeiten.dates, Typ)) 
    for ring in (1:length(ToolBox.config.Detector_PMT_Ringe))
        pmt_values[ring] = Matrix{Float64}(undef, length(ToolBox.config.Detector_PMT_Ringe[ring]), length(Times))
        for pmt in (1:length(ToolBox.config.Detector_PMT_Ringe[ring]))
            pmt_values[ring][pmt,:] = read(file["pmtmean"])[:,ToolBox.config.Detector_PMT_Ringe[ring][pmt]]
        end 
    end
    close(file)
    #provisional_slider
    lsgrid = SliderGrid(fig[2,:],
    (label = "start", range = minimum(Times):1:maximum(Times), startvalue = minimum(Times)),
    (label = "end", range = minimum(Times):1:maximum(Times), startvalue = maximum(Times))
    # formats = [x -> "$(round(x))"],
    )
    Start = lsgrid.sliders[1].value
    End = lsgrid.sliders[2].value
    time_mask = @lift ToolBox.maskTime(Times,($Start,$End))
    obs = Observable(1)
    for i in (1:length(ToolBox.config.Detector_PMT_Ringe[obs[]]))
        scatter!(ax1, Times[time_mask[]], pmt_values[obs[]][i,time_mask[]], color=ToolBox.config.Color[i], alpha=alpha)
    end
    on(obs) do obs
        empty!(ax1)
        for i in (1:length(ToolBox.config.Detector_PMT_Ringe[obs]))
            scatter!(ax1, Times[time_mask[]], pmt_values[obs][i,time_mask[]], color=ToolBox.config.Color[i], alpha=alpha)
        end
    end
    on(button.clicks) do click
        empty!(ax1)
        for i in (1:length(ToolBox.config.Detector_PMT_Ringe[obs[]]))
            scatter!(ax1, Times[time_mask[]], pmt_values[obs[]][i,time_mask[]], color=ToolBox.config.Color[i], alpha=alpha)
        end
    end
    on(ringmenu.selection) do select
        obs[] = select
    end
    return fig
end

function ToolBox.plot_DomData_PMT(DomID::Integer; loadpath::String="../Data/DomData_Doms", alpha::Float64=1.0, slice_length::Integer=6000)
    fig = Figure() # resolution = (3840, 2160))
    ax1 = Axis(fig[1,1],title ="DomData PMTs")
    fig[2, 1] = buttongrid = GridLayout(tellwidth = false)
    buttonlabels = [string(i) for i in (1:ToolBox.PMT_count)]
    buttons = buttongrid[1, 1:11] = [Button(fig, label = l) for l in buttonlabels[1:11]]
    buttons2 = buttongrid[2, 1:11] = [Button(fig, label = l) for l in buttonlabels[12:22]]
    buttons3 = buttongrid[3, 1:9] = [Button(fig, label = l) for l in buttonlabels[23:31]]
    buttons = append!(buttons, buttons2)
    buttons = append!(buttons, buttons3)
    file = h5open(string(loadpath, "/Dom_", DomID,"_",Int32(slice_length/600),".h5"), "r")
    Times = read(file["Time"])
    pmt_values = Vector{Vector{Float64}}(undef, ToolBox.PMT_count)
    Zeiten, Typ = ToolBox.autoscale_time(minimum(Times), maximum(Int64, Times), intervalls=3)
    ax1.xticks[] = (datetime2unix.(Zeiten.dates) , Dates.format.(Zeiten.dates, Typ)) 
    for pmt in (1:ToolBox.PMT_count)
        pmt_values[pmt] = read(file["pmtmean"])[:,pmt]
    end 
    close(file)
    button_mask = zeros(Bool, ToolBox.PMT_count)
    for i in (1:ToolBox.PMT_count)
        on(buttons[i].clicks) do click
            button_mask[i] = !button_mask[i]
            empty!(ax1)
            for pmt in (1:ToolBox.PMT_count)
                button_mask[pmt] && scatter!(ax1, Times, pmt_values[pmt], color=ToolBox.config.Color[pmt%14+1], alpha=alpha, label = string("PMT", pmt))
            end
            [delete!(leg) for leg in fig.content if leg isa Legend]
            count(button_mask)>0 && Legend(fig[1,2], ax1, "active PMTs")#, framevisible = false)
        end
    end
    return fig
end