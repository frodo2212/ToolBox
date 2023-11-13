#das sind die Interaktiven plots, die sind aber nicht so gut...
"""
This was a test of interactive plotting with all the Options in GLMakie available
"""
function ToolBox.plot_DomData_Rings_interactive(DomID::Integer; loadpath::String="../Data/DomData_Doms", alpha::Float64=1.0, slice_length::Integer=6000)
    fig = Figure()
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
    lsgrid = SliderGrid(fig[2,:],
    (label = "start", range = minimum(Times):1:maximum(Times), startvalue = minimum(Times)),
    (label = "end", range = minimum(Times):1:maximum(Times), startvalue = maximum(Times))
    # formats = [x -> "$(round(x))"],
    )
    Start = lsgrid.sliders[1].value
    End = lsgrid.sliders[2].value
    time_mask = lift(ToolBox.maskTime(Times,(Start,End)))
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

"""
given a Dom Id one can activate the PMT numbers that are plottet in the Graph for interactive visibility
"""
function ToolBox.plot_DomData_PMT_interactive(DomID::Integer; loadpath::String="../Data/DomData_Doms", alpha::Float64=1.0, slice_length::Integer=6000, shift_means::Bool=false)
    fig = Figure()
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
    pmt_values = read(file["pmtmean"])
    shift_means && (pmt_values=ToolBox.shift_pmtmeans(pmt_values))
    Zeiten, Typ = ToolBox.autoscale_time(minimum(Times), maximum(Int64, Times), intervalls=3)
    ax1.xticks[] = (datetime2unix.(Zeiten.dates) , Dates.format.(Zeiten.dates, Typ)) 
    close(file)
    button_mask = zeros(Bool, ToolBox.PMT_count)
    for i in (1:ToolBox.PMT_count)
        on(buttons[i].clicks) do click
            button_mask[i] = !button_mask[i]
            empty!(ax1)
            for pmt in (1:ToolBox.PMT_count)
                button_mask[pmt] && scatter!(ax1, Times, pmt_values[:,pmt], color=ToolBox.config.Color[pmt%14+1], alpha=alpha, label = string("PMT", pmt))
            end
            [delete!(leg) for leg in fig.content if leg isa Legend]
            count(button_mask)>0 && Legend(fig[1,2], ax1, "active PMTs", framevisible = false)
        end
    end
    return fig
end

"""
plots all Doms of the whole Detector in a 3D-Object and points out a given reference Dom
one can then adjust all Doms closer than distance r to the referenceDom via Slider
"""
function ToolBox.plot_closeDoms_3D_interactive(Dom_Id::Int32, detector::Detector; range_bound::Tuple{Real,Real}=(0,200), stringnumbers::Bool=false, plot_allDoms::Bool=true)
    #static part
    Dom_positions = ToolBox.pos_Doms(detector)
    aspect=(1, 1, 1)
    perspectiveness=0.5
    fig = Figure(; resolution=(1200, 400))
    ax1 = Axis3(fig[1, 1]; aspect, perspectiveness)
    if plot_allDoms
        all_positions = [v for (k,v) in Dom_positions]
        scatter!(ax1, all_positions; markersize=10)
        scatter!(Dom_positions[Dom_Id]; markersize=18, color=:darkred)
    end
    if stringnumbers
        str_pos = ToolBox.pos_Strings(detector)
        Strings = collect(keys(str_pos))
        string_positions = [(str_pos[i][1],str_pos[i][2],0) for i in Strings]
        scatter!(ax1, string_positions)
        for i in (1:length(Strings))
            text!(ax1, string_positions[i], text = string(Strings[i]), offset = (4, 0))
        end
    end
    #interactive part
    # sl_x = Slider(fig[2, 1], range = range_bound[1]:5:range_bound[2], startvalue = range_bound[1])
    sl_x = SliderGrid(fig[2,1],(label = "range in m", range = range_bound[1]:5:range_bound[2], startvalue = range_bound[1]))
    marked_positions = lift(sl_x.sliders[1].value) do range
        [Dom_positions[i] for i in ToolBox.close_Doms(Dom_Id, range, detector)]
    end
    scatter!(ax1, marked_positions; markersize=15, color=:darkgreen)
    return fig
end