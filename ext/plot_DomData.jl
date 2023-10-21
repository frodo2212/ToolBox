function ToolBox.plot_DomData_Rings(DomID::Int; T_intervall::Tuple{Integer,Integer}=(0,0), loadpath::String="../Data/DomData_Doms", alpha::Float64=1.0)
    file = h5open(string(loadpath, "/Dom_", Dom_object[1],"_",Int32(slice_length/600),".h5"), "r")
    figure = Figure()
    axf1 = Axis(figure[1,1], title="frequencies of PMT Ring A", xlabel="Time", ylabel="frequency in Hz") 
    axf2 = Axis(figure[2,1], title="frequencies of PMT Ring B", xlabel="Time", ylabel="frequency in Hz")
    axf3 = Axis(figure[3,1], title="frequencies of PMT Ring C", xlabel="Time", ylabel="frequency in Hz")
    axf4 = Axis(figure[1,2], title="frequencies of PMT Ring D", xlabel="Time", ylabel="frequency in Hz")
    axf5 = Axis(figure[2,2], title="frequencies of PMT Ring E", xlabel="Time", ylabel="frequency in Hz")
    axf6 = Axis(figure[3,2], title="frequencies of PMT Ring F", xlabel="Time", ylabel="frequency in Hz")
    Ax = [axf1, axf2, axf3, axf4, axf5, axf6]
    Label(figure[0, :], string("Data of Dom ", DomID), fontsize = 30)
    Times = read(file["Time"])
    time_mask = ToolBox.maskTime(Times,T_intervall)
    for ring in (1:length(ToolBox.config.Detector_PMT_Ringe))
        for pmt_number in (1:length(ToolBox.config.Detector_PMT_Ringe[ring]))
            scatter!(Ax[ring], Times[time_mask], read(file["pmtmean"])[:,ToolBox.config.Detector_PMT_Ringe[ring][pmt_number]][time_mask], color=ToolBox.config.Color[pmt_number], alpha=alpha)
        end 
    end
    Zeiten, Typ = ToolBox.autoscale_time(minimum(Times), maximum(Int64, Times), intervalls=3)
    for i in (1:6)
        Ax[i].xticks[] = (datetime2unix.(Zeiten.dates) , Dates.format.(Zeiten.dates, Typ)); 
    end
    close(file)
    return figure
end

function ToolBox.plot_DomData_PMT(Dom_object::Tuple{Integer,Integer}; T_intervall::Tuple{Integer,Integer}=(0,0), loadpath::String="../Data/DomData_Doms", alpha::Float64=1.0, slice_length::Integer=6000, information::Bool=false)
    file = h5open(string(loadpath, "/Dom_", Dom_object[1],"_",Int32(slice_length/600),".h5"), "r")
    figure = Figure()
    axf1 = Axis(figure[1,1], title="frequencies", xlabel="Time", ylabel="frequency in Hz") 
    Label(figure[0, 1], string("Data of Dom ", Dom_object[1], " PMT ", Dom_object[2]), fontsize = 30, tellwidth = false)
    Times = read(file["Time"])
    time_mask = ToolBox.maskTime(Times,T_intervall)
    scatter!(axf1, Times[time_mask], read(file["pmtmean"])[time_mask,Dom_object[2]], alpha=alpha)
    Zeiten, Typ = ToolBox.autoscale_time(minimum(Times), maximum(Int64, Times), intervalls=3)
    axf1.xticks = (datetime2unix.(Zeiten.dates) , Dates.format.(Zeiten.dates, Typ))
    if information
        axf2 = Axis(figure[2,1], title="hrvcount", xlabel="Time", ylabel="hrv count") 
        scatter!(axf2, Times[time_mask], read(file["hrvcount"])[time_mask,Dom_object[2]], alpha=alpha)
        axf2.xticks = (datetime2unix.(Zeiten.dates) , Dates.format.(Zeiten.dates, Typ))
    end
    close(file)
    return figure
end

function ToolBox.plot_DomData_PMT_rel(Dom_Id::Int32; PMTs::Vector{Int64}=Int64[], T_intervall::Tuple{Integer,Integer}=(0,0), loadpath::String="../Data/DomData_Doms", alpha::Float64=0.4, slice_length::Integer=6000, information::Bool=false, activate_shift::Bool=true, mean_goal::Float64=5500.0)
    if PMTs == Int64[]
        PMTs = (1:ToolBox.PMT_count)
    end
    file = h5open(string(loadpath, "/Dom_", Dom_Id,"_",Int32(slice_length/600),".h5"), "r")
    figure = Figure()
    axf1 = Axis(figure[1,1], title="frequencies", xlabel="Time", ylabel="frequency in Hz") 
    Label(figure[0, 1], string("Data of Dom ", Dom_Id), fontsize = 30, tellwidth = false)
    Times = read(file["Time"])
    time_mask = ToolBox.maskTime(Times,T_intervall)
    pmtmean = read(file["pmtmean"])
    activate_shift && (pmtmean = ToolBox.shift_pmtmeans(pmtmean, mean_goal=mean_goal))
    hrv_count = read(file["hrvcount"])
    Zeiten, Typ = ToolBox.autoscale_time(minimum(Times), maximum(Int64, Times), intervalls=3)
    axf1.xticks = (datetime2unix.(Zeiten.dates) , Dates.format.(Zeiten.dates, Typ))
    close(file)
    for pmt in PMTs
        scatter!(axf1, Times[time_mask], pmtmean[time_mask,pmt], alpha=alpha)
        if information
            axf2 = Axis(figure[2,1], title="hrvcount", xlabel="Time", ylabel="hrv count") 
            scatter!(axf2, Times[time_mask], hrv_count[time_mask,pmt], alpha=alpha)
            axf2.xticks = (datetime2unix.(Zeiten.dates) , Dates.format.(Zeiten.dates, Typ))
        end
    end
    return figure
end

function ToolBox.plot_DomData_Floors(Floor::Int64; T_intervall::Tuple{Integer,Integer}=(0,0), loadpath::String="../Data", alpha::Float64=1.0)
    file = h5open(string(loadpath, "/DomDataV3_Floors.h5"), "r")
    figure = Figure()
    axf1 = Axis(figure[1,1], title="frequencies of PMT Ring A", xlabel="Time", ylabel="frequency in Hz") 
    axf2 = Axis(figure[2,1], title="frequencies of PMT Ring B", xlabel="Time", ylabel="frequency in Hz")
    axf3 = Axis(figure[3,1], title="frequencies of PMT Ring C", xlabel="Time", ylabel="frequency in Hz")
    axf4 = Axis(figure[1,2], title="frequencies of PMT Ring D", xlabel="Time", ylabel="frequency in Hz")
    axf5 = Axis(figure[2,2], title="frequencies of PMT Ring E", xlabel="Time", ylabel="frequency in Hz")
    axf6 = Axis(figure[3,2], title="frequencies of PMT Ring F", xlabel="Time", ylabel="frequency in Hz")
    Ax = [axf1, axf2, axf3, axf4, axf5, axf6]
    Label(figure[0, :], string("Data of Floor ", Floor), fontsize = 30)
    Times = read(file["Time"])
    time_mask = ToolBox.maskTime(Times,T_intervall)
    for ring in (1:length(config.Detector_PMT_Ringe))
        for pmt_number in (1:length(config.Detector_PMT_Ringe[ring]))
            scatter!(Ax[ring], Times[time_mask], read(file["pmtmean"])[Floor,:,config.Detector_PMT_Ringe[ring][pmt_number]][time_mask], color=config.Color[pmt_number], alpha=alpha)
        end 
    end
    Zeiten, Typ = ToolBox.autoscale_time(minimum(Times), maximum(Times), intervalls=3)
    for i in (1:6)
        Ax[i].xticks[] = (datetime2unix.(Zeiten.dates) , Dates.format.(Zeiten.dates, Typ)); 
    end
    close(file)
    return figure
end

function ToolBox.plot_DomData_Event(event::Event; loadpath::String="../Data/DomData_Doms", save_picture::Bool=false, time_added::Integer=100)
    file = h5open(string(loadpath, "/Dom_", event.Dom_Id,"_",Int32(event.slice_length/600),".h5"), "r")
    pmtmean = read(file["pmtmean"])[:,event.pmt]
    Times = read(file["Time"])
    figure = Figure()
    axf1 = Axis(figure[1,1:2], title="frequencies of PMT", xlabel="Time", ylabel="frequency in Hz") 
    axf2 = Axis(figure[2,2], title="frequencies of event", xlabel="Time", ylabel="frequency in Hz")
    axf3 = Axis(figure[3,1], title="good_values of event", xlabel="Time")
    axf4 = Axis(figure[3,2], title="hrvcount of event", xlabel="Time", ylabel="hrvcount")
    Ax = [axf2, axf3, axf4]
    #event intervall über die time beschränken
    event_intervall = maskTime(Times, (event.time_start-60*time_added, event.time_end+60*time_added))
    #event_intervall = event.array_start-10:event.array_end+10 
    event_mean = round(mean(read(file["pmtmean"])[event_intervall,event.pmt]),digits=3)
    event_max = round(maximum(read(file["pmtmean"])[event_intervall,event.pmt]),digits=3)
    Label(figure[2,1], string("Event of Dom ", event.Dom_Id, ", PMT ", event.pmt,"\nEventstart: ", unix2datetime(event.time_start), "\nEventlength: ",event.time_length, "min\nmean Frequency=",event_mean,"Hz,\nmax Frequency=",event_max,"Hz"), fontsize = 20, tellheight = false, tellwidth = false, justification=:center)   
    scatter!(axf2, Times[event_intervall], pmtmean[event_intervall])
    scatter!(axf1, Times, pmtmean)
    scatter!(axf1, Times[event_intervall], pmtmean[event_intervall])
    scatter!(axf4, Times[event_intervall], read(file["hrvcount"])[event_intervall,event.pmt])
    scatter!(axf3, Times[event_intervall], read(file["good_values"])[event_intervall])
    close(file)
    Zeiten, Typ = ToolBox.autoscale_time(minimum(Times[event_intervall]), maximum(Times[event_intervall]), intervalls=2)
    for i in (1:3)
        Ax[i].xticks[] = (datetime2unix.(Zeiten.dates) , Dates.format.(Zeiten.dates, Typ)); 
    end
    Zeiten, Typ = ToolBox.autoscale_time(minimum(Times), maximum(Times))
    axf1.xticks[] = (datetime2unix.(Zeiten.dates) , Dates.format.(Zeiten.dates, Typ)); 
    if save_picture
        save(string("temp_pictures/", event.Dom_Id, "_", event.pmt, "_", event.array_start, ".png"), figure)
    end
    return figure
end

function ToolBox.plot_DomData_Event_array(Events::Vector{Event}; loadpath::String="../Data/DomData_Doms", time_added::Integer=100)
    for Event in Events
        ToolBox.plot_DomData_Event(Event, loadpath=loadpath, save_picture=true, time_added=time_added)
    end
end

function ToolBox.plot_DomData_linFit(event::linfitData; loadpath::String="../Data/DomData_Doms", save_picture::Bool=false)
    file = h5open(string(loadpath, "/Dom_", event.Dom_Id,"_",Int32(event.slice_length/600),".h5"), "r")
    pmtmean = read(file["pmtmean"])[:,event.pmt]
    Times = read(file["Time"])
    close(file)
    figure = Figure()
    axf1 = Axis(figure[1,1], title="frequencies of PMT", xlabel="Time", ylabel="frequency in Hz") 
    axf2 = Axis(figure[2,1], title="frequencies of PMT", xlabel="Time", ylabel="frequency in Hz") 
    scatter!(axf1, Times, pmtmean)
    scatter!(axf2, Times, pmtmean)
    scatter!(axf2, Times, [event.params[1]+event.params[2]*i for i in Times], alpha=0.7, marker=:hline)
    if save_picture
        save(string("temp_pictures/lF_", event.Dom_Id, "_", event.pmt, ".png"), figure)
    end
    return figure
end

function ToolBox.plot_DomData_linFit_array(Events::Vector{linfitData}; loadpath::String="../Data/DomData_Doms")
    for Event in Events
        plot_Data_linFit(Event, loadpath=loadpath, save_picture=true)
    end
end

function ToolBox.plot_Strings(detector::Detector)
    ToolBox.plot_Strings(ToolBox.pos_Strings(detector))
end

function ToolBox.plot_Strings(String_positions::Dict{Int64, Tuple{Float64, Float64}}; Strings::Vector{Int64}=Int64[])
    if Strings == Int32[]
        Strings = collect(keys(String_positions))
    else
        Strings = [st for st in Strings if st in collect(keys(String_positions))]
    end
    fig = Figure()
    ax = Axis(fig[1,1], title = "position of the Strings inside the Detector", xlabel = "x-Position in m", ylabel = "y-Position in m")
    elements = []
    for stringi in Strings
        push!(elements,scatter!(ax, String_positions[stringi][1], String_positions[stringi][2]))
        text!(ax, String_positions[stringi], text = string(stringi), offset = (4, 0)) #, align = (:left, :top))
    end
    return fig
end


function ToolBox.plot_Doms(Dom_Ids::Vector{Int32}, detector::Detector; labels::Bool=true)
    positions = ToolBox.pos_Doms(Dom_Ids, detector)
    fig = Figure()
    ax3d = Axis3(fig[1,1], title = "position of the Doms inside the Detector", xlabel = "x-Position in m", ylabel = "y-Position in m", zlabel = "z-Position in m")
    for Dom in Dom_Ids
        scatter!(ax3d, Point3f(positions[Dom]))
        labels && text!(ax3d, Point3f(positions[Dom]), text = string(Dom))
    end
    return fig
end

function ToolBox.plot_Doms(detector::Detector; labels::Bool=false)
    return ToolBox.plot_Doms(ToolBox.optical_DomIds(detector), detector, labels=labels)
end

function ToolBox.plot_Doms(Dom_Ids::Vector{Int32}, values::Vector{T}, detector::Detector; labels::Bool=false, valuetyp::String="", add_markersizes::Bool=false, colored::Bool=true, size_bounds::Tuple{Real,Real}=(5,25)) where T <: Real
    if length(Dom_Ids) == length(values)
        positions = ToolBox.pos_Doms(Dom_Ids, detector)
        fig = Figure()
        ax3d = Axis3(fig[1,1], title = "position of the Doms inside the Detector", xlabel = "x-Position in m", ylabel = "y-Position in m", zlabel = "z-Position in m")
        sizes = values .* ((size_bounds[2]-size_bounds[1])/(maximum(values)-minimum(values)))
        sizes = sizes .+ (size_bounds[1]-minimum(sizes))
        if add_markersizes && colored
            scatter!(ax3d, Point3f(positions), color = values, colormap = :thermal, colorrange = (minimum(values), maximum(values)), markersize=sizes)
        elseif colored
            scatter!(ax3d, Point3f(positions), color = values, colormap = :thermal, colorrange = (minimum(values), maximum(values)))
        elseif add_markersizes
            scatter!(ax3d, Point3f(positions), markersize=sizes)
        else  
            scatter!(ax3d, Point3f(positions))
        end
        labels && text!(ax3d, Point3f(positions), text = string(Dom_Ids))
        colored && Colorbar(fig[1, 2], limits = (minimum(values), maximum(values)), colormap = :thermal, label=valuetyp)
        return fig
    else
        error("wrong Vector lengths")
    end
end

function ToolBox.plot_Doms(DomData_dict::Dict{Int32, T}, detector::Detector; labels::Bool=false, valuetyp::String="", add_markersizes::Bool=false, colored::Bool=true, size_bounds::Tuple{Real,Real}=(5,25)) where T <: Real
    Dom_Ids = collect(keys(DomData_dict))
    values = [DomData_dict[Dom] for Dom in Dom_Ids]
    ToolBox.plot_Doms_colored(Dom_Ids, values, detector, labels=labels, valuetyp=valuetyp, add_markersizes=add_markersizes, colored=colored, size_bounds=size_bounds)
end