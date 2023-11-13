function ToolBox.plot_DomData_PMT(Dom_Id::Int32, PMTs::Vector{Int64}=Int64[]; T_intervall::Tuple{Integer,Integer}=(0,0), loadpath::String="../Data/DomData_Doms", slice_length::Integer=6000, information::Bool=false, activate_shift::Bool=false, add_linFit::Bool=false, savepicture::Bool=true)
    if PMTs == Int64[]
        PMTs = (1:ToolBox.PMT_count)
    end
    file = h5open(string(loadpath, "/Dom_", Dom_Id,"_",Int32(slice_length/600),".h5"), "r")
    figure = Figure(resolution=(900,600),fontsize=30)
    ylabel="detection rates in kHz"
    activate_shift && (ylabel="relative detection rates")
    axf1 = Axis(figure[1,1], xlabel="Time", ylabel=ylabel) #, title="frequencies"
    # Label(figure[0, 1], string("Data of Dom ", Dom_Id), fontsize = 30, tellwidth = false)
    Times = read(file["Time"])
    time_mask = ToolBox.maskTime(Times,T_intervall)
    pmtmean = ToolBox.shift_pmtmeans(read(file["pmtmean"]), activate_shift=activate_shift, shift_to_matrixmean=true)./1000
    hrv_count = read(file["hrvcount"])
    Zeiten, Typ = ToolBox.autoscale_time(minimum(Times[time_mask]), maximum(Times[time_mask]), intervalls=2)
    axf1.xticks = (datetime2unix.(Zeiten.dates) , Dates.format.(Zeiten.dates, Typ))
    close(file)
    for pmt in PMTs
        scatter!(axf1, Times[time_mask], pmtmean[time_mask,pmt], markersize=10, label=string("PMT ", pmt))
        if add_linFit
            funktion = linfit_DomData((Dom_Id, pmt), return_function=true)
            lines!(axf1, Times[time_mask], funktion.(Times[time_mask]))
        end
        if information
            axf2 = Axis(figure[2,1], title="hrvcount", xlabel="Time", ylabel="hrv count") 
            scatter!(axf2, Times[time_mask], hrv_count[time_mask,pmt])
            axf2.xticks = (datetime2unix.(Zeiten.dates) , Dates.format.(Zeiten.dates, Typ))
        end
    end
    columns = ceil(Int64, length(PMTs)/3)
    if length(PMTs) < 13
        axislegend(axf1, nbanks = columns)#, position = :lt)
    end
    if savepicture
        save(string("temp_pictures/PMTs_",Dom_Id,".png"), figure)
    end
    return figure
end


function ToolBox.plot_DomData_Drift(Dom_Id::Int32, PMT::Integer; T_intervall::Tuple{Integer,Integer}=(0,0), loadpath::String="../Data/DomData_Doms", slice_length::Integer=6000, add_linFit::Bool=true, information::Bool=true, savepicture::Bool=true)
    file = h5open(string(loadpath, "/Dom_", Dom_Id,"_",Int32(slice_length/600),".h5"), "r")
    figure = Figure(resolution=(900,600),fontsize=30)
    axf1 = Axis(figure[1,1], xlabel="Time", ylabel="detection rates in kHz") #, title="frequencies"  ,title=string("Dom ",Dom_Id)
    axf1.aspect = 1.5
    # Label(figure[0, 1], string("Dom ", Dom_Id), fontsize = 20, tellwidth = false)
    Times = read(file["Time"])
    time_mask = ToolBox.maskTime(Times,T_intervall)
    pmtmean = read(file["pmtmean"])./1000
    Zeiten, Typ = ToolBox.autoscale_time(minimum(Times[time_mask]), maximum(Int64, Times[time_mask]), intervalls=3)
    axf1.xticks = (datetime2unix.(Zeiten.dates) , Dates.format.(Zeiten.dates, Typ))
    close(file)
    if add_linFit
        y, slope = linfit_DomData((Dom_Id, PMT))
        f(x) = y+slope*x
        fit_values = f.(Times[time_mask])./1000
        slope = round(slope*3600*24, digits=2)
        lines!(axf1, Times[time_mask], fit_values, label=string("linear fit with slope: ", slope,"Hz/Day"), linewidth=7, color=:darkgreen)
    end
    scatter!(axf1, Times[time_mask], pmtmean[time_mask,PMT], markersize=10, label=string("PMT ",PMT))
    if information && add_linFit
        if slope < 0
            axislegend(axf1, position = :lb)
        else 
            axislegend(axf1, position = :lt)
        end
    end
    if savepicture
        save(string("temp_pictures/Drift",Dom_Id,"_",PMT,".png"), figure)
    end
    return figure
end


function ToolBox.plot_DomData(Dom_Ids::Vector{Int32}, values::Vector{T}, detector::Detector; labels::Bool=false, valuetyp::String="", storagename::String=empty, add_markersizes::Bool=false, colored::Bool=true, size_bounds::Tuple{Real,Real}=(5,25), savepicture::Bool=true) where T <: Real
    if length(Dom_Ids) == length(values)
        Dom_positions = ToolBox.pos_Doms(Dom_Ids, detector)
        positions = [v for (k,v) in Dom_positions]
        fig = Figure(resolution=(1000,650),fontsize=30)
        ax3d = Axis3(fig[1,1], xlabel = "x-Position in m", ylabel = "y-Position in m", zlabel = "z-Position in m", xlabeloffset=50, ylabeloffset=50, zlabeloffset=65) #, title = "position of the Doms inside the Detector"
        sizes = values .* ((size_bounds[2]-size_bounds[1])/(maximum(values)-minimum(values)))
        sizes = sizes .+ (size_bounds[1]-minimum(sizes))
        if add_markersizes && colored
            scatter!(ax3d, positions, color = values, colormap = :thermal, colorrange = (minimum(values), maximum(values)), markersize=sizes)
        elseif colored
            scatter!(ax3d, positions, color = values, colormap = :thermal, colorrange = (minimum(values), maximum(values)))
        elseif add_markersizes
            scatter!(ax3d, positions, markersize=sizes)
        else  
            scatter!(ax3d, positions)
        end
        labels && text!(ax3d, Point3f(positions), text = string(Dom_Ids))
        colored && Colorbar(fig[1, 2], limits = (minimum(values), maximum(values)), colormap = :thermal, label=valuetyp)
        # colsize!(fig.layout, 1, Aspect(1, 1))
        if savepicture
            save(string("temp_pictures/DomData_",storagename,".png"), fig)
        end
        return fig
    else
        error("wrong Vector lengths")
    end
end
function ToolBox.plot_DomData(DomData_dict::Dict{Int32, T}, detector::Detector; labels::Bool=false, storagename::String="empty", valuetyp::String="", add_markersizes::Bool=false, colored::Bool=true, size_bounds::Tuple{Real,Real}=(5,25), savepicture::Bool=true) where T <: Real
    Dom_Ids = collect(keys(DomData_dict))
    values = [DomData_dict[Dom] for Dom in Dom_Ids]
    ToolBox.plot_DomData(Dom_Ids, values, detector, labels=labels, valuetyp=valuetyp, storagename=storagename, add_markersizes=add_markersizes, colored=colored, size_bounds=size_bounds, savepicture=savepicture)
end



function ToolBox.plot_lf_Intervalls(Dom_Id::Int32; pmt=0, Name::String="lfData_Intervalls", storagepath::String=".", save_picture::Bool=false) #approx derivative with bad resolution
    timesteps, params = ToolBox.load_linFitData_intervalls(Dom_Id, Name=Name, storagepath=storagepath)
    figure = Figure()
    axf1 = Axis(figure[1,1], title="slopes of PMT", xlabel="Time", ylabel="slope of the linear fit in Hz/h") 
    if pmt isa Number && pmt != 0
        slopes = [params[pmt,i][2] for i in (1:size(params)[2])]*3600
        scatter!(axf1, timesteps, slopes)
    elseif pmt isa Number
        for PMT in (1:ToolBox.PMT_count)
            slopes = [params[PMT,i][2] for i in (1:size(params)[2])]*3600
            scatter!(axf1, timesteps, slopes)
        end
    else 
        for PMT in pmt
            slopes = [params[PMT,i][2] for i in (1:size(params)[2])]*3600
            scatter!(axf1, timesteps, slopes)
        end
    end
    if save_picture
        save(string("temp_pictures/lF_Int_", Dom_Id, ".", pmt, ".png"), figure)
    end
    Zeiten, Typ = ToolBox.autoscale_time(minimum(timesteps), maximum(timesteps))
    axf1.xticks = (datetime2unix.(Zeiten.dates) , Dates.format.(Zeiten.dates, Typ))
    return figure
end

function ToolBox.plot_Doms(Dom_Ids::Vector{Int32}, detector::Detector; labels::Bool=true, plot_allDoms::Bool=false, savepicture::Bool=true, storagename::String="Dom_positions")
    positions = ToolBox.pos_Doms(detector)
    fig = Figure(resolution=(1000,700),fontsize=30)
    ax3d = Axis3(fig[1,1], xlabel = "x-Position in m", ylabel = "y-Position in m", zlabel = "z-Position in m", xlabeloffset=50, ylabeloffset=50, zlabeloffset=65) #, title = "Positions of the Doms inside the Detector"
    if plot_allDoms
        all_positions = [v for (k,v) in positions]
        scatter!(ax3d, all_positions; markersize=9)
    end
    for Dom in Dom_Ids
        scatter!(ax3d, Point3f(positions[Dom]), markersize=13, color=:darkred)
        labels && text!(ax3d, Point3f(positions[Dom]), text = string(Dom))
    end
    if savepicture
        save(string("temp_pictures/Detector_",storagename,".png"), fig)
    end
    return fig
end
function ToolBox.plot_Doms(detector::Detector; labels::Bool=false)
    return ToolBox.plot_Doms(ToolBox.optical_DomIds(detector), detector, labels=labels)
end

function ToolBox.plot_surrounding_Doms(Dom_Id::Int32, range::Real, detector::Detector; T_intervall::Tuple{Integer,Integer}=(0,0), loadpath::String="../Data/DomData_Doms", slice_length::Integer=6000, activate_shift::Bool=true, mainDomInfo::Bool=false, alpha::Float64=1.0, ylim=(nothing,nothing))
    figure = Figure(resolution=(900,600),fontsize=30)
    if mainDomInfo
        axf1 = Axis(figure[1,1], title="detection rates of main Dom", xlabel="Time", ylabel="detection rates in Hz")
        axf2 = Axis(figure[2,1], title="detection rates of his neighbours", xlabel="Time", ylabel="detection rates in Hz") 
        ylimits!(axf2, ylim) 
    else
        axf2 = Axis(figure[1,1], title="detection rates of his neighbours", xlabel="Time", ylabel="detection rates in Hz") 
        ylimits!(axf2, ylim) 
    end
    Label(figure[0, 1], string("Data around Dom ", Dom_Id, " in the range of: ", range, "m"), fontsize = 30, tellwidth = false)
    Doms = ToolBox.close_Doms(Dom_Id, range, detector)
    file = h5open(string(loadpath, "/Dom_", Dom_Id,"_",Int32(slice_length/600),".h5"), "r")
    Times = read(file["Time"])
    main_pmtmean = ToolBox.shift_pmtmeans(read(file["pmtmean"]), activate_shift=activate_shift, shift_to_matrixmean=true)
    close(file)    
    time_mask = ToolBox.maskTime(Times,T_intervall)
    main_Dom_mean = vec(mean(main_pmtmean, dims=2))
    Zeiten, Typ = ToolBox.autoscale_time(minimum(Times[time_mask]), maximum(Int64, Times[time_mask]), intervalls=3)
    if mainDomInfo
        for pmt in (1:ToolBox.PMT_count)
            scatter!(axf1, Times[time_mask], main_pmtmean[time_mask,pmt], alpha=alpha)
        end
        scatter!(axf1, Times[time_mask], main_Dom_mean[time_mask], marker=:cross)
        axf1.xticks = (datetime2unix.(Zeiten.dates) , Dates.format.(Zeiten.dates, Typ))
    end
    axf2.xticks = (datetime2unix.(Zeiten.dates) , Dates.format.(Zeiten.dates, Typ))
    for Dom in Doms
        file = h5open(string(loadpath, "/Dom_", Dom,"_",Int32(slice_length/600),".h5"), "r")
        Times = read(file["Time"])
        pmtmean = ToolBox.shift_pmtmeans(read(file["pmtmean"]), activate_shift=activate_shift, shift_to_matrixmean=true)
        close(file)
        Dom_mean = ToolBox.shift_pmtmeans(vec(mean(pmtmean, dims=2)), activate_shift=activate_shift)
        scatter!(axf2, Times[time_mask], Dom_mean[time_mask], markersize=8, alpha=alpha)
    end
    scatter!(axf2, Times[time_mask], ToolBox.shift_pmtmeans(main_Dom_mean[time_mask], activate_shift=activate_shift), markersize=10, alpha=alpha)
    return figure
end

function ToolBox.plot_multiple_Doms(Doms::Vector{Int32}; T_intervall::Tuple{Integer,Integer}=(0,0), loadpath::String="../Data/DomData_Doms", slice_length::Integer=6000, activate_shift::Bool=true, alpha::Float64=1.0, ylim=(nothing,nothing), legend::Bool=true)
    figure = Figure(resolution=(900,600),fontsize=30)
    axf2 = Axis(figure[1,1], title="detection rates of the Doms", xlabel="Time", ylabel="detection rates in Hz")
    ylimits!(axf2, ylim) 
    Times = Int64[]
    for Dom in Doms
        file = h5open(string(loadpath, "/Dom_", Dom,"_",Int32(slice_length/600),".h5"), "r")
        Times = read(file["Time"])
        time_mask = ToolBox.maskTime(Times,T_intervall)
        pmtmean = ToolBox.shift_pmtmeans(read(file["pmtmean"]), activate_shift=activate_shift, shift_to_matrixmean=true)
        close(file)
        Dom_mean = ToolBox.shift_pmtmeans(vec(mean(pmtmean, dims=2)), activate_shift=activate_shift)
        scatter!(axf2, Times[time_mask], Dom_mean[time_mask], markersize=10, alpha=alpha, label=Dom)
    end
    legend && axislegend(axf2, nbanks=2)
    Zeiten, Typ = ToolBox.autoscale_time(minimum(Times), maximum(Int64, Times), intervalls=3)
    axf2.xticks = (datetime2unix.(Zeiten.dates) , Dates.format.(Zeiten.dates, Typ))
    return figure
end

function ToolBox.plot_DomData_Floors()
    for floor in (1:18)
        fig = ToolBox.plot_DomData_Floors(floor)
        save(string("temp_pictures/floors/floor",floor,".png"), fig)
    end
end
function ToolBox.plot_DomData_Floors(Floor::Int64; T_intervall::Tuple{Integer,Integer}=(0,0), loadpath::String="../Data", alpha::Float64=1.0, slice_length=6000)
    file = h5open(string(loadpath, "/DomData_floors_",Int32(slice_length/600),".h5"), "r")
    figure = Figure(resolution=(900,600),fontsize=30)
    axf1 = Axis(figure[1,1], title="detection rates of PMT Ring A", xlabel="Time", ylabel="detection rates in Hz") 
    axf2 = Axis(figure[2,1], title="detection rates of PMT Ring B", xlabel="Time", ylabel="detection rates in Hz")
    axf3 = Axis(figure[3,1], title="detection rates of PMT Ring C", xlabel="Time", ylabel="detection rates in Hz")
    axf4 = Axis(figure[1,2], title="detection rates of PMT Ring D", xlabel="Time", ylabel="detection rates in Hz")
    axf5 = Axis(figure[2,2], title="detection rates of PMT Ring E", xlabel="Time", ylabel="detection rates in Hz")
    axf6 = Axis(figure[3,2], title="detection rates of PMT Ring F", xlabel="Time", ylabel="detection rates in Hz")
    Ax = [axf1, axf2, axf3, axf4, axf5, axf6]
    Label(figure[0, :], string("Data of Floor ", Floor), fontsize = 30)
    Times = read(file["Time"])
    time_mask = ToolBox.maskTime(Times,T_intervall)
    for ring in (1:length(ToolBox.config.Detector_PMT_Ringe))
        for pmt_number in (1:length(ToolBox.config.Detector_PMT_Ringe[ring]))
            scatter!(Ax[ring], Times[time_mask], read(file["pmtmean"])[Floor,:,ToolBox.config.Detector_PMT_Ringe[ring][pmt_number]][time_mask], color=ToolBox.config.Color[pmt_number], alpha=alpha)
        end 
    end
    Zeiten, Typ = ToolBox.autoscale_time(minimum(Times), maximum(Times), intervalls=3)
    for i in (1:6)
        Ax[i].xticks[] = (datetime2unix.(Zeiten.dates) , Dates.format.(Zeiten.dates, Typ)); 
    end
    close(file)
    return figure
end

function ToolBox.plot_Strings(detector::Detector; Strings::Vector{Int64}=Int64[])
    String_positions = ToolBox.pos_Strings(detector)
    if Strings == Int32[]
        Strings = collect(keys(String_positions))
    else
        Strings = [st for st in Strings if st in collect(keys(String_positions))]
    end
    fig = Figure(resolution=(900,600),fontsize=30)
    ax = Axis(fig[1,1], title = "position of the Strings inside the Detector", xlabel = "x-Position in m", ylabel = "y-Position in m")
    for stringi in Strings
        scatter!(ax, String_positions[stringi][1], String_positions[stringi][2], color=:darkgreen)
        text!(ax, String_positions[stringi], text = string(stringi), offset = (4, 0)) #, align = (:left, :top))
    end
    return fig
end


function ToolBox.plot_Det(detector::Detector; linfit::Bool=false)
    times, frequencies, Det_mean = ToolBox.mean_Detector(detector)
    y, m = ToolBox.linfit_Detector()
    f(x) = y+m*x
    fig = Figure()
    ax = Axis(fig[1,1])
    linfit && scatter!(ax, times, f.(times))
    scatter!(ax, times, frequencies)
    Zeiten, Typ = ToolBox.autoscale_time(minimum(times), maximum(Int64, times))
    ax.xticks = (datetime2unix.(Zeiten.dates) , Dates.format.(Zeiten.dates, Typ))
    return fig
end


function ToolBox.plot_DomData_Event(event::Event; loadpath::String="../Data/DomData_Doms", save_picture::Bool=false, time_added::Integer=100)
    file = h5open(string(loadpath, "/Dom_", event.Dom_Id,"_",Int32(event.slice_length/600),".h5"), "r")
    pmtmean = read(file["pmtmean"])[:,event.pmt]
    Times = read(file["Time"])
    figure = Figure()
    axf1 = Axis(figure[1,1:2], title="detection rates of PMT", xlabel="Time", ylabel="detection rates in Hz") 
    axf2 = Axis(figure[2,2], title="detection rates of event", xlabel="Time", ylabel="detection rates in Hz")
    axf3 = Axis(figure[3,1], title="good_values of event", xlabel="Time")
    axf4 = Axis(figure[3,2], title="hrvcount of event", xlabel="Time", ylabel="hrvcount")
    Ax = [axf2, axf3, axf4]
    #event intervall über die time beschränken
    event_intervall = ToolBox.maskTime(Times, (event.time_start-60*time_added, event.time_end+60*time_added))
    #event_intervall = event.array_start-10:event.array_end+10 
    event_mean = round(mean(read(file["pmtmean"])[event_intervall,event.pmt]),digits=3)
    event_max = round(maximum(read(file["pmtmean"])[event_intervall,event.pmt]),digits=3)
    Label(figure[2,1], string("Event of Dom ", event.Dom_Id, ", PMT ", event.pmt,"\nEventstart: ", unix2datetime(event.time_start), "\nEventlength: ",event.time_length, "min\nmean detection rates=",event_mean,"Hz,\nmax detection rates=",event_max,"Hz"), fontsize = 20, tellheight = false, tellwidth = false, justification=:center)   
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
    axf1 = Axis(figure[1,1], title="detection rates of PMT", xlabel="Time", ylabel="detection rates in Hz") 
    axf2 = Axis(figure[2,1], title="detection rates of PMT", xlabel="Time", ylabel="detection rates in Hz") 
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


function ToolBox.plot_DomData_Rings(DomID::Integer; T_intervall::Tuple{Integer,Integer}=(0,0), loadpath::String="../Data/DomData_Doms", alpha::Float64=1.0, slice_length::Integer=6000, mean_Rings::Bool=true)
    if !mean_Rings
        file = h5open(string(loadpath, "/Dom_", DomID,"_",Int32(slice_length/600),".h5"), "r")
        figure = Figure()
        axf1 = Axis(figure[1,1], title="detection rates of PMT Ring A", xlabel="Time", ylabel="detection rates in Hz") 
        axf2 = Axis(figure[2,1], title="detection rates of PMT Ring B", xlabel="Time", ylabel="detection rates in Hz")
        axf3 = Axis(figure[3,1], title="detection rates of PMT Ring C", xlabel="Time", ylabel="detection rates in Hz")
        axf4 = Axis(figure[1,2], title="detection rates of PMT Ring D", xlabel="Time", ylabel="detection rates in Hz")
        axf5 = Axis(figure[2,2], title="detection rates of PMT Ring E", xlabel="Time", ylabel="detection rates in Hz")
        axf6 = Axis(figure[3,2], title="detection rates of PMT Ring F", xlabel="Time", ylabel="detection rates in Hz")
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
    else
        file = h5open(string(loadpath, "/Dom_", DomID,"_",Int32(slice_length/600),".h5"), "r")
        figure = Figure()
        axf1 = Axis(figure[1,1], title="detection rates of the Rings", xlabel="Time", ylabel="detection rates in Hz") 
        Times = read(file["Time"])
        time_mask = ToolBox.maskTime(Times,T_intervall)
        for ring in (1:length(ToolBox.config.Detector_PMT_Ringe))
            pmt_mask = [ToolBox.config.Detector_PMT_Ring[i]==ring for i in (1:ToolBox.PMT_count)]
            ring_mean = mean(read(file["pmtmean"])[:,pmt_mask], dims=2)
            scatter!(axf1, Times[time_mask], ring_mean[time_mask], color=ToolBox.config.Color[ring], alpha=alpha, label=ring)
        end
        axislegend(axf1, nbanks=3)
        Zeiten, Typ = ToolBox.autoscale_time(minimum(Times), maximum(Int64, Times))
        axf1.xticks = (datetime2unix.(Zeiten.dates) , Dates.format.(Zeiten.dates, Typ)); 
        close(file)
    end
    return figure
end