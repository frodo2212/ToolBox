function ToolBox.plot_DomDataV3_Rings(DomID::Int; T_intervall::Tuple{Integer,Integer}=(0,0), loadpath::String="../Data/DomData_Doms", alpha::Float64=1.0)
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


function ToolBox.plot_DomDataV3_PMT(Dom_object::Tuple{Integer,Integer}; T_intervall::Tuple{Integer,Integer}=(0,0), loadpath::String="../Data/DomData_Doms", alpha::Float64=1.0, slice_length::Integer=6000, information::Bool=false)
    #file = h5open(string(loadpath, "/Dom_", Dom_object[1],"_",Int32(slice_length/600),".h5"), "r")
    file = h5open(string(loadpath, "/Data_", Dom_object[1],".h5"), "r") #das ist alt. muss weg sobald die neuen Daten da sind
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


function ToolBox.plot_DomDataV3_Floor(Floor::Int64; T_intervall::Tuple{Integer,Integer}=(0,0), loadpath::String="../Data", alpha::Float64=1.0)
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



function ToolBox.plot_DataV3_Event_test1(event::Event; loadpath::String="../Data/DomData_Doms", save_picture::Bool=false, time_added::Integer=100)
    file = h5open(string(loadpath, "/Data_", event.Dom_Id,".h5"), "r")
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

function ToolBox.plot_Data_Event_test2(Events::Vector{Event}; loadpath::String="../Data/DomData_Doms", time_added::Integer=100)
    for Event in Events
        ToolBox.plot_DataV3_Event_test1(Event, loadpath=loadpath, save_picture=true, time_added=time_added)
    end
end

#die hier sind für intervalle
function ToolBox.plot_Data_linFit_test1(Event::linfitData; loadpath::String="../Data/DomData_Doms", save_picture::Bool=false)
    file = h5open(string(loadpath, "/Data_", Event.Dom_Id,".h5"), "r")
    pmtmean = read(file["pmtmean"])[:,Event.pmt]
    Times = read(file["Time"])
    event_intervall = ToolBox.maskTime(Times, (Event.time-30*Event.time_intervall, Event.time+30*Event.time_intervall))
    close(file)
    figure = Figure()
    axf1 = Axis(figure[1,1], title="frequencies of PMT", xlabel="Time", ylabel="frequency in Hz") 
    axf2 = Axis(figure[2,1], title="frequencies of PMT", xlabel="Time", ylabel="frequency in Hz") 
    scatter!(axf1, Times, pmtmean)
    scatter!(axf1, Times[event_intervall], pmtmean[event_intervall])
    scatter!(axf2, Times[event_intervall], pmtmean[event_intervall])
    Label(figure[0,1], "bisher ist der y achsenAbschnitt falsch, da linear_fit nicht funktioniert...", fontsize = 20, tellwidth = false, justification=:center)   
    #y_provisorisch = pmtmean[event_intervall][1]-Event.params[2]*Times[event_intervall][1]
    scatter!(axf2, Times[event_intervall], [Event.params[1]+Event.params[2]*i for i in Times[event_intervall]], alpha=0.7)
    if save_picture
        save(string("temp_pictures/", Event.Dom_Id, "_", Event.pmt, "_", round(Event.params[2],digits=2), ".png"), figure)
    end
    return figure
end

function ToolBox.plot_Data_linFit_test2(Events::Vector{linfitData}; loadpath::String="../Data/DomData_Doms")
    for Event in Events
        ToolBox.plot_Data_linFit_test1(Event, loadpath=loadpath, save_picture=true)
    end
end


#doe hier sind noch sehr experimentell,
#primär, weil die Daten, die rein kommen noch nicht gut sind...
function ToolBox.plot_Data_linFit_ganz(Event::linfitData; loadpath::String="../Data/DomData_Doms", save_picture::Bool=false)
    file = h5open(string(loadpath, "/Data_", Event.Dom_Id,".h5"), "r")
    pmtmean = read(file["pmtmean"])[:,Event.pmt]
    Times = read(file["Time"])
    close(file)
    figure = Figure()
    axf1 = Axis(figure[1,1], title="frequencies of PMT", xlabel="Time", ylabel="frequency in Hz") 
    axf2 = Axis(figure[2,1], title="frequencies of PMT", xlabel="Time", ylabel="frequency in Hz") 
    scatter!(axf1, Times, pmtmean)
    Label(figure[0,1], "bisher ist der y achsenAbschnitt falsch, da linear_fit nicht funktioniert...", fontsize = 20, tellwidth = false, justification=:center)   
    y_provisorisch = pmtmean[1]-Event.params[2]*Times[1]
    scatter!(axf2, Times, [y_provisorisch+Event.params[2]*i for i in Times], alpha=0.7)
    #scatter!(axf2, Times, [Event.params[1]+Event.params[2]*i for i in Times], alpha=0.7)
    if save_picture
        save(string("temp_pictures/lF_", Event.Dom_Id, "_", Event.pmt, ".png"), figure)
    end
    return figure
end

function ToolBox.plot_Data_linFit_array(Events::Vector{linfitData}; loadpath::String="../Data/DomData_Doms")
    for Event in Events
        plot_Data_linFit_ganz(Event, loadpath=loadpath, save_picture=true)
    end
end