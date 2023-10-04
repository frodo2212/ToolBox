
function ToolBox.plot_DomDataV3_Rings(DomID::Int; T_intervall::Tuple{Integer,Integer}=(0,0), loadpath::String="../Data/DomData_Doms", alpha::Float64=1.0)
    file = h5open(string(loadpath, "/Data_", DomID,".h5"), "r")
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
    time_mask = maskTime(Times,T_intervall)
    for ring in (1:length(config.Detector_PMT_Ringe))
        for pmt_number in (1:length(config.Detector_PMT_Ringe[ring]))
            scatter!(Ax[ring], Times[time_mask], read(file["pmtmean"])[:,config.Detector_PMT_Ringe[ring][pmt_number]][time_mask], color=config.Color[pmt_number], alpha=alpha)
        end 
    end
    Zeiten, Typ = autoscale_time(minimum(Times), maximum(Int64, Times), intervalls=3)
    for i in (1:6)
        Ax[i].xticks[] = (datetime2unix.(Zeiten.dates) , Dates.format.(Zeiten.dates, Typ)); 
    end
    close(file)
    return figure
end


function ToolBox.plot_DomDataV3_PMT(Dom_object::Tuple{Integer,Integer}; T_intervall::Tuple{Integer,Integer}=(0,0), loadpath::String="../Data/DomData_Doms", alpha::Float64=1.0)
    file = h5open(string(loadpath, "/Data_", Dom_object[1],".h5"), "r")
    figure = Figure()
    axf1 = Axis(figure[1,1], title="frequencies of PMT Ring A", xlabel="Time", ylabel="frequency in Hz") 
    Label(figure[0, :], string("Data of Dom ", Dom_object[1], " PMT ", Dom_object[2]), fontsize = 30)
    Times = read(file["Time"])
    time_mask = maskTime(Times,T_intervall)
    scatter!(axf1, Times[time_mask], read(file["pmtmean"])[time_mask,Dom_object[2]], alpha=alpha)
    Zeiten, Typ = autoscale_time(minimum(Times), maximum(Int64, Times), intervalls=3)
    for i in (1:6)
        axf1.xticks[] = (datetime2unix.(Zeiten.dates) , Dates.format.(Zeiten.dates, Typ)); 
    end
    close(file)
    return figure
end

function ToolBox.plot_DataV3_Events_alt(Event_dict::Dict{Integer,Any}; loadpath::String="../Data/DomData_Doms", temp=true)
    interesting_Doms = collect(keys(Event_dict))
    if temp== true
        pictures = readdir("temp_pictures")
        for file in pictures
            rm(file)
        end
    end
    for Dom in interesting_Doms
        for i in (1:length(Event_dict[Dom]))
            fig = plot_DataV3_Event(Dom, Event_dict[Dom][i], loadpath=loadpath)
            save(string("temp_pictures/HRF_event",Dom,"_",i,".png"), fig)
        end
    end
end

function ToolBox.plot_DataV3_Event_alt(DomID::Integer, event::Tuple{Int64, Int64, Int64}; loadpath::String="../Data/DomData_Doms", save::Bool=false)
    file = h5open(string(loadpath, "/Data_", DomID,".h5"), "r")
    figure = Figure()
    axf1 = Axis(figure[1,1:2], title="frequencies of PMT", xlabel="Time", ylabel="frequency in Hz") 
    axf2 = Axis(figure[2,2], title="frequencies of event", xlabel="Time", ylabel="frequency in Hz")
    axf3 = Axis(figure[3,1], title="good_values of event", xlabel="Time")
    axf4 = Axis(figure[3,2], title="hrvcount of event", xlabel="Time", ylabel="hrvcount")
    Times = read(file["Time"])
    event_intervall = event[2]:event[2]+event[3]+30
    Ax = [axf2, axf3, axf4]
    #paar mehr Daten ausrechnen:
    event_mean = round(mean(read(file["pmtmean"])[event_intervall,event[1]]),digits=3)
    event_max = round(maximum(read(file["pmtmean"])[event_intervall,event[1]]),digits=3)
    Label(figure[2,1], string("Event of Dom ", DomID, ", PMT ", event[1],"\nEventlength: ",event[3]*10, "min\nmean Frequency=",event_mean,"Hz,\nmax Frequency=",event_max,"Hz"), fontsize = 20, tellheight = false, tellwidth = false, justification=:center)
    scatter!(axf2, Times[event_intervall], read(file["pmtmean"])[event_intervall,event[1]])
    #Nur für testzwecke:
    time_mask = [time <= 1664270856 for time in Times]  #die muss wieder raus/anders geschrieben werden
    scatter!(axf1, Times[time_mask], read(file["pmtmean"])[:,event[1]][time_mask])
    scatter!(axf4, Times[event_intervall], read(file["hrvcount"])[event_intervall,event[1]])
    scatter!(axf3, Times[event_intervall], read(file["good_values"])[event_intervall])
    Zeiten, Typ = ToolBox.autoscale_time(minimum(Times[event_intervall]), maximum(Times[event_intervall]), intervalls=3)
    for i in (1:3)
        Ax[i].xticks[] = (datetime2unix.(Zeiten.dates) , Dates.format.(Zeiten.dates, Typ)); 
    end
    Zeiten, Typ = ToolBox.autoscale_time(minimum(Times[time_mask]), maximum(Times[time_mask]))
    axf1.xticks[] = (datetime2unix.(Zeiten.dates) , Dates.format.(Zeiten.dates, Typ)); 
    close(file)
    if save
        save("Bild.png", figure)
    end
    return figure
end

function ToolBox.plot_DomDataV3_Floors(Floor::Int64; T_intervall::Tuple{Integer,Integer}=(0,0), loadpath::String="../Data", alpha::Float64=1.0)
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
    time_mask = maskTime(Times,T_intervall)
    for ring in (1:length(config.Detector_PMT_Ringe))
        for pmt_number in (1:length(config.Detector_PMT_Ringe[ring]))
            scatter!(Ax[ring], Times[time_mask], read(file["pmtmean"])[Floor,:,config.Detector_PMT_Ringe[ring][pmt_number]][time_mask], color=config.Color[pmt_number], alpha=alpha)
        end 
    end
    Zeiten, Typ = autoscale_time(minimum(Times), maximum(Times), intervalls=3)
    for i in (1:6)
        Ax[i].xticks[] = (datetime2unix.(Zeiten.dates) , Dates.format.(Zeiten.dates, Typ)); 
    end
    close(file)
    return figure
end