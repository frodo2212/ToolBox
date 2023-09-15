function plot_DomDataV3_Rings(DomID::Int; T_max::Int64=0, loadpath::String="../Data/DomData_Doms", alpha::Float64=1.0)
    file = h5open(string(loadpath, "/Data_", DomID,".h5"), "r")
    figure = Figure()
    axf1 = Axis(figure[1,1], title="frequencies of PMT Ring A", xlabel="Time", ylabel="frequency in Hz") #vllt positionierung ändern
    axf2 = Axis(figure[2,1], title="frequencies of PMT Ring B", xlabel="Time", ylabel="frequency in Hz")
    axf3 = Axis(figure[3,1], title="frequencies of PMT Ring C", xlabel="Time", ylabel="frequency in Hz")
    axf4 = Axis(figure[1,2], title="frequencies of PMT Ring D", xlabel="Time", ylabel="frequency in Hz")
    axf5 = Axis(figure[2,2], title="frequencies of PMT Ring E", xlabel="Time", ylabel="frequency in Hz")
    axf6 = Axis(figure[3,2], title="frequencies of PMT Ring F", xlabel="Time", ylabel="frequency in Hz")
    Ax = [axf1, axf2, axf3, axf4, axf5, axf6]
    Label(figure[0, :], string("Data of Dom ", DomID), fontsize = 30)
    Times = read(file["Time"])
    if T_max != 0
        time_mask = [time <= T_max for time in Times]
    else 
        time_mask = [time != T_max for time in Times]
    end
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

function plot_DomDataV3_PMT(DomID::Integer, PMT::Integer; T_max::Int64=0, loadpath::String="../Data/DomData_Doms", alpha::Float64=1.0)
    file = h5open(string(loadpath, "/Data_", DomID,".h5"), "r")
    figure = Figure()
    axf1 = Axis(figure[1,1], title="frequencies of PMT Ring A", xlabel="Time", ylabel="frequency in Hz") #vllt positionierung ändern
    Label(figure[0, :], string("Data of Dom ", DomID, " PMT ", PMT), fontsize = 30)
    Times = read(file["Time"])
    if T_max != 0
        time_mask = [time <= T_max for time in Times]
    else 
        time_mask = [time != T_max for time in Times]
    end
    scatter!(axf1, Times[time_mask], read(file["pmtmean"])[time_mask,PMT], alpha=alpha)
    Zeiten, Typ = autoscale_time(minimum(Times), maximum(Int64, Times), intervalls=3)
    for i in (1:6)
        axf1.xticks[] = (datetime2unix.(Zeiten.dates) , Dates.format.(Zeiten.dates, Typ)); 
    end
    close(file)
    return figure
end

function plot_DomDataV3_PMT(Object::Tuple{Integer,Integer}; T_max::Int64=0, loadpath::String="../Data/DomData_Doms", alpha::Float64=1.0)
    file = h5open(string(loadpath, "/Data_", Object[1],".h5"), "r")
    figure = Figure()
    axf1 = Axis(figure[1,1], title="frequencies of PMT Ring A", xlabel="Time", ylabel="frequency in Hz") #vllt positionierung ändern
    Label(figure[0, :], string("Data of Dom ", Object[1], " PMT ", Object[2]), fontsize = 30)
    Times = read(file["Time"])
    if T_max != 0
        time_mask = [time <= T_max for time in Times]
    else 
        time_mask = [time != T_max for time in Times]
    end
    scatter!(axf1, Times[time_mask], read(file["pmtmean"])[time_mask,Object[2]], alpha=alpha)
    Zeiten, Typ = autoscale_time(minimum(Times), maximum(Int64, Times), intervalls=3)
    for i in (1:6)
        axf1.xticks[] = (datetime2unix.(Zeiten.dates) , Dates.format.(Zeiten.dates, Typ)); 
    end
    close(file)
    return figure
end

function plot_DataV3_Events(Event_dict::Dict{Integer,Any}; loadpath::String="../Data/DomData_Doms", alpha::Float64=1.0, temp=true)
    interesting_Doms = collect(keys(Event_dict))
    if temp== true
        pictures = readdir("temp_pictures")
        for file in pictures
            rm(file)
        end
    end
    for Dom in interesting_Doms
        for i in (1:length(Event_dict[Dom]))
            fig = plot_DataV3_Event(Dom, Event_dict[Dom][i], loadpath=loadpath, alpha=alpha)
            save(string("temp_pictures/HRF_event",Dom,"_",i,".png"), fig)
        end
    end
end

#funktioniert, aber bedarf noch verschönerung
function plot_DataV3_Event(DomID::Integer, event::Tuple{Int64, Int64, Int64}; loadpath::String="../Data/DomData_Doms", alpha::Float64=1.0)
    file = h5open(string(loadpath, "/Data_", DomID,".h5"), "r")
    figure = Figure()
    axf1 = Axis(figure[1,1:2], title="frequencies of of Event", xlabel="Time", ylabel="frequency in Hz") #vllt positionierung ändern
    axf2 = Axis(figure[2,1], title="frequencies PMT", xlabel="Time", ylabel="frequency in Hz")
    axf3 = Axis(figure[2,2], title="hrvcount of PMT", xlabel="Time", ylabel="frequency in Hz")
    axf4 = Axis(figure[3,1], title="good_values", xlabel="Time", ylabel="frequency in Hz")
    axf5 = Axis(figure[3,2], title="fifocount of PMT", xlabel="Time", ylabel="frequency in Hz")
    Times = read(file["Time"])
    event_intervall = event[2]:event[2]+event[3]+30
    Ax = [axf2, axf3, axf4, axf5]
    scatter!(axf2, Times[event_intervall], read(file["pmtmean"])[event_intervall,event[1]])
    #Nur für testzwecke:
    time_mask = [time <= 1664270856 for time in Times]  #die muss wieder raus/anders geschrieben werden
    scatter!(axf1, Times[time_mask], read(file["pmtmean"])[:,event[1]][time_mask])
    scatter!(axf3, Times[event_intervall], read(file["hrvcount"])[event_intervall,event[1]])
    scatter!(axf4, Times[event_intervall], read(file["good_values"])[event_intervall])
    scatter!(axf5, Times[event_intervall], read(file["fifocount"])[event_intervall,event[1]])
    Zeiten, Typ = ToolBox.autoscale_time(minimum(Times[event_intervall]), maximum(Times[event_intervall]), intervalls=3)
    for i in (1:4)
        Ax[i].xticks[] = (datetime2unix.(Zeiten.dates) , Dates.format.(Zeiten.dates, Typ)); 
    end
    Zeiten, Typ = ToolBox.autoscale_time(minimum(Times[time_mask]), maximum(Times[time_mask]))
    axf1.xticks[] = (datetime2unix.(Zeiten.dates) , Dates.format.(Zeiten.dates, Typ)); 
    close(file)
    return figure
end

function plot_DomDataV3_Floors(Floor::Int64; T_max::Int64=0, loadpath::String="../Data", alpha::Float64=1.0)
    file = h5open(string(loadpath, "/DomDataV3_Floors.h5"), "r")
    figure = Figure()
    axf1 = Axis(figure[1,1], title="frequencies of PMT Ring A", xlabel="Time", ylabel="frequency in Hz") #vllt positionierung ändern
    axf2 = Axis(figure[2,1], title="frequencies of PMT Ring B", xlabel="Time", ylabel="frequency in Hz")
    axf3 = Axis(figure[3,1], title="frequencies of PMT Ring C", xlabel="Time", ylabel="frequency in Hz")
    axf4 = Axis(figure[1,2], title="frequencies of PMT Ring D", xlabel="Time", ylabel="frequency in Hz")
    axf5 = Axis(figure[2,2], title="frequencies of PMT Ring E", xlabel="Time", ylabel="frequency in Hz")
    axf6 = Axis(figure[3,2], title="frequencies of PMT Ring F", xlabel="Time", ylabel="frequency in Hz")
    Ax = [axf1, axf2, axf3, axf4, axf5, axf6]
    Label(figure[0, :], string("Data of Floor ", Floor), fontsize = 30)
    Times = read(file["Time"])
    if T_max != 0
        time_mask = [time <= T_max for time in Times]
    else 
        time_mask = [time != T_max for time in Times]
    end
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