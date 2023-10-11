"""
takes a lower and upper bound of frequencies and gives back all Doms for which at least one PMT has a specific number of frequencies in between the bounds
this number must lie inside the threshold values 
"""
function Search_DomData(bound::Tuple{Integer, Integer}, threshold::Tuple{Integer, Integer};loadpath::String="../Data/DomData_Doms")
    possible_files = readdir(loadpath)
    possible_files = [files for files in possible_files if files[1:5] == "Data_"]
    interesting_Doms = Int64[]
    for file in possible_files
        datafile = h5open(string(loadpath, "/", file), "r")
        Dom = parse(Int32, file[6:14])
        pmtmeans = read(datafile["pmtmean"])
        high_rate =  [pmtmeans[i,j] >= bound[1] && pmtmeans[i,j] <= bound[2] for j in (1:PMT_count) for i in (1:size(pmtmeans)[1])]
        Signals = count(high_rate)
        if Signals >= threshold[1] && Signals <= threshold[2]
            push!(interesting_Doms, Dom)
        end
        close(datafile)
    end
    return (interesting_Doms)
end

function intensiveSearch_DomData_test2(detector, bound::Tuple{Integer, Integer}; threshold::Integer=6, loadpath::String="../Data/DomData_Doms")
    Dom_ids = optical_DomIds(detector)
    globalEvents = Event[]
    for Dom in Dom_ids
        Events = intensiveSearch_DomData_test1(Dom, bound, threshold=threshold, loadpath=loadpath)
        append!(globalEvents, Events)
    end
    return globalEvents
end

# Die erstellt Events 
#  es müssen mind 'threshold' pmtmeans aufeinanderfolgend innerhalb der bounds sein 
#TODO: erweiterungsvorschläge:
#       Funktion die bounds lokal über mean erstellt
function intensiveSearch_DomData_test1(DomId::Integer, bound::Tuple{Integer, Integer}; threshold::Integer=6, loadpath::String="../Data/DomData_Doms")
    datafile = h5open(string(loadpath, "/Data_", DomId, ".h5"), "r")
    pmtmeans = read(datafile["pmtmean"])
    timestamps = read(datafile["Time"])
    slice_length = 6000 #read(datafile["slice_length"])
    close(datafile)
    Events = Event[]
    for pmt in (1:PMT_count)
        tmp_events = intSearch_PMT(pmtmeans[:,pmt], bound, threshold)
        for tmp_event in tmp_events
            t_s = timestamps[tmp_event[1]]
            t_e = timestamps[tmp_event[2]] 
            missing_timestamps = round(UInt32, (t_e-t_s)/(slice_length/10) - (tmp_event[2] - tmp_event[1]))
            #alternativ alle timestamps durchgehen und wenn der Abstand von 2 >= 600*1.8 ist, missing_timestamps += 1
            push!(Events, Event("iS_t1", DomId, pmt, tmp_event[1], tmp_event[2], tmp_event[3], t_s, t_e, t_e-t_s, missing_timestamps))
        end
    end
    return Events
end

function intSearch_PMT(pmtmeans::Vector{Float64}, bound::Tuple{Integer, Integer}, threshold::Integer)
    len = length(pmtmeans)
    Events = Tuple{UInt32,UInt32, UInt32}[]
    high_rate = [pmtmeans[i] >= bound[1] && pmtmeans[i] <= bound[2] for i in (1:len)]
    absolute = count(high_rate)
    if absolute <= 250  #das ist ein threshold der die die an sich zu hoch sind rausfiltern soll. das kann man aber besser machen
        # das muss auf die länge der pmtmean abgestimmt sein und irgendwie auch den mittelwert beinhalten
        active_event = false
        event_length::UInt32 = 0
        event_start::UInt32 = 0
        for i in (1:len)
            if high_rate[i]
                if active_event   
                    event_length += 1
                else 
                    active_event = true
                    event_start = i
                    event_length = 1
                end
            end
            if !high_rate[i] && active_event
                if !high_rate[i-1]
                    active_event = false
                    if event_length >= threshold 
                        event_end = event_start+event_length-1
                        push!(Events, (event_start,event_end,event_length-1))
                    end
                else
                    event_length += 1
                end
            end
        end
    end
    return Events
end







"""
does a linfit over a single PMT given (DomId,PMT)
returns either the parameters or a finished function  -  return_function=false
one can restrict the time intervall via T_intervall
"""
function linfit_DomData(Dom_object::Tuple{Integer,Integer}; T_intervall::Tuple{Integer,Integer}=(0,0), loadpath::String="../Data/DomData_Doms", return_function::Bool=false)
    datafile = h5open(string(loadpath, "/Data_", Dom_object[1],".h5"), "r")    
    pmtmean = read(datafile["pmtmean"])[:,Dom_object[2]]
    Times = read(datafile["Time"])
    time_mask = ToolBox.maskTime(Times, T_intervall)
    Start = minimum(Times[time_mask])
    fit_Times = Times[time_mask]
    fit_Times = fit_Times .- Start
    fit_pmtmean = pmtmean[time_mask]
    gerade(t, p) = p[1] .+ (p[2] .* t)
    p0 = [fit_pmtmean[1], 0]
    fit = curve_fit(gerade, fit_Times, fit_pmtmean, p0)
    y = fit.param[1]- Start*fit.param[2]
    if return_function
        return f(x)=y+fit.param[2]*x
    end
    return (y,fit.param[2])
end

"""
Takes a Tuple (DomID,PMT) 
possible arguments: Intervall, step (both in minutes)
"""
function linfit_DomData_intervalls(Dom_object::Tuple{Integer,Integer}; T_intervall::Tuple{Integer,Integer}=(0,0), Intervall::Integer=120, step::Integer=60, loadpath::String="../Data/DomData_Doms")
    datafile = h5open(string(loadpath, "/Data_", Dom_object[1],".h5"), "r")    
    pmtmean = read(datafile["pmtmean"])[:,Dom_object[2]]
    Times = read(datafile["Time"])
    slice_length = 6000 #read(datafile["slice_length"])
    close(datafile)
    return inner_linfit_intervalls(pmtmean,Times, Dom_object, slice_length=slice_length,T_intervall=T_intervall, Intervall=Intervall, step=step)
end

function inner_linfit_intervalls(pmtmean,Times, Dom_object::Tuple{Integer,Integer};slice_length=600,T_intervall::Tuple{Integer,Integer}=(0,0), Intervall::Integer=120, step::Integer=60, Typ::String="inner_lf")
    if T_intervall == (0,0)
        Timesteps = collect(range(minimum(Times)+Intervall*30,maximum(Times)-Intervall*30 ,step=step*60))
    else 
        Timesteps = collect(range(T_intervall[1]+Intervall*30,T_intervall[2]-Intervall*30 ,step=step*60))
    end
    points_per_slice = Intervall/(slice_length/600)
    Events = linfitData[]
    for i in Timesteps
        time_mask = ToolBox.maskTime(Times, (i-Intervall*30,i+Intervall*30)) #*30 für die umrechnung in sek eigl *60*0.5
        rel_values = count(time_mask)/points_per_slice
        if count(time_mask) > 0
            Start = minimum(Times[time_mask])
            fit_Times = Times[time_mask]
            fit_Times = fit_Times .- Start
            fit_pmtmean = pmtmean[time_mask]
            gerade(t, p) = p[1] .+ (p[2] .* t)
            p0 = [fit_pmtmean[1], 0]
            fit = curve_fit(gerade, fit_Times, fit_pmtmean, p0)
            y = fit.param[1]- Start*fit.param[2]
            push!(Events, linfitData(Typ, Dom_object[1], Dom_object[2], i, Intervall, (y,fit.param[2]),rel_values))
        end
    end
    return Events
end

#die funktionieren schon, sind auf den momentanen Daten aber viel zu langsam - daher nur ausschnitte machen
#TODO: Rückgabeformat verbessern
# vllt noch f_min,f_max hinzufügen - die stechen aber wahrscheinlich raus, eher f_min,f_max vom gefitteten, da sind die Störugen draußen
function Search_linFit_test2(detector; range=(1,40), threshold=0.15, loadpath="../Data/DomData_Doms")
    Dom_Ids = optical_DomIds(detector)
    Events = linfitData[]
    for Dom in Dom_Ids[range[1]:range[2]]
        if isfile(string(loadpath, "/Data_", Dom,".h5"))
            append!(Events, inner_Search_linFit_test2(Dom, threshold=threshold, loadpath=loadpath, Typ="S_lf_t2"))
        end
    end
    return Events
end
function inner_Search_linFit_test2(Dom_Id; threshold=0.15, loadpath="../Data/DomData_Doms", Typ::String="inner_lf_t2")
    Event = linfitData[]
    for pmt in (1:PMT_count)
        params = linfit_DomData((Dom_Id,pmt),loadpath=loadpath)
        if params[2] >= threshold || params[2] <= -threshold
            push!(Event, linfitData(Typ, Dom_Id, pmt, 0, 0, params,0))
        end
    end
    return Event
end


# vllt auch hier f_mean,f_min,f_max hinzufügen?
function Search_linFit_test1(Dom_Id; threshold=0.15, loadpath="../Data/DomData_Doms", Intervall::Integer=300, step::Integer=100, values_threshold = 0.45, prefilter::Bool=true, T_intervall::Tuple{Integer,Integer}=(0,0))
    Event = linfitData[]
    if isfile(string(loadpath, "/Data_", Dom_Id,".h5"))
        datafile = h5open(string(loadpath, "/Data_", Dom_Id,".h5"), "r")    
        pmtmean = read(datafile["pmtmean"])
        Times = read(datafile["Time"])
        slice_length = 6000 #read(datafile["slice_length"])
        close(datafile)
        for pmt in (1:PMT_count)
            tmp_Events = inner_linfit_intervalls(pmtmean[:,pmt],Times, (Dom_Id, pmt), slice_length=slice_length,T_intervall=T_intervall, Intervall=Intervall, step=step, Typ="lf_t1")
            for i in (1:length(tmp_Events))
                if tmp_Events[i].params[2] >= threshold || tmp_Events[i].params[2] <= -threshold
                    if tmp_Events[i].rel_values >= values_threshold || !prefilter
                        push!(Event, tmp_Events[i])
                    end
                end
            end
        end
    end
    return Event
end