function DomData(Dom_Id::Integer, loadpath::String, storagepath::String; files::Vector{String}=String[], slice_length::Integer=6000)
    loadpath = string(loadpath,"/Runs_sl",Int32(slice_length/600),"Min/")
    #possible_files = readdir(loadpath)
    files1 = glob("*_S.root", loadpath) 
    possible_files = [file[findlast("/", file)[1]+1:length(file)] for file in files1]
    if files == String[]
        files = possible_files
    end
    notfound = String[]
    wrong_sliceTime = String[]
    Times = Int64[]
    good_values = Int64[]
    pmtmean = Matrix{Float64}(undef, 0,31)
    hrvcount = Matrix{Int64}(undef, 0,31)
    fifocount = Matrix{Int64}(undef, 0,31)
    for file in files
        if file in possible_files
            Start, End, Data, used_config = ToolBox.load_Data_run(file,loadpath)
            if used_config[3] == slice_length
                len = length(Data[Dom_Id][1])
                append!(Times, [(Start+Int64(slice_length/10)*(i-1)) for i in (1:len)])
                append!(good_values, Data[Dom_Id][1])
                pmtmean = [pmtmean;(Data[Dom_Id][2])'] 
                hrvcount = [hrvcount;(Data[Dom_Id][3])']
                fifocount = [fifocount;(Data[Dom_Id][4])']
            else 
                push!(wrong_sliceTime, file)
            end   
        else 
            push!(notfound, file)
        end
    end 
    len = length(Times)
    mkpath(storagepath)
    storage = h5open(string(storagepath,"/Dom_",Dom_Id,"_",Int32(slice_length/600),".h5"), "w")
    d_time = create_dataset(storage, "Time", Int64, len)
    d_pmtmean = create_dataset(storage, "pmtmean", Float64, len, 31)
    d_hrvcount = create_dataset(storage, "hrvcount", Int64, len, 31)
    d_fifocount = create_dataset(storage, "fifocount", Int64, len, 31)
    d_good_values = create_dataset(storage, "good_values", Int64, len)
    write(d_time, Times)
    write(d_pmtmean, pmtmean)
    write(d_hrvcount, hrvcount)
    write(d_fifocount, fifocount)
    write(d_good_values, good_values)
    write(storage, "slice_length", slice_length)
    close(storage)
    return (wrong_sliceTime, notfound)
end


# mega unnötig, braucht viel zu lange, wenn ich alle will mach ich das über snakemake
# function allDomDataV3(detector::Detector, loadpath::String, storagepath::String; files::Vector{String}=String[], slice_length::Int64=6000)
#     Doms = optical_DomIds(detector)
#     for Dom in Doms
#         DomDataV3(Dom, loadpath, storagepath, files=files, slice_length=slice_length)
#     end
# end

"""
NOTIEREN WAS SIE GENAU MACHT 
UND WIE - ggf OPTIMIEREN    
"""
function DomData_Floors(detector::Detector, loadpath::String, storagepath::String)
    possible_files = readdir(loadpath)
    floors, Doms_on_floor, max_floor = ToolBox.Floors(detector)
    file = h5open(string(loadpath,"/Data_",Doms_on_floor[1][1],".h5"), "r")  #meaga ranzig so, das muss anders gehen
    Times = read(file["Time"])
    close(file)
    len = length(Times)
    good_values = zeros(Int64, max_floor, len)
    pmtmean = Array{Float64}(undef, max_floor, len, 31)
    hrvcount = Array{Int64}(undef, max_floor, len, 31)
    fifocount = Array{Int64}(undef, max_floor, len, 31)
    for floor in (1:max_floor)
        Dom_count = length(Doms_on_floor[floor])
        prov_pmtmean = zeros(Float64, Dom_count, len, 31)
        for id in (1:Dom_count)
            Dom = Doms_on_floor[floor][id]
            if string("Data_",Dom,".h5") in possible_files
                file = h5open(string(loadpath,"/Data_",Dom,".h5"), "r")
                if Times != read(file["Time"])
                    print("Zeiten passen nicht")
                    break
                else  
                    for i in (1:len)
                        good_values[floor,i] += read(file["good_values"])[i]
                        prov_pmtmean[id,i,:] .= read(file["pmtmean"])[i,:]
                        hrvcount[floor,i,:] += read(file["hrvcount"])[i,:]
                        fifocount[floor,i,:] += read(file["fifocount"])[i,:]
                    end
                end
                close(file)
            end
        end
        for time in (1:len)
            pmtmean[floor, time, :] = [mean(ToolBox.filternan(prov_pmtmean[:,time,i])) for i in (1:31)]
        end
    end    
    mkpath(storagepath)
    storagefile = h5open(string(storagepath,"/DomDataV3_floors.h5"), "w")
    dset_time = create_dataset(storagefile, "Time", Int64, len)
    dset_good_values = create_dataset(storagefile, "good_values", Int64, max_floor, len)
    dset_pmtmean = create_dataset(storagefile, "pmtmean", Float64, max_floor, len, 31)
    dset_hrvcount = create_dataset(storagefile, "hrvcount", Int64, max_floor, len, 31)
    dset_fifocount = create_dataset(storagefile, "fifocount", Int64, max_floor, len, 31)
    write(dset_time, Times)
    write(dset_good_values, good_values)
    write(dset_pmtmean, pmtmean)
    write(dset_hrvcount, hrvcount)
    write(dset_fifocount, fifocount)
    close(storagefile)
end

"""
takes a lower and upper bound of frequencies and gives back all Doms for which at least one PMT has a specific number of frequencies in between the bounds
this number must lie inside the threshold values 
"""
function Search_DomData(bound::Tuple{Integer, Integer}, threshold::Tuple{Integer, Integer};loadpath::String="../Data/DomData_Doms")
    possible_files = readdir(loadpath)
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
    close(datafile)
    Params = linear_fit(Times[time_mask], pmtmean[time_mask])
    if return_function
        return f(x)=Params[1]+Params[2]*x
    end
    return Params
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
        time_mask = ToolBox.maskTime(Times, (i-Intervall*30,i+Intervall*30))
        rel_values = count(time_mask)/points_per_slice
        params = linear_fit(Times[time_mask], pmtmean[time_mask])
        push!(Events, linfitData(Typ, Dom_object[1], Dom_object[2], i, Intervall, params,rel_values))
    end
    return Events
end


function intensiveSearch_DomDataV3_test2(detector, bound::Tuple{Integer, Integer}; threshold::Integer=6, loadpath::String="../Data/DomData_Doms")
    Dom_ids = optical_DomIds(detector)
    globalEvents = Event[]
    for Dom in Dom_ids
        Events = intensiveSearch_DomDataV3_test1(Dom, bound, threshold=threshold, loadpath=loadpath)
        append!(globalEvents, Events)
    end
    return globalEvents
end

# Die erstellt Events 
#  es müssen mind 'threshold' pmtmeans aufeinanderfolgend innerhalb der bounds sein 
#TODO: erweiterungsvorschläge:
#       Funktion die bounds lokal über mean erstellt
function intensiveSearch_DomDataV3_test1(DomId::Integer, bound::Tuple{Integer, Integer}; threshold::Integer=6, loadpath::String="../Data/DomData_Doms")
    datafile = h5open(string(loadpath, "/Data_", DomId, ".h5"), "r")
    pmtmeans = read(datafile["pmtmean"])
    timestamps = read(datafile["Time"])
    slice_length = 6000 #read(datafile["slice_length"])
    close(datafile)
    Events = Event[]
    for pmt in (1:ToolBox.config.PMT_count)
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
end


#die funktionieren schon, sind auf den momentanen Daten aber viel zu langsam - daher nur ausschnitte machen
#TODO: Rückgabeformat verbessern
# vllt noch f_min,f_max hinzufügen - die stechen aber wahrscheinlich raus, eher f_min,f_max vom gefitteten, da sind die Störugen draußen
function Search_linFit_test2(detector; range=(1,40), threshold=0.15, loadpath="../Data/DomData_Doms")
    Dom_Ids = optical_DomIds(detector)
    Events = linfitData[]
    for Dom in Dom_Ids[range[1]:range[2]]
        append!(Events, inner_Search_linFit_test2(Dom, threshold=threshold, loadpath=loadpath, Typ="S_lf_t2"))
    end
    return Events
end
function inner_Search_linFit_test2(Dom_Id; threshold=0.15, loadpath="../Data/DomData_Doms", Typ::String="inner_lf_t2")
    Event = linfitData[]
    for pmt in (1:ToolBox.config.PMT_count)
        params = linfit_DomDataV3((Dom_Id,pmt),loadpath=loadpath)
        if params[2] >= threshold || params[2] <= -threshold
            push!(Event, linfitData(Typ, Dom_Id, pmt, 0, 0, params,0))
            
        end
    end
    return Event
end


# vllt auch hier f_mean,f_min,f_max hinzufügen?
function Search_linFit_test1(Dom_Id; threshold=0.15, loadpath="../Data/DomData_Doms", Intervall::Integer=300, step::Integer=100, values_threshold = 0.45, prefilter::Bool=true)
    Event = linfitData[]
    datafile = h5open(string(loadpath, "/Data_", Dom_object[1],".h5"), "r")    
    pmtmean = read(datafile["pmtmean"])[:,Dom_object[2]]
    Times = read(datafile["Time"])
    slice_length = 6000 #read(datafile["slice_length"])
    close(datafile)
    for pmt in (1:ToolBox.config.PMT_count)
        tmp_Events = inner_linfit_intervalls(pmtmean,Times, (Dom_Id, pmt), slice_length=slice_length,T_intervall=T_intervall, Intervall=Intervall, step=step, Typ="lf_t1")
        for i in (1:length(tmp_Events))
            if tmp_Events[i].params[2] >= threshold || tmp_Events[i].params[2] <= -threshold
                if tmp_Events[i].rel_values >= values_threshold || !prefilter
                    push!(Event, tmp_Events[i])
                end
            end
        end
    end
    return Event
end

struct linfitData
    Typ::String
    Dom_Id::UInt32
    pmt::UInt32
    time::UInt32
    time_intervall::UInt32 
    params::Tuple{Float64,Float64}
    rel_values::Float64
end
