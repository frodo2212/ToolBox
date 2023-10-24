"""
takes a lower and upper bound of frequencies and gives back all Doms for which at least one PMT has a specific number of frequencies in between the bounds
this number must lie inside the threshold values 
"""
function Search_DomData(bound::Tuple{Integer, Integer}, threshold::Tuple{Integer, Integer}; loadpath::String="../Data/DomData_Doms", slice_length::Integer=6000, return_Events::Bool=false)
    possible_files = readdir(loadpath)
    possible_files = [files for files in possible_files if files[1:4] == "Dom_" && files[15:19] == string(Int32(slice_length/600),".h5")]
    if return_Events
        interesting_Doms = Event[]
    else
        interesting_Doms = Tuple{Int32,Int32}[]
    end
    for file in possible_files
        Dom = parse(Int32, file[5:13])
        datafile = h5open(string(loadpath, "/", file), "r")
        pmtmeans = read(datafile["pmtmean"])
        close(datafile)
        for pmt in (1:PMT_count)
            high_rate =  [pmtmeans[i,pmt] >= bound[1] && pmtmeans[i,pmt] <= bound[2] for i in (1:size(pmtmeans)[1])]
            Signals = count(high_rate)
            if Signals >= threshold[1] && Signals <= threshold[2]
                if return_Events
                    push!(interesting_Doms, Event("S_t1", Dom, pmt, 0, 0, 0, 0, 0, 0, 0, slice_length))
                else
                    push!(interesting_Doms, (Dom,pmt))
                end
            end
        end
    end
    return interesting_Doms
end

"""
does an intensive Search for all optical Doms of given Detector
"""
function intSearch_DomData(detector::Detector, bound::Tuple{Integer, Integer}; threshold::Integer=6, loadpath::String="../Data/DomData_Doms", slice_length::Integer=6000)
    Dom_ids = optical_DomIds(detector)
    globalEvents = Event[]
    for Dom in Dom_ids
        Events = intSearch_DomData_Dom(Dom, bound, threshold=threshold, loadpath=loadpath, slice_length=slice_length)
        append!(globalEvents, Events)
    end
    return globalEvents
end

"""
Scans a Dom for anomalies 
    returns a Event for every time mind 'threshold' values are inside the 'bound'
"""
function intSearch_DomData_Dom(DomId::Integer, bound::Tuple{Integer, Integer}; threshold::Integer=6, loadpath::String="../Data/DomData_Doms", slice_length::Integer=6000, filter_empty_events::Bool=true)
    datafile = h5open(string(loadpath, "/Dom_", DomId,"_",Int32(slice_length/600),".h5"), "r") 
    pmtmeans = read(datafile["pmtmean"])
    timestamps = read(datafile["Time"])
    close(datafile)
    Events = Event[]
    for pmt in (1:PMT_count)
        tmp_events = intSearch_PMT(pmtmeans[:,pmt], bound, threshold)
        for tmp_event in tmp_events
            t_s = timestamps[tmp_event[1]]
            t_e = timestamps[tmp_event[2]] 
            missing_timestamps = round(UInt32, (t_e-t_s)/(slice_length/10) - (tmp_event[2] - tmp_event[1]))
            #alternativ alle timestamps durchgehen und wenn der Abstand von 2 >= 600*1.8 ist, missing_timestamps += 1
            if !filter_empty_events || missing_timestamps < tmp_event[3]*2
                push!(Events, Event("iS_t1", DomId, pmt, tmp_event[1], tmp_event[2], tmp_event[3], t_s, t_e, t_e-t_s, missing_timestamps, slice_length))
            end
        end
    end
    return Events
end

"""
primarily a inside function, not for external use
Does the search for anomalies on a given Vector of pmtmeans
"""
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
does a linfit over either a single PMT given (DomId,PMT) or the mean frequencies of a Dom given DomId 
returns either the parameters or a finished function  -  return_function=false
one can restrict the time intervall via T_intervall
"""
function linfit_DomData(Dom_object::Tuple{Integer,Integer}; T_intervall::Tuple{Integer,Integer}=(0,0), loadpath::String="../Data/DomData_Doms", return_function::Bool=false, slice_length::Integer=54000, ignore_highs::Bool=true, additional_Data::Bool=false)
    datafile = h5open(string(loadpath, "/Dom_", Dom_object[1],"_",Int32(slice_length/600),".h5"), "r")    
    pmtmean = read(datafile["pmtmean"])[:,Dom_object[2]]
    Times = read(datafile["Time"])
    close(datafile)
    time_mask = ToolBox.maskTime(Times, T_intervall)
    fit_Times = Times[time_mask] .- Times[time_mask][1]
    fit_pmtmean = pmtmean[time_mask]
    if ignore_highs 
        deleteat!(fit_Times, findall(x->x>=20000,fit_pmtmean))
        deleteat!(fit_pmtmean, findall(x->x>=20000,fit_pmtmean))
    end
    nanmask = maskinfnan(fit_pmtmean) 
    if count(nanmask) > 10 #a threshold how many values i need to count the linfit
        gerade(t, p) = p[1] .+ (p[2] .* t)
        p0 = [fit_pmtmean[nanmask][1], 0]
        fit = curve_fit(gerade, fit_Times[nanmask], fit_pmtmean[nanmask], p0)
        y = fit.param[1]- fit_Times[nanmask][1]*fit.param[2]
        if additional_Data
            return (y,fit.param[2]), (Times[1],Times[length(Times)])
        elseif return_function
            return f(x)=y+fit.param[2]*x
        else
            return (y,fit.param[2])
        end
    elseif additional_Data
        return (0,0), (0,0)
    end
    return (0,0)
end

function linfit_DomData(Dom_Id::Integer; T_intervall::Tuple{Integer,Integer}=(0,0), loadpath::String="../Data/DomData_Doms", return_function::Bool=false, slice_length::Integer=54000, ignore_highs::Bool=true, additional_Data::Bool=false)
    datafile = h5open(string(loadpath, "/Dom_", Dom_Id,"_",Int32(slice_length/600),".h5"), "r")    
    pmtmean = vec(mean(read(datafile["pmtmean"]), dims=2))
    Times = read(datafile["Time"])
    close(datafile)
    time_mask = ToolBox.maskTime(Times, T_intervall)
    fit_Times = Times[time_mask] .- Times[time_mask][1]
    fit_pmtmean = pmtmean[time_mask]
    if ignore_highs 
        deleteat!(fit_Times, findall(x->x>=20000,fit_pmtmean))
        deleteat!(fit_pmtmean, findall(x->x>=20000,fit_pmtmean))
    end
    nanmask = maskinfnan(fit_pmtmean) 
    if count(nanmask) > 10 #a threshold how many values i need to count the linfit
        gerade(t, p) = p[1] .+ (p[2] .* t)
        p0 = [fit_pmtmean[nanmask][1], 0]
        fit = curve_fit(gerade, fit_Times[nanmask], fit_pmtmean[nanmask], p0)
        y = fit.param[1]- fit_Times[nanmask][1]*fit.param[2]
        if additional_Data
            return (y,fit.param[2]), (Times[1],Times[length(Times)])
        elseif return_function
            return f(x)=y+fit.param[2]*x
        else
            return (y,fit.param[2])
        end
    elseif additional_Data
        return (0,0), (0,0)
    end
    return (0,0)
end

"""
does a liner fit over all doms of given detector.
either over the mean frequencies of the Doms or one specific pmt number of all Doms 
"""
function linfit_DomData(detector::Detector; pmt::Integer=0, T_intervall::Tuple{Integer,Integer}=(0,0), loadpath::String="../Data/DomData_Doms", slice_length::Integer=54000, ignore_highs::Bool=true, only_slopes::Bool=false)
    Dom_Ids = optical_DomIds(detector)
    if only_slopes 
        Data_dict = Dict{Int32,Float64}
    else
        Data_dict = Dict{Int32,Tuple{Float64,Float64}}
    end
    for Dom in Dom_Ids
        if pmt == 0 || pmt > 31 #print an error if pmt > 31
            params = linfit_DomData(Dom, T_intervall=T_intervall, loadpath=loadpath, slice_length=slice_length, ignore_highs=ignore_highs)
        else
            params = linfit_DomData((Dom,pmt), T_intervall=T_intervall, loadpath=loadpath, slice_length=slice_length, ignore_highs=ignore_highs)
        end
        if params != (0,0)
            if only_slopes
                Data_dict = merge!(Data_dict, Dict(Dom=>params[2]))
            else
                Data_dict = merge!(Data_dict, Dict(Dom=>params))
            end
        end
    end
    return Data_dict
end
"""
does a fit for all opticalDoms of given Detector and all its pmts over the whole time intervall
returns the Data in a Vector of linfitData structs
"""
function Search_linFitData(detector::Detector, threshold::Float64; mean_pmts::Bool=false, range::Tuple{Integer,Integer}=(1,378), loadpath::String="../Data/DomData_Doms", slice_length::Integer=54000)
    Dom_Ids = optical_DomIds(detector)
    Events = linfitData[]
    for Dom in Dom_Ids[range[1]:range[2]]
        if isfile(string(loadpath, "/Dom_", Dom,"_",Int32(slice_length/600),".h5"))
            append!(Events, inner_Search_linFitData(Dom, threshold, mean_pmts=mean_pmts, loadpath=loadpath, Typ="S_lf_t2", slice_length=slice_length))
        end
    end
    return Events
end

"""
primarily an inside function not for external use
takes a DomId and does a linear fit over the whole time intervall for all pmts
"""
function inner_Search_linFitData(Dom_Id::Int32, threshold::Float64; mean_pmts::Bool=false, loadpath::String="../Data/DomData_Doms", Typ::String="inner_lf_t2", slice_length::Integer=6000)
    Event = linfitData[]
    if mean_pmts
        params, data = linfit_DomData(Dom_Id,loadpath=loadpath, slice_length=slice_length, additional_Data=true)
        if params[2] > threshold || params[2] < -threshold
            time = round(UInt32, (data[2]-data[1])/2)
            time_intervall = UInt32(data[2]-data[1])
            push!(Event, linfitData(Typ, Dom_Id, 0, time, time_intervall, params, 1, slice_length))
        end
    else
        for pmt in (1:PMT_count)
            params, data = linfit_DomData((Dom_Id,pmt),loadpath=loadpath, slice_length=slice_length, additional_Data=true)
            if params[2] > threshold || params[2] < -threshold
                time = round(UInt32, (data[2]-data[1])/2)
                time_intervall = UInt32(data[2]-data[1])
                push!(Event, linfitData(Typ, Dom_Id, pmt, time, time_intervall, params, 1, slice_length))
            end
        end
    end
    return Event
end


"""
given a Dom_Id it does a linear fit over all its pmts, add a specific pmt and you just get the Data for the pmt
possible arguments: Intervall, step (both in minutes)
"""
function linfit_DomData_intervalls(Dom_Id::Int32; pmt::Integer=0, T_intervall::Tuple{Integer,Integer}=(0,0), Intervall::Integer=500, step::Integer=250, loadpath::String="../Data/DomData_Doms", slice_length::Integer=6000, ignore_highs::Bool=true, return_vector::Bool=false)
    datafile = h5open(string(loadpath, "/Dom_", Dom_Id,"_",Int32(slice_length/600),".h5"), "r") 
    pmtmean = read(datafile["pmtmean"])
    Times = read(datafile["Time"])
    close(datafile)
    if T_intervall == (0,0)
        Timesteps = collect(range(minimum(Times)+Intervall*30,maximum(Times)-Intervall*30 ,step=step*60))
    else 
        Timesteps = collect(range(T_intervall[1]+Intervall*30,T_intervall[2]-Intervall*30 ,step=step*60))
    end
    if return_vector
        if pmt == 0 || pmt < PMT_count
            params1 = Matrix{Float64}(undef, PMT_count, length(Timesteps))
            params2 = Matrix{Float64}(undef, PMT_count, length(Timesteps))
            for pmt in (1:PMT_count)
                params1[pmt,:], params2[pmt,:] = inner_linfit_intervalls(pmtmean[:,pmt],Times,Timesteps, (Dom_Id, pmt), slice_length=slice_length,Intervall=Intervall, ignore_highs=ignore_highs, return_vector=true)
            end
            return params1, params2, Timesteps
        end
        return inner_linfit_intervalls(pmtmean[:pmt],Times,Timesteps, (Dom_Id, pmt), slice_length=slice_length,Intervall=Intervall, ignore_highs=ignore_highs, return_vector=true), Timesteps
    end
    Data = linfitData[]
    if pmt == 0 || pmt < PMT_count
        for pmt in (1:PMT_count)
            append!(Data, inner_linfit_intervalls(pmtmean[:,pmt],Times,Timesteps, (Dom_Id, pmt), slice_length=slice_length, Intervall=Intervall, ignore_highs=ignore_highs))
        end
        return Data
    end
    return inner_linfit_intervalls(pmtmean[:pmt],Times,Timesteps, (Dom_Id, pmt), slice_length=slice_length, Intervall=Intervall, ignore_highs=ignore_highs)
end

function store_linfitData_intervalls(detector::Detector; Doms=Int32[], T_intervall::Tuple{Integer,Integer}=(0,0), Intervall::Integer=500, step::Integer=250, Name::String="lfData_Intervalls", storagepath::String=".", loadpath::String="../Data/DomData_Doms", slice_length::Integer=6000, ignore_highs::Bool=true)
    file = h5open(string(storagepath,"/",Name,".h5"), "w")
    if Doms==Int32[]
        Doms = optical_DomIds(detector)
    end
    write(file, "Intervall", Intervall)
    write(file, "step", step)
    write(file, "ignore_highs", ignore_highs)
    write(file, "slice_length", slice_length)
    for Dom_Id in Doms
        fit_params1, fit_params2, Timesteps = linfit_DomData_intervalls(Dom_Id, T_intervall=T_intervall, Intervall=Intervall, step=step, loadpath=loadpath, slice_length=slice_length, ignore_highs=ignore_highs, return_vector=true) 
        Dom = string(Dom_Id)
        create_group(file, Dom)
        len = length(Timesteps)
        Tsteps = create_dataset(file[Dom], "Timesteps", Int64, len)  
        params1 = create_dataset(file[Dom], "params1", Float64, PMT_count, len)  
        params2 = create_dataset(file[Dom], "params2", Float64, PMT_count, len)   
        write(Tsteps, Timesteps)
        write(params1, fit_params1)
        write(params2, fit_params2)
    end
    close(file)
end

function load_linFitData_intervalls(Dom_Id::Int32; pmt::Integer=0, Name::String="lfData_Intervalls", storagepath::String=".", return_Dict::Bool=false, metadata::Bool=false, return_lfD::Bool=false, deletezeros::Bool=true)
    file = h5open(string(storagepath,"/",Name,".h5"), "r")
    possible_Doms = collect(keys(read(file)))
    deleteat!(possible_Doms, findall(x->x in ["Intervall","step","ignore_highs","slice_length"],possible_Doms))
    if Dom_Id in parse.(Int32, possible_Doms)
        params1 = read(file[string(Dom_Id)], "params1")
        params2 = read(file[string(Dom_Id)], "params2")
        Timesteps = read(file[string(Dom_Id)], "Timesteps")
        step = read(file, "step")
        slice_length = read(file, "slice_length")
        Intervall = read(file, "Intervall")
        ignore_highs = read(file, "ignore_highs")
        close(file)
        if pmt > 0 && pmt <=31
            params = [(params1[pmt,j],params2[pmt,j]) for j in (1:length(Timesteps))] 
            if deletezeros 
            deleteat!(Timesteps, findall(x->x==(0,0),params))
            deleteat!(params, findall(x->x==(0,0),params))
            end
        else 
            params = [(params1[i,j],params2[i,j]) for i in (1:PMT_count), j in (1:length(Timesteps))]
            if deletezeros
                zeros = findall(x->x==0,params1[1,:])
                deleteat!(Timesteps, zeros)
                params = params[:, setdiff(1:end, zeros)]
            end
        end
        if return_Dict
            metadata && (return Dict(Timesteps[i]=>params[i] for i in (1:length(Timesteps))), (Intervall, step, ignore_highs,slice_length))
            return Dict(Timesteps[i]=>params[i] for i in (1:length(Timesteps)))
        end
        if return_lfD
            lfData = linfitData[]
            if pmt > 0 && pmt <=31 
                for pmt in (1:PMT_count)
                    for i in (1:length(Timesteps))
                        push!(lfData, linfitData("lf_int", Dom_Id, pmt, Timesteps[i], Intervall, (params1[pmt,i],params2[pmt,i]),0, slice_length))
                    end
                end
            else 
                for i in (1:length(Timesteps))
                    push!(lfData, linfitData("lf_int", Dom_Id, pmt, Timesteps[i], Intervall, (params1[pmt,i],params2[pmt,i]),0, slice_length))
                end
            end
            metadata && (return lfData, (Intervall, step, ignore_highs,slice_length))
            return lfData
        end
        metadata && (return Timesteps, params, (Intervall, step, ignore_highs,slice_length))
        return Timesteps, params
    end
    return 0
end

"""
primarily an inside function not for external use
takes Vectors for pmtmean and Time and does the intervall linear fit on those
"""
function inner_linfit_intervalls(pmtmean::Vector{Float64}, Times::Vector{Int64}, Timesteps::Vector{Int64}, Dom_object::Tuple{Integer,Integer};slice_length::Integer=6000, Intervall::Integer=120, Typ::String="inner_lf", ignore_highs::Bool=true, return_vector::Bool=false)
    points_per_slice = Intervall/(slice_length/600)
    Events = linfitData[]
    params1 = Float64[]
    params2 = Float64[]
    for i in Timesteps
        time_mask = ToolBox.maskTime(Times, (i-Intervall*30,i+Intervall*30)) #*30 für die umrechnung in sek eigl *60*0.5
        rel_values = count(time_mask)/points_per_slice
        if count(time_mask) > 0
            fit_Times = Times[time_mask] .- Times[time_mask][1]
            fit_pmtmean = pmtmean[time_mask]
            if ignore_highs 
                deleteat!(fit_Times, findall(x->x>=20000,fit_pmtmean))
                deleteat!(fit_pmtmean, findall(x->x>=20000,fit_pmtmean))
            end
            nanmask = maskinfnan(fit_pmtmean) #ist das wichtig? scheint ja nichts zu helfen
            if count(nanmask) > 0
                gerade(t, p) = p[1] .+ (p[2] .* t)
                #m0 = (fit_pmtmean[nanmask][length(fit_pmtmean[nanmask])]-fit_pmtmean[nanmask][1])/fit_Times[nanmask][length(fit_Times[nanmask])]
                p0 = [fit_pmtmean[nanmask][1], 0]
                fit = curve_fit(gerade, fit_Times[nanmask], fit_pmtmean[nanmask], p0)
                y = fit.param[1]- fit_Times[nanmask][1]*fit.param[2]
                push!(params1, y)
                push!(params2, fit.param[2])
                push!(Events, linfitData(Typ, Dom_object[1], Dom_object[2], i, Intervall, (y,fit.param[2]),rel_values, slice_length))
            else
                push!(params1, 0)
                push!(params2, 0)
            end
        else
            push!(params1, 0)
            push!(params2, 0)
        end
    end
    if return_vector
        return params1, params2
    end
    return Events
end

# zu rechenintensiv, geht ja erwig 
# da muss man sich erst was überlegen, dass man nur sinnvolle Doms/Pmts sucht und die dann fittet
"""
does a linfit_intervalls over all optical doms of the Detector
"""
function linfit_DomData_intervalls(detector::Detector; T_intervall::Tuple{Integer,Integer}=(0,0), Intervall::Integer=120, step::Integer=60, loadpath::String="../Data/DomData_Doms", slice_length=6000, ignore_highs::Bool=true)
    Doms = optical_DomIds(detector)
    Data = linfitData[]
    for Dom in Doms
        append!(Data, linfit_DomData_intervalls(Dom,T_intervall=T_intervall, Intervall=Intervall, step=step, loadpath=loadpath, slice_length=slice_length, ignore_highs=ignore_highs))
    end
    return Data
end

"""
does a Intervall search of one Dom and checks for every Intervall the slope and gives back the values grater than the threshold
"""
# vllt auch hier f_mean,f_min,f_max hinzufügen?
function Search_linFitData_intervalls(Dom_Id::Int32; threshold::Float64=0.15, loadpath::String="../Data/DomData_Doms", Intervall::Integer=300, step::Integer=100, values_threshold::Float64=0.45, prefilter::Bool=true, T_intervall::Tuple{Integer,Integer}=(0,0), slice_length::Integer=6000)
    Event = linfitData[]
    if isfile(string(loadpath, "/Dom_", Dom_Id,"_",Int32(slice_length/600),".h5"))
        datafile = h5open(string(loadpath, "/Dom_", Dom_Id,"_",Int32(slice_length/600),".h5"), "r") 
        pmtmean = read(datafile["pmtmean"])
        Times = read(datafile["Time"])
        close(datafile)
        if T_intervall == (0,0)
            Timesteps = collect(range(minimum(Times)+Intervall*30,maximum(Times)-Intervall*30 ,step=step*60))
        else 
            Timesteps = collect(range(T_intervall[1]+Intervall*30,T_intervall[2]-Intervall*30 ,step=step*60))
        end
        for pmt in (1:PMT_count)
            tmp_Events = inner_linfit_intervalls(pmtmean[:,pmt],Times,Timesteps, (Dom_Id, pmt), slice_length=slice_length,Intervall=Intervall, Typ="lf_t1")
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



"""
returns a Dictionary of all Doms => their mean frequencies (mean over all Times and PMTs)
    one can additionaly give a Day, then it only means over all PMTs and given Day
"""
function DomData_Dom_mean(detector::Detector; day::Tuple{Integer,Integer,Integer}=(1980,1,1), slice_length::Integer=54000, loadpath::String="../Data/DomData_Doms", filternans::Bool=true)
    Doms = optical_DomIds(detector)
    Dom_means = Dict{Int32, Float64}
    date = datetime2unix(DateTime_autocorrect(day))
    datapoints = Vector{Int64}(undef, length(Doms))
    for Dom in Doms
        datafile = h5open(string(loadpath, "/Dom_", Dom,"_",Int32(slice_length/600),".h5"), "r") 
        if day != (1980,1,1)
            Timemask = maskTime(read(datafile["Time"]),(date,date+86400))
            value = mean(ToolBox.filternan(vec(read(datafile["pmtmean"])[Timemask,:])))
            datapoints[findfirst(x->x==Dom,Doms)] = count(Timemask)
        else 
            value = mean(ToolBox.filternan(vec(read(datafile["pmtmean"]))))
        end
        close(datafile)
        if !filternans || !isnan(value)
            Dom_means = merge!(Dom_means,Dict(Dom=>value))
        end
    end
    if day != (1980,1,1)
        return Dom_means, datapoints
    end
    return Dom_means
end

"""
returns a Dictionary of Doms => Vector of all mean values of each PMT
if given a Day additionaly it only means every PMT over this one day
"""
function DomData_PMT_mean(detector::Detector; day::Tuple{Integer,Integer,Integer}=(1980,1,1), slice_length::Integer=54000, loadpath::String="../Data/DomData_Doms", filternans::Bool=true)
    Doms = optical_DomIds(detector)
    Dom_means = Dict{Int32, Vector{Float64}}
    date = datetime2unix(DateTime_autocorrect(day))
    datapoints = Vector{Int64}(undef, length(Doms))
    for Dom in Doms
        datafile = h5open(string(loadpath, "/Dom_", Dom,"_",Int32(slice_length/600),".h5"), "r") 
        if day != (1980,1,1)
            Timemask = maskTime(read(datafile["Time"]),(date,date+86400))
            means = read(datafile["pmtmean"])[Timemask,:]
            datapoints[findfirst(x->x==Dom,Doms)] = count(Timemask)
        else 
            means = read(datafile["pmtmean"])
        end
        close(datafile)
        value = Float64[]
        for pmt in (1:PMT_count)
            push!(value, mean(ToolBox.filternan(means[:,pmt])))         
        end
        if !filternans || count(isnan.(value))==0
            Dom_means = merge!(Dom_means,Dict(Dom=>value))
        end
    end
    if day != (1980,1,1)
        return Dom_means, datapoints
    end
    return Dom_means
end