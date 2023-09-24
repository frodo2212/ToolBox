function DomDataV3(Dom_Id::Integer, loadpath::String, storagepath::String; files::Vector{String}=String[], slice_length::Int64=6000)
    possible_files = readdir(loadpath) #use glob to check whether the files are useable?
    if files == String[]  #if nothing given use all files
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
            Start, End, Data, used_config = ToolBox.load_DomData_run(file,loadpath)
            if used_config[3] == slice_length
                len = length(Data[Dom_Id][1])
                append!(Times, [(Start+Int64(slice_length/10)*(i-1)) for i in (1:len)])
                append!(good_values, Data[Dom_Id][1])
                pmtmean = [pmtmean;(Data[Dom_Id][2])'] # pmtmean
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
    storage = h5open(string(storagepath,"/Data_",Dom_Id,".h5"), "w")
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
    close(storage)
    # size_time = size(Times)
    # size_pmt = size(pmtmean)
    # size_hrv = size(hrvcount)
    # size_fifo = size(fifocount)
    # size_good = size(good_values)
    # return (size_time, size_good, size_pmt, size_hrv, size_fifo)
    return (wrong_sliceTime, notfound)
end

function allDomDataV3(detector::Detector, loadpath::String, storagepath::String; files::Vector{String}=String[], slice_length::Int64=6000)
    Doms = optical_DomIds(detector)
    for Dom in Doms
        DomDataV3(Dom, loadpath, storagepath, files=files, slice_length=slice_length)
    end
end

function DomDataV3_Floors(detector::Detector, loadpath::String, storagepath::String)
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
                else  #hier müssen additionen über die Doms hin
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

function Search_DomDataV3(bound::Tuple{Integer, Integer}, threshold::Tuple{Integer, Integer};loadpath::String="../Data/DomData_Doms")
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


function linfit_DomDataV3(Dom_object::Tuple{Integer,Integer}; T_intervall::Tuple{Integer,Integer}=(0,0), loadpath::String="../Data/DomData_Doms", return_function::Bool=true)
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

function linfit_DomDataV3_intervalls(Dom_object::Tuple{Integer,Integer}; T_intervall::Tuple{Integer,Integer}=(0,0), Intervall::Integer=120, step::Integer=60, loadpath::String="../Data/DomData_Doms", length=100)
    datafile = h5open(string(loadpath, "/Data_", Dom_object[1],".h5"), "r")    
    pmtmean = read(datafile["pmtmean"])[:,Dom_object[2]]
    Times = read(datafile["Time"])
    close(datafile)
    if T_intervall == (0,0)
        Timesteps = collect(range(minimum(Times)+Intervall*30,maximum(Times)-Intervall*30 ,step=step*60))
    else 
        Timesteps = collect(range(T_intervall[1]+Intervall*30,T_intervall[2]-Intervall*30 ,step=step*60))
    end
    Params = []
    for i in Timesteps
        time_mask = ToolBox.maskTime(Times, (i-Intervall*30,i+Intervall*30))
        push!(Params,linear_fit(Times[time_mask], pmtmean[time_mask]))
    end
    return (Timesteps,Params)
end

#TODO: den Filter verbessern - bisschen willkürlich was ich hier raus filtere
#TODO: wenn der PMT ne Steigung in der Frequenz hat ist frquency_mean nicht optimal... 
function intensiveSearch_DomDataV3(bound::Tuple{Integer, Integer}, threshold::Tuple{Integer, Integer};loadpath::String="../Data/DomData_Doms", length=30)
    int_Doms = ToolBox.Search_DomDataV3(bound, (threshold[1],500), loadpath=loadpath)
    interesting_Doms = Dict{Integer,Any}()
    for Dom in int_Doms
        datafile = h5open(string(loadpath, "/Data_", Dom, ".h5"), "r")
        event_data = Any[]
        pmtmeans = read(datafile["pmtmean"])
        high_rate = zeros(Bool, size(pmtmeans)[1], PMT_count)
        temp_signal = zeros(Int32, size(pmtmeans)[1], PMT_count)
        temp_event = zeros(Bool, size(pmtmeans)[1], PMT_count)
        for i in (1:size(pmtmeans)[1])
            high_rate[i,:] =  [pmtmeans[i,j] >= bound[1] && pmtmeans[i,j] <= bound[2] for j in (1:PMT_count) ]
        end
        for pmt in (1:PMT_count)
            frequency_mean = mean(pmtmeans[:,pmt])
            if bound[1]-frequency_mean >= 2000 #auf was ich die empfindlichkeit setze ist so die frage
                event = 0
                event_length = 0
                event_start = -3
                for time in (1:size(pmtmeans)[1]-length)
                    temp_signal[time] = count(high_rate[time:time+length,pmt])
                    temp_event[time] = temp_signal[time] >= threshold[1] && temp_signal[time] <= threshold[2]
                    if event == 0 && temp_event[time] == 1 
                        event_start = time
                        event = 1
                        event_length = 0
                    elseif temp_event[time] == 1
                        event_length += 1
                    elseif event == 1 && temp_event[time] == 0
                        push!(event_data, (pmt, event_start, event_length))
                        event = 0
                    end
                end
            end
        end
        if event_data != Any[]         
            merge!(interesting_Doms, Dict(Dom=>event_data))
        end
        close(datafile)
    end
    return interesting_Doms
end