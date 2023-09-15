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


