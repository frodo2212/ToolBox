function Data(loadpath::String, Run::Integer, detector::KM3io.Detector, storagepath::String; slice_length=6000, slice_length_threshold=0.7) 
    Datei = ROOTFile(string(loadpath,"/KM3NeT_00000133_000",Run,"_S.root"))
    Slices = Datei.online.summaryslices
    Start = Slices[1].header.t.s
    len = Int64(length(Slices))
    End = Slices[len].header.t.s 
    sections = floor(Int32, len/slice_length)
    last_section_length = len%slice_length
    last_section = last_section_length >= (slice_length*slice_length_threshold)
    max_section = sections+last_section
    Dom_ids = optical_DomIds(detector)
    mkpath(storagepath)
    file = h5open(string(storagepath,"/",Run,"_",Int32(slice_length/600),".h5"), "w")
    create_group(file, "used_config")
    write(file["used_config"], "last_section", last_section)
    write(file["used_config"], "last_section_length", last_section_length)
    write(file["used_config"], "slice_length", slice_length)
    write(file, "start", Start)
    write(file, "end", End)
    close(file)
    for i in (1:ceil(Int32, length(Dom_ids)/30)) #was für schritte hier machen? - im moment 30er
        Start = (i-1)*30+1
        End = minimum([i*30, length(Dom_ids)])
        Doms = Dom_ids[Start:End]
        Dom_count = End-Start+1
        good_values = zeros(Int32, max_section, Dom_count)
        hrvcount = zeros(Int32, Dom_count, PMT_count, max_section)
        fifocount = zeros(Int32, Dom_count, PMT_count, max_section)
        pmtmean = zeros(Float64, Dom_count, PMT_count, max_section)
        for i in (1:sections)
            pmtmean[:,:,i], hrvcount[:,:,i], fifocount[:,:,i], good_values[i,:] = inner_extract_Data(Doms, Slices, slice_length, i, Dom_count)
        end
        if last_section
            pmtmean[:,:,max_section], hrvcount[:,:,max_section], fifocount[:,:,max_section], good_values[max_section,:] = inner_extract_Data(Doms, Slices, last_section_length, max_section, Dom_count)
        end 
        store_Data(Doms, pmtmean, hrvcount, fifocount, good_values, slice_length, Run, storagepath)
    end    
    return 0 
end


function inner_extract_Data(Doms::Vector{Int32}, Slices::KM3io.SummarysliceContainer, slice_length::Integer, i::Integer, Dom_count::Integer)
    inner_hrvcount = zeros(Int32, Dom_count, PMT_count)
    inner_fifocount = zeros(Int32, Dom_count, PMT_count)
    frequencies = zeros(Float64, Dom_count, PMT_count, slice_length)
    good_values = zeros(Int32, Dom_count)
    for j in (1:slice_length)
        for id in (1:length(Slices[(i-1)*slice_length+j].frames))
            frame = Slices[(i-1)*slice_length+j].frames[id]
            if frame.dom_id in Doms
                position = findfirst(x->x==frame.dom_id,Doms)
                good_values[position] += 1
                F = pmtrates(frame)
                !wrstatus(frame) && continue                 
                if hrvstatus(frame) || fifostatus(frame)
                    for pmt in (1:PMT_count)
                        loc_hrv = hrvstatus(frame, pmt-1)
                        loc_fifo = fifostatus(frame, pmt-1)
                        if  loc_hrv ||  loc_fifo
                            F[pmt] = 0
                        end
                        inner_hrvcount[position, pmt] += loc_hrv
                        inner_fifocount[position, pmt] += loc_fifo
                    end
                end
                frequencies[position,:,j] = F
            end
        end 
    end
    pmtmean = [mean(filterzero(frequencies[j,i,:])) for j in (1:Dom_count), i in (1:PMT_count)]
    return pmtmean, inner_hrvcount, inner_fifocount, good_values
end

function store_Data(Doms::Vector{Int32}, pmtmean::Array{Float64}, hrvcount::Array{Int32}, fifocount::Array{Int32}, good_values::Matrix{Int32}, slice_length::Integer, Run::Integer, storagepath::String)
    len = size(good_values[:,1])[1]
    file = h5open(string(storagepath,"/",Run,"_",Int32(slice_length/600),".h5"), "cw")
    for i in (1:length(Doms))
        create_group(file, string(Doms[i]))
        dset = create_dataset(file[string(Doms[i])], "good_values", Int32, (len,))
        dset2 = create_dataset(file[string(Doms[i])], "pmtmean", Float64, (PMT_count,len))
        dset3 = create_dataset(file[string(Doms[i])], "hrvcount", Int32, (PMT_count,len))
        dset4 = create_dataset(file[string(Doms[i])], "fifocount", Int32, (PMT_count,len))
        write(dset, good_values[:,i])
        write(dset2, pmtmean[i,:,:])
        write(dset3, hrvcount[i,:,:])
        write(dset4, fifocount[i,:,:])
    end
    close(file)
end 


function Data(filename::String, detector::KM3io.Detector, storagepath::String; slice_length=6000, slice_length_threshold=0.7)
    Run = parse(Int32, filename[collect(findlast("000", filename))[3]+1:collect(findlast("000", filename))[3]+5])
    loadpath = filename[1:collect(findlast("KM3NeT_00000133", filename))[1]-2]
    Data(loadpath, Run, detector, storagepath, slice_length=slice_length, slice_length_threshold=slice_length_threshold)
end


function inner_loadData(DomId::String, file::HDF5.File)
    good_values = read(file[DomId]["good_values"])
    pmtmean = read(file[DomId]["pmtmean"])
    hrvcount = read(file[DomId]["hrvcount"])
    fifocount = read(file[DomId]["fifocount"])
    return (good_values, pmtmean, hrvcount, fifocount)
end 

function load_Data_run(FileName::String, loadpath::String)
    file = h5open(string(loadpath,"/",FileName), "r")
    Keys = collect(keys(file))
    Keys_Dom_id = [key_ids for key_ids in Keys if key_ids != "start" && key_ids != "end" && key_ids != "used_config"]
    sections = length(read(file[Keys_Dom_id[1]]["good_values"]))
    DataDict = Dict(Id=>(zeros(Int32,sections),zeros(Float64,PMT_count,sections),zeros(Int32,PMT_count,sections),zeros(Int32,PMT_count,sections)) for Id in parse.(Int32, Keys_Dom_id))
    for key in Keys_Dom_id
        DomId = parse(Int32, key)
        DataDict[DomId] = inner_loadData(string(DomId), file)
    end
    Start = Int64(read(file["start"]))
    End = Int64(read(file["end"]))
    last_section = read(file["used_config"]["last_section"])
    last_section_length = read(file["used_config"]["last_section_length"])
    slice_length = read(file["used_config"]["slice_length"])
    close(file)
    return (Start, End, DataDict, (last_section, last_section_length, slice_length))
end 




function DomData(Dom_Id::Integer, storagepath::String; files::Vector{String}=String[], slice_length::Integer=6000, loadpath::String="../Data")
    loadpath = string(loadpath,"/Runs_sl",Int32(slice_length/600),"Min/")
    files1 = glob(string("*",Int32(slice_length/600),".h5"), loadpath) 
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
    return (files, wrong_sliceTime, notfound)
end


"""
NOTIEREN WAS SIE GENAU MACHT 
UND WIE - ggf OPTIMIEREN    
"""
#TODO: durchlaufen lassen - schauen wie lange das geht + was die macht 
function DomData_Floors(detector::Detector, storagepath::String; slice_length=6000, loadpath::String="../Data/DomData_Doms")
    possible_files = glob(string("*",Int32(slice_length/600),".h5"), loadpath) 
    floors, Doms_on_floor, max_floor = ToolBox.Floors(detector)
    file = h5open(string(loadpath,"/Dom_",Doms_on_floor[1][1],"_",string(Int32(slice_length/600)),".h5"), "r")
    Times = read(file["Time"])
    close(file)
    len = length(Times)
    good_values = zeros(Int64, max_floor, len)
    pmtmean = Array{Float64}(undef, max_floor, len, 31)
    shifted_pmtmean = Array{Float64}(undef, max_floor, len, 31)
    hrvcount = Array{Int64}(undef, max_floor, len, 31)
    fifocount = Array{Int64}(undef, max_floor, len, 31)
    for floor in (1:max_floor)
        Dom_count = length(Doms_on_floor[floor])
        prov_pmtmean = zeros(Float64, Dom_count, len, 31)
        sh_prov_pmtmean = zeros(Float64, Dom_count, len, 31)
        for id in (1:Dom_count)
            Dom = Doms_on_floor[floor][id]
            filename = string(loadpath,"\\Dom_",Dom,"_",string(Int32(slice_length/600)),".h5")
            if filename in possible_files
                file = h5open(filename, "r")
                if Times != read(file["Time"])
                    print("Zeiten passen nicht")
                    break
                else  
                    tmp_good_values = read(file["good_values"])
                    tmp_pmtmean = read(file["pmtmean"])
                    sh_tmp_pmtmean = shift_pmtmeans(tmp_pmtmean)
                    tmp_hrvcount = read(file["hrvcount"])
                    tmp_fifocount = read(file["fifocount"])
                    for i in (1:len)
                        good_values[floor,i] += tmp_good_values[i]
                        prov_pmtmean[id,i,:] .= tmp_pmtmean[i,:]
                        sh_prov_pmtmean[id,i,:] .= sh_tmp_pmtmean[i,:]
                        hrvcount[floor,i,:] += tmp_hrvcount[i,:]
                        fifocount[floor,i,:] += tmp_fifocount[i,:]
                    end
                end
                close(file)
            end
        end
        for time in (1:len)
            pmtmean[floor, time, :] = [mean(ToolBox.filternan(prov_pmtmean[:,time,i])) for i in (1:31)]
            shifted_pmtmean[floor, time, :] = [mean(ToolBox.filternan(sh_prov_pmtmean[:,time,i])) for i in (1:31)]
        end
    end    
    mkpath(storagepath)
    storagefile = h5open(string(storagepath,"/DomData_floors_",Int32(slice_length/600),".h5"), "w")
    dset_time = create_dataset(storagefile, "Time", Int64, len)
    dset_good_values = create_dataset(storagefile, "good_values", Int64, max_floor, len)
    dset_pmtmean = create_dataset(storagefile, "pmtmean", Float64, max_floor, len, 31)
    dset_sh_pmtmean = create_dataset(storagefile, "shifted_pmtmean", Float64, max_floor, len, 31)
    dset_hrvcount = create_dataset(storagefile, "hrvcount", Int64, max_floor, len, 31)
    dset_fifocount = create_dataset(storagefile, "fifocount", Int64, max_floor, len, 31)
    write(dset_time, Times)
    write(dset_good_values, good_values)
    write(dset_pmtmean, pmtmean)
    write(dset_sh_pmtmean, shifted_pmtmean)
    write(dset_hrvcount, hrvcount)
    write(dset_fifocount, fifocount)
    close(storagefile)
end




function mean_Detector(det::Detector; loadpath::String="../Data/DomData_Doms", slice_length::Integer=54000)
    Doms = ToolBox.optical_DomIds(det)
    datafile = h5open(string(loadpath, "/Dom_", Doms[1],"_",Int32(slice_length/600),".h5"), "r") 
    Times = read(datafile["Time"])
    len = length(Times)
    close(datafile)
    DomData = Matrix{Float64}(undef, length(Doms), len)
    for i in (1:length(Doms))
        datafile = h5open(string(loadpath, "/Dom_", Doms[i],"_",Int32(slice_length/600),".h5"), "r") 
        DomData[i,:] = mean(read(datafile["pmtmean"]), dims = 2)
        close(datafile)
    end
    Detector_mean = Float64[]
    for i in (1:len)
        push!(Detector_mean, mean(ToolBox.filternan(DomData[:,i])))
    end
    storagefile = h5open("../Data/DetectorData.h5", "w")
    d_time = create_dataset(storagefile, "Time", Int64, len)
    d_freq = create_dataset(storagefile, "frequency", Float64, len)
    write(d_time, Times)
    write(d_freq, Detector_mean)
    close(storagefile)
    return Times, Detector_mean, mean(Detector_mean)
end

function Detectormean(;loadpath::String="../Data", mean_Times::Bool=false)
    recalculate && mean_Detector(det)
    datafile = h5open(string(loadpath,"/DetectorData.h5"), "r")
    frequencies = read(datafile["frequency"])
    close(datafile)
    mean_Times && (return mean(frequencies))
    return frequencies
end

function linfit_Detector(;T_intervall::Tuple{Integer,Integer}=(0,0), loadpath::String="../Data", ignore_highs::Bool=true)
    datafile = h5open(string(loadpath,"/DetectorData.h5"), "r")
    frequencies = read(datafile["frequency"])
    Times = read(datafile["Time"])
    time_mask = ToolBox.maskTime(Times, T_intervall)
    fit_Times = Times[time_mask] .- Times[time_mask][1]    
    fit_freq = frequencies[time_mask]
    if ignore_highs 
        deleteat!(fit_Times, findall(x->x>=20000,fit_freq))
        deleteat!(fit_freq, findall(x->x>=20000,fit_freq))
    end
    gerade(t, p) = p[1] .+ (p[2] .* t)
    p0 = [fit_freq[1], 0]
    fit = curve_fit(gerade, fit_Times, fit_freq, p0)
    y = fit.param[1] - Times[1]*fit.param[2]
    return (y,fit.param[2])
end