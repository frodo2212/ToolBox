function extract_Data(Slices::KM3io.SummarysliceContainer, detector::KM3io.Detector; slice_length=6000, slice_length_threshold=0.7) 
    Start = Slices[1].header.t.s
    len = Int64(length(Slices))
    End = Slices[len].header.t.s 
    sections = floor(Int32, len/slice_length)
    last_section_length = len%slice_length
    last_section = last_section_length >= (slice_length*slice_length_threshold)
    Dom_ids = optical_DomIds(detector)
    Dom_count = length(Dom_ids) 
    array_pos = Dict(Dom_ids[i]=>(1:Dom_count)[i] for i in (1:Dom_count))
    Data = Dict{Int32, Tuple{Vector{Int32},Matrix{Float64},Matrix{Int32},Matrix{Int32}}}(Id=>(zeros(Int32,sections+last_section),zeros(Float64,PMT_count,sections+last_section),zeros(Int32,PMT_count,sections+last_section),zeros(Int32,PMT_count,sections+last_section)) for Id in Dom_ids)
    for i in (1:sections)
        inner_extract_Data(Slices, slice_length, array_pos, i, Dom_count, Dom_ids, Data)
    end
    if last_section
        inner_extract_Data(Slices, last_section_length, array_pos, sections+1, Dom_count, Dom_ids, Data)
    end 
    return (Start, End, Data, (last_section, last_section_length, slice_length))
end

function inner_extract_Data(Slices::KM3io.SummarysliceContainer, slice_length::Integer, array_pos, i::Integer, Dom_count::Integer, Dom_ids, Data)
    frequencies, inner_hrvcount, inner_fifocount, good_values = extract_loops(Slices, slice_length, array_pos, i, Dom_count)
    for allDoms in Dom_ids
        Data[allDoms][1][i] = good_values[array_pos[allDoms]]
        Data[allDoms][2][:,i] = [mean(filterzero(frequencies[array_pos[allDoms],i,:])) for i in (1:PMT_count)]
        Data[allDoms][3][:,i] = inner_hrvcount[array_pos[allDoms],:]
        Data[allDoms][4][:,i] = inner_fifocount[array_pos[allDoms],:]
    end
end

function extract_loops(Slices::KM3io.SummarysliceContainer, slice_length::Integer, array_pos, i::Integer, Dom_count::Integer)
    inner_hrvcount = zeros(Int32, Dom_count, PMT_count)
    inner_fifocount = zeros(Int32, Dom_count, PMT_count)
    frequencies = zeros(Float64, Dom_count, PMT_count, slice_length)
    good_values = zeros(Int32, Dom_count)
    for j in (1:slice_length)
        for id in (1:length(Slices[(i-1)*slice_length+j].frames))
            frame = Slices[(i-1)*slice_length+j].frames[id]
            position = array_pos[frame.dom_id]
            good_values[position] +=1
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
    return (frequencies, inner_hrvcount, inner_fifocount, good_values)
end

function store_Data(Data::Tuple{UInt32, UInt32, Dict{Int32, Tuple{Vector{Int32}, Matrix{Float64}, Matrix{Int32}, Matrix{Int32}}}, Tuple{Bool, Integer, Integer}}, Run::Int32, storagepath::String)
    Dom_ids = collect(keys(Data[3]))
    len = length(Data[3][Dom_ids[1]][1])
    mkpath(string(storagepath,"/Run_",Run,"_",Int32(Data[4][3]/600),".h5"))
    file = h5open(string(storagepath,"/Run_",Run,"_",Int32(Data[4][3]/600),".h5"), "w")
    for Dom in Dom_ids
        create_group(file, string(Dom))
        dset = create_dataset(file[string(Dom)], "good_values", Int32, (len,))
        dset2 = create_dataset(file[string(Dom)], "pmtmean", Float64, (PMT_count,len))
        dset3 = create_dataset(file[string(Dom)], "hrvcount", Int32, (PMT_count,len))
        dset4 = create_dataset(file[string(Dom)], "fifocount", Int32, (PMT_count,len))
        write(dset, Data[3][Dom][1])
        write(dset2, Data[3][Dom][2])
        write(dset3, Data[3][Dom][3])
        write(dset4, Data[3][Dom][4])
    end
    write(file, "start", Data[1])
    write(file, "end", Data[2])
    create_group(file, "used_config")
    write(file["used_config"], "last_section", Data[4][1])
    write(file["used_config"], "last_section_length", Data[4][2])
    write(file["used_config"], "slice_length", Data[4][3])
    close(file)
end 


function Data(filename::String, detector::KM3io.Detector, storagepath::String; slice_length=6000)
    Datei = ROOTFile(filename)
    Data = extract_Data(Datei.online.summaryslices, detector, slice_length=slice_length)
    Runnumber = parse(Int32, filename[collect(findlast("000", filename))[3]+1:collect(findlast("000", filename))[3]+5])
    store_Data(Data, Runnumber, string(storagepath,"/Runs_sl",Int32(slice_length/600),"Min"))
end

function Data(Runnumber::Integer, detector::KM3io.Detector, loadpath::String, storagepath::String; slice_length=6000)
    Datei = ROOTFile(string(loadpath,"/KMeNeT_00000133_000",Runnumber,"_S.root"))
    Data = extract_Data(Datei.online.summaryslices, detector, slice_length=slice_length)
    store_Data(Data, Runnumber, string(storagepath,"/Runs_sl",Int32(slice_length/600),"Min"))
end


function inner_loadData(DomId::Int32, file::HDF5.File)
    good_values = read(file[string(DomId)]["good_values"])
    pmtmean = read(file[string(DomId)]["pmtmean"])
    hrvcount = read(file[string(DomId)]["hrvcount"])
    fifocount = read(file[string(DomId)]["fifocount"])
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
        DataDict[DomId] = inner_loadData(DomId, file)
    end
    Start = Int64(read(file["start"]))
    End = Int64(read(file["end"]))
    last_section = read(file["used_config"]["last_section"])
    last_section_length = read(file["used_config"]["last_section_length"])
    slice_length = read(file["used_config"]["slice_length"])
    close(file)
    return (Start, End, DataDict, (last_section, last_section_length, slice_length))
end 


#Funktionen ab hier sind eigentlich sinnlos

# mega unnötig, wenn ich alle haben will, dann nehme ich snakemake
# function Data_folder(datapath::String, detectorpath::String, storagepath::String)
#     files = readdir(datapath)
#     det = Detector(detectorpath)
#     for file in files
#         Data(string(datapath, "/",file), det, storagepath)
#     end
# end
