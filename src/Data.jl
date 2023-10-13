function Data(loadpath::String, Run::Integer, detector::KM3io.Detector, storagepath::String; slice_length=6000, slice_length_threshold=0.7) 
    Datei = ROOTFile(string(loadpath,"/KM3NeT_00000133_000",Run,"_S.root"))
    Slices = Datei.online.summaryslices
    Start = Slices[1].header.t.s
    len = Int64(length(Slices))
    End = Slices[len].header.t.s 
    sections = floor(Int32, len/slice_length)
    last_section_length = len%slice_length
    last_section = last_section_length >= (slice_length*slice_length_threshold)
    Dom_ids = optical_DomIds(detector)
    mkpath(storagepath)
    file = h5open(string(storagepath,"/",Run,"_",Int32(slice_length/600),".h5"), "w")
    create_group(file, "used_config")
    write(file["used_config"], "last_section", last_section)
    write(file["used_config"], "last_section_length", last_section_length)
    write(file["used_config"], "slice_length", slice_length)
    close(file)
    for i in (1:ceil(Int32, length(Dom_ids)/30)) #was für schritte hier machen? - im moment 30er
        Start = (i-1)*30+1
        End = minimum([i*30, length(Dom_ids)])
        Doms = Dom_ids[Start:End]
        Dom_count = End-Start+1
        good_values = zeros(Int32, sections+last_section, Dom_count)
        hrvcount = zeros(Int32, Dom_count, PMT_count, sections+last_section)
        fifocount = zeros(Int32, Dom_count, PMT_count, sections+last_section)
        pmtmean = zeros(Float64, Dom_count, PMT_count, sections+last_section)
        for i in (1:sections)
            pmtmean[:,:,i], hrvcount[:,:,i], fifocount[:,:,i], good_values[i,:] = inner_extract_Data(Doms, Slices, slice_length, i, Dom_count)
        end
        if last_section
            pmtmean[:,:,i], hrvcount[:,:,i], fifocount[:,:,i], good_values[i,:] = inner_extract_Data(Doms, Slices, last_section_length, sections+1, Dom_count)
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


function Data(filename::String, detector::KM3io.Detector, loadpath::String, storagepath::String; slice_length=6000, slice_length_threshold=0.7)
    Run = parse(Int32, filename[collect(findlast("000", filename))[3]+1:collect(findlast("000", filename))[3]+5])
    loadpath = filename[1:collect(findfirst(string(Run), filename))[1]]
    Data(loadpath, Run, detector, storagepath, slice_length=slice_length, slice_length_threshold=slice_length_threshold)
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
