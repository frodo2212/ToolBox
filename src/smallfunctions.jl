#was davon verwende ich überhaupt noch aktiv?
"""
returns the Direction and Ring of a given PMT in a Standard_Dom
if given Dom_Id it returns the exact values of given PMT in Dom
"""
function PMT_Direction(det::Detector, PMT::Integer, Dom_Id::Integer=0)
    PMT_Ring = Int32[6,5,5,5,6,6,5,6,6,6,5,5,4,3,2,4,4,3,2,2,3,3,1,4,2,2,2,4,3,3,4]
    if Dom_Id != 0 && Dom_Id in collect(keys(det.modules))
        dir_x = det.modules[Dom_Id].pmts[PMT].dir.x
        dir_y = det.modules[Dom_Id].pmts[PMT].dir.y
        dir_z = det.modules[Dom_Id].pmts[PMT].dir.z
    else 
        dir_x = 0 #das hier noch anders gestalten
        dir_y = 0 #gibts den standard DOM schon, aus dem ich die geometrie raus lesen kann?
        standard_dirz = [-1,-0.850,-0.555,-0.295,0.295,0.555]
        dir_z = standard_dirz[config.Detector_PMT_Ring[PMT]]
    end
    return ((dir_x, dir_y, dir_z), PMT_Ring[PMT])
end
"""
give it start and end as unix time and it returns a Dataframe of Times to scale the Axis with
possible Argument: intervalls [default=4] sets the Maximum of generated timestamps
"""
function autoscale_time(time_start::Real, time_end::Real;intervalls::Integer=4)
    time_start = round(Int64, time_start)
    time_end = round(Int64, time_end)
    intervall_minutes = (time_end-time_start)/(60*intervalls)
    if intervall_minutes <= 60
        Zeiten = DataFrame(dates=unix2datetime(time_start) : Minute(ceil(intervall_minutes)) : unix2datetime(time_end))
        Typ = "mm/dd/yyyyTHH:MM"
    elseif (intervall_minutes/60) <=24
        Zeiten = DataFrame(dates=unix2datetime(time_start) : Hour(ceil(intervall_minutes/60)) : unix2datetime(time_end))
        Typ = "mm/dd/yyyyTHH"
    else 
        Zeiten = DataFrame(dates=unix2datetime(time_start) : Day(ceil(intervall_minutes/1440)) : unix2datetime(time_end))
        Typ = "mm/dd/yyyy"
    end
    #print(Zeiten)
    return (Zeiten, Typ)
end
"""
takes a Detector, returns:
floors = Vector of all the possible floors = keys of Doms_on_floor
Doms_on_dloor = Dict{Floor_Number::Int => DomIds_on_this_floor::Vector{Int64}}
floor_count = Number of possible floors
floor_size = Number of Doms on a floor
"""
function Floors(detector::Detector)
    #das muss doch schlauer gehen!
    max_floor = maximum([collect(keys(detector.locations))[i][2] for i in (1:length(keys(detector.locations)))])
    min_floor = minimum([collect(keys(detector.locations))[i][2] for i in (1:length(keys(detector.locations)))])
    floors = collect(min_floor:max_floor)
    Doms_on_Floor = Dict(i=>Int32[] for i in floors)
    for i in floors
        for j in detector.strings
            push!(Doms_on_Floor[i],detector.locations[j,i].id)
        end
    end 
    floor_size = length(Doms_on_Floor[1])
    return (floors, Doms_on_Floor, max_floor, floor_size)
end

function pos_Strings(detector::Detector)
    floors, = Floors(detector)
    String_positions = Dict(i=>(detector.modules[detector.locations[i,floors[1]].id].pos.x,detector.modules[detector.locations[i,floors[1]].id].pos.y) for i in detector.strings)
    return String_positions
end

function pos_Dom(Dom_Id::Integer, detector::Detector)
    return detector.modules[Dom_Id].pos.x, detector.modules[Dom_Id].pos.y, detector.modules[Dom_Id].pos.z
end

function pos_Doms(Dom_Ids::Vector{Int32}, detector::Detector)
    positions = Dict{Int32, Tuple{Float64,Float64,Float64}}
    for Dom in Dom_Ids
        positions = merge!(positions, Dict(Dom=>pos_Dom(Dom, detector)))
    end
    return positions
end

function pos_Doms(detector::Detector)
    return pos_Doms(optical_DomIds(detector), detector)
end    

function close_Doms(Dom_Id::Integer, detector::Detector, range::Real)
    Doms = optical_DomIds(detector)
    deleteat!(Doms, findfirst(x->x==Dom_Id,Doms))
    positions = ToolBox.pos_Doms(detector)
    close_Doms = Int32[]
    for Dom in Doms
        distance(positions[Dom_Id], positions[Dom]) <= range && push!(close_Doms, Dom)
    end
    return close_Doms
end

function distance(PointA::Tuple{Float64,Float64,Float64}, PointB::Tuple{Float64,Float64,Float64})
    x = PointA[1]-PointB[1]
    y = PointA[2]-PointB[2]
    z = PointA[3]-PointB[3]
    return sqrt(x*x+y*y+z*z)
end

function optical_DomIds(detector::Detector)    
    floors = collect(1:maximum([collect(keys(detector.locations))[i][2] for i in (1:length(keys(detector.locations)))]))
    optical_Dom_Ids = Int32[]
    for i in floors
        for j in detector.strings
            push!(optical_Dom_Ids,detector.locations[j,i].id)
        end
    end
    return optical_Dom_Ids
end

"""
works just like DateTime for Dates, but doesnt care about argument out of range errors
"""
function inner_DateTime_autocorrect((Year, Month, Day)::Tuple{Integer,Integer,Integer}; Datetime::Bool=false)
    if Month >= 13
        Year += floor(Int64, Month/12)
        Month = Month%12
    end
    while Day > daysinmonth(Date(Year, Month))
        Day -= daysinmonth(Date(Year, Month))
        Month += 1
        if Month == 13
            Month = 1
            Year += 1
        end
    end
    if Datetime == true
        return DateTime(Year, Month, Day)
    else
        return (Year, Month, Day)
    end
end

#TODO: is this function necessary
function DateTime_autocorrect(Year::Integer, Month::Integer, Day::Integer,Hour::Integer, Minute::Integer; Datetime::Bool=true)
    return DateTime_autocorrect((Year,Month,Day),Time=(Hour,Minute), Datetime=Datetime)
end
"""
works just like DateTime but instead of giving bound errors it calculates a corresponding date
"""
function DateTime_autocorrect((Year, Month, Day)::Tuple{Integer,Integer,Integer}; Time=(0,0), Datetime::Bool=true)
    Day_from_time = 0
    Hour = Time[1]
    Minute = Time[2]
    if Minute >= 60
        Hour += floor(Int64, Minute/60)
        Minute = Minute%60
    end
    if Hour >= 24
        Day_from_time += floor(Int64, Hour/24)
        Hour = Hour%24
    end
    year, month, day = inner_DateTime_autocorrect((Year, Month, Day+Day_from_time), Datetime=false)
    if Datetime
        return DateTime(year, month, day, Hour, Minute)
    else 
        return year, month, day, Hour, Minute
    end
end

#TODO: find a better Name and maybe a better argument structure
function T_intervall2(Start::Tuple{Integer,Integer,Integer}=(0,0,0), Intervall::Tuple{Integer,Integer,Integer}=(0,0,0))
    start_DT = inner_DateTime_autocorrect(Start, Datetime=true)
    end_DT = inner_DateTime_autocorrect((Start[1]+Intervall[1],Start[2]+Intervall[2],Start[3]+Intervall[3]), Datetime=true) 
    return (Int64(datetime2unix(start_DT)), Int64(datetime2unix(end_DT)))
end

function T_intervall(Start_Date::Tuple{Integer,Integer,Integer}, End_Date::Tuple{Integer,Integer,Integer}; Start_Time=(0,0), End_Time=(0,0))
    Start = DateTime_autocorrect(Start_Date, Time=Start_Time)
    End = DateTime_autocorrect(End_Date, Time=End_Time)
    return (Int64(datetime2unix(Start)), Int64(datetime2unix(End)))
end

"""
Takes the Array of Timestamps and two bounds and masks every Timestamp inside the bound
"""
function maskTime(Times, T_intervall::Tuple{Real,Real})
    len = length(Times)
    if T_intervall[2] != 0
        T_mask = [(Times[i] >= T_intervall[1] && Times[i] <= T_intervall[2]) for i in (1:len)]
    elseif T_intervall[1] != 0
        T_mask = [Times[i] >= T_intervall[1] for i in (1:len)]
    else 
        T_mask = ones(Bool, len)
    end
    return T_mask
end
  

"""
takes a Vector and returns a mask with true for all elements, that aren't zero
"""
function maskzero(Vektor::Vector{T}) where T <: Real
    mask = Vector{Bool}(undef, length(Vektor))
    for i in (1:length(Vektor))
        mask[i] = true
        if Vektor[i] == 0
            mask[i] = false
        end
    end
    return mask
end
"""
takes a Vector of Vectors and removes all Elements that are Zero 
"""
function filterzero(array::Vector{Vector{T}}) where T<:Real
    neu = Vector{Vector{T}}(undef, length(array))
    for i in (1:length(array))
        neu[i] = filterzero(array[i])
    end
    return neu
end
"""
takes a Vector and removes all Elements that are Zero 
"""
function filterzero(array::Vector{T}) where T<:Real
    len = length(array)-count(j->(j==0), array)
    neu = Vector{T}(undef, len)
    j=0
    for i in (1:length(array))
        if array[i] == 0
            j = j+1
        else
            neu[i-j] = array[i]
        end
    end
    return neu
end
"""
takes a Vector and removes all Elements that are NaN 
"""
function filternan(array::Vector{T}) where T<:Real
    len = length(array)-count(isnan.(array)) 
    neu = Vector{T}(undef, len)
    j=0
    for i in (1:length(array))
        if isnan(array[i]) 
            j = j+1
        else
            neu[i-j] = array[i]
        end
    end
    return neu
end

function masknan(Vektor::Vector{T}) where T <: Real
    mask = Vector{Bool}(undef, length(Vektor))
    for i in (1:length(Vektor))
        mask[i] = true
        if isnan(Vektor[i])
            mask[i] = false
        end
    end
    return mask
end

function maskinfnan(Vektor::Vector{T}) where T <: Real
    mask = Vector{Bool}(undef, length(Vektor))
    for i in (1:length(Vektor))
        mask[i] = true
        if isnan(Vektor[i]) || isinf(Vektor[i])
            mask[i] = false
        end
    end
    return mask
end


function save_Events(Events::Vector{Event}, Name::String; storagepath::String="./")
    file = h5open(string(storagepath,"/",Name,".h5"), "w")
    for i in (1:length(Events))
        create_group(file, string(i))
        write(file[string(i)], "Typ", Events[i].Typ)
        write(file[string(i)], "Dom_Id", Events[i].Dom_Id)
        write(file[string(i)], "pmt", Events[i].pmt)
        write(file[string(i)], "array_start", Events[i].array_start)
        write(file[string(i)], "array_end", Events[i].array_end)
        write(file[string(i)], "array_length", Events[i].array_length)
        write(file[string(i)], "time_start", Events[i].time_start)
        write(file[string(i)], "time_end", Events[i].time_end)
        write(file[string(i)], "time_length", Events[i].time_length)
        write(file[string(i)], "missing_timestamps", Events[i].missing_timestamps)
        write(file[string(i)], "slice_length", Events[i].slice_length)
    end
    close(file)
end
function load_Events(Name::String; loadpath::String="./")
    Events = Event[]
    file = h5open(string(loadpath,"/",Name,".h5"), "r")
    for i in keys(file)
        Typ = read(file[string(i)], "Typ")
        Dom_Id = read(file[string(i)], "Dom_Id")
        pmt = read(file[string(i)], "pmt")
        array_start = read(file[string(i)], "array_start")
        array_end =read(file[string(i)], "array_end")
        array_length = read(file[string(i)], "array_length")
        time_start = read(file[string(i)], "time_start")
        time_end = read(file[string(i)], "time_end")
        time_length = read(file[string(i)], "time_length")
        missing_timestamps = read(file[string(i)], "missing_timestamps")
        slice_length = read(file[string(i)], "slice_length")
        push!(Events, Event(Typ, Dom_Id, pmt, array_start, array_end, array_length, time_start, time_end, time_length, missing_timestamps, slice_length))
    end
    close(file)
    return Events
end

function save_linFitData(Events::Vector{linfitData}, Name::String; storagepath::String="./")
    file = h5open(string(storagepath,"/",Name,".h5"), "w")
    for i in (1:length(Events))
        create_group(file, string(i))
        write(file[string(i)], "Typ", Events[i].Typ)
        write(file[string(i)], "Dom_Id", Events[i].Dom_Id)
        write(file[string(i)], "pmt", Events[i].pmt)
        write(file[string(i)], "time", Events[i].time)
        write(file[string(i)], "time_intervall", Events[i].time_intervall)
        write(file[string(i)], "params1", Events[i].params[1])
        write(file[string(i)], "params2", Events[i].params[2])
        write(file[string(i)], "rel_values", Events[i].rel_values)
        write(file[string(i)], "slice_length", Events[i].slice_length)
    end
    close(file)
end
function load_linfitData(Name::String; loadpath::String="./")
    Events = linfitData[]
    file = h5open(string(loadpath,"/",Name,".h5"), "r")
    for i in keys(file)
        Typ = read(file[string(i)], "Typ")
        Dom_Id = read(file[string(i)], "Dom_Id")
        pmt = read(file[string(i)], "pmt")
        time = read(file[string(i)], "time")
        time_intervall = read(file[string(i)], "time_intervall")
        params = (read(file[string(i)], "params1"),read(file[string(i)], "params2"))
        rel_values = read(file[string(i)], "rel_values")
        slice_length = read(file[string(i)], "slice_length")
        push!(Events, linfitData(Typ, Dom_Id, pmt, time, time_intervall, params, rel_values, slice_length))
    end
    close(file)
    return Events
end