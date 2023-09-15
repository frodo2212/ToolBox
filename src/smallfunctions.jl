#was davon verwende ich überhaupt noch aktiv?
"""
returns the Direction and Ring of a given PMT in a Standard_Dom
if given Dom_Id it returns the exact values of given PMT in Dom
"""
function PMT_Position(det::Detector, PMT::Integer, Dom_Id::Integer=0)
    PMT_Ring = Int32[6,5,5,5,6,6,5,6,6,6,5,5,4,3,2,4,4,3,2,2,3,3,1,4,2,2,2,4,3,3,4]
    if Dom_Id != 0 && Dom_Id in collect(keys(det.modules))
        dir_x = det.modules[Dom_Id].pmts[PMT].dir.x
        dir_y = det.modules[Dom_Id].pmts[PMT].dir.y
        dir_z = det.modules[Dom_Id].pmts[PMT].dir.z
    else 
        dir_x = 0 #das hier noch anders gestalten
        dir_y = 0 #gibts den standard DOM schon, aus dem ich die geometrie raus lesen kann?
        standard_dirz = [-1,-0.850,-0.555,-0.295,0.295,0.555]
        dir_z = standard_dirz[PMT_Ring[PMT]]
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
    #println((time_end-time_start))
    intervall_minutes = (time_end-time_start)/(60*intervalls)
    #println(intervall_minutes)
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
    floors = collect(0:maximum([collect(keys(detector.locations))[i][2] for i in (1:length(keys(detector.locations)))]))
    max_floor = maximum(floors)
    Doms_on_Floor = Dict(i=>Int32[] for i in floors)
    for i in floors
        for j in detector.strings
            push!(Doms_on_Floor[i],detector.locations[j,i].id)
        end
    end 
    floor_size = length(Doms_on_Floor[1])
    return (floors, Doms_on_Floor, max_floor, floor_size)
end

function pos_Strings(detector::Detector; plotten::Bool=false)
    Groundmodules = [detector.locations[j,0].id for j in detector.strings]
    len = length(Groundmodules)
    pos_x = [detector.modules[Groundmodules[i]].pos.x for i in (1:len)]
    pos_y = [detector.modules[Groundmodules[i]].pos.y for i in (1:len)]
    if plotten 
        fig = Figure()
        ax = Axis(fig[1,1], title="Position der Bodenmodule der Strings")
        elements = Vector{Any}(undef, len)
        label= Vector{String}(undef, len)
        for i in (1:len)
            elements[i] = scatter!(ax, pos_x[i], pos_y[i])
            label[i] = string("String", i)
            
        end
        Legend(fig[1,2], elements, label)
        return fig, (pos_x, pos_y)
    else 
        return (pos_x, pos_y)
    end
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
takes a Vector of Vectors and returns a 1d Vector with all the elements
"""
function vectorize(array::Vector{Vector{T}}) where T<:Real
    len = length(array)
    inner_length = Vector{Int32}(undef,len)
    for i in (1:len)
        inner_length[i] = length(array[i])
    end
    newvector = Vector{typeof(array[1][1])}(undef, sum(inner_length))
    for i in (1:len)
        for j in (1:inner_length[i])
            newvector[sum(inner_length[1:i-1])+j] = array[i][j]
        end
    end
    return newvector
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
    neu = Vector{}(undef, len)
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

