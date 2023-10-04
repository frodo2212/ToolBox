#Die funktion zählt in einem intervall die Frequenzwerte innerhalb der bound
#sind genug frequenzwerte drin, erhält man ein Event - ob diese am stück drüber sind ist irrelevant
#solange genug drüber bleiben ist alles das selbe Event - eventlength wächst
#ein neues Event startet, wenn es zwischendurch mal Timesteps gab für die die Bedingung nicht erfüllt war 
#
#TODO: den Filter verbessern - bisschen willkürlich was ich hier raus filtere
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
            #TODO: wenn der PMT ne Steigung in der Frequenz hat ist frquency_mean nicht optimal... 
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
            # Rückgabewerte sind bisschen random, in neuer funktion verbessert 
            merge!(interesting_Doms, Dict(Dom=>event_data))
        end
        close(datafile)
    end
    return interesting_Doms
end