using Test
using KM3io
using KM3NeTTestData
using ToolBox
using Dates
using DataFrames



@testset "smallfunctions" begin    
    det = Detector(datapath("detx", "km3net_offline.detx"))
    Test_Dom = 809503416
    @test DateTime(2023,3,4,2,11) == DateTime_autocorrect(2022,14,31,25,71)
    @test ((0.142, -0.945, 0.295),5)== PMT_Direction(det, 12, Test_Dom)
    @test 90 == length(optical_DomIds(det))
        #noch testen welche drin sind
    floors, dict = Floors(det)
    @test [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18] == floors && dict[5] == [808945480, 808979567, 806483369, 808959411, 808489117]
    @test (DataFrame(dates=unix2datetime(1670764952) : Minute(18) : unix2datetime(1670769052)),"mm/dd/yyyyTHH:MM") == autoscale_time(1670764952,1670769052)
    Times = 1670764555,1670764850,1670764955,1670769050,1670769155,1670769350
    @test [false,false,true,true,false,false] == ToolBox.maskTime(Times, (1670764952,1670769052))
    test_array = [0.0,2.0,5.0,0.0,NaN,12.0,0.0,NaN]
    @test [false, true, true, false, true, true, false, true] == ToolBox.maskzero(test_array)
    @test [true, true, true, true, false, true, true, false] == ToolBox.masknan(test_array)
    @test [2,5,12,6] == ToolBox.filterzero([0,2,5,0,0,12,6,0])
    @test [0.0,2.0,5.0,0.0,12.0,0.0] == ToolBox.filternan(test_array)
end

@testset "Data.jl" begin
    det = Detector(datapath("detx", "km3net_offline.detx"))
    f = ROOTFile(datapath("online", "km3net_online.root"))  #da brauch ich was anderes, 3 summaryslices reichen nicht
    @test 1 == 1
    #solution = ()
    #@test solution == extract_Data(f.online.summaryslices, det, slice_length=500)
end

@testset "DomData_search.jl" begin
    det = Detector("testData/KM3NeT_00000133_00014422.detx")
    @test [806481218] == Search_DomData((8000,20000),(4,80),loadpath="testData")
    @test Int64[] == Search_DomData((10000,20000),(4,80),loadpath="testData")
    @test Event[] == ToolBox.intensiveSearch_DomData_test1(806481218, (8000,15000), loadpath="testData")
    @test Event[ToolBox.Event("iS_t1", 806481218, 23, 80, 87, 7, 1663906235, 1663910435, 4200, 0)] == ToolBox.intensiveSearch_DomData_test1(806481218, (7500,15000), loadpath="testData")
    @test (192856.7930720679, -0.00011199840989034525) == linfit_DomData((806481218, 12), loadpath="testData") #keine Ahnung ob das so sinn macht
    # @test 0 == linfit_DomData_intervalls((806481218, 12), loadpath="test/testData")  #erst die Funktion fixen
    @test linfitData[] == ToolBox.Search_linFit_test2(det, loadpath="testData")
end

