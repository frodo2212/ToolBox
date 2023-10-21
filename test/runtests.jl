using Test
using KM3io
using KM3NeTTestData
using ToolBox
using Dates
using DataFrames


#TODO:
@testset "smallfunctions" begin
    test_array = [0.0,2.0,Inf,5.0,0.0,NaN,12.0,Inf,0.0,NaN]
    @test [false, true, true, true, false, true, true, true, false, true] == ToolBox.maskzero(test_array)
    @test [true, true, true, true, true, false, true, true, true, false] == ToolBox.masknan(test_array)
    @test isequal([2.0, Inf, 5.0, NaN, 12.0, Inf, NaN],ToolBox.filterzero(test_array))
    @test [0.0,2.0,Inf,5.0,0.0,12.0,Inf,0.0] == ToolBox.filternan(test_array)
    @test [true, true, false, true, true, false, true, false, true, false] == ToolBox.maskinfnan(test_array)
    test_array = [1,3,0,5,8,0,6,0]
    @test [1,3,5,8,6] == ToolBox.filterzero(test_array)
    test_array = [[1,5,2,9,0],[12,0,0,0,3,9],[1,2,3,4]]
    @test [[1,5,2,9],[12,3,9],[1,2,3,4]] == ToolBox.filterzero(test_array)

    det = Detector("testData/TestData14422.detx")
    det2 = Detector(datapath("detx", "km3net_offline.detx"))
    Test_Dom = 809503416
    @test DateTime(2023,3,4,2,11) == DateTime_autocorrect(2022,14,31,25,71)
    @test (2023,3,4,2,11) == DateTime_autocorrect((2022,14,31),Time=(25,71),Datetime=false)
    @test ((0.142, -0.945, 0.295),5) == ToolBox.PMT_Direction(det2, 12, Test_Dom) #hier noch nen test mit det schreiben?
    @test 378 == length(optical_DomIds(det))
    @test 90 == length(optical_DomIds(det2))
    @test 808972605 in optical_DomIds(det)
    floors, dict, max_floor, floor_size = ToolBox.Floors(det2)
    @test (18, 5) == (max_floor, floor_size)
    @test [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18] == floors 
    @test dict[5] == Int32[808945480, 808979567, 806483369, 808959411, 808489117]
    @test (DataFrame(dates=unix2datetime(1670764952) : Minute(18) : unix2datetime(1670769052)),"mm/dd/yyyyTHH:MM") == ToolBox.autoscale_time(1670764952,1670769052)
    Times = [1670764555,1670764850,1670764955,1670769050,1670769155,1670769350]
    @test [false,false,true,true,false,false] == ToolBox.maskTime(Times, (1670764952,1670769052))
    String_dict = ToolBox.pos_Strings(det)
    @test String_dict[5] == (-165.18, 91.84)
    @test length(keys(String_dict)) == 21
    @test (1702166400, 1704412800) == ToolBox.T_intervall((2023,11,40),(2023,13,5))
    @test (1702166400, 1714608000) == ToolBox.T_intervall2((2023,11,40),(0,3,52))
end

@testset "Data.jl" begin
    # det = Detector(datapath("detx", "km3net_offline.detx"))
    # f = ROOTFile(datapath("online", "km3net_online.root"))  #da brauch ich was anderes, 3 summaryslices reichen nicht
    @test 1 == 1
    #solution = ()
    #@test solution == extract_Data(f.online.summaryslices, det, slice_length=500)
end

@testset "DomData_search.jl" begin
    det = Detector("testData/TestData14422.detx")
    @test [806481218] == Search_DomData((8000,20000),(4,80),loadpath="testData")
    @test Int64[] == Search_DomData((10000,20000),(4,80),loadpath="testData")
    @test Event[] == ToolBox.intSearch_DomData_Dom(806481218, (8000,15000), loadpath="testData")
    @test Event[ToolBox.Event("iS_t1", 806481218, 23, 80, 87, 7, 1663906235, 1663910435, 4200, 0)] == ToolBox.intSearch_DomData_Dom(806481218, (7500,15000), loadpath="testData")
    @test (6507.236166251622, -0.00011199840989034525) == linfit_DomData((806481218, 12), loadpath="testData") #keine Ahnung ob das so sinn macht
    # @test 0 == linfit_DomData_intervalls((806481218, 12), loadpath="test/testData")  #erst die Funktion fixen
    @test linfitData[] == ToolBox.Search_linFitData(det, loadpath="testData")
end

@testset "DomData.jl" begin
    @test 1 == 1
end

