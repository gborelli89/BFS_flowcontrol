# -----------------------------------------------------------------------
# Test set for critical functions
# -----------------------------------------------------------------------
using Test
include("GA.jl")

# Test set for stepper movements
function coordsteps_test()
    @testset "Coordinates and steps" begin
        @testset "One cylinder" begin
            p1 = OneCylType(0.0,0.0)
            p2 = OneCylType(10.0,0.0)
            p3 = OneCylType(20.0,0.0)
            p2_cor = steps2coord(StepperMotionType([41,16],[0]),xorig=[-8.8],d=2.4, Rh=80.0, Re=10.0)
            p3_cor = steps2coord(StepperMotionType([82,-10],[0]),xorig=[-8.8],d=2.4, Rh=80.0, Re=10.0)
            @test coord2steps(p1,p2,xorig=[-8.8],d=2.4, Rh=80.0, Re=10.0)[1].step ≈ [41,16] atol=1
            @test coord2steps(p1,p3,xorig=[-8.8],d=2.4, Rh=80.0, Re=10.0)[1].step ≈ [82,-10] atol=1
            @test coord2steps(p2,p3,xorig=[-8.8],d=2.4, Rh=80.0, Re=10.0)[1].step ≈ [41,-25] atol=1
            @test coord2steps(p1,p3,xorig=[-8.8],d=2.4, Rh=80.0, Re=10.0)[1].strategy == [1,2]
            @test coord2steps(p3,p2,xorig=[-8.8],d=2.4, Rh=80.0, Re=10.0)[1].strategy == [1,2]
            @test convert(Array, p2_cor) ≈ convert(Array, p2) atol=0.2
            @test convert(Array, p3_cor) ≈ convert(Array, p3) atol=0.2
        end
        @testset "Two cylinders" begin
            p1 = TwoCylType(0.0,0.0,2.4,0.0)
            p2 = TwoCylType(10.0,0.0,2.4,0.0)
            p3 = TwoCylType(20.0,0.0,2.4,0.0)
            p2_cor = steps2coord(StepperMotionType([41,16,41,47],[0]),xorig=[-8.8,-18.8],d=2.4, Rh=80.0, Re=10.0)
            p3_cor = steps2coord(StepperMotionType([82,-10,82,53],[0]),xorig=[-8.8,-18.8],d=2.4, Rh=80.0, Re=10.0)
            @test coord2steps(p1,p2,xorig=[-8.8,-18.8],d=2.4, Rh=80.0, Re=10.0)[1].step ≈ [41,16,41,47] atol=1
            @test coord2steps(p1,p3,xorig=[-8.8,-18.8],d=2.4, Rh=80.0, Re=10.0)[1].step ≈ [82,-10,82,53] atol=1
            @test coord2steps(p2,p3,xorig=[-8.8,-18.8],d=2.4, Rh=80.0, Re=10.0)[1].step ≈ [41,-25,41,6] atol=1
            @test coord2steps(p1,p3,xorig=[-8.8,-18.8],d=2.4, Rh=80.0, Re=10.0)[1].strategy == [3,4,1,2]
            @test coord2steps(p3,p2,xorig=[-8.8,-18.8],d=2.4, Rh=80.0, Re=10.0)[1].strategy == [1,2,3,4]
            @test convert(Array, p2_cor) ≈ convert(Array, p2) atol=0.2  
            @test convert(Array, p3_cor) ≈ convert(Array, p3) atol=0.2
        end
        @testset "Three cyliders" begin
            p1 = ThreeCylType(0.0,0.0,2.4,0.0,2.4,0.0)
            p2 = ThreeCylType(10.0,0.0,2.4,0.0,2.4,0.0)
            p3 = ThreeCylType(20.0,0.0,2.4,0.0,2.4,0.0)
            p2_cor = steps2coord(StepperMotionType([41,16,41,47,42,80],[0]),xorig=[-8.8,-18.8,-28.8],d=2.4, Rh=80.0, Re=10.0)
            p3_cor = steps2coord(StepperMotionType([82,-10,82,53,83,117],[0]),xorig=[-8.8,-18.8,-28.8],d=2.4,Rh=80.0, Re=10.0)
            @test coord2steps(p1,p2,xorig=[-8.8,-18.8,-28.8],d=2.4, Rh=80.0, Re=10.0)[1].step ≈ [41,16,41,47,42,80] atol=1
            @test coord2steps(p1,p3,xorig=[-8.8,-18.8,-28.8],d=2.4, Rh=80.0, Re=10.0)[1].step ≈ [82,-10,82,53,83,117] atol=1
            @test coord2steps(p2,p3,xorig=[-8.8,-18.8,-28.8],d=2.4, Rh=80.0, Re=10.0)[1].step ≈ [41,-25,41,6,41,37] atol=1
            @test coord2steps(p1,p3,xorig=[-8.8,-18.8,-28.8],d=2.4, Rh=80.0, Re=10.0)[1].strategy ≈ [5,6,3,4,1,2] 
            @test coord2steps(p3,p2,xorig=[-8.8,-18.8,-28.8],d=2.4, Rh=80.0, Re=10.0)[1].strategy ≈ [1,2,3,4,5,6]
            @test convert(Array, p2_cor) ≈ convert(Array, p2) atol=0.2
            @test convert(Array, p3_cor) ≈ convert(Array, p3) atol=0.2
        end
    end
end

function bitsconversion_tests()
    @testset "Bits and coordinates transformations" begin
        
        c11 = OneCylType(1.0, 2.0)
        c11_bit = cylpos2bits(c11, size=[6,5])
        @test bits2cylpos(c11_bit) == c11
        c12 = OneCylType(1.0, -2.0)
        c12_bit = cylpos2bits(c12, size=[6,5])
        @test bits2cylpos(c12_bit) == c12

        c21 = TwoCylType(2.0,1.0,3.0,4.0)
        c21_bit = cylpos2bits(c21, size=[6,5,4,4])
        @test bits2cylpos(c21_bit) == c21
        c22 = TwoCylType(2.0,-1.0,3.0,-4.0)
        c22_bit = cylpos2bits(c22, size=[6,5,4,4])
        @test bits2cylpos(c22_bit) == c22

        c31 = ThreeCylType(10.0,2.0,2.0,3.0,2.0,1.0)
        c31_bit = cylpos2bits(c31, size=[6,5,4,4,4,4])
        @test bits2cylpos(c31_bit) == c31
        c32 = ThreeCylType(10.0,-2.0,2.0,-3.0,2.0,-1.0)
        c32_bit = cylpos2bits(c32, size=[6,5,4,4,4,4])
        @test bits2cylpos(c32_bit) == c32
        c33 = ThreeCylType(0,0,0,0,0,0)
        c33_bit = cylpos2bits(c33, size=[6,5,4,4,4,4])
        @test bits2cylpos(c33_bit) == c33

    end
end


# TEST ALL
function testall()
    coordsteps_test()
    bitsconversion_tests()
end
