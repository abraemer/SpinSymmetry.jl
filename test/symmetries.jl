@testset "symmetries.jl" begin

    function order_test(symm, N)
        indices = 0:2^N-1 #indices
        current = indices
        order = SpinSymmetry._order(symm)
        for i in 1:order-1
            current = symm.(current)
            if current == indices
                println("Order too long! Reported: $order, brute-forced: $i. Symmetry: $symm")
                return false
            end
        end
        ok = symm.(current) == indices
        !ok && println("Order too short! Reported: $order, Symmetry: $symm")
        return ok
    end

    @testset "Shift" begin
        # test default arg
        @test Shift(10) == Shift(10,1)

        # small correctness test
        @test Shift(2).(0:3) == [0,2,1,3]

        # order test
        for N in 2:10
            for amount in 1:div(N,2)
                @test order_test(Shift(N, amount), N)
            end
        end
    end

    @testset "Flip" begin
        # test default arg
        @test Flip(10) == Flip(collect(1:10))

        # small correctness test
        @test Flip(2).(0:3) == [3,2,1,0]
        @test Flip([1,2]).(0:7) == vcat([3,2,1,0], 4 .+ [3,2,1,0])

        # order test
        for N in 2:10
            @test order_test(Flip(N), N)
            @test order_test(Flip(div(N,2)), N)
        end
    end

    @testset "Swap" begin
        # small correctness test
        @test Swap(1,2).(0:3) == [0,2,1,3]
        @test Swap(2,2).(0:31) == 0:31 # identity

        # order test
        for N in 3:2:16
            pos1 = round(Int, N/4)
            pos2 = round(Int, 3*N/4)
            @test order_test(Swap(pos1,pos2), N)
            @test order_test(Swap(1,N), N)
            @test Swap(pos1,pos2).(0:2^N-1) == Swap(pos2,pos1).(0:2^N-1)
        end
    end

    # rather pointless
    # GenericSymmetry is as correct as the values it contains...
    @testset "GenericSymmetry" begin
        symm = GenericSymmetry(identity, 1)

        @test symm.(0:255) == 0:255
        @test SpinSymmetry._order(symm) == 1
    end

end
