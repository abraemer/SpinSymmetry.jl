@testset "basis.jl" begin
    @testset "ZBasis" begin
        @testset "FullZBasis" begin
            @test_throws ArgumentError zbasis(-1)
            @test zbasis(5) isa FullZBasis
            @test SpinSymmetry._indices(zbasis(5)) == 1:2^5

            @test zbasis(5) == FullZBasis(5)
            @test zbasis(6) != FullZBasis(5)
            @test hash(zbasis(5)) == hash(FullZBasis(5))
            @test hash(zbasis(6)) != hash(FullZBasis(5))
        end

        @testset "ZBlockBasis" begin
            @test_throws ArgumentError zbasis(-1, 0)
            @test_throws ArgumentError zbasis(2,-1)
            @test_throws ArgumentError zbasis(2,3)

            @test zbasis(2,1) isa ZBlockBasis

            @test zbasis(5,2) == ZBlockBasis(5,2)
            @test zbasis(6,2) != ZBlockBasis(5,2)
            @test zbasis(6,3) != ZBlockBasis(6,2)
            @test hash(zbasis(5,2)) == hash(ZBlockBasis(5,2))
            @test hash(zbasis(6,2)) != hash(ZBlockBasis(5,2))
            @test hash(zbasis(6,3)) != hash(ZBlockBasis(6,2))

            # small correctness check
            @test SpinSymmetry._indices(zbasis(2,0)) == [1]
            @test SpinSymmetry._indices(zbasis(2,1)) == [2,3]
            @test SpinSymmetry._indices(zbasis(2,2)) == [4]


            digitsum(k) = sum(parse.(Int, [string(k-1; base=2)...]))
            for k in 0:10
                zblock = zbasis(10, k)
                @test all(digitsum.(SpinSymmetry._indices(zblock)) .== k)
            end
        end
    end

    @testset "symmetrized_basis" begin
        # different argument forms
        @test symmetrized_basis(zbasis(6), Flip(6), 0) == symmetrized_basis(6, Flip(6), 0)
        @test symmetrized_basis(zbasis(6, 2), Flip(6), 0) == symmetrized_basis(6, 2, Flip(6), 0)

        # equality
        @test symmetrized_basis(6, Flip(6), 0, Shift(6), 1) == symmetrized_basis(6, Shift(6), 1, Flip(6), 0)
        @test symmetrized_basis(6, Flip(6), 0, Shift(6, 5), 1) == symmetrized_basis(6, Shift(6, -1), 1, Flip(6), 0)
        @test symmetrized_basis(6, Flip(6), 0, Shift(6), 0) != symmetrized_basis(6, Shift(6), 1, Flip(6), 0)
        @test symmetrized_basis(6, Flip(6), 1, Shift(6), 0) != symmetrized_basis(6, Shift(6), 1, Flip(6), 0)

        let b1 = symmetrized_basis(6, Flip(6), 0, Shift(6), 1),
            b2 = symmetrized_basis(6, Shift(6), 1, Flip(6), 0)
            @test SpinSymmetry._phase_factors(1:2^6, b1.symmetries, b1.sectors)  ==
                SpinSymmetry._phase_factors(1:2^6, b2.symmetries, b2.sectors)
        end
    end

    @testset "basissize" begin
        let basis = zbasis(10)
            @test basissize(basis) == length(SpinSymmetry._indices(basis))
        end

        for k in 0:10
            let basis = zbasis(10, k)
                @test basissize(basis) == length(SpinSymmetry._indices(basis))
            end
        end

        @test basissize(symmetrized_basis(10, Flip(10), 0)) == 2^9
        @test basissize(symmetrized_basis(10, 5,  Flip(10), 0)) == binomial(10, 5)/2
        @test basissize(symmetrized_basis(10, 4,  Flip(10), 0)) == binomial(10, 4)
    end

    @testset "symmetrize stuff" begin
        using LinearAlgebra, Random
        Random.seed!(5)

        N = 5

        state = normalize!(rand(2^N))
        operator = rand(2^N, 2^N)
        base = symmetrized_basis(N, Flip(N), 0)
        mid = 2^(N-1)
        front = 1:mid
        back = 2^N:-1:mid+1

        function compare_bases(b1, b2)
            return all(SpinSymmetry._indices(b1.basis) .== SpinSymmetry._indices(b2.basis)) &&
                   all(b1.symmetries .== b2. symmetries) &&
                   all(b1.sectors .== b2.sectors)
        end

        @test compare_bases(base, symmetrized_basis(zbasis(N), Flip(N), 0))
        @test compare_bases(symmetrized_basis(zbasis(N,1), Flip(N), 1), symmetrized_basis(N, 1, Flip(N), 1))

        symm_state = symmetrize_state(state, base)
        @test symm_state ≈ symmetrize_state(state, N, Flip(N), 0)
        @test symm_state ≈ 1/√2 .* (state[front] .+ state[back])

        asymm_state = symmetrize_state(state, N, Flip(N), 1)
        @test asymm_state ≈ 1/√2 .* (state[front] .- state[back])

        @test_throws ArgumentError symmetrize_state(zeros(2^N+1), base)


        symm_operator = symmetrize_operator(operator, base)
        @test symm_operator ≈ symmetrize_operator(operator, N, Flip(N), 0)
        @test symm_operator ≈ 1/2 .* (operator[front, front] .+ operator[back, back] .+ operator[front, back] .+ operator[back, front])

        asymm_operator = symmetrize_operator(operator, N, Flip(N), 1)
        @test asymm_operator ≈ 1/2 .* (operator[front, front] .+ operator[back, back] .- operator[front, back] .- operator[back, front])

        @test_throws ArgumentError symmetrize_operator(zeros(2^N+1, 2^N+1), base)
        @test_throws ArgumentError symmetrize_operator(zeros(2^(N-1), 2^(N+1)), base)
    end
end
