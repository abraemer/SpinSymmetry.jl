@testset "basis.jl" begin
    @testset "ZBasis" begin
        @test_throws ArgumentError zbasis(-1)
        @test SpinSymmetry._indices(zbasis(5)) == 1:2^5

        @test_throws ArgumentError zbasis(-1, 0)
        @test_throws ArgumentError zbasis(2,-1)
        @test_throws ArgumentError zbasis(2,3)

        @test SpinSymmetry._indices(zbasis(2,0)) == [1]
        @test SpinSymmetry._indices(zbasis(2,1)) == [2,3]
        @test SpinSymmetry._indices(zbasis(2,2)) == [4]
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
