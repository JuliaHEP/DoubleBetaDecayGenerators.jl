using Test
using DoubleBetaDecayGenerators
using Random
using LinearAlgebra
using Logging

@testset "Test control of random seed" begin
    rng = Xoshiro()

    Random.seed!(rng, 999)

    dat_2nu = TwoNuDBDData(:Ge76)
    E1, E2, cosθ12 = rand(rng, dat_2nu)
    @test E1 ≈ 374.2355832686694 atol=1e-15
    @test E2 ≈ 1251.059139039826 atol=1e-15
    @test cosθ12 ≈ 0.9180839309783714 atol=1e-15

    Random.seed!(rng, 999)

    dat_0nu = ZeroNuDBDData(:Ge76)
    E1, E2, cosθ12 = rand(rng, dat_0nu)
    @test E1 ≈ 2016.2355832686694 atol=1e-15
    @test E2 ≈ 22.764416731330584 atol=1e-15
    @test cosθ12 ≈ -0.9113939723175647 atol=1e-15
end

@testset "acute_inv_cdf" begin
    inv_cdf = DoubleBetaDecayGenerators.acute_inv_cdf
    # q = 0 is the uniform distribution on [-1, 1], not NaN
    @test inv_cdf(0.0, 0) == -1
    @test inv_cdf(0.5, 0) == 0
    @test inv_cdf(1.0, 0) == 1
    # endpoints map to [-1, 1] for any slope
    for q in (-0.9, -0.3, 0.4, 1.0)
        @test inv_cdf(0.0, q) ≈ -1
        @test inv_cdf(1.0, q) ≈ 1
    end
end

@testset "Debug logging in constructors" begin
    @test_logs (:debug, r"Reading raw data") min_level=Logging.Debug match_mode=:any TwoNuDBDData(:Ge76)
    @test_logs (:debug, r"Reading raw data") min_level=Logging.Debug match_mode=:any ZeroNuDBDData(:Ge76)
end

@testset "Legacy.rand_2nbb" begin
    Legacy = DoubleBetaDecayGenerators.Legacy
    data = TwoNuDBDData(:Ge76)

    # momenta carry the sampled energies and opening angle
    E1, E2, cosθ12 = rand(Xoshiro(7), data)
    p1, p2 = Legacy.rand_2nbb(Xoshiro(7), data)
    @test norm(p1) ≈ Legacy.momentum(E1)
    @test norm(p2) ≈ Legacy.momentum(E2)
    @test dot(p1, p2) / (norm(p1) * norm(p2)) ≈ cosθ12

    # the RNG controls every random draw
    @test Legacy.rand_2nbb(Xoshiro(1), data) == Legacy.rand_2nbb(Xoshiro(1), data)
end

@testset "Legacy.dk0gen" begin
    mktempdir() do dir
        output = joinpath(dir, "ge76.dk0")
        DoubleBetaDecayGenerators.Legacy.dk0gen(3; output, isotope=:Ge76, rng=Xoshiro(1))
        lines = readlines(output)
        @test any(contains("event type: Ge76"), lines)
        @test count(startswith("  3 "), lines) == 6
    end
end
