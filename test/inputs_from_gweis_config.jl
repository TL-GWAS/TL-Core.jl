module TestGweisEstimands

using Test
using SnpArrays
using TargeneCore
using Arrow
using DataFrames
using Serialization
using TMLE
using CSV

TESTDIR = joinpath(pkgdir(TargeneCore), "test")

include(joinpath(TESTDIR, "testutils.jl"))

function get_summary_stats(estimands)
    outcomes = [TargeneCore.get_outcome(Ψ) for Ψ in estimands]
    results = DataFrame(ESTIMAND = estimands, OUTCOME = outcomes)
    return sort(combine(groupby(results, :OUTCOME), nrow), :OUTCOME)
end

function check_estimands_levels_interactions(estimands)
    for Ψ in estimands
        # If the two components are present, the first is the 0 -> 1 and the second is the 1 -> 2
        # The variant should always be the last key
        variant = last(collect(keys(Ψ.args[1].treatment_values)))
        if length(Ψ.args) == 2
            @test Ψ.args[1].treatment_values[variant] == (control = 0x00, case = 0x01)
            @test Ψ.args[2].treatment_values[variant] == (control = 0x01, case = 0x02)
        else
            # Otherwise we check they are one or the other
            arg = only(Ψ.args)
            @test arg.treatment_values[variant]==(control = 0x00, case = 0x01) ||
            arg.treatment_values[variant]==( control = 0x01, case = 0x02)
        end
   end
end

@testset "Test inputs_from_config gweis: no positivity constraint" begin
    tmpdir = mktempdir()
    copy!(ARGS, [
        "estimation-inputs",
        joinpath(TESTDIR, "data", "config_gweis_first_order.yaml"),
        string("--traits-file=", joinpath(TESTDIR, "data", "ukbb_traits.csv")),
        string("--pcs-file=", joinpath(TESTDIR, "data", "ukbb_pcs.csv")),
        string("--genotypes-prefix=", joinpath(TESTDIR, "data", "ukbb", "genotypes" , "ukbb_1.")),
        string("--outprefix=", joinpath(tmpdir, "final")), 
        "--batchsize=5",
        "--verbosity=0",
        "--positivity-constraint=0"
    ])
    TargeneCore.julia_main()
    # Check dataset
    dataset = DataFrame(Arrow.Table(joinpath(tmpdir, "final.data.arrow")))
    @test size(dataset) == (1940, 886)

    # Check estimands
    estimands = []
    for file in readdir(tmpdir, join=true)
        if endswith(file, "jls")
            append!(estimands, deserialize(file).estimands)
        end
    end
    @test all(e isa JointEstimand for e in estimands)

    # There are 875 variants in the dataset
    summary_stats = get_summary_stats(estimands)
    @test summary_stats == DataFrame(
        OUTCOME = [:BINARY_1, :BINARY_2, :CONTINUOUS_1, :CONTINUOUS_2, :TREAT_1], 
        nrow = repeat([875], 5)
    )

    check_estimands_levels_interactions(estimands)
end


@testset "Test inputs_from_config gweis: positivity constraint" begin
    tmpdir = mktempdir()
    copy!(ARGS, [
        "estimation-inputs",
        joinpath(TESTDIR, "data", "config_gweis_first_order.yaml"),
        string("--traits-file=", joinpath(TESTDIR, "data", "ukbb_traits.csv")),
        string("--pcs-file=", joinpath(TESTDIR, "data", "ukbb_pcs.csv")),
        string("--genotypes-prefix=", joinpath(TESTDIR, "data", "ukbb", "genotypes" , "ukbb_1.")),
        string("--outprefix=", joinpath(tmpdir, "final")), 
        "--batchsize=5",
        "--verbosity=0",
        "--positivity-constraint=0.2"
    ])
    TargeneCore.julia_main()
    # Check dataset
    dataset = DataFrame(Arrow.Table(joinpath(tmpdir, "final.data.arrow")))
    @test size(dataset) == (1940, 886)
    # Check estimands
    estimands = []
    for file in readdir(tmpdir, join=true)
        if endswith(file, "jls")
            append!(estimands, deserialize(file).estimands)
        end
    end
    # The positivity constraint reduces the number of variants
    @test all(e isa JointEstimand for e in estimands)
    summary_stats = get_summary_stats(estimands)
    @test summary_stats == DataFrame(
        OUTCOME = [:BINARY_1, :BINARY_2, :CONTINUOUS_1, :CONTINUOUS_2, :TREAT_1], 
        nrow = repeat([142], 5)
    )
   
    check_estimands_levels_interactions(estimands)
end

@testset "Test inputs_from_config gweis: no positivity constraint and four-point interaction" begin
    tmpdir = mktempdir()
    copy!(ARGS, [
        "estimation-inputs",
        joinpath(TESTDIR, "data", "config_gweis_higher_order.yaml"),
        string("--traits-file=", joinpath(TESTDIR, "data", "ukbb_traits.csv")),
        string("--pcs-file=", joinpath(TESTDIR, "data", "ukbb_pcs.csv")),
        string("--genotypes-prefix=", joinpath(TESTDIR, "data", "ukbb", "genotypes" , "ukbb_1.")),
        string("--outprefix=", joinpath(tmpdir, "final")), 
        "--batchsize=5",
        "--verbosity=0",
        "--positivity-constraint=0"
    ])
    TargeneCore.julia_main()
    # Check dataset
    dataset = DataFrame(Arrow.Table(joinpath(tmpdir, "final.data.arrow")))
    @test size(dataset) == (1940, 886)

    # Check estimands
    estimands = []
    for file in readdir(tmpdir, join=true)
        if endswith(file, "jls")
            append!(estimands, deserialize(file).estimands)
        end
    end
    @test all(e isa JointEstimand for e in estimands)

    # There are 875 variants in the dataset
    summary_stats = get_summary_stats(estimands)
    @test summary_stats == DataFrame(
        OUTCOME = [:CONTINUOUS_1, :CONTINUOUS_2, :TREAT_1], 
        nrow = repeat([875], 3)
    )

    check_estimands_levels_interactions(estimands)
end

@testset "Test inputs_from_config gweis: positivity constraint and four-point interaction" begin
    tmpdir = mktempdir()
    copy!(ARGS, [
        "estimation-inputs",
        joinpath(TESTDIR, "data", "config_gweis_higher_order.yaml"),
        string("--traits-file=", joinpath(TESTDIR, "data", "ukbb_traits.csv")),
        string("--pcs-file=", joinpath(TESTDIR, "data", "ukbb_pcs.csv")),
        string("--genotypes-prefix=", joinpath(TESTDIR, "data", "ukbb", "genotypes" , "ukbb_1.")),
        string("--outprefix=", joinpath(tmpdir, "final")), 
        "--batchsize=5",
        "--verbosity=0",
        "--positivity-constraint=0.02"
    ])
    TargeneCore.julia_main()
    # Check dataset
    dataset = DataFrame(Arrow.Table(joinpath(tmpdir, "final.data.arrow")))
    @test size(dataset) == (1940, 886)

    # Check estimands
    estimands = []
    for file in readdir(tmpdir, join=true)
        if endswith(file, "jls")
            append!(estimands, deserialize(file).estimands)
        end
    end
    @test all(e isa JointEstimand for e in estimands)

    # There are 784 treatments in the dataset after positivity_constraint
    summary_stats = get_summary_stats(estimands)
    @test summary_stats == DataFrame(
        OUTCOME = [:CONTINUOUS_1, :CONTINUOUS_2, :TREAT_1], 
        nrow = repeat([784], 3)
    )

    check_estimands_levels_interactions(estimands)
end

end
true