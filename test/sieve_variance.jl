module TestSievePlateau

using TargeneCore
using Test
using DataFrames
using CSV 
using JLD2
using TMLE
using CategoricalArrays
using TMLECLI
using StableRNGs
using Distributions
using LogExpFunctions

TESTDIR = joinpath(pkgdir(TargeneCore), "test")

function statistical_estimands_only_config()
    configuration = Configuration(
        estimands=[
            AIE(
                outcome = Symbol("CONTINUOUS, OUTCOME"), 
                treatment_values = (
                    T1 = (case = true, control = false), 
                    T2 = (case = true, control = false)), 
                treatment_confounders = (T1 = (:W1, :W2), T2 = (:W1, :W2)), 
                outcome_extra_covariates = (:C1,)
            ),
            ATE(
                outcome = Symbol("CONTINUOUS, OUTCOME"), 
                treatment_values = (T1 = (case = true, control = false),), 
                treatment_confounders = (T1 = (:W1, :W2),), 
                outcome_extra_covariates = ()
            ),
            AIE(
                outcome = Symbol("CONTINUOUS, OUTCOME"), 
                treatment_values = (
                    T1 = (case = true, control = false), 
                    T2 = (case = false, control = true)
                ), 
                treatment_confounders = (T1 = (:W1, :W2), T2 = (:W1, :W2)), 
                outcome_extra_covariates = ()
            ),
            AIE(
                outcome = Symbol("BINARY/OUTCOME"), 
                treatment_values = (
                    T1 = (case = true, control = false), 
                    T2 = (case = false, control = true)
                ), 
                treatment_confounders = (T1 = (:W1, :W2), T2 = (:W1, :W2)), 
                outcome_extra_covariates = (:C1,)
            ),
            AIE(
                outcome = Symbol("BINARY/OUTCOME"), 
                treatment_values = (
                    T1 = (case = true, control = false), 
                    T2 = (case = true, control = false)), 
                treatment_confounders = (T1 = (:W1, :W2), T2 = (:W1, :W2)), 
                outcome_extra_covariates = (:C1,)
            ),
            CM(
                outcome = Symbol("COUNT_OUTCOME"), 
                treatment_values = (
                    T1 = true, 
                    T2 = false), 
                treatment_confounders = (T1 = (:W1, :W2), T2 = (:W1, :W2)),
                outcome_extra_covariates = (:C1,)
            )
        ]
    )
    return configuration
end

function causal_and_joint_estimands_config()
    ATE₁ = ATE(
        outcome = Symbol("CONTINUOUS, OUTCOME"), 
        treatment_values = (T1 = (case = true, control = false),), 
    )
    ATE₂ = ATE(
        outcome = Symbol("CONTINUOUS, OUTCOME"), 
        treatment_values = (T1 = (case = false, control = true),), 
    )
    joint = JointEstimand(ATE₁, ATE₂)
    scm = StaticSCM(
        outcomes = ["CONTINUOUS, OUTCOME"],
        treatments = ["T1"],
        confounders = [:W1, :W2]
    )
    configuration = Configuration(
        estimands = [ATE₁, ATE₂, joint],
        scm       = scm
    )
    return configuration
end

function write_sieve_dataset(datafile, sample_ids)
    rng = StableRNG(123)
    n = size(sample_ids, 1)
    # Confounders
    W₁ = rand(rng, Uniform(), n)
    W₂ = rand(rng, Uniform(), n)
    # Covariates
    C₁ = rand(rng, n)
    # Treatment | Confounders
    T₁ = rand(rng, Uniform(), n) .< logistic.(0.5sin.(W₁) .- 1.5W₂)
    T₂ = rand(rng, Uniform(), n) .< logistic.(-3W₁ - 1.5W₂)
    # target | Confounders, Covariates, Treatments
    μ = 1 .+ 2W₁ .+ 3W₂ .- 4C₁.*T₁ .+ T₁ + T₂.*W₂.*T₁
    y₁ = μ .+ rand(rng, Normal(0, 0.01), n)
    y₂ = rand(rng, Uniform(), n) .< logistic.(μ)
    # Add some missingness
    y₂ = vcat(missing, y₂[2:end])

    dataset = DataFrame(
        SAMPLE_ID = string.(sample_ids),
        T1 = categorical(T₁),
        T2 = categorical(T₂),
        W1 = W₁, 
        W2 = W₂,
        C1 = C₁,
    )

    dataset[!, "CONTINUOUS, OUTCOME"] = y₁
    dataset[!, "BINARY/OUTCOME"] = categorical(y₂)
    dataset[!, "COUNT_OUTCOME"] = rand(rng, [1, 2, 3, 4], n)

    CSV.write(datafile, dataset)
end

"""
This function mimics the outputs generated by multiple estimation runs.
"""
function build_tmle_output_file(workdir, configuration, sample_ids_df; pvalue_threshold=1.)
    n_div_2 = length(configuration.estimands) ÷ 2
    # Write Dataset
    dataset_file = joinpath(workdir, "dataset.csv")
    write_sieve_dataset(dataset_file, sample_ids_df.SAMPLE_ID)
    # Run estimation 1
    tmle_output_1 = joinpath(workdir, "tmle_output_1.hdf5")
    outputs_1 = TMLECLI.Outputs(hdf5=tmle_output_1)
    estimands_file_1 = joinpath(workdir, "configuration_1.json")
    TMLE.write_json(
        estimands_file_1, 
        TMLE.Configuration(estimands=configuration.estimands[1:n_div_2], scm=configuration.scm)
    )
    tmle(dataset_file; 
        estimands=estimands_file_1, 
        estimators="wtmle-ose--glm", 
        outputs=outputs_1,
        pvalue_threshold=pvalue_threshold,
        save_sample_ids=true
    )
    # Run estimation 2
    tmle_output_2 = joinpath(workdir, "tmle_output_2.hdf5")
    outputs_2 = TMLECLI.Outputs(hdf5=tmle_output_2)
    estimands_file_2 = joinpath(workdir, "configuration_2.json")
    TMLE.write_json(
        estimands_file_2, 
        TMLE.Configuration(estimands=configuration.estimands[n_div_2+1:end], scm=configuration.scm)
    )
    tmle(dataset_file; 
        estimands=estimands_file_2, 
        estimators="wtmle-ose--glm", 
        outputs=outputs_2,
        pvalue_threshold=pvalue_threshold,
        save_sample_ids=true
    )
end

function basic_variance_implementation(matrix_distance, influence_curve, n_obs)
    variance = 0.f0
    n_samples = size(influence_curve, 1)
    for i in 1:n_samples
        for j in 1:n_samples
            variance += matrix_distance[i, j]*influence_curve[i]* influence_curve[j]
        end
    end
    variance/n_obs
end

function distance_vector_to_matrix!(matrix_distance, vector_distance, n_samples)
    index = 1
    for i in 1:n_samples
        for j in 1:i
            # enforce indicator = 1 when i =j 
            if i == j
                matrix_distance[i, j] = 1
            else
                matrix_distance[i, j] = vector_distance[index]
                matrix_distance[j, i] = vector_distance[index]
            end
            index += 1
        end
    end
end

function test_initial_output(output, expected_output)
    # Metadata columns
    for col in [:PARAMETER_TYPE, :TREATMENTS, :CASE, :CONTROL, :OUTCOME, :CONFOUNDERS, :COVARIATES]
        for index in eachindex(output[!, col])
            if expected_output[index, col] === missing
                @test expected_output[index, col] === output[index, col]
            else
                @test expected_output[index, col] == output[index, col]
            end
        end
    end
end

@testset "Test readGRM" begin
    prefix = joinpath(TESTDIR, "data", "grm", "test.grm")
    GRM, ids = TargeneCore.readGRM(prefix)
    @test eltype(ids.SAMPLE_ID) == String
    @test size(GRM, 1) == 18915
    @test size(ids, 1) == 194
end

@testset "Test build_work_list" begin
    sample_ids = TargeneCore.GRMIDs(joinpath(TESTDIR, "data", "grm", "test.grm.id"))
    configuration = statistical_estimands_only_config()

    # CASE_1: pval = 1.
    tmpdir = mktempdir()
    tmle_prefix = joinpath(tmpdir, "tmle_output")
    build_tmle_output_file(tmpdir, configuration, sample_ids; pvalue_threshold=1.)
    results, influence_curves, n_obs = TargeneCore.build_work_list(tmle_prefix, sample_ids)
    # Check n_obs
    @test n_obs == [194, 194, 194, 193, 193, 194]
    # Check influence curves
    expected_influence_curves = [size(r.IC, 1) == 194 ? r.IC : vcat(0, r.IC) for r in results]
    for rowindex in 1:6
        @test convert(Vector{Float32}, expected_influence_curves[rowindex]) == influence_curves[rowindex, :]
    end
    # Check results
    all(x isa TMLE.TMLEstimate for x in results)
    all(size(x.IC, 1) > 0 for x in results)

    # CASE_2: pval = 0.1
    tmpdir = mktempdir()
    tmle_prefix = joinpath(tmpdir, "tmle_output")
    build_tmle_output_file(tmpdir, configuration, sample_ids; pvalue_threshold= 0.1)

    results, influence_curves, n_obs = TargeneCore.build_work_list(tmle_prefix, sample_ids)
    # Check n_obs
    @test n_obs == [194, 194, 193, 193, 194]
    # Check influence curves
    expected_influence_curves = [size(r.IC, 1) == 194 ? r.IC : vcat(0, r.IC) for r in results]
    for rowindex in 1:4
        @test convert(Vector{Float32}, expected_influence_curves[rowindex]) == influence_curves[rowindex, :]
    end
    # Check results
    @test all(x isa TMLE.TMLEstimate for x in results)
    @test all(size(x.IC, 1) > 0 for x in results)
end

@testset "Test bit_distance" begin
    sample_grm = Float32[-0.6, -0.8, -0.25, -0.3, -0.1, 0.1, 0.7, 0.5, 0.2, 1.]
    nτs = 6
    τs = TargeneCore.default_τs(nτs, max_τ=0.75)
    @test τs == Float32[0.0, 0.15, 0.3, 0.45, 0.6, 0.75]
    τs = TargeneCore.default_τs(nτs)
    @test τs == Float32[0., 0.4, 0.8, 1.2, 1.6, 2.0]
    d = TargeneCore.bit_distances(sample_grm, τs)
    @test d == [0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0
                0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  1.0
                0.0  0.0  0.0  0.0  0.0  0.0  1.0  1.0  1.0  1.0
                0.0  0.0  0.0  0.0  1.0  1.0  1.0  1.0  1.0  1.0
                1.0  0.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
                1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0]
end

@testset "Test aggregate_variances" begin
    # 2 influence curves containing 5 individuals
    influence_curves = [1. 2. 3. 4. 5.
                        6. 7. 8. 9. 10.]
    # distance indicator with 3 τs and corresponding to row 4
    indicator = [1. 0. 0. 0.2
                 0. 0. 1. 1.
                 1. 0. 1. 1.]
    sample = 4
    var_ = TargeneCore.aggregate_variances(influence_curves, indicator, sample)
    @test var_ == [24.0  189.0
                   40.0  225.0
                   48.0  333.0]
end

@testset "Test normalize!" begin
    # 2 τs and 3 curves
    n_obs = [10, 10, 100]
    variances = [1. 2. 3.
                 4. 5. 6.]
    TargeneCore.normalize!(variances, n_obs)
    @test variances == [0.1 0.2 0.03
                        0.4 0.5 0.06]
end

@testset "Test compute_variances" begin
    n_curves = 3
    n_samples = 5
    nτs = 5
    n_obs = [3, 4, 4]
    τs = TargeneCore.default_τs(nτs)
    # The GRM has 15 lower triangular elements
    grm = Float32[0.4, 0.1, 0.5, 0.2, -0.2, 0.6, 0.3, -0.6, 
                  0.4, 0.3, 0.6, 0.3, 0.7, 0.3, 0.1]
    influence_curves = Float32[0.1 0. 0.1 0.3 0.
                               0.1 0.2 0.1 0.0 0.2
                               0.0 0. 0.1 0.3 0.2]
                  
    
    variances = TargeneCore.compute_variances(influence_curves, grm, τs, n_obs)
    @test size(variances) == (nτs, n_curves)

    # when τ=2, all elements are used
    for curve_id in 1:n_curves
        s = sum(influence_curves[curve_id, :])
        var = sum(s*influence_curves[curve_id, i] for i in 1:n_samples)/n_obs[curve_id]
        @test variances[end, curve_id] ≈ var
    end

    # Decreasing variances with τ as all inf curves are positives
    for nτ in 1:nτs-1
        @test all(variances[nτ, :] .<= variances[nτ+1, :])
    end

    # Check against basic_variance_implementation
    matrix_distance = zeros(Float32, n_samples, n_samples)
    for τ_id in 1:nτs
        vector_distance = TargeneCore.bit_distances(grm, [τs[τ_id]])
        distance_vector_to_matrix!(matrix_distance, vector_distance, n_samples)
        for curve_id in 1:n_curves
            influence_curve = influence_curves[curve_id, :]
            var_ = basic_variance_implementation(matrix_distance, influence_curve, n_obs[curve_id])
            @test variances[τ_id, curve_id] ≈ var_
        end
    end

    # Check by hand for a single τ=0.5
    @test variances[2, :] ≈ Float32[0.03666667, 0.045, 0.045]
end

@testset "Test grm_rows_bounds" begin
    n_samples = 5
    grm_bounds = TargeneCore.grm_rows_bounds(n_samples)
    @test grm_bounds == [1 => 1
                         2 => 3
                         4 => 6
                         7 => 10
                         11 => 15]
end

@testset "Test corrected_stderrors" begin
    variances = [
        1. 2. 6.
        4. 5. 3.
    ]
    stderrors = TargeneCore.corrected_stderrors(variances)
    # sanity check
    @test stderrors == sqrt.([4., 5., 6.])
end

@testset "Test SVP" begin
    # Generate data
    sample_ids = TargeneCore.GRMIDs(joinpath(TESTDIR, "data", "grm", "test.grm.id"))
    configuration = statistical_estimands_only_config()

    tmpdir = mktempdir()
    tmle_prefix = joinpath(tmpdir, "tmle_output")

    build_tmle_output_file(tmpdir, configuration, sample_ids; pvalue_threshold=0.1)


    # Using the main command
    output = joinpath(tmpdir, "svp.hdf5")
    main([
        "svp", 
        joinpath(tmpdir, "tmle_output"),
        "--out", output,
        "--grm-prefix", joinpath(TESTDIR, "data", "grm", "test.grm"), 
        "--max-tau", "0.75"
    ])

    io = jldopen(output)
    # Check τs
    @test io["taus"] == TargeneCore.default_τs(10; max_τ=0.75)
    # Check variances
    @test size(io["variances"]) == (10, 5)
    # Check results
    svp_results = io["results"]
    @test size(svp_results, 1) == 5

    tmleout1 = jldopen(x -> x["Batch_1"], string(tmle_prefix, "_1.hdf5"))
    tmleout2 = jldopen(x -> x["Batch_1"], string(tmle_prefix, "_2.hdf5"))
    src_results = [tmleout1..., tmleout2...]

    for svp_result in svp_results
        src_result_index = findall(x.WTMLE_GLM_GLM.estimand == svp_result.SVP.estimand for x in src_results)
        src_result = src_results[only(src_result_index)]
        @test src_result.WTMLE_GLM_GLM.std != svp_result.SVP.std
        @test src_result.WTMLE_GLM_GLM.estimate == svp_result.SVP.estimate
        @test src_result.WTMLE_GLM_GLM.n == svp_result.SVP.n
        @test svp_result.SVP.IC == []
    end

    close(io)
end

@testset "Test SVP: causal and composed estimands" begin
    # Generate data
    sample_ids = TargeneCore.GRMIDs(joinpath(TESTDIR, "data", "grm", "test.grm.id"))
    tmpdir = mktempdir()
    tmle_prefix = joinpath(tmpdir, "tmle_output")

    configuration = causal_and_joint_estimands_config()
    build_tmle_output_file(
        tmpdir,
        configuration,
        sample_ids, 
    )

    # Using the main command
    output = joinpath(tmpdir, "svp.hdf5")
    main([
        "svp", 
        joinpath(tmpdir, "tmle_output"),
        "--out", output,
        "--grm-prefix", joinpath(TESTDIR, "data", "grm", "test.grm"), 
        "--max-tau", "0.75",
        "--estimator-key", "OSE_GLM_GLM"
    ])

    # The JointEstimate std is not updated but each component is.
    svp_results = jldopen(io -> io["results"], output)
    from_one_dimensional_estimates = svp_results[1:2]
    from_joint_estimates = svp_results[3:4]
    @test from_one_dimensional_estimates[1].SVP.estimand == from_joint_estimates[1].SVP.estimand
    @test from_one_dimensional_estimates[2].SVP.estimand == from_joint_estimates[2].SVP.estimand

    tmleout1 = jldopen(io -> io["Batch_1"], string(tmle_prefix, "_1.hdf5"))
    tmleout2 = jldopen(io -> io["Batch_1"], string(tmle_prefix, "_2.hdf5"))
    src_results = [tmleout1..., tmleout2...]

    # Check std has been updated
    for i in 1:2
        @test from_one_dimensional_estimates[i].SVP.estimand == src_results[i].OSE_GLM_GLM.estimand
        @test from_one_dimensional_estimates[i].SVP.estimate == src_results[i].OSE_GLM_GLM.estimate
        @test from_one_dimensional_estimates[i].SVP.std != src_results[i].OSE_GLM_GLM.std
    end
end

end

true
