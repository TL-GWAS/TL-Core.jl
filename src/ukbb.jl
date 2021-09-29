function parse_queries(queryfile::String)
    config = TOML.parsefile(queryfile)
    queries = Dict()
    for (queryname, querydict) in config
        if lowercase(queryname) ∉ ("threshold", "snps")
            rsids = collect(keys(querydict))
            vals = [split(filter(x->!isspace(x), querydict[rsid]), "->") for rsid in rsids]
            rsids_symbols = Tuple(Symbol(x) for x in rsids)
            queries[queryname] = NamedTuple{rsids_symbols}(vals)
        end
    end
    return queries
end


function read_bgen(bgen_file::String)
    kwargs = Dict{Symbol, Any}(:sample_path => nothing, :idx_path => nothing)
    if bgen_file[end-3:end] == "bgen"
        base = bgen_file[1:end-4]

        samplefile = base * "sample"
        isfile(samplefile) ? kwargs[:sample_path] = samplefile : nothing

        bgifile = bgen_file * ".bgi"
        isfile(bgifile) ? kwargs[:idx_path] = bgifile : nothing
    end
    return Bgen(bgen_file; kwargs...)
end


function samples_genotype(probabilities, variant_genotypes, threshold=0.9)
    n = size(probabilities)[2]
    # The default value is missing
    t = Vector{Union{String, Missing}}(missing, n)
    for i in 1:n
        # If no allele has been annotated with sufficient confidence
        # the sample is declared as missing for this variant
        sample_gen_index = findfirst(x -> x >= threshold, probabilities[:, i])
        sample_gen_index isa Nothing || (t[i] = variant_genotypes[sample_gen_index])
    end
    return t
end


"""

A heterozygous genotype can be specified as (ALLELE₁, ALLELE₂) or (ALLELE₂, ALLELE₁).
Here we align this heterozygous specification on the query and default to 
(ALLELE₁, ALLELE₂) provided in the BGEN file if nothing is specified in the query.
"""
function variant_genotypes(variant::Variant, query::NamedTuple)
    all₁, all₂ = alleles(variant)
    # Either (ALLELE₂, ALLELE₁) is provided in the query
    # and we return it as the heterozygous genotype. 
    # Or the other specification will do in all other cases.
    if all₂*all₁ in query[Symbol(variant.rsid)]
        return [all₁*all₁, all₂*all₁, all₂*all₂]
    end
    return [all₁*all₁, all₁*all₂, all₂*all₂]
end


function UKBBGenotypes(queryfile, query)
    config = TOML.parsefile(queryfile)
    snps = config["SNPS"]
    threshold = config["threshold"]
    # Let's load the variants by the files they are in
    bgen_groups = Dict()
    for (rsid, path) in snps
        haskey(bgen_groups, path) ? push!(bgen_groups[path], rsid) : bgen_groups[path] = [rsid]
    end

    genotypes = nothing
    for (path, rsids) in bgen_groups
        b = GenesInteraction.read_bgen(path)
        chr_genotypes = DataFrame(SAMPLE_ID=b.samples)

        # Iterate over variants in this chromosome
        for rsid in rsids
            v = variant_by_rsid(b, rsid)
            variant_gens = variant_genotypes(v, query)
            probabilities = probabilities!(b, v)
            chr_genotypes[!, rsid] = samples_genotype(probabilities, variant_gens, threshold)
        end
        # I think concatenating should suffice but I still join as a safety
        genotypes isa Nothing ? genotypes = chr_genotypes :
            genotypes = innerjoin(genotypes, chr_genotypes, on=:SAMPLE_ID)
    end
    return genotypes
end


function filter_data(T, W, y)
    # Combine all elements together
    fulldata = hcat(T, W)
    fulldata[!, "Y"] = y

    # Filter based on missingness
    filtered_data = dropmissing(fulldata)

    return filtered_data[!, names(T)], filtered_data[!, names(W)], filtered_data[!, "Y"]
end


function TMLEEpistasisUKBB(phenotypefile, 
    confoundersfile, 
    queryfile, 
    estimatorfile,
    outfile;
    verbosity=1)
    
    # Build tmle
    tmle = tmle_from_toml(TOML.parsefile(estimatorfile))

    # Parse queries
    queries = parse_queries(queryfile)

    # Build Genotypes
    T = UKBBGenotypes(queryfile, queries)

    # Read Confounders
    W = CSV.File(confoundersfile) |> DataFrame

    # Read Target
    y = DataFrame(CSV.File(phenotypefile;header=false))[:, 3]

    # Filter data based on missingness
    T, W, y = filter_data(T, W, y)

    # Run TMLE over potential epistatic SNPS
    mach = machine(tmle, T, W, y)
    fit!(mach, verbosity)
    println(mach.fitresult.estimate)
    println(mach.fitresult.stderror)
end