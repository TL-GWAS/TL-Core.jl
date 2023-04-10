const CHR_REG = r"chr[1-9]+"

"""
    read_data(filepath)

The SAMPLE_ID column should be read as a String.
"""
read_data(filepath) = CSV.read(filepath, DataFrame, types=Dict(:SAMPLE_ID => String))

yaml_out_path(prefix, index) = string(prefix, ".param_$index.yaml")

function write_tmle_inputs(outprefix, final_dataset, param_files)
    # Write final_dataset
    CSV.write(string(outprefix, ".data.csv"), final_dataset)
    # Write param_files
    for (index, param_file) in enumerate(param_files)
        YAML.write_file(yaml_out_path(outprefix, index), param_file)
    end
end


function call_genotypes(probabilities::AbstractArray, variant_genotypes::AbstractVector{T}, threshold::Real) where T
    n = size(probabilities, 2)
    t = Vector{Union{T, Missing}}(missing, n)
    for i in 1:n
        # If no allele has been annotated with sufficient confidence
        # the sample is declared as missing for this variant
        genotype_index = findfirst(x -> x >= threshold, probabilities[:, i])
        genotype_index isa Nothing || (t[i] = variant_genotypes[genotype_index])
    end
    return t
end

function is_numbered_chromosome_file(filename, prefix)
    if occursin(prefix, filename) && endswith(filename, "bgen")
        regexp_match = match(CHR_REG, filename)
        if regexp_match !== nothing
            return true
        end
    end
    return false
end

function read_bgen(filepath)
    sample_filepath = string(filepath[1:end-4], "sample")
    idx_filepath = string(filepath, ".bgi")
    return Bgen(filepath, sample_path=sample_filepath, idx_path=idx_filepath)
end

all_snps_called(found_snps::Set, snp_list) = Set(snp_list) == found_snps

"""
    genotypes_encoding(variant; asint=true)

If asint is true then the number of minor alleles is reported, otherwise string genotypes are reported.
"""
function genotypes_encoding(variant; asint=true)
    minor = minor_allele(variant)
    all₁, all₂ = alleles(variant)
    if asint
        if all₁ == minor
            return [2, 1, 0]
        else
            return [0, 1, 2]
        end
    else
        return [all₁*all₁, all₁*all₂, all₂*all₂]
    end
end

NotAllVariantsFoundError(found_snps, snp_list) = 
    ArgumentError(string("Some variants were not found in the genotype files: ", join(setdiff(snp_list, found_snps), ", ")))

"""
    bgen_files(snps, bgen_prefix)

This function assumes the UK-Biobank structure
"""
function call_genotypes(bgen_prefix::String, snp_list, threshold::Real; asint=true)
    chr_dir_, prefix_ = splitdir(bgen_prefix)
    chr_dir = chr_dir_ == "" ? "." : chr_dir_
    genotypes = nothing
    found_snps = Set()
    for filename in readdir(chr_dir)
        all_snps_called(found_snps, snp_list) ? break : nothing
        if is_numbered_chromosome_file(filename, prefix_)
            bgenfile = read_bgen(joinpath(chr_dir_, filename))
            chr_genotypes = DataFrame(SAMPLE_ID=bgenfile.samples)
            for variant in BGEN.iterator(bgenfile)
                rsid_ = rsid(variant)
                if rsid_ ∈ snp_list
                    push!(found_snps, rsid_)
                    if n_alleles(variant) != 2
                        @warn("Skipping $rsid_, not bi-allelic")
                        continue
                    end
                    minor_allele_dosage!(bgenfile, variant)
                    variant_genotypes = genotypes_encoding(variant; asint=asint)
                    probabilities = probabilities!(bgenfile, variant)
                    chr_genotypes[!, rsid_] = call_genotypes(probabilities, variant_genotypes, threshold)
                end
            end
            genotypes = genotypes isa Nothing ? chr_genotypes :
                    innerjoin(genotypes, chr_genotypes, on=:SAMPLE_ID)
        end
    end
    all_snps_called(found_snps, snp_list) || throw(NotAllVariantsFoundError(found_snps, snp_list))
    return genotypes
end

function satisfies_positivity(interaction_setting, freqs; positivity_constraint=0.01)
    for base_setting in Iterators.product(interaction_setting...)
        if !haskey(freqs, base_setting) || freqs[base_setting] < positivity_constraint
            return false
        end
    end
    return true
end

function frequency_table(treatments, treatment_tuple::AbstractVector)
    freqs = Dict()
    N = nrow(treatments)
    for (key, group) in pairs(groupby(treatments, treatment_tuple; skipmissing=true))
        freqs[values(key)] = nrow(group) / N
    end
    return freqs
end

read_txt_file(path::Nothing) = nothing
read_txt_file(path) = CSV.read(path, DataFrame, header=false)[!, 1]

pcnames(pcs) = filter(!=("SAMPLE_ID"), names(pcs))

all_confounders(pcs, extraW::Nothing) = pcnames(pcs)
all_confounders(pcs, extraW) = vcat(pcnames(pcs), extraW)

targets_from_traits(traits, non_targets) = filter(x -> x ∉ non_targets, names(traits))


function add_batchified_param_files!(new_param_files, param_file, variables, batch_size)
    param_files = batched_param_files(param_file, variables, batch_size)
    append!(new_param_files, param_files)
end

function batched_param_files(param_file, phenotypes_list, batch_size::Nothing)
    param_file = copy(param_file)
    param_file["Y"] = phenotypes_list
    return [param_file]
end

function batched_param_files(param_file, phenotypes_list, batch_size::Int)
    new_param_files = []
    for batch in Iterators.partition(phenotypes_list, batch_size)
        batched_param_file = copy(param_file)
        batched_param_file["Y"] = batch
        push!(new_param_files, batched_param_file)
    end
    return new_param_files
end

function merge(traits, pcs, genotypes)
    return innerjoin(
        innerjoin(traits, pcs, on="SAMPLE_ID"),
        genotypes,
        on="SAMPLE_ID"
    )
end


"""
    tmle_inputs(parsed_args)

Support for the generation of parameters according to 2 strategies:
- from-actors
- from-param-files
"""
function tmle_inputs(parsed_args)
    if parsed_args["%COMMAND%"] == "from-actors"
        tmle_inputs_from_actors(parsed_args)
    elseif parsed_args["%COMMAND%"] == "from-param-files"
        tmle_inputs_from_param_files(parsed_args)
    else
        throw(ArgumentError("Unrecognized command."))
    end
end


