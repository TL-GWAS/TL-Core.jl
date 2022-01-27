using ArgParse
using TMLEEpistasis


function parse_commandline()
    s = ArgParseSettings(
        description = "This program computes estimates of Genetic Epistatis from the UK-Biobank using TMLE."*
                      "Here is a list of the arguments that should be provided, you can also have a look at the "*
                      "test/data and test/config folders to see some examples.",
        commands_are_required = false,
        version = "0.2",
        add_version = true)


    @add_arg_table s begin
        "phenotypes"
            help = "A file (.csv format). The first row contains the column names with `eid` the sample ID"*
                   " and the rest of the columns are phenotypes of interest."
            required = true
        "confounders"
            help = "A file (.csv format) containing the confounding variables values and the sample ids associated"*
                   " with the participants. The first line of the file should contain the columns names and the sample ids "*
                   " column name should be: `SAMPLE_ID`."
            required = true
        "queries"
            help = "A file (.toml format) see: config/sample_query.toml for more information"
            required = true
        "estimator"
            help = "A file (.toml format) describing the tmle estimator to use, see config/sample_estimator.toml"*
                   " for a basic example."
            required = true
        "output"
            help = "A path where the results will be serialized (.bin format). One entry will be saved"*
                    " for each phenotype. see: --savefull."
            required = true
        "--phenotypes-list", "-p"
            help = "A file, one line for each phenotype, containing a restrictions of the phenotypes "*
                   "to consider for the analysis."
            required = false
            arg_type = String
        "--adaptive-cv", "-a"
            help = "Adaptively selects the number of folds used in cross validation and overrides the default used in the estimator file."
            default = true
            required = false
            arg_type = Bool
        "--savefull", "-f"
            help = "To save the full machine for each phenotype. Otherwise only QueryReport(s) are saved."
            default = false
            required = false
            arg_type = Bool
        "--verbosity", "-v"
            help = "Verbosity level"
            arg_type = Int
            default = 1
    end

    return parse_args(s)
end

parsed_args = parse_commandline()
println(parsed_args)
# UKBBVariantRun(parsed_args)


