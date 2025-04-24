##for installing packages
using Pkg

##for paralleling
using Distributed

addprocs(2)

##for sending all the packages and dependencies to workers
@everywhere using Pkg
@everywhere Pkg.activate("~/.julia/environments/v1.11")
@everywhere Pkg.instantiate()
@everywhere using PANDA
@everywhere using JLD2

# List of tree files of amphibians
tree_files = ["niche_position_age/results/data/metadata/extinction_dynamics$(i).tre" for i in [1:2]]


# Load necessary packages on all workers
@everywhere using JLD2, PANDA.ClaDS

# Function to process a single tree
@everywhere function process_tree(file)
    tree = PANDA.ClaDS.load_tree(file)
    result = PANDA.ClaDS.infer_ClaDS(tree, 100)
    
    filename = split(file, "/")[end]  # Extract file name
    output_file = "spe_rates/extinction_dynamics/claDS_results_$filename.jld2"
    
    JLD2.save(output_file, "result", result)  # Save results
    return output_file  # Return saved filename
end

# Step 2: Run in parallel using `pmap` (One tree per worker)
results = pmap(process_tree, tree_files)


##when the CladsOutput are generated you should convert them to .RData format for further evaluation

#calling packages needed
using JLD2, PANDA

##naming the JLD2 (CladsOutputs) files
files = filter(x -> endswith(x, ".jld2"), readdir("spe_rates", join=true))


##naming the output directory
output_dir = "Extinciton_dyn_iberic_RData"
mkpath(output_dir)  # Ensure output directory exists

## loop for loading files, convert them to RData, and save them 
for file in files
    data = load(file)  # Load the JLD2 file
    
    if haskey(data, "result") && isa(data["result"], CladsOutput)
        clads_result = data["result"]  # Extract the correct object
        filename = splitext(basename(file))[1] * ".RData"  # Change extension
        output_path = joinpath(output_dir, filename)
        
        save_ClaDS_in_R(clads_result, output_path)  # Save as RData
    else
        println("Skipping file: ", file, " (CladsOutput not found)")
    end
end