##########ClaDS for the plant phylogeny (collapsed to genus level) ############

###########################################################
###### we ran this Julia code in USP server,  #############
######        which operates under Linux       ############
##########################################################

#open an independent julia session
tmux new -s julia_session

#open Julia
julia
using Pkg
Pkg.activate("~/.julia/environments/v1.11")
Pkg.instantiate()

using PANDA
using JLD2
using PANDA.ClaDS

# Single tree file
tree_file = "plant_genus_phylo.tre"

# Function to process the tree
function process_tree(file)
    tree = PANDA.ClaDS.load_tree(file)
    result = PANDA.ClaDS.infer_ClaDS(tree, 100)
    
    filename = split(file, "/")[end]  # Extract file name
    output_file = "claDS_results_$filename.jld2"
    
    JLD2.save(output_file, "result", result)
    
    return output_file
end

# Run for a single tree
result_file = process_tree(tree_file)

####### after running ClaDs ########

#transforming jld2 file into RData

# Input file
file = "claDS_results_plant_genus_phylo.tre.jld2"

# Load the JLD2 file
data = load(file)

#assigning only the result part of the data 
clads_result = data["result"]

# Save directly in the workspace
save_ClaDS_in_R(clads_result, "clads_plant_phylo.RData")





