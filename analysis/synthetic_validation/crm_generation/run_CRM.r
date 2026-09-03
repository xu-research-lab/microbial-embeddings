# --- 0. Load packages ---

library(miaSim)
library(SummarizedExperiment) # Provides assay().
library(gtools)
library(parallel)

cat("Starting microbial abundance simulations...\n\n")

# --- Global parameters ---
n_global_species <- 1000
n_global_resources <- 10
n_simulations <- 1000000

global_species_names <- paste0("GlobalSpecies", sprintf("%04d", 1:n_global_species))
global_resource_names <- paste0("Resource", sprintf("%02d", 1:n_global_resources))

set.seed(123)

# --- 1. Define the global species-resource interaction matrix ---
cat("Step 1: Generating the global species-resource interaction matrix (E_global)...\n")
E_global <- miaSim::randomE(
    n_species = n_global_species,
    n_resources = n_global_resources,
    names_species = global_species_names,
    names_resources = global_resource_names,
    mean_production = 1
)
cat("E_global dimensions:", dim(E_global), "\n")
E_global_filename <- "data/E_global.csv"

cat("Step 1.5: Generating the global species-resource Monod constant matrix (monod_constant_global)...\n")
# Match the function default by assuming a maximum resource amount of 100.
max_resource_assumed <- 100 
monod_constant_global <- matrix(rgamma(n = n_global_species * n_global_resources,
                                       shape = 50 * max_resource_assumed,
                                       rate = 1),
                                nrow = n_global_species,
                                # Preserve names for later indexing.
                                dimnames = list(global_species_names, global_resource_names))
cat("monod_constant_global dimensions:", dim(monod_constant_global), "\n\n")

# --- 2. Define prevalence for each microbe ---
cat("Step 2: Defining prevalence for each microbe...\n")
prevalences <- rbeta(n_global_species, 1, 10) 
names(prevalences) <- global_species_names
cat("Prevalence summary:\n")
print(summary(prevalences))
cat("Prevalence variance:", var(prevalences), "\n\n")


# --- 3. Run simulations in parallel and record abundances ---
run_single_simulation <- function(sim_idx, E_global_matrix, monod_constant_global_matrix, global_species_names_vec, prevalences_vec, n_global_resources_val, global_resource_names_vec, min_species_val, simulation_t_end_val) {
    # Sample the species present in this community from their prevalences.
    present_species_logical <- runif(length(global_species_names_vec)) < prevalences_vec
    subset_species_names <- global_species_names_vec[present_species_logical]
    n_subset_species <- length(subset_species_names)

    current_simulation_abundances <- rep(NA, length(global_species_names_vec))
    names(current_simulation_abundances) <- global_species_names_vec
    
    simulation_status <- "skipped_too_few_species" 

    if (n_subset_species < min_species_val) {
        return(list(abundances = current_simulation_abundances, status = simulation_status, sim_idx = sim_idx, n_subset = n_subset_species))
    }

    # Subset the global interaction and Monod constant matrices.
    E_current_sample <- E_global_matrix[subset_species_names, , drop = FALSE]
    monod_constant_current_sample <- monod_constant_global_matrix[subset_species_names, , drop = FALSE]
    
    current_x0 <- rep(10, n_subset_species)
    names(current_x0) <- subset_species_names

    # NULL lets simulateConsumerResource generate random initial resources.
    current_resources <- NULL 

    tse_result <- NULL
    simulation_successful_flag <- FALSE
    error_message <- ""

    tryCatch({
        tse_result <- miaSim::simulateConsumerResource(
            n_species = n_subset_species,
            n_resources = n_global_resources_val,
            names_species = subset_species_names,
            names_resources = global_resource_names_vec,
            E = E_current_sample,
            x0 = current_x0,
            resources = current_resources,
            monod_constant = monod_constant_current_sample,
            t_end = simulation_t_end_val,
            t_step = 1,
            growth_rates = 10,
            inflow_rate = 10,
            outflow_rate = 10,
            stochastic = FALSE,
            error_variance = 0,
        )
        simulation_successful_flag <- TRUE
        simulation_status <- "success"
    }, error = function(e) {
        error_message <- conditionMessage(e)
        simulation_status <- paste("error:", error_message)
    })

    if (simulation_successful_flag && !is.null(tse_result)) {
        if ("counts" %in% SummarizedExperiment::assayNames(tse_result)) {
            abundance_timeseries <- SummarizedExperiment::assay(tse_result, "counts")
            if (ncol(abundance_timeseries) > 0) {
                abundance_at_t_end <- abundance_timeseries[, 501, drop = FALSE]
                common_names_in_result <- intersect(subset_species_names, rownames(abundance_at_t_end))
                current_simulation_abundances[common_names_in_result] <- abundance_at_t_end[common_names_in_result, 1]
            } else {
                simulation_status <- "no_timeseries_data"
            }
        } else {
            simulation_status <- "no_counts_assay"
        }
    }
    
    return(list(abundances = current_simulation_abundances, status = simulation_status, sim_idx = sim_idx, n_subset = n_subset_species))
}

# Set up the parallel cluster.
cat("Step 3: Starting ", n_simulations, " simulations in parallel...\n")
min_species_for_simulation <- 10

num_cores <- detectCores() - 5
if (num_cores < 1) num_cores <- 1 
cat("Detected logical CPU cores:", detectCores(), "using:", num_cores, "cores for parallel processing.\n")

cl <- makeCluster(num_cores)

# Export shared inputs to worker processes.
clusterExport(cl, varlist = c("E_global", "monod_constant_global", "global_species_names", "prevalences", 
                                "n_global_resources", "global_resource_names", 
                                "min_species_for_simulation", "run_single_simulation",
                                "n_global_species")) 

clusterEvalQ(cl, {
    library(miaSim)
    library(SummarizedExperiment)
})

original_rngkind <- RNGkind() 
RNGkind("L'Ecuyer-CMRG")
clusterSetRNGStream(cl, iseed = 4567) 

cat("Starting parallel simulations...\n")
start_time <- Sys.time()

results_list <- parLapply(cl, 1:n_simulations, function(idx) {
    run_single_simulation(
        sim_idx = idx,
        E_global_matrix = E_global, 
        monod_constant_global_matrix = monod_constant_global,
        global_species_names_vec = global_species_names,
        prevalences_vec = prevalences,
        n_global_resources_val = n_global_resources,
        global_resource_names_vec = global_resource_names,
        min_species_val = min_species_for_simulation,
        simulation_t_end_val = 1000 # t_end
    )
})

end_time <- Sys.time()
cat("Parallel simulations completed. Elapsed time:", format(end_time - start_time), "\n")

stopCluster(cl)
RNGkind(original_rngkind[1], original_rngkind[2], original_rngkind[3]) 

# --- Process parallel results ---
cat("\nProcessing parallel simulation results...\n")
final_abundance_table <- matrix(NA, nrow = n_simulations, ncol = n_global_species,
                                  dimnames = list(paste0("Sample", 1:n_simulations), global_species_names))

successful_sim_count <- 0
failed_sim_count <- 0
skipped_count <- 0

for (i in 1:length(results_list)) {
    result_item <- results_list[[i]]
     if (is.null(result_item) || !is.list(result_item) || is.null(result_item$status)) { 
        failed_sim_count <- failed_sim_count + 1
        next
    }
    if (result_item$status == "success") {
        if(length(result_item$abundances) == n_global_species) {
            final_abundance_table[result_item$sim_idx, ] <- result_item$abundances
             successful_sim_count <- successful_sim_count + 1
        } else {
            failed_sim_count <- failed_sim_count + 1
        }
    } else if (startsWith(result_item$status, "skipped")) {
        skipped_count <- skipped_count + 1
    } else { 
        failed_sim_count <- failed_sim_count + 1
    }
    if (i %% (max(1, n_simulations/10)) == 0 || i == n_simulations || i == 1 ) { 
        cat("Processed results:", i, "/", n_simulations,
            " (successful:", successful_sim_count,
            "failed:", failed_sim_count,
            "skipped:", skipped_count, ")\n")
    }
}

cat("\nAll simulation results have been processed.\n")
cat("Successful simulations:", successful_sim_count, "\n")
cat("Failed simulations (errors or data extraction issues):", failed_sim_count, "\n")
cat("Skipped simulations (too few species):", skipped_count, "\n")

if(!file.exists(E_global_filename)){
    write.csv(E_global, E_global_filename, row.names = TRUE)
    cat("Global interaction matrix E_global saved to", E_global_filename, "\n")
}

abundance_table_filename <- "data/final_abundance_table.csv"
write.csv(final_abundance_table, abundance_table_filename, row.names = TRUE, na = "NA")
cat("Final abundance table saved to", abundance_table_filename, "\n")

cat("\nScript completed.\n")
