# --- 0. 安装和加载必要的包 ---
# (安装命令保持不变，确保已执行)

library(miaSim)
library(SummarizedExperiment) # assay() 函数在这里
library(gtools)
library(parallel) # 用于并行计算

cat("开始执行微生物组丰度数据生成脚本 (并行版本)...\n\n")

# --- 全局参数定义 ---
n_global_species <- 1000
n_global_resources <- 10 # 你在脚本中设置的是10
n_simulations <- 1000000    # 你在脚本中设置的是1000

global_species_names <- paste0("GlobalSpecies", sprintf("%04d", 1:n_global_species))
global_resource_names <- paste0("Resource", sprintf("%02d", 1:n_global_resources))

set.seed(123)

# --- 1. 定义全局物种资源交互矩阵 (E_global) ---
cat("步骤 1: 使用 randomE 生成全局物种-资源交互矩阵 (E_global)...\n")
E_global <- miaSim::randomE(
    n_species = n_global_species,
    n_resources = n_global_resources,
    names_species = global_species_names,
    names_resources = global_resource_names,
    mean_production = 1
)
cat("E_global 维度:", dim(E_global), "\n")
E_global_filename <- "E_global_v46.csv"

cat("步骤 1.5: 生成全局物种-资源 Monod 常数矩阵 (monod_constant_global)...\n")
# 假设资源量的最大值为 100 来设定 shape 参数，这与函数默认行为一致
max_resource_assumed <- 100 
monod_constant_global <- matrix(rgamma(n = n_global_species * n_global_resources,
                                       shape = 50 * max_resource_assumed,
                                       rate = 1),
                                nrow = n_global_species,
                                # 确保矩阵有行列名称，便于后续按名称索引
                                dimnames = list(global_species_names, global_resource_names))
cat("monod_constant_global 维度:", dim(monod_constant_global), "\n\n")

# --- 2. 为每个微生物定义 prevalence ---
# (这部分代码与你提供的脚本一致)
cat("步骤 2: 为每个微生物定义 prevalence...\n")
prevalences <- rbeta(n_global_species, 1, 10) 
names(prevalences) <- global_species_names
cat("Prevalence 分布摘要:\n")
print(summary(prevalences))
cat("Prevalence 方差:", var(prevalences), "\n\n")


# --- 3. 执行迭代模拟并记录丰度数据 (并行版本) ---
# [修改代码] 增加 monod_constant_global_matrix 参数
run_single_simulation <- function(sim_idx, E_global_matrix, monod_constant_global_matrix, global_species_names_vec, prevalences_vec, n_global_resources_val, global_resource_names_vec, min_species_val, simulation_t_end_val) {
    # 3.1 根据 prevalence 采样当前样本的物种组成
    present_species_logical <- runif(length(global_species_names_vec)) < prevalences_vec
    subset_species_names <- global_species_names_vec[present_species_logical]
    n_subset_species <- length(subset_species_names)

    current_simulation_abundances <- rep(NA, length(global_species_names_vec))
    names(current_simulation_abundances) <- global_species_names_vec
    
    simulation_status <- "skipped_too_few_species" 

    if (n_subset_species < min_species_val) {
        return(list(abundances = current_simulation_abundances, status = simulation_status, sim_idx = sim_idx, n_subset = n_subset_species))
    }

    # 从全局矩阵中提取当前样本的子集
    E_current_sample <- E_global_matrix[subset_species_names, , drop = FALSE]
    # [新增代码] 从全局 Monod 矩阵中提取当前样本的子集
    monod_constant_current_sample <- monod_constant_global_matrix[subset_species_names, , drop = FALSE]
    
    current_x0 <- rep(10, n_subset_species)
    names(current_x0) <- subset_species_names

    # 注意：这里的 resources 为 NULL，函数会默认生成随机资源。
    # 如果希望所有样本的初始资源也相同，需要在这里手动设定一个固定的向量。
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
            # [修改代码] 传入固定的 Monod 常数子集
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

# 设置并行计算集群
cat("步骤 3: 开始 ", n_simulations, " 次迭代模拟 (使用并行计算)...\n")
min_species_for_simulation <- 10

num_cores <- detectCores() - 5
if (num_cores < 1) num_cores <- 1 
cat("检测到CPU核心数 (逻辑):", detectCores(), "将使用:", num_cores, "核心进行并行计算。\n")

cl <- makeCluster(num_cores)

# [修改代码] 将新的全局 Monod 矩阵导出到工作进程
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

cat("开始并行执行模拟任务...\n")
start_time <- Sys.time()

# [修改代码] 在调用时传递新的全局 Monod 矩阵
results_list <- parLapply(cl, 1:n_simulations, function(idx) {
    run_single_simulation(
        sim_idx = idx,
        E_global_matrix = E_global, 
        monod_constant_global_matrix = monod_constant_global, # <-- 传递新矩阵
        global_species_names_vec = global_species_names,
        prevalences_vec = prevalences,
        n_global_resources_val = n_global_resources,
        global_resource_names_vec = global_resource_names,
        min_species_val = min_species_for_simulation,
        simulation_t_end_val = 1000 # t_end
    )
})

end_time <- Sys.time()
cat("并行模拟任务完成。耗时:", format(end_time - start_time), "\n")

stopCluster(cl)
RNGkind(original_rngkind[1], original_rngkind[2], original_rngkind[3]) 

# --- 处理并行运算的结果 ---
# (这部分代码与你脚本中一致，无需改动)
cat("\n正在处理并行模拟结果...\n")
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
        cat("已处理结果:", i, "/", n_simulations, 
            " (成功:", successful_sim_count, 
            "失败:", failed_sim_count, 
            "跳过:", skipped_count, ")\n")
    }
}

cat("\n所有模拟数据处理完毕。\n")
cat("最终成功模拟次数:", successful_sim_count, "\n")
cat("最终失败模拟次数 (因错误或数据提取问题):", failed_sim_count, "\n")
cat("最终跳过模拟次数 (因物种太少):", skipped_count, "\n")

if(!file.exists(E_global_filename)){
    write.csv(E_global, E_global_filename, row.names = TRUE)
    cat("全局交互矩阵 E_global 已保存到", E_global_filename, "\n")
}

abundance_table_filename <- "final_abundance_table_v46.csv"
write.csv(final_abundance_table, abundance_table_filename, row.names = TRUE, na = "NA")
cat("最终丰度表已保存到", abundance_table_filename, "\n")

cat("\n脚本执行完毕。\n")