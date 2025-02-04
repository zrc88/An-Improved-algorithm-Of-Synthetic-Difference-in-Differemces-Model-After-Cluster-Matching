#' 使用 future.apply 进行并行计算的 synthdid 聚类和 ATT 估计
#'
#' @param Y Outcome 矩阵
#' @param N0 控制组个数（假设 Y 的前 N0 行为控制组）
#' @param T0 治疗前时间期数
#' @return 一个列表，包含两个部分：
#'         - cluster_results: 每个簇的 tau_hat 与 se 的列表
#'         - total_ATT: 加权后的整体 ATT
#'         - total_se: 加权后的整体标准误
#' @export
synthdid_cluster <- function(Y, N0, T0) {
  # 检查并加载 future.apply 包
  if (!requireNamespace("future.apply", quietly = TRUE)) {
    stop("需要安装 'future.apply' 包，请先运行 install.packages('future.apply')")
  }
  library(future.apply)

  # 设置并行计划。Windows 系统推荐使用 multisession。
  future::plan(multisession)

  # 1. 通过 elbow 方法确定聚类个数。
  # 计算不同聚类个数下的总簇内平方误差（WSS）
  wss <- sapply(1:10, function(k) {
    kmeans(Y, centers = k, nstart = 20)$tot.withinss
  })

  # 计算一阶差分（变化率）与二阶差分（加速度）
  wss_diff  <- diff(wss)
  wss_diff2 <- diff(wss_diff)

  # 找到下降幅度最大的点，至少保证 2 个聚类
  optimal_k <- max(2, which.max(abs(wss_diff2)) + 1)

  # 2. 使用确定的 optimal_k 进行 k-means 聚类
  kmeans_result <- kmeans(Y, centers = optimal_k, nstart = 20)

  # 按照聚类结果分组行索引
  cluster_groups <- split(1:nrow(Y), kmeans_result$cluster)

  # 3. 使用 future_lapply 并行计算每个聚类内的 synthdid 估计
  synthdid_results <- future.apply::future_lapply(names(cluster_groups), function(cluster_id) {
    # 获取该簇的行索引和对应的子矩阵
    selected_rows <- cluster_groups[[cluster_id]]
    cluster_matrix <- Y[selected_rows, , drop = FALSE]

    # 假定前 N0 行为控制组，因此在该簇中控制组个数为：
    control_count <- sum(selected_rows <= N0)

    # 计算 synthdid 估计与标准误
    tau.hat <- synthdid_estimate(Y = cluster_matrix, N0 = control_count, T0 = T0)
    se      <- synthdid_se(tau.hat, method = "bootstrap", replications = 200)

    list(tau_hat = as.numeric(tau.hat), se = as.numeric(se))
  })

  # 为结果命名（名称与聚类 id 对应）
  names(synthdid_results) <- names(cluster_groups)

  # 4. 根据每个聚类的单位数量计算权重
  total_units <- nrow(Y)
  group_sizes <- sapply(cluster_groups, length)
  weights <- group_sizes / total_units

  # 5. 计算整体的 ATT 和标准误（加权汇总各聚类结果）
  total_ATT <- sum(sapply(seq_along(synthdid_results), function(i) {
    synthdid_results[[i]]$tau_hat * weights[i]
  }))

  total_se <- sum(sapply(seq_along(synthdid_results), function(i) {
    synthdid_results[[i]]$se * weights[i]^2
  }))

  # 返回包含每个聚类结果和整体 ATT 与 se 的列表
  return(list(
    cluster_results = synthdid_results,
    total_ATT = total_ATT,
    total_se = total_se
  ))
}
