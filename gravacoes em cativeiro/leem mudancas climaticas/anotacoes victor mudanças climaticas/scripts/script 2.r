###############################################################################
# SCRIPT 02 — CLASSIFICACAO NAO SUPERVISIONADA RF + PAM
###############################################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(123)

NTREE_FINAL <- 2000
NTREE_BOOT <- 500
N_BOOT <- 30
K_MAX <- 8
PROPORCAO_BOOT <- 0.80
SEED_BASE <- 123

pacotes <- c(
  "randomForest", "cluster", "mclust",
  "dplyr", "tidyr", "ggplot2", "scales", "patchwork"
)
instalar <- pacotes[!pacotes %in% rownames(installed.packages())]
if (length(instalar) > 0) install.packages(instalar, dependencies = TRUE)

suppressPackageStartupMessages({
  library(randomForest)
  library(cluster)
  library(mclust)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

base <- "C:/Users/abiel/Downloads/Patrick mudancas climaticas"
raiz <- file.path(base, "Pipeline_V2_Podocnemis")

pasta_pre <- file.path(
  raiz,
  "10_Pipeline_RF_PAM_FINAL",
  "01_Pre_processamento"
)
pasta_saida <- file.path(
  raiz,
  "10_Pipeline_RF_PAM_FINAL",
  "02_Classificacao_RF_PAM"
)
pasta_figuras <- file.path(pasta_saida, "Figuras_relatorio")
pasta_diagnosticos <- file.path(pasta_saida, "Diagnosticos")

dir.create(pasta_figuras, recursive = TRUE, showWarnings = FALSE)
dir.create(pasta_diagnosticos, recursive = TRUE, showWarnings = FALSE)

arquivos_necessarios <- c(
  file.path(pasta_pre, "dados_analise_final.rds"),
  file.path(pasta_pre, "matriz_descritores_padronizada.rds"),
  file.path(pasta_pre, "descritores_finais.rds"),
  file.path(pasta_pre, "PCA_exploratoria.rds")
)

if (!all(file.exists(arquivos_necessarios))) {
  stop("Arquivos do Script 01 nao encontrados.", call. = FALSE)
}

dados <- readRDS(arquivos_necessarios[1])
matriz_padronizada <- as.matrix(readRDS(arquivos_necessarios[2]))
descritores_finais <- readRDS(arquivos_necessarios[3])
pca_exploratoria <- readRDS(arquivos_necessarios[4])

ajustar_urf <- function(matriz, ntree, seed) {
  set.seed(seed)
  randomForest::randomForest(
    x = matriz,
    y = NULL,
    ntree = ntree,
    proximity = TRUE,
    oob.prox = TRUE,
    keep.forest = TRUE
  )
}

proximidade_para_dissimilaridade <- function(proximidade) {
  proximidade <- as.matrix(proximidade)
  proximidade[!is.finite(proximidade)] <- 0
  diag(proximidade) <- 1
  
  dissimilaridade <- 1 - proximidade
  dissimilaridade[dissimilaridade < 0] <- 0
  dissimilaridade[dissimilaridade > 1] <- 1
  diag(dissimilaridade) <- 0
  
  stats::as.dist(dissimilaridade)
}

avaliar_k_pam <- function(dissimilaridade, k_max) {
  n <- attr(dissimilaridade, "Size")
  k_limite <- min(k_max, floor(sqrt(n)), n - 1)
  
  dplyr::bind_rows(
    lapply(
      2:k_limite,
      function(k) {
        ajuste <- cluster::pam(
          dissimilaridade,
          k = k,
          diss = TRUE
        )
        data.frame(
          k = k,
          silhouette_media = ajuste$silinfo$avg.width
        )
      }
    )
  )
}

selecionar_numero_pcs <- function(pca, proporcao = 0.80) {
  acumulada <- cumsum(pca$sdev^2 / sum(pca$sdev^2))
  indice <- which(acumulada >= proporcao)[1]
  max(2, min(indice, 20, ncol(pca$x)))
}

message("Ajustando Random Forest nao supervisionado...")

urf_final <- ajustar_urf(
  matriz = matriz_padronizada,
  ntree = NTREE_FINAL,
  seed = SEED_BASE
)

proximidade_rf <- urf_final$proximity
dissimilaridade_rf <- proximidade_para_dissimilaridade(proximidade_rf)

avaliacao_k_rf <- avaliar_k_pam(
  dissimilaridade = dissimilaridade_rf,
  k_max = K_MAX
)

melhor_k <- avaliacao_k_rf$k[
  which.max(avaliacao_k_rf$silhouette_media)
]

pam_rf_final <- cluster::pam(
  dissimilaridade_rf,
  k = melhor_k,
  diss = TRUE
)

grupo_rf <- factor(pam_rf_final$clustering)

n_pcs <- selecionar_numero_pcs(pca_exploratoria, proporcao = 0.80)
scores_pca_80 <- pca_exploratoria$x[, seq_len(n_pcs), drop = FALSE]
dissimilaridade_pca <- stats::dist(scores_pca_80)

avaliacao_k_pca <- avaliar_k_pam(
  dissimilaridade = dissimilaridade_pca,
  k_max = K_MAX
)

melhor_k_pca <- avaliacao_k_pca$k[
  which.max(avaliacao_k_pca$silhouette_media)
]

pam_pca_final <- cluster::pam(
  dissimilaridade_pca,
  k = melhor_k_pca,
  diss = TRUE
)

grupo_pca <- factor(pam_pca_final$clustering)

ari_rf_pca <- mclust::adjustedRandIndex(grupo_rf, grupo_pca)

message("Executando reamostragens de estabilidade...")

estabilidade <- dplyr::bind_rows(
  lapply(
    seq_len(N_BOOT),
    function(repeticao) {
      set.seed(SEED_BASE + 1000 + repeticao)
      
      indices <- sort(
        sample(
          seq_len(nrow(matriz_padronizada)),
          size = floor(PROPORCAO_BOOT * nrow(matriz_padronizada)),
          replace = FALSE
        )
      )
      
      matriz_sub <- matriz_padronizada[indices, , drop = FALSE]
      
      urf_sub <- ajustar_urf(
        matriz = matriz_sub,
        ntree = NTREE_BOOT,
        seed = SEED_BASE + 2000 + repeticao
      )
      
      diss_sub <- proximidade_para_dissimilaridade(urf_sub$proximity)
      
      pam_sub <- cluster::pam(
        diss_sub,
        k = melhor_k,
        diss = TRUE
      )
      
      data.frame(
        repeticao = repeticao,
        ARI = mclust::adjustedRandIndex(
          grupo_rf[indices],
          pam_sub$clustering
        )
      )
    }
  )
)

mds_rf <- stats::cmdscale(
  dissimilaridade_rf,
  k = 2,
  eig = TRUE,
  add = TRUE
)

scores_rf <- as.data.frame(mds_rf$points)
names(scores_rf) <- c("Eixo1", "Eixo2")

dados_cluster <- dados %>%
  dplyr::select(
    id_vocalizacao, tipo, sala, ambiente, data, unidade
  ) %>%
  dplyr::mutate(
    grupo_computacional = grupo_rf
  )

scores_rf <- scores_rf %>%
  dplyr::bind_cols(dados_cluster)

matriz_descritores <- as.data.frame(matriz_padronizada)
matriz_descritores$grupo_computacional <- grupo_rf

efeito_descritores <- dplyr::bind_rows(
  lapply(
    descritores_finais,
    function(descritor) {
      x <- matriz_descritores[[descritor]]
      grupo <- matriz_descritores$grupo_computacional
      teste <- stats::kruskal.test(x ~ grupo)
      
      n <- length(x)
      k <- nlevels(grupo)
      epsilon2 <- max(
        0,
        (unname(teste$statistic) - k + 1) / (n - k)
      )
      
      medianas <- tapply(x, grupo, stats::median, na.rm = TRUE)
      
      data.frame(
        descritor = descritor,
        H = unname(teste$statistic),
        p = teste$p.value,
        epsilon2 = epsilon2,
        mediana_grupo_1 = unname(medianas[1]),
        mediana_grupo_2 = if (length(medianas) >= 2) {
          unname(medianas[2])
        } else {
          NA_real_
        }
      )
    }
  )
) %>%
  dplyr::mutate(
    p_ajustado_BH = stats::p.adjust(p, method = "BH")
  ) %>%
  dplyr::arrange(dplyr::desc(epsilon2))

top_descritores <- efeito_descritores %>%
  dplyr::slice_head(n = min(12, nrow(efeito_descritores)))

cores_grupos <- c(
  "1" = "#D73027", "2" = "#4575B4", "3" = "#1A9850",
  "4" = "#984EA3", "5" = "#FF7F00", "6" = "#A65628",
  "7" = "#F781BF", "8" = "#666666"
)
cores_tipo <- c(
  tipo1 = "#D73027",
  tipo2 = "#4575B4",
  tipo3 = "#1A9850"
)
rotulos_tipo <- c(
  tipo1 = "Tipo 1",
  tipo2 = "Tipo 2",
  tipo3 = "Tipo 3"
)

tema <- ggplot2::theme_bw(base_size = 14) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 17),
    plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 12),
    axis.title = ggplot2::element_text(face = "bold"),
    strip.text = ggplot2::element_text(face = "bold"),
    panel.grid.minor = ggplot2::element_blank()
  )
ggplot2::theme_set(tema)

avaliacao_k_comparada <- dplyr::bind_rows(
  avaliacao_k_rf %>%
    dplyr::mutate(metodo = "RF nao supervisionado + PAM"),
  avaliacao_k_pca %>%
    dplyr::mutate(metodo = "PCA 80% + PAM")
)

figura_k <- ggplot2::ggplot(
  avaliacao_k_comparada,
  ggplot2::aes(
    x = k,
    y = silhouette_media,
    color = metodo,
    group = metodo
  )
) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_vline(xintercept = melhor_k, linetype = 2) +
  ggplot2::scale_x_continuous(
    breaks = sort(unique(avaliacao_k_comparada$k))
  ) +
  ggplot2::labs(
    title = "Escolha do numero de agrupamentos acusticos",
    subtitle = paste0(
      "Modelo principal: RF nao supervisionado + PAM; k = ",
      melhor_k
    ),
    x = "Numero de grupos (k)",
    y = "Silhouette media",
    color = "Metodo"
  )

ggplot2::ggsave(
  file.path(pasta_figuras, "Figura_03_escolha_k_RF_PAM.png"),
  figura_k,
  width = 11,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

figura_grupos <- ggplot2::ggplot(
  scores_rf,
  ggplot2::aes(
    x = Eixo1,
    y = Eixo2,
    color = grupo_computacional
  )
) +
  ggplot2::geom_point(alpha = 0.55, size = 1.9) +
  ggplot2::scale_color_manual(
    values = cores_grupos[
      levels(scores_rf$grupo_computacional)
    ]
  ) +
  ggplot2::labs(
    title = "Agrupamentos acusticos identificados pelo modelo nao supervisionado",
    subtitle = "Ordenacao da dissimilaridade derivada das proximidades do Random Forest",
    x = "Eixo 1 da ordenacao",
    y = "Eixo 2 da ordenacao",
    color = "Grupo acustico"
  )

ggplot2::ggsave(
  file.path(pasta_figuras, "Figura_04_ordenacao_RF_grupos.png"),
  figura_grupos,
  width = 11,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)

painel_manual <- ggplot2::ggplot(
  scores_rf,
  ggplot2::aes(x = Eixo1, y = Eixo2, color = tipo)
) +
  ggplot2::geom_point(alpha = 0.52, size = 1.8) +
  ggplot2::scale_color_manual(values = cores_tipo, labels = rotulos_tipo) +
  ggplot2::labs(
    title = "Tipos classificados manualmente",
    x = "Eixo 1",
    y = "Eixo 2",
    color = "Tipo manual"
  )

painel_computacional <- ggplot2::ggplot(
  scores_rf,
  ggplot2::aes(x = Eixo1, y = Eixo2, color = grupo_computacional)
) +
  ggplot2::geom_point(alpha = 0.52, size = 1.8) +
  ggplot2::scale_color_manual(
    values = cores_grupos[
      levels(scores_rf$grupo_computacional)
    ]
  ) +
  ggplot2::labs(
    title = "Grupos identificados computacionalmente",
    x = "Eixo 1",
    y = "Eixo 2",
    color = "Grupo acustico"
  )

figura_lado_a_lado <- (
  painel_manual | painel_computacional
) +
  patchwork::plot_annotation(
    title = "Comparacao entre classificacao manual e agrupamento nao supervisionado",
    subtitle = "As mesmas vocalizacoes sao exibidas no mesmo espaco de ordenacao"
  )

ggplot2::ggsave(
  file.path(
    pasta_figuras,
    "Figura_05_manual_computacional_lado_a_lado.png"
  ),
  figura_lado_a_lado,
  width = 16,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

figura_estabilidade <- ggplot2::ggplot(
  estabilidade,
  ggplot2::aes(x = ARI)
) +
  ggplot2::geom_histogram(
    bins = 15,
    color = "white",
    fill = "#555555"
  ) +
  ggplot2::geom_vline(
    xintercept = stats::median(estabilidade$ARI, na.rm = TRUE),
    linetype = 2
  ) +
  ggplot2::labs(
    title = "Estabilidade dos agrupamentos por reamostragem",
    subtitle = paste0(
      "Mediana do ARI = ",
      round(stats::median(estabilidade$ARI, na.rm = TRUE), 3)
    ),
    x = "Indice de Rand ajustado (ARI)",
    y = "Numero de reamostragens"
  )

ggplot2::ggsave(
  file.path(pasta_diagnosticos, "Diagnostico_estabilidade_clusters.png"),
  figura_estabilidade,
  width = 9,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

figura_descritores <- ggplot2::ggplot(
  top_descritores %>%
    dplyr::mutate(
      descritor = factor(descritor, levels = rev(descritor))
    ),
  ggplot2::aes(x = epsilon2, y = descritor)
) +
  ggplot2::geom_col(fill = "#4575B4") +
  ggplot2::labs(
    title = "Descritores que mais diferenciam os grupos acusticos",
    subtitle = "Tamanho de efeito do teste de Kruskal-Wallis",
    x = "Tamanho de efeito (epsilon quadrado)",
    y = "Descritor acustico"
  )

ggplot2::ggsave(
  file.path(pasta_figuras, "Figura_06_descritores_discriminantes.png"),
  figura_descritores,
  width = 10,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)

medias_grupo <- matriz_descritores %>%
  dplyr::group_by(grupo_computacional) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(top_descritores$descritor),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(
    cols = -grupo_computacional,
    names_to = "descritor",
    values_to = "media_padronizada"
  )

figura_heatmap <- ggplot2::ggplot(
  medias_grupo,
  ggplot2::aes(
    x = grupo_computacional,
    y = descritor,
    fill = media_padronizada
  )
) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.2f", media_padronizada)),
    size = 3.5
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#4575B4",
    mid = "white",
    high = "#D73027",
    midpoint = 0
  ) +
  ggplot2::labs(
    title = "Perfil acustico medio dos grupos computacionais",
    subtitle = "Valores padronizados dos descritores com maior tamanho de efeito",
    x = "Grupo acustico",
    y = "Descritor",
    fill = "Media\npadronizada"
  )

ggplot2::ggsave(
  file.path(pasta_figuras, "Figura_07_heatmap_perfil_acustico_grupos.png"),
  figura_heatmap,
  width = 10,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)

utils::write.csv(
  dados_cluster,
  file.path(pasta_saida, "Grupos_computacionais_por_vocalizacao.csv"),
  row.names = FALSE
)
utils::write.csv(
  avaliacao_k_comparada,
  file.path(pasta_saida, "Silhouette_comparacao_metodos.csv"),
  row.names = FALSE
)
utils::write.csv(
  estabilidade,
  file.path(pasta_saida, "Estabilidade_reamostragem.csv"),
  row.names = FALSE
)
utils::write.csv(
  efeito_descritores,
  file.path(pasta_saida, "Descritores_diferencas_entre_grupos.csv"),
  row.names = FALSE
)
utils::write.csv(
  scores_rf,
  file.path(pasta_saida, "Scores_ordenacao_RF.csv"),
  row.names = FALSE
)

saveRDS(urf_final, file.path(pasta_saida, "Modelo_URF_final.rds"))
saveRDS(pam_rf_final, file.path(pasta_saida, "Modelo_PAM_RF_final.rds"))
saveRDS(dissimilaridade_rf, file.path(pasta_saida, "Dissimilaridade_RF.rds"))

writeLines(
  c(
    paste("Vocalizacoes analisadas:", nrow(dados_cluster)),
    paste("Descritores usados:", ncol(matriz_padronizada)),
    paste("Numero de arvores RF:", NTREE_FINAL),
    paste("Numero de grupos RF + PAM:", melhor_k),
    paste(
      "Silhouette RF + PAM:",
      round(pam_rf_final$silinfo$avg.width, 4)
    ),
    paste("Numero de grupos PCA + PAM:", melhor_k_pca),
    paste(
      "Silhouette PCA + PAM:",
      round(pam_pca_final$silinfo$avg.width, 4)
    ),
    paste(
      "ARI entre RF + PAM e PCA + PAM:",
      round(ari_rf_pca, 4)
    ),
    paste(
      "Mediana ARI da estabilidade RF + PAM:",
      round(stats::median(estabilidade$ARI, na.rm = TRUE), 4)
    ),
    "",
    "Os tipos manuais nao foram usados para formar os grupos.",
    paste(
      "O Random Forest foi utilizado em modo nao supervisionado",
      "para produzir proximidades entre as vocalizacoes."
    )
  ),
  file.path(pasta_saida, "Resumo_classificacao_RF_PAM.txt")
)

message(
  "\nSCRIPT 02 CONCLUIDO.\nNumero de grupos selecionado: ",
  melhor_k,
  "\nResultados em:\n",
  pasta_saida
)

