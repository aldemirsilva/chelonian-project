###############################################################################
# SCRIPT 01 — PRE-PROCESSAMENTO DEFINITIVO DOS DESCRITORES
###############################################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(123)

pacotes <- c("dplyr", "ggplot2", "tidyr", "scales")
instalar <- pacotes[!pacotes %in% rownames(installed.packages())]
if (length(instalar) > 0) install.packages(instalar, dependencies = TRUE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(scales)
})

base <- "C:/Users/abiel/Downloads/Patrick mudancas climaticas"
raiz <- file.path(base, "Pipeline_V2_Podocnemis")
pasta_base <- file.path(raiz, "01_Base_processada")
pasta_saida <- file.path(raiz, "10_Pipeline_RF_PAM_FINAL", "01_Pre_processamento")
pasta_figuras <- file.path(pasta_saida, "Figuras")
dir.create(pasta_figuras, recursive = TRUE, showWarnings = FALSE)

arquivo_dados <- file.path(pasta_base, "dados_modelo.rds")
arquivo_descritores <- file.path(pasta_base, "descritores_validos.rds")

if (!file.exists(arquivo_dados) || !file.exists(arquivo_descritores)) {
  stop("Arquivos do Script 01 anterior nao encontrados.", call. = FALSE)
}

dados <- readRDS(arquivo_dados)
descritores_iniciais <- readRDS(arquivo_descritores)

colunas_obrigatorias <- c(
  "id_vocalizacao", "tipo", "sala", "ambiente", "data", "unidade"
)
faltantes <- setdiff(colunas_obrigatorias, names(dados))
if (length(faltantes) > 0) {
  stop(
    paste("Colunas obrigatorias ausentes:", paste(faltantes, collapse = ", ")),
    call. = FALSE
  )
}

dados <- dados %>%
  dplyr::mutate(
    tipo = factor(tipo, levels = c("tipo1", "tipo2", "tipo3")),
    sala = factor(sala, levels = c("Sala1", "Sala2", "Sala3", "Sala4")),
    ambiente = factor(ambiente, levels = c("Externo", "Interno"))
  )

descritores_presentes <- intersect(descritores_iniciais, names(dados))
matriz_bruta <- dados[, descritores_presentes, drop = FALSE]

descritores_numericos <- names(matriz_bruta)[
  vapply(matriz_bruta, is.numeric, logical(1))
]
matriz_bruta <- matriz_bruta[, descritores_numericos, drop = FALSE]

if (ncol(matriz_bruta) < 2) {
  stop("Menos de dois descritores numericos encontrados.", call. = FALSE)
}

proporcao_problematicos <- vapply(
  matriz_bruta,
  function(x) mean(is.na(x) | !is.finite(x)),
  numeric(1)
)

descritores_apos_ausentes <- names(
  proporcao_problematicos[proporcao_problematicos <= 0.20]
)
matriz_bruta <- matriz_bruta[, descritores_apos_ausentes, drop = FALSE]

linhas_validas <- stats::complete.cases(matriz_bruta) &
  apply(matriz_bruta, 1, function(x) all(is.finite(x)))

dados_analise <- dados[linhas_validas, , drop = FALSE]
matriz_bruta <- matriz_bruta[linhas_validas, , drop = FALSE]

variancias <- vapply(matriz_bruta, stats::var, numeric(1), na.rm = TRUE)
descritores_variaveis <- names(
  variancias[is.finite(variancias) & variancias > 1e-10]
)
matriz_bruta <- matriz_bruta[, descritores_variaveis, drop = FALSE]

remover_correlacionados <- function(x, limite = 0.90) {
  correlacao <- stats::cor(
    x,
    method = "spearman",
    use = "pairwise.complete.obs"
  )
  remover <- character(0)
  
  if (ncol(correlacao) >= 2) {
    for (i in seq_len(ncol(correlacao) - 1)) {
      for (j in (i + 1):ncol(correlacao)) {
        nome_i <- colnames(correlacao)[i]
        nome_j <- colnames(correlacao)[j]
        valor <- correlacao[i, j]
        
        if (
          is.finite(valor) &&
          abs(valor) >= limite &&
          !nome_i %in% remover &&
          !nome_j %in% remover
        ) {
          remover <- c(remover, nome_j)
        }
      }
    }
  }
  
  list(
    manter = setdiff(colnames(x), unique(remover)),
    remover = unique(remover),
    correlacao = correlacao
  )
}

selecao_correlacao <- remover_correlacionados(matriz_bruta, limite = 0.90)
descritores_finais <- selecao_correlacao$manter
descritores_removidos_correlacao <- selecao_correlacao$remover
matriz_final <- matriz_bruta[, descritores_finais, drop = FALSE]

if (ncol(matriz_final) < 2) {
  stop("Menos de dois descritores permaneceram.", call. = FALSE)
}

matriz_padronizada <- scale(matriz_final)

pca <- stats::prcomp(
  matriz_padronizada,
  center = FALSE,
  scale. = FALSE
)

variancia_pca <- 100 * pca$sdev^2 / sum(pca$sdev^2)

scores_pca <- as.data.frame(pca$x[, 1:2, drop = FALSE]) %>%
  dplyr::bind_cols(
    dados_analise %>%
      dplyr::select(id_vocalizacao, tipo, sala, ambiente, data, unidade)
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
    panel.grid.minor = ggplot2::element_blank()
  )
ggplot2::theme_set(tema)

figura_pca_manual <- ggplot2::ggplot(
  scores_pca,
  ggplot2::aes(x = PC1, y = PC2, color = tipo)
) +
  ggplot2::geom_point(alpha = 0.50, size = 1.8) +
  ggplot2::scale_color_manual(values = cores_tipo, labels = rotulos_tipo) +
  ggplot2::labs(
    title = "Espaco acustico das vocalizacoes classificadas manualmente",
    subtitle = "A classificacao manual nao participa da formacao dos grupos computacionais",
    x = paste0("PC1 (", round(variancia_pca[1], 1), "%)"),
    y = paste0("PC2 (", round(variancia_pca[2], 1), "%)"),
    color = "Tipo vocal"
  )

ggplot2::ggsave(
  file.path(pasta_figuras, "Figura_01_PCA_tipos_manuais_sem_elipses.png"),
  figura_pca_manual,
  width = 11,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)

resumo_sala_ambiente <- dados_analise %>%
  dplyr::count(sala, ambiente, name = "n_vocalizacoes")

resumo_tipo_sala <- dados_analise %>%
  dplyr::count(sala, tipo, name = "n_vocalizacoes")

resumo_tipo_sala_ambiente <- dados_analise %>%
  dplyr::count(sala, ambiente, tipo, name = "n_vocalizacoes")

utils::write.csv(
  resumo_sala_ambiente,
  file.path(pasta_saida, "Resumo_vocalizacoes_sala_ambiente.csv"),
  row.names = FALSE
)
utils::write.csv(
  resumo_tipo_sala,
  file.path(pasta_saida, "Resumo_tipos_manuais_por_sala.csv"),
  row.names = FALSE
)
utils::write.csv(
  resumo_tipo_sala_ambiente,
  file.path(pasta_saida, "Resumo_tipos_manuais_sala_ambiente.csv"),
  row.names = FALSE
)

figura_quantidade <- ggplot2::ggplot(
  resumo_sala_ambiente,
  ggplot2::aes(x = sala, y = n_vocalizacoes, fill = ambiente)
) +
  ggplot2::geom_col(
    position = ggplot2::position_dodge(width = 0.80),
    color = "grey25"
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = n_vocalizacoes),
    position = ggplot2::position_dodge(width = 0.80),
    vjust = -0.35,
    size = 4
  ) +
  ggplot2::scale_fill_manual(
    values = c(Externo = "#555555", Interno = "#E66101")
  ) +
  ggplot2::scale_x_discrete(
    labels = c(
      Sala1 = "Sala 1",
      Sala2 = "Sala 2",
      Sala3 = "Sala 3",
      Sala4 = "Sala 4"
    )
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0.12))
  ) +
  ggplot2::labs(
    title = "Numero de vocalizacoes analisadas por sala e condicao de gravacao",
    subtitle = "Interno e externo representam condicoes metodologicas, nao tratamentos climaticos",
    x = "Sala experimental",
    y = "Numero de vocalizacoes",
    fill = "Condicao de gravacao"
  )

ggplot2::ggsave(
  file.path(pasta_figuras, "Figura_02_numero_vocalizacoes_sala_condicao.png"),
  figura_quantidade,
  width = 11,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

saveRDS(dados_analise, file.path(pasta_saida, "dados_analise_final.rds"))
saveRDS(matriz_final, file.path(pasta_saida, "matriz_descritores_final.rds"))
saveRDS(
  matriz_padronizada,
  file.path(pasta_saida, "matriz_descritores_padronizada.rds")
)
saveRDS(descritores_finais, file.path(pasta_saida, "descritores_finais.rds"))
saveRDS(pca, file.path(pasta_saida, "PCA_exploratoria.rds"))

utils::write.csv(
  data.frame(
    descritor = descritores_iniciais,
    presente = descritores_iniciais %in% descritores_presentes,
    selecionado_final = descritores_iniciais %in% descritores_finais,
    removido_correlacao = descritores_iniciais %in%
      descritores_removidos_correlacao
  ),
  file.path(pasta_saida, "Selecao_descritores.csv"),
  row.names = FALSE
)

writeLines(
  c(
    paste("Vocalizacoes originais:", nrow(dados)),
    paste("Vocalizacoes analisadas:", nrow(dados_analise)),
    paste("Descritores iniciais:", length(descritores_iniciais)),
    paste("Descritores finais:", length(descritores_finais)),
    paste(
      "Descritores removidos por correlacao:",
      length(descritores_removidos_correlacao)
    ),
    "",
    "As quatro salas representam os cenarios climaticos experimentais.",
    "Interno e externo representam somente condicoes metodologicas de gravacao."
  ),
  file.path(pasta_saida, "Resumo_pre_processamento.txt")
)

message("\nSCRIPT 01 CONCLUIDO.\nResultados em:\n", pasta_saida)

table(dados_analise$sala)
unique(dados_analise$sala)

