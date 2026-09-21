###############################################################################
# SCRIPT 03 — COMPARACAO MANUAL × COMPUTACIONAL
###############################################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

pacotes <- c("dplyr", "tidyr", "ggplot2", "scales", "ggalluvial")
instalar <- pacotes[!pacotes %in% rownames(installed.packages())]
if (length(instalar) > 0) install.packages(instalar, dependencies = TRUE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(ggalluvial)
})

base <- "C:/Users/abiel/Downloads/Patrick mudancas climaticas"
raiz <- file.path(base, "Pipeline_V2_Podocnemis")

pasta_classificacao <- file.path(
  raiz,
  "10_Pipeline_RF_PAM_FINAL",
  "02_Classificacao_RF_PAM"
)
pasta_saida <- file.path(
  raiz,
  "10_Pipeline_RF_PAM_FINAL",
  "03_Comparacao_manual_computacional"
)
pasta_figuras <- file.path(pasta_saida, "Figuras_relatorio")
dir.create(pasta_figuras, recursive = TRUE, showWarnings = FALSE)

arquivo_grupos <- file.path(
  pasta_classificacao,
  "Grupos_computacionais_por_vocalizacao.csv"
)

if (!file.exists(arquivo_grupos)) {
  stop("Arquivo do Script 02 nao encontrado.", call. = FALSE)
}

dados <- utils::read.csv(
  arquivo_grupos,
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    tipo = factor(tipo, levels = c("tipo1", "tipo2", "tipo3")),
    grupo_computacional = factor(grupo_computacional)
  )

correspondencia <- dados %>%
  dplyr::count(
    tipo,
    grupo_computacional,
    name = "n_vocalizacoes"
  ) %>%
  dplyr::group_by(tipo) %>%
  dplyr::mutate(
    proporcao_dentro_tipo = n_vocalizacoes / sum(n_vocalizacoes)
  ) %>%
  dplyr::ungroup()

utils::write.csv(
  correspondencia,
  file.path(pasta_saida, "Correspondencia_manual_computacional.csv"),
  row.names = FALSE
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
cores_grupos <- c(
  "1" = "#D73027", "2" = "#4575B4", "3" = "#1A9850",
  "4" = "#984EA3", "5" = "#FF7F00", "6" = "#A65628",
  "7" = "#F781BF", "8" = "#666666"
)

tema <- ggplot2::theme_bw(base_size = 14) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 17),
    plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 12),
    axis.title = ggplot2::element_text(face = "bold"),
    panel.grid.minor = ggplot2::element_blank()
  )
ggplot2::theme_set(tema)

figura_porcentagem <- ggplot2::ggplot(
  correspondencia,
  ggplot2::aes(
    x = tipo,
    y = proporcao_dentro_tipo,
    fill = grupo_computacional
  )
) +
  ggplot2::geom_col(color = "grey25") +
  ggplot2::geom_text(
    ggplot2::aes(
      label = scales::percent(
        proporcao_dentro_tipo,
        accuracy = 0.1
      )
    ),
    position = ggplot2::position_stack(vjust = 0.5),
    size = 4
  ) +
  ggplot2::scale_x_discrete(labels = rotulos_tipo) +
  ggplot2::scale_y_continuous(labels = scales::percent_format()) +
  ggplot2::scale_fill_manual(
    values = cores_grupos[
      levels(dados$grupo_computacional)
    ]
  ) +
  ggplot2::labs(
    title = "Correspondencia proporcional entre tipos manuais e grupos computacionais",
    subtitle = "Porcentagem de cada tipo manual atribuida a cada grupo acustico",
    x = "Tipo vocal manual",
    y = "Proporcao das vocalizacoes",
    fill = "Grupo acustico"
  )

ggplot2::ggsave(
  file.path(pasta_figuras, "Figura_08_correspondencia_percentual.png"),
  figura_porcentagem,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

dados_alluvial <- correspondencia %>%
  dplyr::mutate(
    tipo_rotulo = dplyr::recode(
      as.character(tipo),
      tipo1 = "Tipo 1",
      tipo2 = "Tipo 2",
      tipo3 = "Tipo 3"
    ),
    grupo_rotulo = paste("Grupo", grupo_computacional)
  )

figura_alluvial <- ggplot2::ggplot(
  dados_alluvial,
  ggplot2::aes(
    axis1 = tipo_rotulo,
    axis2 = grupo_rotulo,
    y = n_vocalizacoes
  )
) +
  ggalluvial::geom_alluvium(
    ggplot2::aes(fill = tipo),
    width = 0.15,
    alpha = 0.75
  ) +
  ggalluvial::geom_stratum(
    width = 0.15,
    fill = "grey90",
    color = "grey30"
  ) +
  ggplot2::geom_text(
    stat = "stratum",
    ggplot2::aes(label = after_stat(stratum)),
    size = 4
  ) +
  ggplot2::scale_x_discrete(
    limits = c(
      "Classificacao manual",
      "Agrupamento computacional"
    ),
    expand = c(0.08, 0.08)
  ) +
  ggplot2::scale_fill_manual(
    values = cores_tipo,
    labels = rotulos_tipo
  ) +
  ggplot2::labs(
    title = "Fluxo das vocalizacoes entre tipos manuais e grupos computacionais",
    subtitle = "A largura dos fluxos representa o numero de vocalizacoes",
    x = NULL,
    y = "Numero de vocalizacoes",
    fill = "Tipo manual"
  ) +
  ggplot2::theme(panel.grid = ggplot2::element_blank())

ggplot2::ggsave(
  file.path(
    pasta_figuras,
    "Figura_09_alluvial_manual_computacional.png"
  ),
  figura_alluvial,
  width = 11,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)

figura_heatmap <- ggplot2::ggplot(
  correspondencia,
  ggplot2::aes(
    x = grupo_computacional,
    y = tipo,
    fill = n_vocalizacoes
  )
) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::geom_text(
    ggplot2::aes(
      label = paste0(
        n_vocalizacoes,
        "\n",
        scales::percent(
          proporcao_dentro_tipo,
          accuracy = 0.1
        )
      )
    ),
    size = 4
  ) +
  ggplot2::scale_y_discrete(labels = rotulos_tipo) +
  ggplot2::scale_fill_gradient(
    low = "white",
    high = "#4575B4"
  ) +
  ggplot2::labs(
    title = "Matriz de correspondencia manual × computacional",
    subtitle = "Contagem e porcentagem dentro de cada tipo manual",
    x = "Grupo acustico computacional",
    y = "Tipo vocal manual",
    fill = "Numero de\nvocalizacoes"
  )

ggplot2::ggsave(
  file.path(pasta_figuras, "Figura_10_heatmap_correspondencia.png"),
  figura_heatmap,
  width = 9,
  height = 6,
  units = "in",
  dpi = 300,
  bg = "white"
)

message("\nSCRIPT 03 CONCLUIDO.\nResultados em:\n", pasta_saida)

