###############################################################################
# SCRIPT 04 — EFEITOS DOS TRATAMENTOS CLIMATICOS
###############################################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(123)

pacotes <- c("vegan", "dplyr", "tidyr", "ggplot2", "scales")
instalar <- pacotes[!pacotes %in% rownames(installed.packages())]
if (length(instalar) > 0) install.packages(instalar, dependencies = TRUE)

suppressPackageStartupMessages({
  library(vegan)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
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
  "04_Efeito_tratamentos"
)
pasta_figuras <- file.path(pasta_saida, "Figuras_relatorio")
pasta_diagnosticos <- file.path(
  pasta_saida,
  "Diagnosticos_metodologicos"
)

dir.create(pasta_figuras, recursive = TRUE, showWarnings = FALSE)
dir.create(pasta_diagnosticos, recursive = TRUE, showWarnings = FALSE)

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
)

dados <- dados %>%
  dplyr::mutate(
    tipo = factor(
      tipo,
      levels = c("tipo1", "tipo2", "tipo3")
    ),
    sala = factor(
      sala,
      levels = c("Sala1", "Sala2", "Sala3", "Sala4")
    ),
    ambiente = factor(
      ambiente,
      levels = c("Externo", "Interno")
    ),
    grupo_computacional = factor(
      grupo_computacional
    )
  )
arquivo_dados_analise <- file.path(
  raiz,
  "10_Pipeline_RF_PAM_FINAL",
  "01_Pre_processamento",
  "dados_analise_final.rds"
)

dados_descritores <- readRDS(
  arquivo_dados_analise
)


dados_completo <- dados %>%
  left_join(
    dados_descritores,
    by = "id_vocalizacao",
    suffix = c(".rf", ".desc")
  )

dados_completo <- dados_completo %>%
  mutate(
    sala = factor(
      sala.rf,
      levels = c("Sala1","Sala2","Sala3","Sala4")
    ),
    ambiente = factor(
      ambiente.rf,
      levels = c("Externo","Interno")
    ),
    grupo_computacional = factor(grupo_computacional)
  )

criar_composicao_com_zeros <- function(dados_entrada) {
  metadados <- dados_entrada %>%
    dplyr::distinct(unidade, sala, ambiente, data)
  
  grupos <- sort(
    unique(as.character(dados_entrada$grupo_computacional))
  )
  
  grade <- tidyr::crossing(
    metadados,
    grupo_computacional = grupos
  )
  
  contagens <- dados_entrada %>%
    dplyr::count(
      unidade,
      sala,
      ambiente,
      data,
      grupo_computacional,
      name = "n"
    ) %>%
    dplyr::mutate(
      grupo_computacional = as.character(grupo_computacional)
    )
  
  grade %>%
    dplyr::left_join(
      contagens,
      by = c(
        "unidade",
        "sala",
        "ambiente",
        "data",
        "grupo_computacional"
      )
    ) %>%
    dplyr::mutate(
      n = tidyr::replace_na(n, 0L)
    ) %>%
    dplyr::group_by(unidade, sala, ambiente, data) %>%
    dplyr::mutate(
      total = sum(n),
      proporcao = dplyr::if_else(
        total > 0,
        n / total,
        0
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      grupo_computacional = factor(
        grupo_computacional,
        levels = grupos
      )
    )
}

composicao <- criar_composicao_com_zeros(dados)

composicao_larga <- composicao %>%
  dplyr::select(
    unidade,
    sala,
    ambiente,
    data,
    grupo_computacional,
    proporcao
  ) %>%
  tidyr::pivot_wider(
    names_from = grupo_computacional,
    names_prefix = "grupo_",
    values_from = proporcao,
    values_fill = 0
  )

colunas_grupos <- grep(
  "^grupo_",
  names(composicao_larga),
  value = TRUE
)

matriz_composicao <- as.matrix(
  composicao_larga[, colunas_grupos, drop = FALSE]
)

matriz_hellinger <- vegan::decostand(
  matriz_composicao,
  method = "hellinger"
)

# Teste principal:
# sala = tratamento biologico;
# ambiente = controle metodologico da condicao de gravacao.
permanova_principal <- vegan::adonis2(
  matriz_hellinger ~ ambiente + sala,
  data = composicao_larga,
  permutations = 9999,
  method = "euclidean",
  by = "margin"
)

dispersao_sala <- vegan::betadisper(
  stats::dist(matriz_hellinger),
  composicao_larga$sala
)

permdisp_sala <- vegan::permutest(
  dispersao_sala,
  permutations = 9999
)

analisar_sala_em_ambiente <- function(ambiente_atual) {
  dados_sub <- composicao_larga %>%
    dplyr::filter(ambiente == ambiente_atual) %>%
    droplevels()
  
  matriz_sub <- as.matrix(
    dados_sub[, colunas_grupos, drop = FALSE]
  )
  
  matriz_sub <- vegan::decostand(
    matriz_sub,
    method = "hellinger"
  )
  
  teste <- vegan::adonis2(
    matriz_sub ~ sala,
    data = dados_sub,
    permutations = 9999,
    method = "euclidean"
  )
  
  data.frame(
    analise = paste("Sala - somente", ambiente_atual),
    termo = "sala",
    F = teste$F[1],
    R2 = teste$R2[1],
    p = teste$`Pr(>F)`[1]
  )
}

resultado_interno <- analisar_sala_em_ambiente("Interno")
resultado_externo <- analisar_sala_em_ambiente("Externo")

formatar_p <- function(p) {
  ifelse(
    is.na(p),
    NA_character_,
    ifelse(
      p < 0.001,
      "< 0,001",
      sub("\\.", ",", sprintf("%.3f", p))
    )
  )
}

resultado_principal <- data.frame(
  analise = c(
    "Condicao de gravacao",
    "Sala controlando condicao de gravacao"
  ),
  termo = c("ambiente", "sala"),
  F = c(
    permanova_principal$F[
      rownames(permanova_principal) == "ambiente"
    ],
    permanova_principal$F[
      rownames(permanova_principal) == "sala"
    ]
  ),
  R2 = c(
    permanova_principal$R2[
      rownames(permanova_principal) == "ambiente"
    ],
    permanova_principal$R2[
      rownames(permanova_principal) == "sala"
    ]
  ),
  p = c(
    permanova_principal$`Pr(>F)`[
      rownames(permanova_principal) == "ambiente"
    ],
    permanova_principal$`Pr(>F)`[
      rownames(permanova_principal) == "sala"
    ]
  )
)

resultados_permanova <- dplyr::bind_rows(
  resultado_principal,
  resultado_interno,
  resultado_externo
) %>%
  dplyr::mutate(
    p_formatado = formatar_p(p),
    resultado = dplyr::if_else(
      p < 0.05,
      "Significativo",
      "Nao significativo"
    )
  )

resultado_permdisp <- data.frame(
  analise = "PERMDISP entre salas",
  F = permdisp_sala$tab$F[1],
  p = permdisp_sala$tab$`Pr(>F)`[1],
  p_formatado = formatar_p(
    permdisp_sala$tab$`Pr(>F)`[1]
  )
)

utils::write.csv(
  resultados_permanova,
  file.path(pasta_saida, "Resultados_PERMANOVA.csv"),
  row.names = FALSE
)
utils::write.csv(
  resultado_permdisp,
  file.path(pasta_saida, "Resultado_PERMDISP.csv"),
  row.names = FALSE
)
utils::write.csv(
  composicao,
  file.path(pasta_saida, "Composicao_grupos_por_sessao.csv"),
  row.names = FALSE
)

resultados_kruskal <- dplyr::bind_rows(
  lapply(
    levels(composicao$grupo_computacional),
    function(grupo_atual) {
      dados_grupo <- composicao %>%
        dplyr::filter(grupo_computacional == grupo_atual)
      
      teste <- stats::kruskal.test(
        proporcao ~ sala,
        data = dados_grupo
      )
      
      data.frame(
        grupo_computacional = grupo_atual,
        H = unname(teste$statistic),
        gl = unname(teste$parameter),
        p = teste$p.value
      )
    }
  )
) %>%
  dplyr::mutate(
    p_ajustado_BH = stats::p.adjust(p, method = "BH"),
    p_formatado = formatar_p(p_ajustado_BH)
  )

utils::write.csv(
  resultados_kruskal,
  file.path(pasta_saida, "Resultados_Kruskal_por_grupo.csv"),
  row.names = FALSE
)

quantidade_sala <- dados %>%
  dplyr::count(
    sala,
    grupo_computacional,
    name = "n_vocalizacoes"
  )

quantidade_sala_ambiente <- dados %>%
  dplyr::count(
    sala,
    ambiente,
    grupo_computacional,
    name = "n_vocalizacoes"
  )

totais_sessao <- dados %>%
  dplyr::count(
    unidade,
    sala,
    ambiente,
    data,
    name = "n_vocalizacoes"
  )

pcoa <- stats::cmdscale(
  stats::dist(matriz_hellinger),
  k = 2,
  eig = TRUE,
  add = TRUE
)

scores_pcoa <- dplyr::bind_cols(
  composicao_larga %>%
    dplyr::select(unidade, sala, ambiente, data),
  setNames(
    as.data.frame(pcoa$points),
    c("PCoA1", "PCoA2")
  )
)

set.seed(123)
nmds <- vegan::metaMDS(
  matriz_hellinger,
  distance = "euclidean",
  k = 2,
  trymax = 100,
  autotransform = FALSE,
  trace = FALSE
)

scores_nmds <- as.data.frame(
  vegan::scores(nmds, display = "sites")
) %>%
  dplyr::bind_cols(
    composicao_larga %>%
      dplyr::select(unidade, sala, ambiente, data)
  )

cores_sala <- c(
  Sala1 = "#222222",
  Sala2 = "#E66101",
  Sala3 = "#5E3C99",
  Sala4 = "#1B9E77"
)

cores_grupos <- c(
  "1" = "#D73027", "2" = "#4575B4", "3" = "#1A9850",
  "4" = "#984EA3", "5" = "#FF7F00", "6" = "#A65628",
  "7" = "#F781BF", "8" = "#666666"
)

rotulos_sala <- c(
  Sala1 = "Sala 1",
  Sala2 = "Sala 2",
  Sala3 = "Sala 3",
  Sala4 = "Sala 4"
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

figura_quantidade_sala <- ggplot2::ggplot(
  quantidade_sala,
  ggplot2::aes(
    x = sala,
    y = n_vocalizacoes,
    fill = grupo_computacional
  )
) +
  ggplot2::geom_col(color = "grey25") +
  ggplot2::geom_text(
    ggplot2::aes(
      label = ifelse(
        n_vocalizacoes > 0,
        n_vocalizacoes,
        ""
      )
    ),
    position = ggplot2::position_stack(vjust = 0.5),
    size = 4
  ) +
  ggplot2::scale_x_discrete(labels = rotulos_sala) +
  ggplot2::scale_fill_manual(
    values = cores_grupos[
      levels(dados$grupo_computacional)
    ]
  ) +
  ggplot2::labs(
    title = "Quantidade dos grupos acusticos por cenario climatico",
    subtitle = "Numeros absolutos de vocalizacoes atribuidas a cada grupo",
    x = "Cenario climatico experimental / sala",
    y = "Numero de vocalizacoes",
    fill = "Grupo acustico"
  )

ggplot2::ggsave(
  file.path(pasta_figuras, "Figura_11_quantidade_grupos_por_sala.png"),
  figura_quantidade_sala,
  width = 11,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

figura_sessoes <- ggplot2::ggplot(
  composicao,
  ggplot2::aes(
    x = sala,
    y = proporcao,
    fill = sala
  )
) +
  ggplot2::geom_boxplot(
    alpha = 0.82,
    outlier.shape = NA
  ) +
  ggplot2::geom_point(
    position = ggplot2::position_jitter(width = 0.06),
    size = 2.4
  ) +
  ggplot2::facet_wrap(
    ~ grupo_computacional,
    nrow = 1,
    labeller = ggplot2::labeller(
      grupo_computacional = function(x) paste("Grupo", x)
    )
  ) +
  ggplot2::scale_fill_manual(values = cores_sala) +
  ggplot2::scale_x_discrete(labels = rotulos_sala) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format()
  ) +
  ggplot2::guides(fill = "none") +
  ggplot2::labs(
    title = "Composicao dos grupos acusticos entre os cenarios climaticos",
    subtitle = "Cada ponto representa uma sessao de gravacao",
    x = "Cenario climatico experimental / sala",
    y = "Proporcao do grupo na sessao"
  )

ggplot2::ggsave(
  file.path(pasta_figuras, "Figura_12_composicao_grupos_por_sessao.png"),
  figura_sessoes,
  width = 14,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

figura_pcoa <- ggplot2::ggplot(
  scores_pcoa,
  ggplot2::aes(
    x = PCoA1,
    y = PCoA2,
    color = sala,
    shape = ambiente
  )
) +
  ggplot2::geom_point(size = 4) +
  ggplot2::scale_color_manual(
    values = cores_sala,
    labels = rotulos_sala
  ) +
  ggplot2::labs(
    title = "Ordenacao da composicao dos grupos por sessao",
    subtitle = "Cores representam os cenarios climaticos; formas representam a condicao de gravacao",
    x = "PCoA1",
    y = "PCoA2",
    color = "Cenario / sala",
    shape = "Condicao de gravacao"
  )

ggplot2::ggsave(
  file.path(pasta_figuras, "Figura_13_PCoA_sessoes_tratamentos.png"),
  figura_pcoa,
  width = 11,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)

figura_nmds <- ggplot2::ggplot(
  scores_nmds,
  ggplot2::aes(
    x = NMDS1,
    y = NMDS2,
    color = sala,
    shape = ambiente
  )
) +
  ggplot2::geom_point(size = 4) +
  ggplot2::scale_color_manual(
    values = cores_sala,
    labels = rotulos_sala
  ) +
  ggplot2::labs(
    title = "NMDS da composicao dos grupos acusticos",
    subtitle = paste0(
      "Visualizacao complementar; stress = ",
      round(nmds$stress, 3)
    ),
    x = "NMDS1",
    y = "NMDS2",
    color = "Cenario / sala",
    shape = "Condicao de gravacao"
  )

ggplot2::ggsave(
  file.path(pasta_diagnosticos, "Complementar_NMDS_sessoes.png"),
  figura_nmds,
  width = 11,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)

perfil_sala <- composicao %>%
  dplyr::group_by(sala, grupo_computacional) %>%
  dplyr::summarise(
    media = mean(proporcao, na.rm = TRUE),
    mediana = stats::median(proporcao, na.rm = TRUE),
    .groups = "drop"
  )

figura_perfil <- ggplot2::ggplot(
  perfil_sala,
  ggplot2::aes(
    x = sala,
    y = media,
    group = grupo_computacional,
    color = grupo_computacional
  )
) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::geom_point(size = 3) +
  ggplot2::scale_color_manual(
    values = cores_grupos[
      levels(composicao$grupo_computacional)
    ]
  ) +
  ggplot2::scale_x_discrete(labels = rotulos_sala) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format()
  ) +
  ggplot2::labs(
    title = "Perfil medio da composicao acustica entre os cenarios",
    subtitle = "Media das proporcoes observadas nas sessoes de cada sala",
    x = "Cenario climatico experimental / sala",
    y = "Proporcao media",
    color = "Grupo acustico"
  )

ggplot2::ggsave(
  file.path(pasta_figuras, "Figura_14_perfil_medio_grupos_por_sala.png"),
  figura_perfil,
  width = 11,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

heatmap_sala <- perfil_sala %>%
  dplyr::mutate(
    sala_rotulo = dplyr::recode(
      as.character(sala),
      Sala1 = "Sala 1",
      Sala2 = "Sala 2",
      Sala3 = "Sala 3",
      Sala4 = "Sala 4"
    )
  )

figura_heatmap_sala <- ggplot2::ggplot(
  heatmap_sala,
  ggplot2::aes(
    x = grupo_computacional,
    y = sala_rotulo,
    fill = media
  )
) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::geom_text(
    ggplot2::aes(
      label = scales::percent(media, accuracy = 0.1)
    ),
    size = 4
  ) +
  ggplot2::scale_fill_gradient(
    low = "white",
    high = "#4575B4",
    labels = scales::percent_format()
  ) +
  ggplot2::labs(
    title = "Composicao media dos grupos em cada cenario climatico",
    subtitle = "Porcentagem media calculada a partir das sessoes de gravacao",
    x = "Grupo acustico",
    y = "Cenario / sala",
    fill = "Proporcao\nmedia"
  )

ggplot2::ggsave(
  file.path(pasta_figuras, "Figura_15_heatmap_sala_grupo.png"),
  figura_heatmap_sala,
  width = 9,
  height = 6,
  units = "in",
  dpi = 300,
  bg = "white"
)

resultados_plot <- resultados_permanova %>%
  dplyr::mutate(
    analise = factor(analise, levels = rev(analise)),
    rotulo = paste0(
      "R2 = ",
      sprintf("%.3f", R2),
      "; p = ",
      p_formatado
    ),
    classe = dplyr::case_when(
      p < 0.05 ~ "Significativo",
      p < 0.10 ~ "Evidencia limitada",
      TRUE ~ "Nao significativo"
    )
  )

figura_efeitos <- ggplot2::ggplot(
  resultados_plot,
  ggplot2::aes(
    x = R2,
    y = analise,
    fill = classe
  )
) +
  ggplot2::geom_col(width = 0.65) +
  ggplot2::geom_text(
    ggplot2::aes(label = rotulo),
    hjust = -0.04,
    size = 4
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Significativo" = "#1B9E77",
      "Evidencia limitada" = "#E66101",
      "Nao significativo" = "#777777"
    )
  ) +
  ggplot2::scale_x_continuous(
    limits = c(
      0,
      max(resultados_plot$R2, na.rm = TRUE) * 1.55
    ),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  ) +
  ggplot2::labs(
    title = "Variacao explicada nas comparacoes entre os cenarios climaticos",
    subtitle = "Sala representa o tratamento; a condicao de gravacao e um controle metodologico",
    x = "Variacao explicada (R2)",
    y = NULL,
    fill = "Resultado"
  )

ggplot2::ggsave(
  file.path(pasta_figuras, "Figura_16_resumo_PERMANOVA_tratamentos.png"),
  figura_efeitos,
  width = 12,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)

figura_ambiente <- ggplot2::ggplot(
  quantidade_sala_ambiente,
  ggplot2::aes(
    x = sala,
    y = n_vocalizacoes,
    fill = grupo_computacional
  )
) +
  ggplot2::geom_col(color = "grey25") +
  ggplot2::facet_wrap(~ ambiente, nrow = 1) +
  ggplot2::scale_fill_manual(
    values = cores_grupos[
      levels(dados$grupo_computacional)
    ]
  ) +
  ggplot2::scale_x_discrete(labels = rotulos_sala) +
  ggplot2::labs(
    title = "Diagnostico da condicao de gravacao",
    subtitle = "Interno e externo avaliam possivel influencia dos ruidos das salas",
    x = "Sala experimental",
    y = "Numero de vocalizacoes",
    fill = "Grupo acustico"
  )

ggplot2::ggsave(
  file.path(pasta_diagnosticos, "Diagnostico_condicao_gravacao.png"),
  figura_ambiente,
  width = 14,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

linha_sala <- resultados_permanova %>%
  dplyr::filter(
    analise == "Sala controlando condicao de gravacao"
  ) %>%
  dplyr::slice_head(n = 1)

writeLines(
  c(
    paste("Vocalizacoes:", nrow(dados)),
    paste("Sessoes:", dplyr::n_distinct(dados$unidade)),
    paste(
      "Grupos acusticos:",
      nlevels(dados$grupo_computacional)
    ),
    "",
    paste(
      "Resultado principal para sala, controlando condicao de gravacao:",
      "R2 =",
      round(linha_sala$R2, 4),
      "; p =",
      formatar_p(linha_sala$p)
    ),
    paste(
      "PERMDISP entre salas: p =",
      formatar_p(resultado_permdisp$p)
    ),
    "",
    "Sala responde aos cenarios climaticos experimentais.",
    "Interno e externo sao apenas condicoes metodologicas de gravacao.",
    paste(
      "Como cada cenario esta representado por uma unica sala,",
      "os resultados devem ser interpretados como exploratorios."
    )
  ),
  file.path(pasta_saida, "Resumo_efeito_tratamentos.txt")
)

figura_17 <- ggplot(
  perfil_sala,
  aes(
    x = grupo_computacional,
    y = media,
    fill = sala
  )
) +
  geom_col(
    position = "dodge"
  ) +
  scale_fill_manual(values = cores_sala) +
  scale_y_continuous(
    labels = scales::percent_format()
  ) +
  scale_x_discrete(
    labels = function(x) paste("Grupo", x)
  ) +
  labs(
    title = "Frequência relativa dos grupos acústicos por tratamento",
    subtitle = "Proporção média observada em cada sala experimental",
    x = "Grupo acústico",
    y = "Proporção média",
    fill = "Tratamento"
  )

ggsave(
  file.path(
    pasta_figuras,
    "Figura_17_frequencia_relativa_grupos.png"
  ),
  figura_17,
  width = 11,
  height = 7,
  dpi = 300
)

controle <- perfil_sala %>%
  filter(sala == "Sala1") %>%
  select(
    grupo_computacional,
    media_controle = media
  )

comparacao <- perfil_sala %>%
  left_join(
    controle,
    by = "grupo_computacional"
  ) %>%
  mutate(
    diferenca = media - media_controle
  )

figura_18 <- ggplot(
  comparacao %>% filter(sala != "Sala1"),
  aes(
    x = grupo_computacional,
    y = sala,
    fill = diferenca
  )
) +
  geom_tile(color = "white") +
  geom_text(
    aes(
      label = round(diferenca, 3)
    ),
    size = 4
  ) +
  scale_fill_gradient2(
    low = "#D73027",
    mid = "white",
    high = "#1A9850",
    midpoint = 0
  ) +
  scale_y_discrete(labels = rotulos_sala) +
  labs(
    title = "Diferença em relação ao controle",
    subtitle = "Valores positivos indicam aumento relativo do grupo acústico",
    x = "Grupo acústico",
    y = "Tratamento",
    fill = "Diferença"
  )

ggsave(
  file.path(
    pasta_figuras,
    "Figura_18_heatmap_diferenca_controle.png"
  ),
  figura_18,
  width = 10,
  height = 6,
  dpi = 300
)

figura_19 <- ggplot(
  resultados_plot,
  aes(
    x = reorder(analise, R2),
    y = R2,
    fill = classe
  )
) +
  geom_col() +
  coord_flip() +
  geom_text(
    aes(
      label = paste0(
        "p = ",
        p_formatado
      )
    ),
    hjust = -0.15,
    size = 4
  ) +
  scale_fill_manual(
    values = c(
      "Significativo" = "#1B9E77",
      "Evidencia limitada" = "#E66101",
      "Nao significativo" = "#777777"
    )
  ) +
  expand_limits(y = max(resultados_plot$R2) * 1.25) +
  labs(
    title = "Magnitude dos efeitos detectados",
    subtitle = "Comparação entre os testes realizados",
    x = NULL,
    y = "R²",
    fill = "Resultado"
  )

ggsave(
  file.path(
    pasta_figuras,
    "Figura_19_ranking_efeitos.png"
  ),
  figura_19,
  width = 10,
  height = 6,
  dpi = 300
)

###############################################################################
# COMPLEMENTO — ESTRUTURA ACUSTICA ENTRE OS TRATAMENTOS
###############################################################################

# Recuperar exatamente os 17 descritores selecionados no Script 01

arquivo_descritores_finais <- file.path(
  raiz,
  "10_Pipeline_RF_PAM_FINAL",
  "01_Pre_processamento",
  "descritores_finais.rds"
)

descritores <- readRDS(
  arquivo_descritores_finais
)

# Conferencia
message(
  "Descritores usados nesta etapa: ",
  paste(descritores, collapse = ", ")
)

###############################################################################
# ANALISE DOS DESCRITORES ENTRE AS SALAS
###############################################################################

# Para esta etapa usamos diretamente dados_descritores,
# que contem os descritores acusticos e os metadados originais.

dados_acusticos <- dados_descritores %>%
  dplyr::mutate(
    sala = factor(
      sala,
      levels = c("Sala1", "Sala2", "Sala3", "Sala4")
    ),
    ambiente = factor(
      ambiente,
      levels = c("Externo", "Interno")
    ),
    tipo = factor(
      tipo,
      levels = c("tipo1", "tipo2", "tipo3")
    )
  )

# Garantia de que somente descritores numericos serao utilizados

descritores <- intersect(
  descritores,
  names(dados_acusticos)
)

descritores <- descritores[
  vapply(
    dados_acusticos[, descritores, drop = FALSE],
    is.numeric,
    logical(1)
  )
]

message(
  "Numero de descritores acusticos analisados: ",
  length(descritores)
)

###############################################################################
# KRUSKAL-WALLIS — DESCRITOR × SALA
###############################################################################

resultados_descritores_sala <- dplyr::bind_rows(
  lapply(
    descritores,
    function(descritor_atual) {
      
      formula_teste <- stats::as.formula(
        paste0("`", descritor_atual, "` ~ sala")
      )
      
      teste <- stats::kruskal.test(
        formula_teste,
        data = dados_acusticos
      )
      
      n <- nrow(dados_acusticos)
      k <- nlevels(dados_acusticos$sala)
      
      epsilon2 <- max(
        0,
        (
          unname(teste$statistic) - k + 1
        ) /
          (n - k)
      )
      
      data.frame(
        descritor = descritor_atual,
        H = unname(teste$statistic),
        gl = unname(teste$parameter),
        p = teste$p.value,
        epsilon2 = epsilon2
      )
    }
  )
) %>%
  dplyr::mutate(
    p_ajustado_BH = stats::p.adjust(
      p,
      method = "BH"
    ),
    p_formatado = formatar_p(
      p_ajustado_BH
    )
  ) %>%
  dplyr::arrange(
    dplyr::desc(epsilon2)
  )

utils::write.csv(
  resultados_descritores_sala,
  file.path(
    pasta_saida,
    "Resultados_Kruskal_Descritores_Sala.csv"
  ),
  row.names = FALSE
)

###############################################################################
# FIGURA 20 — HEATMAP DOS DESCRITORES ENTRE AS SALAS
###############################################################################

# Primeiro padronizamos cada vocalizacao.
# Isso permite comparar descritores medidos em escalas diferentes.

matriz_z <- scale(
  dados_acusticos[, descritores, drop = FALSE]
)

dados_z <- as.data.frame(
  matriz_z
)

dados_z$sala <- dados_acusticos$sala

perfil_descritores <- dados_z %>%
  dplyr::group_by(sala) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(descritores),
      \(x) mean(x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(
    cols = -sala,
    names_to = "descritor",
    values_to = "media_padronizada"
  )

# Ordenar descritores pelo tamanho de efeito

ordem_descritores <- resultados_descritores_sala$descritor

perfil_descritores$descritor <- factor(
  perfil_descritores$descritor,
  levels = rev(ordem_descritores)
)

figura_20 <- ggplot2::ggplot(
  perfil_descritores,
  ggplot2::aes(
    x = sala,
    y = descritor,
    fill = media_padronizada
  )
) +
  ggplot2::geom_tile(
    color = "white",
    linewidth = 0.6
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = sprintf(
        "%.2f",
        media_padronizada
      )
    ),
    size = 3.6
  ) +
  ggplot2::scale_x_discrete(
    labels = rotulos_sala
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#4575B4",
    mid = "white",
    high = "#D73027",
    midpoint = 0
  ) +
  ggplot2::labs(
    title = "Perfil acústico das vocalizações entre os tratamentos",
    subtitle = "Média padronizada dos descritores acústicos em cada sala experimental",
    x = "Tratamento",
    y = "Descritor acústico",
    fill = "Média\npadronizada"
  ) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank()
  )

ggplot2::ggsave(
  file.path(
    pasta_figuras,
    "Figura_20_heatmap_descritores_por_sala.png"
  ),
  figura_20,
  width = 11,
  height = 9,
  units = "in",
  dpi = 300,
  bg = "white"
)

###############################################################################
# FIGURA 21 — PCA DOS DESCRITORES POR TRATAMENTO
###############################################################################

pca_tratamentos <- stats::prcomp(
  matriz_z,
  center = FALSE,
  scale. = FALSE
)

variancia_pca_tratamentos <- (
  100 * pca_tratamentos$sdev^2 /
    sum(pca_tratamentos$sdev^2)
)

scores_tratamentos <- as.data.frame(
  pca_tratamentos$x[, 1:2, drop = FALSE]
)

scores_tratamentos$sala <- dados_acusticos$sala

figura_21 <- ggplot2::ggplot(
  scores_tratamentos,
  ggplot2::aes(
    x = PC1,
    y = PC2,
    color = sala
  )
) +
  ggplot2::geom_point(
    alpha = 0.50,
    size = 1.9
  ) +
  ggplot2::scale_color_manual(
    values = cores_sala,
    labels = rotulos_sala
  ) +
  ggplot2::labs(
    title = "Espaço acústico das vocalizações entre os tratamentos",
    subtitle = "PCA calculada a partir dos descritores acústicos padronizados",
    x = paste0(
      "PC1 (",
      round(
        variancia_pca_tratamentos[1],
        1
      ),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(
        variancia_pca_tratamentos[2],
        1
      ),
      "%)"
    ),
    color = "Tratamento"
  )

ggplot2::ggsave(
  file.path(
    pasta_figuras,
    "Figura_21_PCA_descritores_por_sala.png"
  ),
  figura_21,
  width = 11,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)

###############################################################################
# FIGURA 22 — DESCRITORES COM MAIOR DIFERENCA ENTRE SALAS
###############################################################################

top_descritores_sala <- resultados_descritores_sala %>%
  dplyr::slice_head(
    n = min(
      8,
      nrow(resultados_descritores_sala)
    )
  )

figura_22 <- ggplot2::ggplot(
  top_descritores_sala %>%
    dplyr::mutate(
      descritor = factor(
        descritor,
        levels = rev(descritor)
      )
    ),
  ggplot2::aes(
    x = epsilon2,
    y = descritor
  )
) +
  ggplot2::geom_col(
    fill = "#4575B4",
    width = 0.72
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = paste0(
        "p = ",
        p_formatado
      )
    ),
    hjust = -0.08,
    size = 3.8
  ) +
  ggplot2::scale_x_continuous(
    limits = c(
      0,
      max(
        top_descritores_sala$epsilon2,
        na.rm = TRUE
      ) * 1.35
    ),
    expand = ggplot2::expansion(
      mult = c(0, 0.02)
    )
  ) +
  ggplot2::labs(
    title = "Descritores acústicos que mais variam entre os tratamentos",
    subtitle = "Tamanho de efeito do teste de Kruskal-Wallis",
    x = "Tamanho de efeito (epsilon quadrado)",
    y = "Descritor acústico"
  )

ggplot2::ggsave(
  file.path(
    pasta_figuras,
    "Figura_22_descritores_maior_efeito_sala.png"
  ),
  figura_22,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

###############################################################################
# FIGURA 23 — DISTRIBUICAO DOS PRINCIPAIS DESCRITORES
###############################################################################

descritores_box <- top_descritores_sala$descritor[
  seq_len(
    min(
      6,
      nrow(top_descritores_sala)
    )
  )
]

dados_box <- dados_acusticos %>%
  dplyr::select(
    sala,
    dplyr::all_of(descritores_box)
  ) %>%
  tidyr::pivot_longer(
    cols = -sala,
    names_to = "descritor",
    values_to = "valor"
  )

dados_box$descritor <- factor(
  dados_box$descritor,
  levels = descritores_box
)

figura_23 <- ggplot2::ggplot(
  dados_box,
  ggplot2::aes(
    x = sala,
    y = valor,
    fill = sala
  )
) +
  ggplot2::geom_boxplot(
    alpha = 0.78,
    outlier.alpha = 0.20,
    linewidth = 0.5
  ) +
  ggplot2::facet_wrap(
    ~ descritor,
    scales = "free_y",
    ncol = 3
  ) +
  ggplot2::scale_fill_manual(
    values = cores_sala
  ) +
  ggplot2::scale_x_discrete(
    labels = rotulos_sala
  ) +
  ggplot2::guides(
    fill = "none"
  ) +
  ggplot2::labs(
    title = "Distribuição dos descritores acústicos com maior variação",
    subtitle = "Comparação das vocalizações entre as quatro salas experimentais",
    x = "Tratamento",
    y = "Valor do descritor"
  )

ggplot2::ggsave(
  file.path(
    pasta_figuras,
    "Figura_23_boxplots_descritores_por_sala.png"
  ),
  figura_23,
  width = 13,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)

message(
  "\nANALISE DOS DESCRITORES ENTRE SALAS CONCLUIDA.\n",
  "Resultados em:\n",
  pasta_saida
)

###############################################################################
# VALIDACAO EXPLORATORIA DOS DESCRITORES POR SESSAO
#
# Objetivo:
# verificar se os padroes observados entre as 854 vocalizacoes permanecem
# quando cada sessao de gravacao e tratada como unidade analitica.
###############################################################################

###############################################################################
# 1. RESUMO DOS DESCRITORES POR SESSAO
###############################################################################

# A mediana e usada porque varios descritores possuem distribuicoes
# assimetricas e valores extremos.

dados_sessao_descritores <- dados_acusticos %>%
  dplyr::group_by(
    unidade,
    sala,
    ambiente,
    data
  ) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(descritores),
      \(x) stats::median(x, na.rm = TRUE)
    ),
    n_vocalizacoes = dplyr::n(),
    .groups = "drop"
  )

message(
  "Numero de sessoes usadas na validacao: ",
  nrow(dados_sessao_descritores)
)

print(
  table(dados_sessao_descritores$sala)
)

utils::write.csv(
  dados_sessao_descritores,
  file.path(
    pasta_saida,
    "Descritores_resumidos_por_sessao.csv"
  ),
  row.names = FALSE
)

###############################################################################
# 2. KRUSKAL-WALLIS ENTRE SALAS — AGORA POR SESSAO
###############################################################################

resultados_descritores_sessao <- dplyr::bind_rows(
  lapply(
    descritores,
    function(descritor_atual) {
      
      formula_teste <- stats::as.formula(
        paste0(
          "`",
          descritor_atual,
          "` ~ sala"
        )
      )
      
      teste <- stats::kruskal.test(
        formula_teste,
        data = dados_sessao_descritores
      )
      
      n <- nrow(dados_sessao_descritores)
      k <- nlevels(
        droplevels(
          dados_sessao_descritores$sala
        )
      )
      
      epsilon2 <- max(
        0,
        (
          unname(teste$statistic) - k + 1
        ) /
          (n - k)
      )
      
      data.frame(
        descritor = descritor_atual,
        H = unname(teste$statistic),
        gl = unname(teste$parameter),
        p = teste$p.value,
        epsilon2 = epsilon2
      )
    }
  )
) %>%
  dplyr::mutate(
    p_ajustado_BH = stats::p.adjust(
      p,
      method = "BH"
    ),
    p_formatado = formatar_p(
      p_ajustado_BH
    ),
    resultado = dplyr::case_when(
      p_ajustado_BH < 0.05 ~ "Significativo",
      p_ajustado_BH < 0.10 ~ "Evidencia limitada",
      TRUE ~ "Nao significativo"
    )
  ) %>%
  dplyr::arrange(
    dplyr::desc(epsilon2)
  )

utils::write.csv(
  resultados_descritores_sessao,
  file.path(
    pasta_saida,
    "Resultados_Kruskal_Descritores_por_Sessao.csv"
  ),
  row.names = FALSE
)

print(
  resultados_descritores_sessao
)

###############################################################################
# 3. COMPARACOES PAR A PAR ENTRE AS SALAS
#
# O Wilcoxon e usado de forma exploratoria.
# Com poucas sessoes por sala, o poder estatistico sera baixo.
###############################################################################

comparacoes_pairwise <- dplyr::bind_rows(
  lapply(
    descritores,
    function(descritor_atual) {
      
      x <- dados_sessao_descritores[
        ,
        c("sala", descritor_atual)
      ]
      
      names(x)[2] <- "valor"
      
      teste <- stats::pairwise.wilcox.test(
        x$valor,
        x$sala,
        p.adjust.method = "BH",
        exact = FALSE
      )
      
      matriz_p <- teste$p.value
      
      if (is.null(matriz_p)) {
        return(NULL)
      }
      
      resultado <- as.data.frame(
        as.table(matriz_p),
        stringsAsFactors = FALSE
      )
      
      names(resultado) <- c(
        "sala_1",
        "sala_2",
        "p_ajustado"
      )
      
      resultado %>%
        dplyr::filter(
          !is.na(p_ajustado)
        ) %>%
        dplyr::mutate(
          descritor = descritor_atual,
          p_formatado = formatar_p(
            p_ajustado
          )
        ) %>%
        dplyr::select(
          descritor,
          sala_1,
          sala_2,
          p_ajustado,
          p_formatado
        )
    }
  )
)

utils::write.csv(
  comparacoes_pairwise,
  file.path(
    pasta_saida,
    "Comparacoes_pairwise_Descritores_por_Sessao.csv"
  ),
  row.names = FALSE
)

###############################################################################
# 4. COMPARAR ANALISE POR VOCALIZACAO × ANALISE POR SESSAO
###############################################################################

comparacao_niveis <- resultados_descritores_sala %>%
  dplyr::select(
    descritor,
    epsilon2_vocalizacao = epsilon2,
    p_vocalizacao = p_ajustado_BH
  ) %>%
  dplyr::left_join(
    resultados_descritores_sessao %>%
      dplyr::select(
        descritor,
        epsilon2_sessao = epsilon2,
        p_sessao = p_ajustado_BH
      ),
    by = "descritor"
  ) %>%
  dplyr::mutate(
    permanece_significativo =
      p_vocalizacao < 0.05 &
      p_sessao < 0.05
  )

utils::write.csv(
  comparacao_niveis,
  file.path(
    pasta_saida,
    "Comparacao_Vocalizacao_vs_Sessao.csv"
  ),
  row.names = FALSE
)

###############################################################################
# 5. FIGURA 24 — TAMANHO DE EFEITO POR SESSAO
###############################################################################

top_sessao <- resultados_descritores_sessao %>%
  dplyr::slice_head(
    n = min(
      10,
      nrow(resultados_descritores_sessao)
    )
  )

figura_24 <- ggplot2::ggplot(
  top_sessao %>%
    dplyr::mutate(
      descritor = factor(
        descritor,
        levels = rev(descritor)
      )
    ),
  ggplot2::aes(
    x = epsilon2,
    y = descritor,
    fill = resultado
  )
) +
  ggplot2::geom_col(
    width = 0.72
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = paste0(
        "p = ",
        p_formatado
      )
    ),
    hjust = -0.08,
    size = 3.7
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Significativo" = "#1B9E77",
      "Evidencia limitada" = "#E66101",
      "Nao significativo" = "#777777"
    )
  ) +
  ggplot2::scale_x_continuous(
    limits = c(
      0,
      max(
        top_sessao$epsilon2,
        na.rm = TRUE
      ) * 1.40
    ),
    expand = ggplot2::expansion(
      mult = c(0, 0.02)
    )
  ) +
  ggplot2::labs(
    title = "Variação dos descritores acústicos entre tratamentos",
    subtitle = "Análise exploratória utilizando as sessões de gravação como unidades analíticas",
    x = "Tamanho de efeito (epsilon quadrado)",
    y = "Descritor acústico",
    fill = "Resultado"
  )

ggplot2::ggsave(
  file.path(
    pasta_figuras,
    "Figura_24_descritores_por_sessao_efeito.png"
  ),
  figura_24,
  width = 11,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)

###############################################################################
# 6. FIGURA 25 — DISTRIBUICAO DAS SESSOES ENTRE AS SALAS
#
# Cada ponto agora representa uma sessao, e nao uma vocalizacao.
###############################################################################

descritores_top_sessao <- top_sessao$descritor[
  seq_len(
    min(
      6,
      nrow(top_sessao)
    )
  )
]

dados_sessao_long <- dados_sessao_descritores %>%
  dplyr::select(
    unidade,
    sala,
    dplyr::all_of(
      descritores_top_sessao
    )
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(
      descritores_top_sessao
    ),
    names_to = "descritor",
    values_to = "valor"
  )

dados_sessao_long$descritor <- factor(
  dados_sessao_long$descritor,
  levels = descritores_top_sessao
)

set.seed(123)

figura_25 <- ggplot2::ggplot(
  dados_sessao_long,
  ggplot2::aes(
    x = sala,
    y = valor,
    fill = sala
  )
) +
  ggplot2::geom_boxplot(
    alpha = 0.60,
    outlier.shape = NA,
    linewidth = 0.55
  ) +
  ggplot2::geom_point(
    position = ggplot2::position_jitter(
      width = 0.07,
      height = 0
    ),
    size = 2.8,
    alpha = 0.85
  ) +
  ggplot2::facet_wrap(
    ~ descritor,
    scales = "free_y",
    ncol = 3
  ) +
  ggplot2::scale_fill_manual(
    values = cores_sala
  ) +
  ggplot2::scale_x_discrete(
    labels = rotulos_sala
  ) +
  ggplot2::guides(
    fill = "none"
  ) +
  ggplot2::labs(
    title = "Descritores acústicos entre os tratamentos",
    subtitle = "Cada ponto representa a mediana de uma sessão de gravação",
    x = "Tratamento",
    y = "Mediana da sessão"
  )

ggplot2::ggsave(
  file.path(
    pasta_figuras,
    "Figura_25_boxplots_descritores_por_sessao.png"
  ),
  figura_25,
  width = 13,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)

###############################################################################
# 7. FIGURA 26 — PERFIL DAS SALAS EM RELACAO AO CONTROLE
###############################################################################

# Padronizar os descritores no nivel das sessoes

matriz_sessao_z <- scale(
  dados_sessao_descritores[
    ,
    descritores,
    drop = FALSE
  ]
)

sessao_z <- as.data.frame(
  matriz_sessao_z
)

sessao_z$sala <- dados_sessao_descritores$sala

perfil_sessao_sala <- sessao_z %>%
  dplyr::group_by(sala) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(descritores),
      \(x) mean(x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(
    cols = -sala,
    names_to = "descritor",
    values_to = "media_z"
  )

perfil_controle <- perfil_sessao_sala %>%
  dplyr::filter(
    sala == "Sala1"
  ) %>%
  dplyr::select(
    descritor,
    controle = media_z
  )

diferenca_controle_descritores <- perfil_sessao_sala %>%
  dplyr::left_join(
    perfil_controle,
    by = "descritor"
  ) %>%
  dplyr::mutate(
    diferenca_controle =
      media_z - controle
  ) %>%
  dplyr::filter(
    sala != "Sala1"
  )

# Ordenacao conforme efeito por sessao

diferenca_controle_descritores$descritor <- factor(
  diferenca_controle_descritores$descritor,
  levels = rev(
    resultados_descritores_sessao$descritor
  )
)

figura_26 <- ggplot2::ggplot(
  diferenca_controle_descritores,
  ggplot2::aes(
    x = sala,
    y = descritor,
    fill = diferenca_controle
  )
) +
  ggplot2::geom_tile(
    color = "white",
    linewidth = 0.6
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = sprintf(
        "%.2f",
        diferenca_controle
      )
    ),
    size = 3.5
  ) +
  ggplot2::scale_x_discrete(
    labels = rotulos_sala
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#4575B4",
    mid = "white",
    high = "#D73027",
    midpoint = 0
  ) +
  ggplot2::labs(
    title = "Diferenças acústicas em relação ao tratamento controle",
    subtitle = "Sala 1 utilizada como referência; valores calculados a partir das sessões",
    x = "Tratamento",
    y = "Descritor acústico",
    fill = "Diferença\npadronizada"
  ) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank()
  )

ggplot2::ggsave(
  file.path(
    pasta_figuras,
    "Figura_26_descritores_diferenca_controle_sessao.png"
  ),
  figura_26,
  width = 10,
  height = 9,
  units = "in",
  dpi = 300,
  bg = "white"
)

###############################################################################
# 8. RESUMO AUTOMATICO
###############################################################################

n_sig_vocalizacao <- sum(
  resultados_descritores_sala$p_ajustado_BH < 0.05,
  na.rm = TRUE
)

n_sig_sessao <- sum(
  resultados_descritores_sessao$p_ajustado_BH < 0.05,
  na.rm = TRUE
)

descritores_persistentes <- comparacao_niveis %>%
  dplyr::filter(
    permanece_significativo
  ) %>%
  dplyr::pull(
    descritor
  )

writeLines(
  c(
    paste(
      "Numero de vocalizacoes:",
      nrow(dados_acusticos)
    ),
    paste(
      "Numero de sessoes:",
      nrow(dados_sessao_descritores)
    ),
    "",
    paste(
      "Descritores significativos na analise por vocalizacao:",
      n_sig_vocalizacao
    ),
    paste(
      "Descritores significativos na analise por sessao:",
      n_sig_sessao
    ),
    "",
    paste(
      "Descritores significativos nos dois niveis:",
      ifelse(
        length(descritores_persistentes) == 0,
        "Nenhum",
        paste(
          descritores_persistentes,
          collapse = ", "
        )
      )
    ),
    "",
    paste(
      "A analise por sessao reduz a pseudorreplicacao",
      "associada ao uso de cada vocalizacao como observacao independente."
    ),
    paste(
      "Entretanto, cada tratamento continua representado",
      "por uma unica sala experimental."
    ),
    paste(
      "Os resultados devem permanecer interpretados como exploratorios."
    )
  ),
  file.path(
    pasta_saida,
    "Resumo_validacao_descritores_por_sessao.txt"
  )
)

message(
  "\nVALIDACAO DOS DESCRITORES POR SESSAO CONCLUIDA.\n",
  "Resultados em:\n",
  pasta_saida
)

###############################################################################
# ANALISE DOS TIPOS MANUAIS ENTRE OS TRATAMENTOS
###############################################################################

###############################################################################
# 1. COMPOSICAO DOS TIPOS MANUAIS POR SESSAO
###############################################################################

tipos_sessao <- dados_acusticos %>%
  dplyr::count(
    unidade,
    sala,
    ambiente,
    data,
    tipo,
    name = "n"
  ) %>%
  tidyr::complete(
    unidade,
    tipo,
    fill = list(n = 0)
  ) %>%
  dplyr::group_by(unidade) %>%
  dplyr::mutate(
    total = sum(n),
    proporcao = dplyr::if_else(
      total > 0,
      n / total,
      0
    )
  ) %>%
  dplyr::ungroup()

# Recuperar metadados que podem ter sido perdidos pelo complete()
metadados_sessao <- dados_acusticos %>%
  dplyr::distinct(
    unidade,
    sala,
    ambiente,
    data
  )

tipos_sessao <- tipos_sessao %>%
  dplyr::select(
    unidade,
    tipo,
    n,
    total,
    proporcao
  ) %>%
  dplyr::left_join(
    metadados_sessao,
    by = "unidade"
  ) %>%
  dplyr::mutate(
    sala = factor(
      sala,
      levels = c(
        "Sala1",
        "Sala2",
        "Sala3",
        "Sala4"
      )
    ),
    tipo = factor(
      tipo,
      levels = c(
        "tipo1",
        "tipo2",
        "tipo3"
      )
    )
  )

utils::write.csv(
  tipos_sessao,
  file.path(
    pasta_saida,
    "Tipos_manuais_composicao_por_sessao.csv"
  ),
  row.names = FALSE
)

###############################################################################
# 2. KRUSKAL-WALLIS — FREQUENCIA RELATIVA DOS TIPOS ENTRE SALAS
###############################################################################

resultado_tipos_frequencia <- dplyr::bind_rows(
  lapply(
    levels(tipos_sessao$tipo),
    function(tipo_atual) {
      
      dados_tipo <- tipos_sessao %>%
        dplyr::filter(
          tipo == tipo_atual
        )
      
      teste <- stats::kruskal.test(
        proporcao ~ sala,
        data = dados_tipo
      )
      
      n <- nrow(dados_tipo)
      k <- nlevels(
        droplevels(
          dados_tipo$sala
        )
      )
      
      epsilon2 <- max(
        0,
        (
          unname(teste$statistic) - k + 1
        ) /
          (n - k)
      )
      
      data.frame(
        tipo = tipo_atual,
        H = unname(teste$statistic),
        gl = unname(teste$parameter),
        p = teste$p.value,
        epsilon2 = epsilon2
      )
    }
  )
) %>%
  dplyr::mutate(
    p_ajustado_BH = stats::p.adjust(
      p,
      method = "BH"
    ),
    p_formatado = formatar_p(
      p_ajustado_BH
    ),
    resultado = dplyr::case_when(
      p_ajustado_BH < 0.05 ~ "Significativo",
      p_ajustado_BH < 0.10 ~ "Evidencia limitada",
      TRUE ~ "Nao significativo"
    )
  )

utils::write.csv(
  resultado_tipos_frequencia,
  file.path(
    pasta_saida,
    "Resultados_Tipos_Manuais_Frequencia_por_Sala.csv"
  ),
  row.names = FALSE
)

print(
  resultado_tipos_frequencia
)

###############################################################################
# 3. FIGURA 27 — COMPOSICAO MEDIA DOS TIPOS MANUAIS
###############################################################################

perfil_tipos_sala <- tipos_sessao %>%
  dplyr::group_by(
    sala,
    tipo
  ) %>%
  dplyr::summarise(
    media = mean(
      proporcao,
      na.rm = TRUE
    ),
    .groups = "drop"
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

figura_27 <- ggplot2::ggplot(
  perfil_tipos_sala,
  ggplot2::aes(
    x = sala,
    y = media,
    group = tipo,
    color = tipo
  )
) +
  ggplot2::geom_line(
    linewidth = 1
  ) +
  ggplot2::geom_point(
    size = 3.5
  ) +
  ggplot2::scale_color_manual(
    values = cores_tipo,
    labels = rotulos_tipo
  ) +
  ggplot2::scale_x_discrete(
    labels = rotulos_sala
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format()
  ) +
  ggplot2::labs(
    title = "Composição dos tipos vocais manuais entre os tratamentos",
    subtitle = "Proporção média observada nas sessões de cada sala",
    x = "Tratamento",
    y = "Proporção média",
    color = "Tipo vocal"
  )

ggplot2::ggsave(
  file.path(
    pasta_figuras,
    "Figura_27_tipos_manuais_composicao_salas.png"
  ),
  figura_27,
  width = 11,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

###############################################################################
# 4. FIGURA 28 — DISTRIBUICAO DOS TIPOS POR SESSAO
###############################################################################

figura_28 <- ggplot2::ggplot(
  tipos_sessao,
  ggplot2::aes(
    x = sala,
    y = proporcao,
    fill = sala
  )
) +
  ggplot2::geom_boxplot(
    alpha = 0.65,
    outlier.shape = NA
  ) +
  ggplot2::geom_point(
    position = ggplot2::position_jitter(
      width = 0.07,
      height = 0
    ),
    size = 2.8
  ) +
  ggplot2::facet_wrap(
    ~ tipo,
    nrow = 1,
    labeller = ggplot2::labeller(
      tipo = rotulos_tipo
    )
  ) +
  ggplot2::scale_fill_manual(
    values = cores_sala
  ) +
  ggplot2::scale_x_discrete(
    labels = rotulos_sala
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format()
  ) +
  ggplot2::guides(
    fill = "none"
  ) +
  ggplot2::labs(
    title = "Frequência relativa dos tipos vocais entre os tratamentos",
    subtitle = "Cada ponto representa uma sessão de gravação",
    x = "Tratamento",
    y = "Proporção do tipo na sessão"
  )

ggplot2::ggsave(
  file.path(
    pasta_figuras,
    "Figura_28_tipos_manuais_por_sessao.png"
  ),
  figura_28,
  width = 13,
  height = 6,
  units = "in",
  dpi = 300,
  bg = "white"
)

###############################################################################
# ESTRUTURA ACUSTICA DENTRO DE CADA TIPO MANUAL — CORRIGIDO
###############################################################################

resultados_descritores_por_tipo <- list()
dados_descritores_tipo_sessao <- list()

for (tipo_atual in levels(dados_acusticos$tipo)) {
  
  dados_tipo <- dados_acusticos %>%
    dplyr::filter(tipo == tipo_atual)
  
  resumo_tipo_sessao <- dados_tipo %>%
    dplyr::group_by(
      unidade,
      sala,
      ambiente,
      data
    ) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(descritores),
        \(x) stats::median(x, na.rm = TRUE)
      ),
      n_vocalizacoes_tipo = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      tipo = tipo_atual,
      sala = factor(
        sala,
        levels = c("Sala1", "Sala2", "Sala3", "Sala4")
      )
    )
  
  dados_descritores_tipo_sessao[[tipo_atual]] <- resumo_tipo_sessao
  
  resultado_tipo <- dplyr::bind_rows(
    lapply(
      descritores,
      function(descritor_atual) {
        
        # CORRECAO PRINCIPAL
        valores <- resumo_tipo_sessao[[descritor_atual]]
        
        # Remove valores nao finitos para avaliacao
        valores_validos <- valores[is.finite(valores)]
        
        # Nao testar se houver pouca informacao ou ausencia de variacao
        if (
          length(valores_validos) < 4 ||
          length(unique(valores_validos)) < 2
        ) {
          
          return(
            data.frame(
              tipo = tipo_atual,
              descritor = descritor_atual,
              H = NA_real_,
              gl = NA_real_,
              p = NA_real_,
              epsilon2 = NA_real_,
              n_sessoes = nrow(resumo_tipo_sessao)
            )
          )
        }
        
        dados_teste <- resumo_tipo_sessao %>%
          dplyr::filter(
            is.finite(.data[[descritor_atual]])
          ) %>%
          droplevels()
        
        # Precisamos de pelo menos duas salas representadas
        if (dplyr::n_distinct(dados_teste$sala) < 2) {
          
          return(
            data.frame(
              tipo = tipo_atual,
              descritor = descritor_atual,
              H = NA_real_,
              gl = NA_real_,
              p = NA_real_,
              epsilon2 = NA_real_,
              n_sessoes = nrow(dados_teste)
            )
          )
        }
        
        teste <- tryCatch(
          stats::kruskal.test(
            x = dados_teste[[descritor_atual]],
            g = dados_teste$sala
          ),
          error = function(e) NULL
        )
        
        if (is.null(teste)) {
          
          return(
            data.frame(
              tipo = tipo_atual,
              descritor = descritor_atual,
              H = NA_real_,
              gl = NA_real_,
              p = NA_real_,
              epsilon2 = NA_real_,
              n_sessoes = nrow(dados_teste)
            )
          )
        }
        
        n <- nrow(dados_teste)
        k <- dplyr::n_distinct(dados_teste$sala)
        
        epsilon2 <- if (n > k) {
          max(
            0,
            (
              unname(teste$statistic) - k + 1
            ) /
              (n - k)
          )
        } else {
          NA_real_
        }
        
        data.frame(
          tipo = tipo_atual,
          descritor = descritor_atual,
          H = unname(teste$statistic),
          gl = unname(teste$parameter),
          p = teste$p.value,
          epsilon2 = epsilon2,
          n_sessoes = n
        )
      }
    )
  )
  
  resultados_descritores_por_tipo[[tipo_atual]] <- resultado_tipo
}

###############################################################################
# JUNTAR RESULTADOS
###############################################################################

dados_descritores_tipo_sessao <- dplyr::bind_rows(
  dados_descritores_tipo_sessao
)

resultados_descritores_por_tipo <- dplyr::bind_rows(
  resultados_descritores_por_tipo
) %>%
  dplyr::group_by(tipo) %>%
  dplyr::mutate(
    p_ajustado_BH = stats::p.adjust(
      p,
      method = "BH"
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    p_formatado = formatar_p(p_ajustado_BH),
    resultado = dplyr::case_when(
      is.na(p_ajustado_BH) ~ "Nao testado",
      p_ajustado_BH < 0.05 ~ "Significativo",
      p_ajustado_BH < 0.10 ~ "Evidencia limitada",
      TRUE ~ "Nao significativo"
    )
  ) %>%
  dplyr::arrange(
    tipo,
    dplyr::desc(epsilon2)
  )

print(resultados_descritores_por_tipo)

utils::write.csv(
  dados_descritores_tipo_sessao,
  file.path(
    pasta_saida,
    "Descritores_Tipos_Manuais_por_Sessao.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  resultados_descritores_por_tipo,
  file.path(
    pasta_saida,
    "Resultados_Descritores_Tipos_Manuais_entre_Salas.csv"
  ),
  row.names = FALSE
)

###############################################################################
# FIGURA 29 — MAIORES EFEITOS DENTRO DOS TIPOS MANUAIS
###############################################################################

top_tipo_descritores <- resultados_descritores_por_tipo %>%
  dplyr::filter(
    !is.na(epsilon2)
  ) %>%
  dplyr::group_by(tipo) %>%
  dplyr::slice_max(
    order_by = epsilon2,
    n = 5,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    descritor_tipo = paste(
      descritor,
      tipo,
      sep = "___"
    )
  )

ordem_descritores <- top_tipo_descritores %>%
  dplyr::arrange(epsilon2) %>%
  dplyr::pull(descritor_tipo)

top_tipo_descritores$descritor_tipo <- factor(
  top_tipo_descritores$descritor_tipo,
  levels = ordem_descritores
)

figura_29 <- ggplot2::ggplot(
  top_tipo_descritores,
  ggplot2::aes(
    x = epsilon2,
    y = descritor_tipo,
    fill = tipo
  )
) +
  ggplot2::geom_col(
    width = 0.70
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = paste0(
        "p = ",
        p_formatado
      )
    ),
    hjust = -0.05,
    size = 3.7
  ) +
  ggplot2::facet_wrap(
    ~ tipo,
    scales = "free_y",
    ncol = 1,
    labeller = ggplot2::labeller(
      tipo = rotulos_tipo
    )
  ) +
  ggplot2::scale_fill_manual(
    values = cores_tipo,
    labels = rotulos_tipo
  ) +
  ggplot2::scale_y_discrete(
    labels = function(x) {
      sub("___.*$", "", x)
    }
  ) +
  ggplot2::scale_x_continuous(
    expand = ggplot2::expansion(
      mult = c(0, 0.18)
    )
  ) +
  ggplot2::labs(
    title = "Variação da estrutura acústica dentro dos tipos vocais",
    subtitle = "Maiores tamanhos de efeito entre as salas; cada sessão é uma unidade analítica",
    x = "Tamanho de efeito (epsilon quadrado)",
    y = "Descritor acústico",
    fill = "Tipo vocal"
  ) +
  ggplot2::guides(
    fill = "none"
  )

ggplot2::ggsave(
  file.path(
    pasta_figuras,
    "Figura_29_descritores_tipos_manuais_efeitos.png"
  ),
  figura_29,
  width = 11,
  height = 10,
  units = "in",
  dpi = 300,
  bg = "white"
)

###############################################################################
# RESUMO
###############################################################################

n_descritores_tipo_sig <- sum(
  resultados_descritores_por_tipo$p_ajustado_BH < 0.05,
  na.rm = TRUE
)

# ============================================================
# POS-TESTE DO TIPO 1 ENTRE AS SALAS
# Dunn post-hoc + correcao BH
# Unidade analitica = sessao de gravacao
# ============================================================


# ------------------------------------------------------------
# 1. Conferir se dados_acusticos existe
# ------------------------------------------------------------

if (!exists("dados_acusticos")) {
  
  arquivo_dados_acusticos <- file.path(
    raiz,
    "10_Pipeline_RF_PAM_FINAL",
    "01_Pre_processamento",
    "dados_analise_final.rds"
  )
  
  if (!file.exists(arquivo_dados_acusticos)) {
    stop(
      "Nao foi possivel encontrar dados_acusticos nem dados_analise_final.rds.",
      call. = FALSE
    )
  }
  
  dados_acusticos <- readRDS(
    arquivo_dados_acusticos
  )
}


# ------------------------------------------------------------
# 2. Garantir fatores corretos
# ------------------------------------------------------------

dados_acusticos <- dados_acusticos %>%
  dplyr::mutate(
    
    sala = factor(
      sala,
      levels = c(
        "Sala1",
        "Sala2",
        "Sala3",
        "Sala4"
      )
    ),
    
    ambiente = factor(
      ambiente,
      levels = c(
        "Externo",
        "Interno"
      )
    ),
    
    tipo = factor(
      tipo,
      levels = c(
        "tipo1",
        "tipo2",
        "tipo3"
      )
    )
  )


# ------------------------------------------------------------
# 3. Calcular proporcao do Tipo 1 EM CADA SESSAO
#
# Importante:
# uma sessao sem Tipo 1 recebe proporcao = 0.
# ------------------------------------------------------------

dados_tipo1_posthoc <- dados_acusticos %>%
  
  dplyr::group_by(
    unidade,
    sala,
    ambiente,
    data
  ) %>%
  
  dplyr::summarise(
    
    n_total = dplyr::n(),
    
    n_tipo1 = sum(
      tipo == "tipo1",
      na.rm = TRUE
    ),
    
    proporcao = n_tipo1 / n_total,
    
    .groups = "drop"
  ) %>%
  
  dplyr::mutate(
    
    sala = factor(
      sala,
      levels = c(
        "Sala1",
        "Sala2",
        "Sala3",
        "Sala4"
      )
    )
  )


# ------------------------------------------------------------
# 4. Conferencia
# ------------------------------------------------------------

print(dados_tipo1_posthoc)

table(dados_tipo1_posthoc$sala)


# ------------------------------------------------------------
# 5. Resumo por sala
# ------------------------------------------------------------

resumo_tipo1_posthoc <- dados_tipo1_posthoc %>%
  
  dplyr::group_by(sala) %>%
  
  dplyr::summarise(
    
    n_sessoes = dplyr::n(),
    
    media = mean(
      proporcao,
      na.rm = TRUE
    ),
    
    mediana = stats::median(
      proporcao,
      na.rm = TRUE
    ),
    
    Q1 = stats::quantile(
      proporcao,
      probs = 0.25,
      na.rm = TRUE
    ),
    
    Q3 = stats::quantile(
      proporcao,
      probs = 0.75,
      na.rm = TRUE
    ),
    
    minimo = min(
      proporcao,
      na.rm = TRUE
    ),
    
    maximo = max(
      proporcao,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

print(resumo_tipo1_posthoc)


# ------------------------------------------------------------
# 6. Repetir Kruskal-Wallis global do Tipo 1
# ------------------------------------------------------------

teste_kw_tipo1 <- stats::kruskal.test(
  proporcao ~ sala,
  data = dados_tipo1_posthoc
)

print(teste_kw_tipo1)


# ------------------------------------------------------------
# 7. Pos-teste de Dunn
# ------------------------------------------------------------

if (!requireNamespace("FSA", quietly = TRUE)) {
  install.packages("FSA")
}

teste_dunn_tipo1 <- FSA::dunnTest(
  proporcao ~ sala,
  data = dados_tipo1_posthoc,
  method = "bh"
)


resultado_dunn_tipo1 <- teste_dunn_tipo1$res %>%
  
  dplyr::transmute(
    
    comparacao = Comparison,
    
    Z = Z,
    
    p_nao_ajustado = P.unadj,
    
    p_ajustado_BH = P.adj,
    
    p_formatado = ifelse(
      P.adj < 0.001,
      "< 0,001",
      sub(
        "\\.",
        ",",
        sprintf("%.3f", P.adj)
      )
    ),
    
    resultado = dplyr::case_when(
      
      P.adj < 0.05 ~
        "Significativo",
      
      P.adj < 0.10 ~
        "Evidencia limitada",
      
      TRUE ~
        "Nao significativo"
    )
  ) %>%
  
  dplyr::arrange(
    p_ajustado_BH
  )


print(resultado_dunn_tipo1)


# ------------------------------------------------------------
# 8. Comparacoes significativas
# ------------------------------------------------------------

comparacoes_significativas_tipo1 <-
  resultado_dunn_tipo1 %>%
  dplyr::filter(
    p_ajustado_BH < 0.05
  )

print(comparacoes_significativas_tipo1)


# ------------------------------------------------------------
# 9. Salvar resultados
# ------------------------------------------------------------

utils::write.csv(
  
  resumo_tipo1_posthoc,
  
  file.path(
    pasta_saida,
    "Resumo_Tipo1_por_Sala.csv"
  ),
  
  row.names = FALSE
)


utils::write.csv(
  
  resultado_dunn_tipo1,
  
  file.path(
    pasta_saida,
    "Pos_teste_Dunn_Tipo1_entre_Salas.csv"
  ),
  
  row.names = FALSE
)

print(resumo_tipo1_posthoc)
print(teste_kw_tipo1)
print(resultado_dunn_tipo1)

# ------------------------------------------------------------
# 10. Final
# ------------------------------------------------------------

###############################################################################
# FIGURAS ADICIONAIS PARA A APRESENTACAO
# CONTINUACAO DO SCRIPT 04
###############################################################################

# Garantir fatores corretos
dados <- dados %>%
  dplyr::mutate(
    sala = factor(
      sala,
      levels = c("Sala1", "Sala2", "Sala3", "Sala4")
    ),
    tipo = factor(
      tipo,
      levels = c("tipo1", "tipo2", "tipo3")
    ),
    grupo_computacional = factor(
      grupo_computacional
    )
  )

cores_sala <- c(
  Sala1 = "#222222",
  Sala2 = "#E66101",
  Sala3 = "#5E3C99",
  Sala4 = "#1B9E77"
)

rotulos_sala <- c(
  Sala1 = "Sala 1",
  Sala2 = "Sala 2",
  Sala3 = "Sala 3",
  Sala4 = "Sala 4"
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
  "1" = "#D73027",
  "2" = "#4575B4",
  "3" = "#1A9850",
  "4" = "#984EA3",
  "5" = "#FF7F00",
  "6" = "#A65628",
  "7" = "#F781BF",
  "8" = "#666666"
)


###############################################################################
# FIGURA 30
# NUMERO TOTAL DE VOCALIZACOES POR TRATAMENTO
###############################################################################

quantidade_total_sala <- dados %>%
  dplyr::count(
    sala,
    name = "n_vocalizacoes"
  )

figura_30 <- ggplot2::ggplot(
  quantidade_total_sala,
  ggplot2::aes(
    x = sala,
    y = n_vocalizacoes,
    fill = sala
  )
) +
  ggplot2::geom_col(
    width = 0.70,
    color = "grey20"
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = n_vocalizacoes
    ),
    vjust = -0.45,
    size = 5,
    fontface = "bold"
  ) +
  ggplot2::scale_fill_manual(
    values = cores_sala
  ) +
  ggplot2::scale_x_discrete(
    labels = rotulos_sala
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(0, 0.12)
    )
  ) +
  ggplot2::guides(
    fill = "none"
  ) +
  ggplot2::labs(
    title = "Número de vocalizações registradas entre os tratamentos",
    subtitle = "Número total de vocalizações incluídas na análise acústica",
    x = "Tratamento",
    y = "Número de vocalizações"
  )

ggplot2::ggsave(
  file.path(
    pasta_figuras,
    "Figura_30_numero_vocalizacoes_tratamentos.png"
  ),
  figura_30,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)


###############################################################################
# FIGURA 31
# COMPOSICAO DOS TIPOS MANUAIS POR TRATAMENTO
###############################################################################

composicao_manual <- dados %>%
  dplyr::count(
    sala,
    tipo,
    name = "n"
  ) %>%
  dplyr::group_by(sala) %>%
  dplyr::mutate(
    proporcao = n / sum(n)
  ) %>%
  dplyr::ungroup()

figura_31 <- ggplot2::ggplot(
  composicao_manual,
  ggplot2::aes(
    x = sala,
    y = proporcao,
    fill = tipo
  )
) +
  ggplot2::geom_col(
    width = 0.72,
    color = "white",
    linewidth = 0.5
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = ifelse(
        proporcao >= 0.04,
        scales::percent(
          proporcao,
          accuracy = 0.1
        ),
        ""
      )
    ),
    position = ggplot2::position_stack(
      vjust = 0.5
    ),
    size = 4.5
  ) +
  ggplot2::scale_fill_manual(
    values = cores_tipo,
    labels = rotulos_tipo
  ) +
  ggplot2::scale_x_discrete(
    labels = rotulos_sala
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(
      accuracy = 1
    ),
    expand = c(0, 0)
  ) +
  ggplot2::labs(
    title = "Composição dos tipos vocais manuais entre os tratamentos",
    subtitle = "Proporção das vocalizações classificadas visual e auditivamente",
    x = "Tratamento",
    y = "Proporção das vocalizações",
    fill = "Tipo vocal"
  )

ggplot2::ggsave(
  file.path(
    pasta_figuras,
    "Figura_31_composicao_tipos_manuais_tratamentos.png"
  ),
  figura_31,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)


###############################################################################
# FIGURA 32
# COMPOSICAO DOS GRUPOS COMPUTACIONAIS POR TRATAMENTO
###############################################################################

composicao_computacional <- dados %>%
  dplyr::count(
    sala,
    grupo_computacional,
    name = "n"
  ) %>%
  dplyr::group_by(sala) %>%
  dplyr::mutate(
    proporcao = n / sum(n)
  ) %>%
  dplyr::ungroup()

figura_32 <- ggplot2::ggplot(
  composicao_computacional,
  ggplot2::aes(
    x = sala,
    y = proporcao,
    fill = grupo_computacional
  )
) +
  ggplot2::geom_col(
    width = 0.72,
    color = "white",
    linewidth = 0.5
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = ifelse(
        proporcao >= 0.055,
        scales::percent(
          proporcao,
          accuracy = 0.1
        ),
        ""
      )
    ),
    position = ggplot2::position_stack(
      vjust = 0.5
    ),
    size = 4
  ) +
  ggplot2::scale_fill_manual(
    values = cores_grupos[
      levels(dados$grupo_computacional)
    ]
  ) +
  ggplot2::scale_x_discrete(
    labels = rotulos_sala
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(
      accuracy = 1
    ),
    expand = c(0, 0)
  ) +
  ggplot2::labs(
    title = "Composição dos grupos acústicos computacionais entre os tratamentos",
    subtitle = "Classificação exploratória obtida por Random Forest não supervisionado + PAM",
    x = "Tratamento",
    y = "Proporção das vocalizações",
    fill = "Grupo acústico"
  )

ggplot2::ggsave(
  file.path(
    pasta_figuras,
    "Figura_32_composicao_grupos_computacionais_tratamentos.png"
  ),
  figura_32,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)


###############################################################################
# FIGURA 33
# PRINCIPAIS DESCRITORES ACUSTICOS ENTRE OS TRATAMENTOS
#
# Usa os descritores do Script 01 e resume por SESSAO,
# evitando tratar cada vocalizacao como replica independente.
###############################################################################

arquivo_descritores_figuras <- file.path(
  raiz,
  "10_Pipeline_RF_PAM_FINAL",
  "01_Pre_processamento",
  "dados_analise_final.rds"
)

dados_descritores_figuras <- readRDS(
  arquivo_descritores_figuras
)

# Mantem apenas metadados necessários + descritores
descritores_figuras <- c(
  "startdom",
  "peakt",
  "maxdom",
  "modindx",
  "sp.ent",
  "skew"
)

descritores_figuras <- intersect(
  descritores_figuras,
  names(dados_descritores_figuras)
)

dados_sessao_descritores <- dados_descritores_figuras %>%
  dplyr::mutate(
    sala = factor(
      sala,
      levels = c(
        "Sala1",
        "Sala2",
        "Sala3",
        "Sala4"
      )
    )
  ) %>%
  dplyr::group_by(
    unidade,
    sala,
    ambiente,
    data
  ) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(
        descritores_figuras
      ),
      \(x) stats::median(
        x,
        na.rm = TRUE
      )
    ),
    .groups = "drop"
  )

dados_long_descritores <- dados_sessao_descritores %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(
      descritores_figuras
    ),
    names_to = "descritor",
    values_to = "valor"
  ) %>%
  dplyr::filter(
    is.finite(valor)
  )

figura_33 <- ggplot2::ggplot(
  dados_long_descritores,
  ggplot2::aes(
    x = sala,
    y = valor,
    fill = sala
  )
) +
  ggplot2::geom_boxplot(
    width = 0.62,
    alpha = 0.76,
    outlier.shape = NA
  ) +
  ggplot2::geom_point(
    position = ggplot2::position_jitter(
      width = 0.07,
      height = 0
    ),
    size = 2.8,
    alpha = 0.85,
    color = "black"
  ) +
  ggplot2::facet_wrap(
    ~ descritor,
    scales = "free_y",
    ncol = 3
  ) +
  ggplot2::scale_fill_manual(
    values = cores_sala
  ) +
  ggplot2::scale_x_discrete(
    labels = rotulos_sala
  ) +
  ggplot2::guides(
    fill = "none"
  ) +
  ggplot2::labs(
    title = "Descritores acústicos com maior variação entre os tratamentos",
    subtitle = "Cada ponto representa a mediana de uma sessão de gravação",
    x = "Tratamento",
    y = "Mediana da sessão"
  )

ggplot2::ggsave(
  file.path(
    pasta_figuras,
    "Figura_33_descritores_acusticos_tratamentos.png"
  ),
  figura_33,
  width = 14,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)


###############################################################################
# FIGURA 34
# DIFERENCAS PADRONIZADAS DOS DESCRITORES EM RELACAO AO CONTROLE
###############################################################################

perfil_descritores_sala <- dados_sessao_descritores %>%
  dplyr::group_by(sala) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(
        descritores_figuras
      ),
      \(x) mean(
        x,
        na.rm = TRUE
      )
    ),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(
    cols = -sala,
    names_to = "descritor",
    values_to = "media"
  )

# Padroniza cada descritor entre as quatro salas
perfil_descritores_padronizado <-
  perfil_descritores_sala %>%
  dplyr::group_by(descritor) %>%
  dplyr::mutate(
    media_z = as.numeric(
      scale(media)
    )
  ) %>%
  dplyr::ungroup()

controle_descritores <-
  perfil_descritores_padronizado %>%
  dplyr::filter(
    sala == "Sala1"
  ) %>%
  dplyr::select(
    descritor,
    controle_z = media_z
  )

diferenca_descritores <-
  perfil_descritores_padronizado %>%
  dplyr::left_join(
    controle_descritores,
    by = "descritor"
  ) %>%
  dplyr::mutate(
    diferenca_z = media_z - controle_z
  ) %>%
  dplyr::filter(
    sala != "Sala1"
  )

figura_34 <- ggplot2::ggplot(
  diferenca_descritores,
  ggplot2::aes(
    x = sala,
    y = descritor,
    fill = diferenca_z
  )
) +
  ggplot2::geom_tile(
    color = "white",
    linewidth = 1
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = sprintf(
        "%.2f",
        diferenca_z
      )
    ),
    size = 4.5
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#4575B4",
    mid = "white",
    high = "#D73027",
    midpoint = 0
  ) +
  ggplot2::scale_x_discrete(
    labels = rotulos_sala
  ) +
  ggplot2::labs(
    title = "Mudanças dos descritores acústicos em relação ao tratamento controle",
    subtitle = "Sala 1 utilizada como referência; diferenças padronizadas entre as médias das sessões",
    x = "Tratamento",
    y = "Descritor acústico",
    fill = "Diferença\npadronizada"
  )

ggplot2::ggsave(
  file.path(
    pasta_figuras,
    "Figura_34_descritores_diferenca_controle.png"
  ),
  figura_34,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

###############################################################################
# FIGURA 35
# ATIVIDADE VOCAL POR SESSAO
###############################################################################

atividade_sessao <- dados %>%
  dplyr::count(
    unidade,
    sala,
    ambiente,
    data,
    name = "n_vocalizacoes"
  ) %>%
  dplyr::mutate(
    sala = factor(
      sala,
      levels = c("Sala1", "Sala2", "Sala3", "Sala4")
    )
  )

figura_35 <- ggplot2::ggplot(
  atividade_sessao,
  ggplot2::aes(
    x = sala,
    y = n_vocalizacoes,
    fill = sala
  )
) +
  ggplot2::geom_boxplot(
    width = 0.60,
    alpha = 0.75,
    outlier.shape = NA
  ) +
  ggplot2::geom_point(
    position = ggplot2::position_jitter(
      width = 0.07
    ),
    size = 3.2,
    color = "black"
  ) +
  ggplot2::scale_fill_manual(
    values = cores_sala
  ) +
  ggplot2::scale_x_discrete(
    labels = rotulos_sala
  ) +
  ggplot2::guides(
    fill = "none"
  ) +
  ggplot2::labs(
    title = "Atividade vocal entre os tratamentos",
    subtitle = "Cada ponto representa uma sessão de gravação",
    x = "Tratamento",
    y = "Número de vocalizações por sessão"
  )

ggplot2::ggsave(
  file.path(
    pasta_figuras,
    "Figura_35_atividade_vocal_por_sessao.png"
  ),
  figura_35,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

###############################################################################
# FIGURA 36
# PERFIL ACUSTICO DOS TRATAMENTOS
###############################################################################

descritores_perfil <- c(
  "skew",
  "maxdom",
  "modindx",
  "sp.ent",
  "freq.Q25",
  "meandom"
)

descritores_perfil <- intersect(
  descritores_perfil,
  names(dados_descritores_figuras)
)

perfil_acustico <- dados_descritores_figuras %>%
  dplyr::select(
    sala,
    dplyr::all_of(descritores_perfil)
  ) %>%
  dplyr::group_by(sala) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(descritores_perfil),
      \(x) mean(x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(
    -sala,
    names_to = "descritor",
    values_to = "media"
  ) %>%
  dplyr::group_by(descritor) %>%
  dplyr::mutate(
    z = as.numeric(scale(media))
  ) %>%
  dplyr::ungroup()

figura_36 <- ggplot2::ggplot(
  perfil_acustico,
  ggplot2::aes(
    x = descritor,
    y = z,
    group = sala,
    color = sala
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = 2,
    color = "grey60"
  ) +
  ggplot2::geom_line(
    linewidth = 1
  ) +
  ggplot2::geom_point(
    size = 3
  ) +
  ggplot2::scale_color_manual(
    values = cores_sala,
    labels = rotulos_sala
  ) +
  ggplot2::labs(
    title = "Perfil acústico dos tratamentos",
    subtitle = "Médias padronizadas dos descritores com maior variação",
    x = "Descritor acústico",
    y = "Média padronizada",
    color = "Tratamento"
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 35,
      hjust = 1
    )
  )

ggplot2::ggsave(
  file.path(
    pasta_figuras,
    "Figura_36_perfil_acustico_tratamentos.png"
  ),
  figura_36,
  width = 11,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

###############################################################################
# FIGURA 37
###############################################################################
# FIGURA 38
# SINTese DOS DESCRITORES QUE MAIS RESPONDEM AOS TRATAMENTOS
###############################################################################

# Usa o objeto resultados_descritores_sala se ele já existir.
# Caso não exista, calcula novamente usando as sessoes.

if (!exists("resultados_descritores_sala")) {
  
  resultados_descritores_sala <- dplyr::bind_rows(
    lapply(
      descritores_figuras,
      function(d) {
        
        dados_teste <- dados_sessao_descritores %>%
          dplyr::select(
            sala,
            dplyr::all_of(d)
          ) %>%
          dplyr::filter(
            is.finite(.data[[d]])
          )
        
        teste <- stats::kruskal.test(
          dados_teste[[d]],
          dados_teste$sala
        )
        
        n <- nrow(dados_teste)
        k <- dplyr::n_distinct(
          dados_teste$sala
        )
        
        epsilon2 <- max(
          0,
          (
            unname(teste$statistic) -
              k +
              1
          ) /
            (
              n -
                k
            )
        )
        
        data.frame(
          descritor = d,
          H = unname(
            teste$statistic
          ),
          p = teste$p.value,
          epsilon2 = epsilon2
        )
      }
    )
  ) %>%
    dplyr::mutate(
      p_BH = stats::p.adjust(
        p,
        method = "BH"
      ),
      classe = dplyr::case_when(
        p_BH < 0.05 ~ "Significativo",
        p_BH < 0.10 ~ "Evidência limitada",
        TRUE ~ "Não significativo"
      )
    )
}

top_descritores_tratamento <-
  resultados_descritores_sala %>%
  dplyr::arrange(
    dplyr::desc(
      epsilon2
    )
  ) %>%
  dplyr::slice_head(
    n = min(
      8,
      nrow(
        resultados_descritores_sala
      )
    )
  )

figura_38 <- ggplot2::ggplot(
  top_descritores_tratamento %>%
    dplyr::mutate(
      descritor = factor(
        descritor,
        levels = rev(
          descritor
        )
      )
    ),
  ggplot2::aes(
    x = epsilon2,
    y = descritor,
    fill = classe
  )
) +
  ggplot2::geom_col(
    width = 0.70
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = paste0(
        "p ajust. = ",
        ifelse(
          p_BH < 0.001,
          "<0,001",
          sub(
            "\\.",
            ",",
            sprintf(
              "%.3f",
              p_BH
            )
          )
        )
      )
    ),
    hjust = -0.04,
    size = 4
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Significativo" = "#1B9E77",
      "Evidência limitada" = "#E66101",
      "Não significativo" = "#777777"
    )
  ) +
  ggplot2::scale_x_continuous(
    expand = ggplot2::expansion(
      mult = c(
        0,
        0.30
      )
    )
  ) +
  ggplot2::labs(
    title = "Descritores acústicos com maior resposta aos tratamentos",
    subtitle = "Tamanho de efeito calculado entre sessões de gravação",
    x = "Tamanho de efeito (epsilon quadrado)",
    y = "Descritor acústico",
    fill = "Resultado"
  )

ggplot2::ggsave(
  file.path(
    pasta_figuras,
    "Figura_38_descritores_maior_resposta_tratamentos.png"
  ),
  figura_38,
  width = 11,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

###############################################################################
# FIGURA 39
# DIFERENCAS ACUSTICAS DENTRO DOS TIPOS MANUAIS
###############################################################################

figura_tipo_manual <- resultados_descritores_por_tipo %>%
  dplyr::filter(
    !is.na(epsilon2),
    epsilon2 > 0
  ) %>%
  dplyr::group_by(tipo) %>%
  dplyr::slice_max(
    epsilon2,
    n = 5,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    descritor = factor(
      descritor,
      levels = unique(
        descritor[
          order(epsilon2)
        ]
      )
    )
  )

figura_39 <- ggplot2::ggplot(
  figura_tipo_manual,
  ggplot2::aes(
    x = epsilon2,
    y = descritor,
    fill = tipo
  )
) +
  ggplot2::geom_col(
    width = 0.68
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = paste0(
        "p = ",
        p_formatado
      )
    ),
    hjust = -0.05,
    size = 4
  ) +
  ggplot2::facet_wrap(
    ~ tipo,
    ncol = 1,
    scales = "free_y",
    labeller = ggplot2::labeller(
      tipo = rotulos_tipo
    )
  ) +
  ggplot2::scale_fill_manual(
    values = cores_tipo
  ) +
  ggplot2::scale_x_continuous(
    expand = ggplot2::expansion(
      mult = c(0, 0.20)
    )
  ) +
  ggplot2::guides(
    fill = "none"
  ) +
  ggplot2::labs(
    title = "Variação da estrutura acústica dos tipos vocais manuais",
    subtitle = "Descritores com maiores tamanhos de efeito entre os tratamentos",
    x = "Tamanho de efeito (epsilon quadrado)",
    y = "Descritor acústico"
  )

ggplot2::ggsave(
  file.path(
    pasta_figuras,
    "Figura_39_diferencas_acusticas_tipos_manuais.png"
  ),
  figura_39,
  width = 10,
  height = 9,
  units = "in",
  dpi = 300,
  bg = "white"
)

###############################################################################
# CORRECAO DA FIGURA 38
# ESTRUTURA ACUSTICA GERAL ENTRE OS TRATAMENTOS
###############################################################################

# ------------------------------------------------------------
# 1. Recalcular se o objeto nao existir
# ------------------------------------------------------------

if (!exists("resultados_descritores_sala")) {
  
  resultados_descritores_sala <- dplyr::bind_rows(
    lapply(
      descritores_figuras,
      function(d) {
        
        dados_teste <- dados_sessao_descritores %>%
          dplyr::select(
            sala,
            dplyr::all_of(d)
          ) %>%
          dplyr::filter(
            is.finite(.data[[d]])
          )
        
        teste <- stats::kruskal.test(
          dados_teste[[d]],
          dados_teste$sala
        )
        
        n <- nrow(dados_teste)
        
        k <- dplyr::n_distinct(
          dados_teste$sala
        )
        
        epsilon2 <- if (n > k) {
          
          max(
            0,
            (
              unname(teste$statistic) -
                k +
                1
            ) /
              (
                n -
                  k
              )
          )
          
        } else {
          
          NA_real_
          
        }
        
        data.frame(
          descritor = d,
          H = unname(teste$statistic),
          p = teste$p.value,
          epsilon2 = epsilon2
        )
      }
    )
  )
}


# ------------------------------------------------------------
# 2. GARANTIR que p_BH e classe existam
#
# Esta parte fica FORA do if()
# ------------------------------------------------------------

resultados_descritores_sala <-
  resultados_descritores_sala %>%
  dplyr::mutate(
    
    p_BH = stats::p.adjust(
      p,
      method = "BH"
    ),
    
    classe = dplyr::case_when(
      
      is.na(p_BH) ~ "Nao testado",
      
      p_BH < 0.05 ~
        "Significativo",
      
      p_BH < 0.10 ~
        "Evidencia limitada",
      
      TRUE ~
        "Nao significativo"
    ),
    
    p_formatado = dplyr::case_when(
      
      is.na(p_BH) ~
        "NA",
      
      p_BH < 0.001 ~
        "< 0,001",
      
      TRUE ~
        sub(
          "\\.",
          ",",
          sprintf("%.3f", p_BH)
        )
    )
  )


# Conferencia
print(resultados_descritores_sala)

names(resultados_descritores_sala)


# ------------------------------------------------------------
# 3. Selecionar descritores com maiores efeitos
# ------------------------------------------------------------

top_descritores_tratamento <-
  resultados_descritores_sala %>%
  dplyr::filter(
    !is.na(epsilon2)
  ) %>%
  dplyr::arrange(
    dplyr::desc(epsilon2)
  ) %>%
  dplyr::slice_head(
    n = min(
      8,
      nrow(resultados_descritores_sala)
    )
  )


# ------------------------------------------------------------
# 4. Ordenacao correta para o grafico
# ------------------------------------------------------------

top_descritores_tratamento <-
  top_descritores_tratamento %>%
  dplyr::mutate(
    descritor = factor(
      descritor,
      levels = rev(descritor)
    )
  )


# ------------------------------------------------------------
# 5. FIGURA 38
# ------------------------------------------------------------

figura_38 <- ggplot2::ggplot(
  top_descritores_tratamento,
  ggplot2::aes(
    x = epsilon2,
    y = descritor,
    fill = classe
  )
) +
  
  ggplot2::geom_col(
    width = 0.70
  ) +
  
  ggplot2::geom_text(
    ggplot2::aes(
      label = paste0(
        "p ajust. = ",
        p_formatado
      )
    ),
    hjust = -0.04,
    size = 4
  ) +
  
  ggplot2::scale_fill_manual(
    values = c(
      "Significativo" = "#1B9E77",
      "Evidencia limitada" = "#E66101",
      "Nao significativo" = "#777777",
      "Nao testado" = "#BBBBBB"
    )
  ) +
  
  ggplot2::scale_x_continuous(
    expand = ggplot2::expansion(
      mult = c(
        0,
        0.35
      )
    )
  ) +
  
  ggplot2::labs(
    title = "Variação da estrutura acústica geral entre os tratamentos",
    subtitle = "Descritores das vocalizações resumidos por sessão de gravação",
    x = "Tamanho de efeito (epsilon quadrado)",
    y = "Descritor acústico",
    fill = "Resultado"
  )


# ------------------------------------------------------------
# 6. SALVAR
# ------------------------------------------------------------

ggplot2::ggsave(
  file.path(
    pasta_figuras,
    "Figura_38_descritores_maior_resposta_tratamentos.png"
  ),
  figura_38,
  width = 11,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

names(resultados_descritores_sala)

###############################################################################
# FIGURA 39
# VARIACAO DA ESTRUTURA ACUSTICA DENTRO DOS TIPOS MANUAIS
###############################################################################

# Garantir nomes padronizados dos tipos
rotulos_tipo <- c(
  tipo1 = "Tipo 1",
  tipo2 = "Tipo 2",
  tipo3 = "Tipo 3"
)

cores_tipo <- c(
  tipo1 = "#D73027",
  tipo2 = "#4575B4",
  tipo3 = "#1A9850"
)

# Selecionar os descritores com maiores tamanhos de efeito dentro de cada tipo
dados_figura_39 <- resultados_descritores_por_tipo %>%
  dplyr::filter(
    !is.na(epsilon2)
  ) %>%
  dplyr::group_by(tipo) %>%
  dplyr::slice_max(
    order_by = epsilon2,
    n = 5,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    classe = dplyr::case_when(
      is.na(p_ajustado_BH) ~ "Nao testado",
      p_ajustado_BH < 0.05 ~ "Significativo",
      p_ajustado_BH < 0.10 ~ "Evidencia limitada",
      TRUE ~ "Nao significativo"
    ),
    rotulo_descritor = paste(
      tipo,
      descritor,
      sep = "___"
    )
  )

# Ordenacao independente dentro de cada painel
ordem_39 <- dados_figura_39 %>%
  dplyr::arrange(
    tipo,
    epsilon2
  ) %>%
  dplyr::pull(rotulo_descritor)

dados_figura_39$rotulo_descritor <- factor(
  dados_figura_39$rotulo_descritor,
  levels = ordem_39
)

figura_39 <- ggplot2::ggplot(
  dados_figura_39,
  ggplot2::aes(
    x = epsilon2,
    y = rotulo_descritor,
    fill = classe
  )
) +
  ggplot2::geom_col(
    width = 0.70
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = paste0(
        "p ajust. = ",
        p_formatado
      )
    ),
    hjust = -0.04,
    size = 3.7
  ) +
  ggplot2::facet_wrap(
    ~ tipo,
    ncol = 1,
    scales = "free_y",
    labeller = ggplot2::labeller(
      tipo = rotulos_tipo
    )
  ) +
  ggplot2::scale_y_discrete(
    labels = function(x) {
      sub("^.*___", "", x)
    }
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Significativo" = "#1B9E77",
      "Evidencia limitada" = "#E66101",
      "Nao significativo" = "#777777",
      "Nao testado" = "#BBBBBB"
    )
  ) +
  ggplot2::scale_x_continuous(
    expand = ggplot2::expansion(
      mult = c(0, 0.25)
    )
  ) +
  ggplot2::labs(
    title = "Variação da estrutura acústica dentro dos tipos vocais manuais",
    subtitle = "Descritores com maiores tamanhos de efeito entre os tratamentos",
    x = "Tamanho de efeito (epsilon quadrado)",
    y = "Descritor acústico",
    fill = "Resultado"
  )

ggplot2::ggsave(
  file.path(
    pasta_figuras,
    "Figura_39_estrutura_acustica_tipos_manuais.png"
  ),
  figura_39,
  width = 11,
  height = 10,
  units = "in",
  dpi = 300,
  bg = "white"
)


###############################################################################
# FIGURA 40
# VARIACAO DA ESTRUTURA ACUSTICA DENTRO DOS GRUPOS COMPUTACIONAIS
###############################################################################

# ------------------------------------------------------------
# 1. Montar base unindo descritores e grupos computacionais
# ------------------------------------------------------------

dados_grupos_descritores <- dados %>%
  dplyr::select(
    id_vocalizacao,
    grupo_computacional
  ) %>%
  dplyr::left_join(
    dados_descritores,
    by = "id_vocalizacao"
  ) %>%
  dplyr::mutate(
    sala = factor(
      sala,
      levels = c(
        "Sala1",
        "Sala2",
        "Sala3",
        "Sala4"
      )
    ),
    grupo_computacional = factor(
      grupo_computacional
    )
  )


# ------------------------------------------------------------
# 2. Definir apenas descritores numericos usados nas analises
# ------------------------------------------------------------

descritores_grupos <- intersect(
  descritores,
  names(dados_grupos_descritores)
)

descritores_grupos <- descritores_grupos[
  vapply(
    dados_grupos_descritores[
      descritores_grupos
    ],
    is.numeric,
    logical(1)
  )
]


# ------------------------------------------------------------
# 3. Resumir cada grupo por sessao
# ------------------------------------------------------------

dados_descritores_grupo_sessao <- dados_grupos_descritores %>%
  dplyr::group_by(
    unidade,
    sala,
    ambiente,
    data,
    grupo_computacional
  ) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(
        descritores_grupos
      ),
      \(x) stats::median(
        x,
        na.rm = TRUE
      )
    ),
    n_vocalizacoes_grupo = dplyr::n(),
    .groups = "drop"
  )


# ------------------------------------------------------------
# 4. Kruskal-Wallis por descritor dentro de cada grupo
# ------------------------------------------------------------

resultados_descritores_por_grupo <- list()

for (
  grupo_atual in levels(
    dados_grupos_descritores$grupo_computacional
  )
) {
  
  dados_grupo_atual <- dados_descritores_grupo_sessao %>%
    dplyr::filter(
      grupo_computacional == grupo_atual
    ) %>%
    droplevels()
  
  resultado_grupo_atual <- dplyr::bind_rows(
    lapply(
      descritores_grupos,
      function(descritor_atual) {
        
        dados_teste <- dados_grupo_atual %>%
          dplyr::filter(
            is.finite(
              .data[[descritor_atual]]
            )
          ) %>%
          droplevels()
        
        # Condicoes minimas para testar
        if (
          nrow(dados_teste) < 4 ||
          dplyr::n_distinct(
            dados_teste$sala
          ) < 2 ||
          length(
            unique(
              dados_teste[
                [descritor_atual]
              ]
            )
          ) < 2
        ) {
          
          return(
            data.frame(
              grupo_computacional = grupo_atual,
              descritor = descritor_atual,
              H = NA_real_,
              gl = NA_real_,
              p = NA_real_,
              epsilon2 = NA_real_,
              n_sessoes = nrow(
                dados_teste
              )
            )
          )
        }
        
        teste <- tryCatch(
          stats::kruskal.test(
            x = dados_teste[
              [descritor_atual]
            ],
            g = dados_teste$sala
          ),
          error = function(e) NULL
        )
        
        if (is.null(teste)) {
          
          return(
            data.frame(
              grupo_computacional = grupo_atual,
              descritor = descritor_atual,
              H = NA_real_,
              gl = NA_real_,
              p = NA_real_,
              epsilon2 = NA_real_,
              n_sessoes = nrow(
                dados_teste
              )
            )
          )
        }
        
        n <- nrow(dados_teste)
        
        k <- dplyr::n_distinct(
          dados_teste$sala
        )
        
        epsilon2 <- if (
          n > k
        ) {
          
          max(
            0,
            (
              unname(
                teste$statistic
              ) -
                k +
                1
            ) /
              (
                n -
                  k
              )
          )
          
        } else {
          
          NA_real_
          
        }
        
        data.frame(
          grupo_computacional = grupo_atual,
          descritor = descritor_atual,
          H = unname(
            teste$statistic
          ),
          gl = unname(
            teste$parameter
          ),
          p = teste$p.value,
          epsilon2 = epsilon2,
          n_sessoes = n
        )
      }
    )
  )
  
  resultados_descritores_por_grupo[
    [grupo_atual]
  ] <- resultado_grupo_atual
}


# ------------------------------------------------------------
# 5. Juntar e corrigir p por BH dentro de cada grupo
# ------------------------------------------------------------

resultados_descritores_por_grupo <-
  dplyr::bind_rows(
    resultados_descritores_por_grupo
  ) %>%
  dplyr::group_by(
    grupo_computacional
  ) %>%
  dplyr::mutate(
    p_ajustado_BH = stats::p.adjust(
      p,
      method = "BH"
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    p_formatado = dplyr::case_when(
      is.na(
        p_ajustado_BH
      ) ~ "NA",
      p_ajustado_BH < 0.001 ~
        "< 0,001",
      TRUE ~
        sub(
          "\\.",
          ",",
          sprintf(
            "%.3f",
            p_ajustado_BH
          )
        )
    ),
    classe = dplyr::case_when(
      is.na(
        p_ajustado_BH
      ) ~ "Nao testado",
      p_ajustado_BH < 0.05 ~
        "Significativo",
      p_ajustado_BH < 0.10 ~
        "Evidencia limitada",
      TRUE ~
        "Nao significativo"
    )
  ) %>%
  dplyr::arrange(
    grupo_computacional,
    dplyr::desc(
      epsilon2
    )
  )


# Conferencia
print(
  resultados_descritores_por_grupo
)


# ------------------------------------------------------------
# 6. Salvar resultados
# ------------------------------------------------------------

utils::write.csv(
  resultados_descritores_por_grupo,
  file.path(
    pasta_saida,
    "Resultados_Descritores_Grupos_Computacionais_entre_Salas.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 7. Selecionar top 5 de cada grupo
# ------------------------------------------------------------

dados_figura_40 <-
  resultados_descritores_por_grupo %>%
  dplyr::filter(
    !is.na(epsilon2)
  ) %>%
  dplyr::group_by(
    grupo_computacional
  ) %>%
  dplyr::slice_max(
    order_by = epsilon2,
    n = 5,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    grupo_rotulo = paste(
      "Grupo",
      grupo_computacional
    ),
    descritor_grupo = paste(
      grupo_computacional,
      descritor,
      sep = "___"
    )
  )


ordem_40 <- dados_figura_40 %>%
  dplyr::arrange(
    grupo_computacional,
    epsilon2
  ) %>%
  dplyr::pull(
    descritor_grupo
  )

dados_figura_40$descritor_grupo <- factor(
  dados_figura_40$descritor_grupo,
  levels = ordem_40
)


# ------------------------------------------------------------
# 8. FIGURA 40
# ------------------------------------------------------------

figura_40 <- ggplot2::ggplot(
  dados_figura_40,
  ggplot2::aes(
    x = epsilon2,
    y = descritor_grupo,
    fill = classe
  )
) +
  ggplot2::geom_col(
    width = 0.70
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = paste0(
        "p ajust. = ",
        p_formatado
      )
    ),
    hjust = -0.04,
    size = 3.5
  ) +
  ggplot2::facet_wrap(
    ~ grupo_rotulo,
    scales = "free_y",
    ncol = 2
  ) +
  ggplot2::scale_y_discrete(
    labels = function(x) {
      sub("^.*___", "", x)
    }
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Significativo" = "#1B9E77",
      "Evidencia limitada" = "#E66101",
      "Nao significativo" = "#777777",
      "Nao testado" = "#BBBBBB"
    )
  ) +
  ggplot2::scale_x_continuous(
    expand = ggplot2::expansion(
      mult = c(
        0,
        0.30
      )
    )
  ) +
  ggplot2::labs(
    title = "Variação da estrutura acústica dentro dos grupos computacionais",
    subtitle = "Descritores com maiores tamanhos de efeito entre os tratamentos",
    x = "Tamanho de efeito (epsilon quadrado)",
    y = "Descritor acústico",
    fill = "Resultado"
  )

ggplot2::ggsave(
  file.path(
    pasta_figuras,
    "Figura_40_estrutura_acustica_grupos_computacionais.png"
  ),
  figura_40,
  width = 13,
  height = 10,
  units = "in",
  dpi = 300,
  bg = "white"
)

###############################################################################
# CORRECAO DA FIGURA 40
# ESTRUTURA ACUSTICA DOS GRUPOS COMPUTACIONAIS ENTRE OS TRATAMENTOS
###############################################################################

# Reinicia apenas o objeto que ficou vazio por causa do erro anterior
resultados_descritores_por_grupo <- list()


# ---------------------------------------------------------------------------
# 1. TESTAR OS DESCRITORES DENTRO DE CADA GRUPO COMPUTACIONAL
# ---------------------------------------------------------------------------

for (
  grupo_atual in levels(
    dados_grupos_descritores$grupo_computacional
  )
) {
  
  dados_grupo_atual <- dados_descritores_grupo_sessao %>%
    dplyr::filter(
      grupo_computacional == grupo_atual
    ) %>%
    droplevels()
  
  resultado_grupo_atual <- dplyr::bind_rows(
    
    lapply(
      descritores_grupos,
      
      function(descritor_atual) {
        
        dados_teste <- dados_grupo_atual %>%
          dplyr::filter(
            is.finite(
              .data[[descritor_atual]]
            )
          ) %>%
          droplevels()
        
        
        # ---------------------------------------------------------------
        # Verificacoes minimas
        # ---------------------------------------------------------------
        
        valores <- dados_teste[[descritor_atual]]
        
        if (
          nrow(dados_teste) < 4 ||
          dplyr::n_distinct(dados_teste$sala) < 2 ||
          length(unique(valores)) < 2
        ) {
          
          return(
            data.frame(
              grupo_computacional = grupo_atual,
              descritor = descritor_atual,
              H = NA_real_,
              gl = NA_real_,
              p = NA_real_,
              epsilon2 = NA_real_,
              n_sessoes = nrow(dados_teste)
            )
          )
        }
        
        
        # ---------------------------------------------------------------
        # Kruskal-Wallis
        # ---------------------------------------------------------------
        
        teste <- tryCatch(
          
          stats::kruskal.test(
            x = dados_teste[[descritor_atual]],
            g = dados_teste$sala
          ),
          
          error = function(e) NULL
        )
        
        
        if (is.null(teste)) {
          
          return(
            data.frame(
              grupo_computacional = grupo_atual,
              descritor = descritor_atual,
              H = NA_real_,
              gl = NA_real_,
              p = NA_real_,
              epsilon2 = NA_real_,
              n_sessoes = nrow(dados_teste)
            )
          )
        }
        
        
        # ---------------------------------------------------------------
        # Tamanho de efeito
        # ---------------------------------------------------------------
        
        n <- nrow(dados_teste)
        
        k <- dplyr::n_distinct(
          dados_teste$sala
        )
        
        epsilon2 <- if (n > k) {
          
          max(
            0,
            (
              unname(teste$statistic) -
                k +
                1
            ) /
              (
                n -
                  k
              )
          )
          
        } else {
          
          NA_real_
        }
        
        
        # ---------------------------------------------------------------
        # Resultado
        # ---------------------------------------------------------------
        
        data.frame(
          grupo_computacional = grupo_atual,
          descritor = descritor_atual,
          H = unname(teste$statistic),
          gl = unname(teste$parameter),
          p = teste$p.value,
          epsilon2 = epsilon2,
          n_sessoes = n
        )
      }
    )
  )
  
  
  resultados_descritores_por_grupo[[grupo_atual]] <-
    resultado_grupo_atual
}

resultados_descritores_por_grupo <-
  dplyr::bind_rows(
    resultados_descritores_por_grupo
  )

print(resultados_descritores_por_grupo)

dim(resultados_descritores_por_grupo)

# ---------------------------------------------------------------------------
# 2. JUNTAR TODOS OS GRUPOS
# ---------------------------------------------------------------------------

resultados_descritores_por_grupo <-
  dplyr::bind_rows(
    resultados_descritores_por_grupo
  ) %>%
  dplyr::mutate(
    grupo_computacional = factor(
      grupo_computacional,
      levels = levels(
        dados_grupos_descritores$grupo_computacional
      )
    )
  )


# Conferencia importante
print(resultados_descritores_por_grupo)

dim(resultados_descritores_por_grupo)


# ---------------------------------------------------------------------------
# 3. CORRECAO BH SEPARADAMENTE DENTRO DE CADA GRUPO
# ---------------------------------------------------------------------------

resultados_descritores_por_grupo <-
  resultados_descritores_por_grupo %>%
  
  dplyr::group_by(
    grupo_computacional
  ) %>%
  
  dplyr::mutate(
    
    p_ajustado_BH = stats::p.adjust(
      p,
      method = "BH"
    )
    
  ) %>%
  
  dplyr::ungroup() %>%
  
  dplyr::mutate(
    
    p_formatado = dplyr::case_when(
      
      is.na(p_ajustado_BH) ~
        "NA",
      
      p_ajustado_BH < 0.001 ~
        "< 0,001",
      
      TRUE ~
        sub(
          "\\.",
          ",",
          sprintf(
            "%.3f",
            p_ajustado_BH
          )
        )
    ),
    
    classe = dplyr::case_when(
      
      is.na(p_ajustado_BH) ~
        "Nao testado",
      
      p_ajustado_BH < 0.05 ~
        "Significativo",
      
      p_ajustado_BH < 0.10 ~
        "Evidencia limitada",
      
      TRUE ~
        "Nao significativo"
    )
    
  ) %>%
  
  dplyr::arrange(
    grupo_computacional,
    dplyr::desc(epsilon2)
  )


# ---------------------------------------------------------------------------
# 4. CONFERENCIA DOS RESULTADOS
# ---------------------------------------------------------------------------

print(
  resultados_descritores_por_grupo,
  n = Inf
)


# Quantidade de resultados significativos
n_significativos_grupos <- sum(
  resultados_descritores_por_grupo$p_ajustado_BH < 0.05,
  na.rm = TRUE
)

message(
  "\nComparacoes significativas descritor x grupo apos BH: ",
  n_significativos_grupos,
  "\n"
)


# ---------------------------------------------------------------------------
# 5. SALVAR RESULTADOS
# ---------------------------------------------------------------------------

utils::write.csv(
  
  resultados_descritores_por_grupo,
  
  file.path(
    pasta_saida,
    "Resultados_Descritores_Grupos_Computacionais_entre_Salas.csv"
  ),
  
  row.names = FALSE
)


###############################################################################
# FIGURA 40
# DESCRITORES COM MAIOR VARIACAO EM CADA GRUPO COMPUTACIONAL
###############################################################################

dados_figura_40 <-
  resultados_descritores_por_grupo %>%
  
  dplyr::filter(
    !is.na(epsilon2)
  ) %>%
  
  dplyr::group_by(
    grupo_computacional
  ) %>%
  
  dplyr::slice_max(
    order_by = epsilon2,
    n = 5,
    with_ties = FALSE
  ) %>%
  
  dplyr::ungroup() %>%
  
  dplyr::mutate(
    
    grupo_rotulo = paste(
      "Grupo",
      grupo_computacional
    ),
    
    descritor_grupo = paste(
      grupo_computacional,
      descritor,
      sep = "___"
    )
  )


# ---------------------------------------------------------------------------
# 6. Ordenacao dos descritores
# ---------------------------------------------------------------------------

ordem_40 <- dados_figura_40 %>%
  
  dplyr::arrange(
    grupo_computacional,
    epsilon2
  ) %>%
  
  dplyr::pull(
    descritor_grupo
  )


# Remove eventual duplicacao de níveis
ordem_40 <- unique(ordem_40)


dados_figura_40$descritor_grupo <- factor(
  
  dados_figura_40$descritor_grupo,
  
  levels = ordem_40
)


# ---------------------------------------------------------------------------
# 7. Construir figura
# ---------------------------------------------------------------------------

figura_40 <- ggplot2::ggplot(
  
  dados_figura_40,
  
  ggplot2::aes(
    x = epsilon2,
    y = descritor_grupo,
    fill = classe
  )
  
) +
  
  ggplot2::geom_col(
    width = 0.70
  ) +
  
  ggplot2::geom_text(
    
    ggplot2::aes(
      label = paste0(
        "p ajust. = ",
        p_formatado
      )
    ),
    
    hjust = -0.04,
    size = 3.5
  ) +
  
  ggplot2::facet_wrap(
    
    ~ grupo_rotulo,
    
    scales = "free_y",
    
    ncol = 2
    
  ) +
  
  ggplot2::scale_y_discrete(
    
    labels = function(x) {
      sub(
        "^.*___",
        "",
        x
      )
    }
    
  ) +
  
  ggplot2::scale_fill_manual(
    
    values = c(
      
      "Significativo" =
        "#1B9E77",
      
      "Evidencia limitada" =
        "#E66101",
      
      "Nao significativo" =
        "#777777",
      
      "Nao testado" =
        "#BBBBBB"
    )
    
  ) +
  
  ggplot2::scale_x_continuous(
    
    expand = ggplot2::expansion(
      mult = c(
        0,
        0.35
      )
    )
    
  ) +
  
  ggplot2::labs(
    
    title =
      "Variação da estrutura acústica dentro dos grupos computacionais",
    
    subtitle =
      "Descritores com maiores tamanhos de efeito entre os tratamentos",
    
    x =
      "Tamanho de efeito (epsilon quadrado)",
    
    y =
      "Descritor acústico",
    
    fill =
      "Resultado"
  )


# ---------------------------------------------------------------------------
# 8. Salvar Figura 40
# ---------------------------------------------------------------------------

ggplot2::ggsave(
  
  file.path(
    pasta_figuras,
    "Figura_40_estrutura_acustica_grupos_computacionais.png"
  ),
  
  figura_40,
  
  width = 13,
  height = 10,
  units = "in",
  dpi = 300,
  bg = "white"
)

cat("\nNúmero total de linhas:\n")
nrow(dados_sessao_descritores)

cat("\nNúmero de unidades/sessões únicas:\n")
dplyr::n_distinct(dados_sessao_descritores$unidade)

cat("\nLinhas por sala:\n")
dados_sessao_descritores %>%
  dplyr::count(sala)

cat("\nSessões únicas por sala:\n")
dados_sessao_descritores %>%
  dplyr::summarise(
    n_unidades = dplyr::n_distinct(unidade),
    .by = sala
  )

cat("\nDuplicação de unidades:\n")
dados_sessao_descritores %>%
  dplyr::count(unidade, sala) %>%
  dplyr::filter(n > 1)

cat("\nEstrutura do objeto:\n")
str(dados_sessao_descritores)

message(
  "\nFIGURA 40 CONCLUIDA COM SUCESSO.\n"
)

message(
  "\nFIGURAS 39 E 40 CONCLUIDAS.\n"
)

message(
  "\nFIGURA 38 GERADA COM SUCESSO.\n"
)

message(
  "\nFIGURAS ADICIONAIS PARA A APRESENTACAO CONCLUIDAS.\n",
  "Figuras 30 a 34 salvas em:\n",
  pasta_figuras
)

message(
  "\nPOS-TESTE DO TIPO 1 CONCLUIDO.\n",
  "Comparacoes significativas apos BH: ",
  nrow(comparacoes_significativas_tipo1),
  "\n"
)

message(
  "\nANALISE DA ESTRUTURA DOS TIPOS CONCLUIDA.\n",
  "Comparacoes significativas descritor x tipo: ",
  n_descritores_tipo_sig,
  "\n"
)
message("\nSCRIPT 04 CONCLUIDO.\nResultados em:\n", pasta_saida
        
        resultado_tipos_frequencia
        
        perfil_tipos_sala

        resultados_descritores_por_tipo %>%
          dplyr::group_by(tipo) %>%
          dplyr::slice_head(n = 5)        
        
        
        resultados_descritores_sala %>%
          dplyr::arrange(p) %>%
          dplyr::select(
            descritor,
            H,
            gl,
            p,
            epsilon2,
            p_ajustado_BH,
            p_BH
          ) %>%
          print(n = Inf)
        
        descritores_auditoria <- c(
          "startdom",
          "peakt",
          "maxdom",
          "modindx",
          "sp.ent",
          "skew"
        )
        
        auditoria_38 <- dplyr::bind_rows(
          lapply(
            descritores_auditoria,
            function(d) {
              
              x <- dados_sessao_descritores[[d]]
              g <- dados_sessao_descritores$sala
              
              ok <- is.finite(x) & !is.na(g)
              
              x <- x[ok]
              g <- droplevels(g[ok])
              
              teste <- stats::kruskal.test(
                x = x,
                g = g
              )
              
              n <- length(x)
              k <- nlevels(g)
              
              eps2 <- max(
                0,
                (
                  unname(teste$statistic) -
                    k + 1
                ) /
                  (n - k)
              )
              
              data.frame(
                descritor = d,
                n = n,
                H = unname(teste$statistic),
                gl = unname(teste$parameter),
                p = teste$p.value,
                epsilon2 = eps2
              )
            }
          )
        ) %>%
          dplyr::mutate(
            p_BH = stats::p.adjust(
              p,
              method = "BH"
            )
          ) %>%
          dplyr::arrange(p_BH)
        
        print(auditoria_38)
        
        ###############################################################################
        # FIGURA 38 CORRIGIDA
        # ESTRUTURA ACUSTICA GERAL ENTRE OS TRATAMENTOS
        # UNIDADE ANALITICA = SESSAO DE GRAVACAO
        ###############################################################################
        
        resultados_gerais_corretos <- auditoria_38 %>%
          dplyr::mutate(
            
            classe = dplyr::case_when(
              p_BH < 0.05 ~ "Significativo",
              p_BH < 0.10 ~ "Evidencia limitada",
              TRUE ~ "Nao significativo"
            ),
            
            p_formatado = dplyr::case_when(
              p_BH < 0.001 ~ "< 0,001",
              TRUE ~ sub(
                "\\.",
                ",",
                sprintf("%.3f", p_BH)
              )
            )
          ) %>%
          dplyr::arrange(
            dplyr::desc(epsilon2)
          )
        
        
        resultados_gerais_corretos$descritor <- factor(
          resultados_gerais_corretos$descritor,
          levels = rev(
            resultados_gerais_corretos$descritor
          )
        )
        
        
        figura_38_corrigida <- ggplot2::ggplot(
          resultados_gerais_corretos,
          ggplot2::aes(
            x = epsilon2,
            y = descritor,
            fill = classe
          )
        ) +
          
          ggplot2::geom_col(
            width = 0.70
          ) +
          
          ggplot2::geom_text(
            ggplot2::aes(
              label = paste0(
                "p ajust. = ",
                p_formatado
              )
            ),
            hjust = -0.05,
            size = 4
          ) +
          
          ggplot2::scale_fill_manual(
            values = c(
              "Significativo" = "#1B9E77",
              "Evidencia limitada" = "#E66101",
              "Nao significativo" = "#777777"
            )
          ) +
          
          ggplot2::scale_x_continuous(
            expand = ggplot2::expansion(
              mult = c(0, 0.25)
            )
          ) +
          
          ggplot2::labs(
            title = "Variação da estrutura acústica geral entre os tratamentos",
            subtitle = "Análise exploratória com as sessões de gravação como unidades independentes",
            x = "Tamanho de efeito (epsilon quadrado)",
            y = "Descritor acústico",
            fill = "Resultado"
          )
        
        
        ggplot2::ggsave(
          file.path(
            pasta_figuras,
            "Figura_38_CORRIGIDA_estrutura_acustica_geral.png"
          ),
          figura_38_corrigida,
          width = 11,
          height = 7,
          units = "in",
          dpi = 300,
          bg = "white"
        )

        ###############################################################################
        # AUDITORIA DAS ANALISES PARA A APRESENTACAO
        # Objetivo:
        # verificar quais objetos existem e qual e a unidade amostral utilizada
        ###############################################################################
        
        cat("\n============================================================\n")
        cat("AUDITORIA DAS ANALISES DA APRESENTACAO\n")
        cat("============================================================\n\n")
        
        
        # ---------------------------------------------------------------------------
        # 1. OBJETOS IMPORTANTES DISPONIVEIS
        # ---------------------------------------------------------------------------
        
        objetos_auditoria <- c(
          "dados",
          "dados_descritores",
          "dados_sessao_descritores",
          "composicao_tipos_sessao",
          "perfil_tipos_sala",
          "resultados_descritores_por_tipo",
          "dados_grupos_descritores",
          "dados_descritores_grupo_sessao",
          "resultados_descritores_por_grupo"
        )
        
        cat("OBJETOS DISPONIVEIS:\n\n")
        
        for (obj in objetos_auditoria) {
          
          cat(
            obj,
            ": ",
            ifelse(
              exists(obj),
              "SIM",
              "NAO"
            ),
            "\n",
            sep = ""
          )
        }
        
        
        # ---------------------------------------------------------------------------
        # 2. BASE GERAL POR SESSAO
        # ---------------------------------------------------------------------------
        
        if (exists("dados_sessao_descritores")) {
          
          cat("\n\n============================================================\n")
          cat("1. ESTRUTURA ACUSTICA GERAL\n")
          cat("============================================================\n")
          
          cat(
            "\nNumero de linhas: ",
            nrow(dados_sessao_descritores),
            "\n",
            sep = ""
          )
          
          cat(
            "Numero de sessoes unicas: ",
            dplyr::n_distinct(
              dados_sessao_descritores$unidade
            ),
            "\n",
            sep = ""
          )
          
          cat("\nLinhas por sala:\n")
          
          print(
            dados_sessao_descritores %>%
              dplyr::count(sala)
          )
          
          cat("\nSessoes unicas por sala:\n")
          
          print(
            dados_sessao_descritores %>%
              dplyr::summarise(
                n_sessoes =
                  dplyr::n_distinct(unidade),
                .by = sala
              )
          )
          
          cat("\nDuplicacoes unidade x sala:\n")
          
          print(
            dados_sessao_descritores %>%
              dplyr::count(
                unidade,
                sala
              ) %>%
              dplyr::filter(
                n > 1
              )
          )
        }
        
        
        # ---------------------------------------------------------------------------
        # 3. TIPOS VOCAIS MANUAIS POR SESSAO
        # ---------------------------------------------------------------------------
        
        if (exists("composicao_tipos_sessao")) {
          
          cat("\n\n============================================================\n")
          cat("2. TIPOS VOCAIS MANUAIS\n")
          cat("============================================================\n")
          
          cat(
            "\nNumero de linhas: ",
            nrow(composicao_tipos_sessao),
            "\n",
            sep = ""
          )
          
          cat(
            "Numero de sessoes unicas: ",
            dplyr::n_distinct(
              composicao_tipos_sessao$unidade
            ),
            "\n",
            sep = ""
          )
          
          cat("\nSessoes por sala:\n")
          
          print(
            composicao_tipos_sessao %>%
              dplyr::summarise(
                n_sessoes =
                  dplyr::n_distinct(unidade),
                .by = sala
              )
          )
          
          cat("\nLinhas por tipo vocal:\n")
          
          print(
            composicao_tipos_sessao %>%
              dplyr::count(
                tipo
              )
          )
          
          cat("\nSessoes por tipo vocal e sala:\n")
          
          print(
            composicao_tipos_sessao %>%
              dplyr::summarise(
                n_sessoes =
                  dplyr::n_distinct(unidade),
                .by = c(
                  tipo,
                  sala
                )
              )
          )
          
          cat("\nDuplicacoes unidade x sala x tipo:\n")
          
          print(
            composicao_tipos_sessao %>%
              dplyr::count(
                unidade,
                sala,
                tipo
              ) %>%
              dplyr::filter(
                n > 1
              )
          )
        }
        
        
        # ---------------------------------------------------------------------------
        # 4. ESTRUTURA ACUSTICA DENTRO DOS TIPOS MANUAIS
        # ---------------------------------------------------------------------------
        
        if (exists("resultados_descritores_por_tipo")) {
          
          cat("\n\n============================================================\n")
          cat("3. ESTRUTURA DENTRO DOS TIPOS MANUAIS\n")
          cat("============================================================\n")
          
          cat("\nNumero de resultados por tipo:\n")
          
          print(
            resultados_descritores_por_tipo %>%
              dplyr::count(
                tipo
              )
          )
          
          if (
            "n_sessoes" %in%
            names(resultados_descritores_por_tipo)
          ) {
            
            cat("\nNumero de sessoes utilizado nos testes:\n")
            
            print(
              resultados_descritores_por_tipo %>%
                dplyr::summarise(
                  minimo = min(
                    n_sessoes,
                    na.rm = TRUE
                  ),
                  maximo = max(
                    n_sessoes,
                    na.rm = TRUE
                  ),
                  .by = tipo
                )
            )
          }
          
          cat("\nMenores p ajustados por tipo:\n")
          
          print(
            resultados_descritores_por_tipo %>%
              dplyr::filter(
                !is.na(p_ajustado_BH)
              ) %>%
              dplyr::arrange(
                tipo,
                p_ajustado_BH
              ) %>%
              dplyr::group_by(
                tipo
              ) %>%
              dplyr::slice_head(
                n = 3
              ) %>%
              dplyr::ungroup() %>%
              dplyr::select(
                dplyr::any_of(
                  c(
                    "tipo",
                    "descritor",
                    "n_sessoes",
                    "H",
                    "gl",
                    "p",
                    "epsilon2",
                    "p_ajustado_BH"
                  )
                )
              )
          )
        }
        
        
        # ---------------------------------------------------------------------------
        # 5. GRUPOS COMPUTACIONAIS POR SESSAO
        # ---------------------------------------------------------------------------
        
        if (exists("dados_descritores_grupo_sessao")) {
          
          cat("\n\n============================================================\n")
          cat("4. GRUPOS COMPUTACIONAIS\n")
          cat("============================================================\n")
          
          cat(
            "\nNumero de linhas: ",
            nrow(dados_descritores_grupo_sessao),
            "\n",
            sep = ""
          )
          
          cat(
            "Numero de sessoes unicas: ",
            dplyr::n_distinct(
              dados_descritores_grupo_sessao$unidade
            ),
            "\n",
            sep = ""
          )
          
          cat("\nSessoes por grupo computacional:\n")
          
          print(
            dados_descritores_grupo_sessao %>%
              dplyr::summarise(
                n_sessoes =
                  dplyr::n_distinct(unidade),
                .by =
                  grupo_computacional
              )
          )
          
          cat("\nSessoes por grupo e sala:\n")
          
          print(
            dados_descritores_grupo_sessao %>%
              dplyr::summarise(
                n_sessoes =
                  dplyr::n_distinct(unidade),
                .by = c(
                  grupo_computacional,
                  sala
                )
              )
          )
          
          cat("\nDuplicacoes unidade x sala x grupo:\n")
          
          print(
            dados_descritores_grupo_sessao %>%
              dplyr::count(
                unidade,
                sala,
                grupo_computacional
              ) %>%
              dplyr::filter(
                n > 1
              )
          )
        }
        
        
        # ---------------------------------------------------------------------------
        # 6. RESULTADOS DOS GRUPOS COMPUTACIONAIS
        # ---------------------------------------------------------------------------
        
        if (
          exists("resultados_descritores_por_grupo") &&
          is.data.frame(
            resultados_descritores_por_grupo
          )
        ) {
          
          cat("\n\n============================================================\n")
          cat("5. RESULTADOS DOS GRUPOS COMPUTACIONAIS\n")
          cat("============================================================\n")
          
          cat("\nNumero de resultados por grupo:\n")
          
          print(
            resultados_descritores_por_grupo %>%
              dplyr::count(
                grupo_computacional
              )
          )
          
          if (
            "n_sessoes" %in%
            names(resultados_descritores_por_grupo)
          ) {
            
            cat("\nNumero de sessoes utilizado:\n")
            
            print(
              resultados_descritores_por_grupo %>%
                dplyr::summarise(
                  minimo = min(
                    n_sessoes,
                    na.rm = TRUE
                  ),
                  maximo = max(
                    n_sessoes,
                    na.rm = TRUE
                  ),
                  .by =
                    grupo_computacional
                )
            )
          }
          
          cat("\nMenores p ajustados por grupo:\n")
          
          print(
            resultados_descritores_por_grupo %>%
              dplyr::filter(
                !is.na(p_ajustado_BH)
              ) %>%
              dplyr::arrange(
                grupo_computacional,
                p_ajustado_BH
              ) %>%
              dplyr::group_by(
                grupo_computacional
              ) %>%
              dplyr::slice_head(
                n = 3
              ) %>%
              dplyr::ungroup() %>%
              dplyr::select(
                dplyr::any_of(
                  c(
                    "grupo_computacional",
                    "descritor",
                    "n_sessoes",
                    "H",
                    "gl",
                    "p",
                    "epsilon2",
                    "p_ajustado_BH"
                  )
                )
              )
          )
        }
        
        
        # ---------------------------------------------------------------------------
        # 7. BASE ORIGINAL
        # ---------------------------------------------------------------------------
        
        if (exists("dados")) {
          
          cat("\n\n============================================================\n")
          cat("6. BASE ORIGINAL DE VOCALIZACOES\n")
          cat("============================================================\n")
          
          cat(
            "\nNumero total de vocalizacoes: ",
            nrow(dados),
            "\n",
            sep = ""
          )
          
          cat("\nColunas disponíveis:\n")
          
          print(
            names(dados)
          )
        }
        
        
        cat("\n\n============================================================\n")
        cat("AUDITORIA CONCLUIDA\n")
        cat("============================================================\n")
        
        ###############################################################################
        # AUDITORIA FINAL — ATIVIDADE VOCAL E COMPOSICAO
        ###############################################################################
        
        cat("\n============================================================\n")
        cat("AUDITORIA DE ATIVIDADE E COMPOSICAO\n")
        cat("============================================================\n")
        
        
        # ---------------------------------------------------------------------------
        # 1. VOCALIZACOES POR SESSAO
        # ---------------------------------------------------------------------------
        
        atividade_sessao_auditoria <- dados %>%
          dplyr::count(
            unidade,
            sala,
            ambiente,
            data,
            name = "n_vocalizacoes"
          )
        
        cat("\n1. ATIVIDADE VOCAL POR SESSAO\n")
        
        cat(
          "\nNumero de sessoes:",
          dplyr::n_distinct(
            atividade_sessao_auditoria$unidade
          ),
          "\n"
        )
        
        print(
          atividade_sessao_auditoria %>%
            dplyr::summarise(
              n_sessoes = dplyr::n(),
              minimo = min(n_vocalizacoes),
              mediana = median(n_vocalizacoes),
              maximo = max(n_vocalizacoes),
              .by = sala
            )
        )
        
        
        # ---------------------------------------------------------------------------
        # 2. TIPOS VOCAIS — CONTAGEM POR SESSAO
        # ---------------------------------------------------------------------------
        
        tipos_sessao_auditoria <- dados %>%
          dplyr::count(
            unidade,
            sala,
            ambiente,
            data,
            tipo,
            name = "n_tipo"
          )
        
        cat("\n\n2. TIPOS VOCAIS POR SESSAO\n")
        
        cat(
          "\nNumero de sessoes unicas:",
          dplyr::n_distinct(
            tipos_sessao_auditoria$unidade
          ),
          "\n"
        )
        
        cat("\nNumero de sessoes contendo cada tipo:\n")
        
        print(
          tipos_sessao_auditoria %>%
            dplyr::summarise(
              n_sessoes =
                dplyr::n_distinct(unidade),
              .by = tipo
            )
        )
        
        cat("\nNumero de sessoes contendo cada tipo por sala:\n")
        
        print(
          tipos_sessao_auditoria %>%
            dplyr::summarise(
              n_sessoes =
                dplyr::n_distinct(unidade),
              .by = c(
                tipo,
                sala
              )
            )
        )
        
        
        # ---------------------------------------------------------------------------
        # 3. PROPORCAO DE CADA TIPO DENTRO DE CADA SESSAO
        # ---------------------------------------------------------------------------
        
        composicao_sessao_auditoria <- dados %>%
          dplyr::count(
            unidade,
            sala,
            ambiente,
            data,
            tipo,
            name = "n_tipo"
          ) %>%
          dplyr::group_by(
            unidade,
            sala,
            ambiente,
            data
          ) %>%
          dplyr::mutate(
            total_sessao = sum(n_tipo),
            proporcao = n_tipo / total_sessao
          ) %>%
          dplyr::ungroup()
        
        cat("\n\n3. COMPOSICAO POR SESSAO\n")
        
        cat(
          "\nSessoes representadas:",
          dplyr::n_distinct(
            composicao_sessao_auditoria$unidade
          ),
          "\n"
        )
        
        cat("\nSoma das proporcoes dentro das sessoes:\n")
        
        print(
          composicao_sessao_auditoria %>%
            dplyr::summarise(
              soma = sum(proporcao),
              .by = unidade
            )
        )
        
        
        # ---------------------------------------------------------------------------
        # 4. GRUPOS COMPUTACIONAIS — CONTAGEM POR SESSAO
        # ---------------------------------------------------------------------------
        
        grupos_sessao_auditoria <- dados %>%
          dplyr::count(
            unidade,
            sala,
            ambiente,
            data,
            grupo_computacional,
            name = "n_grupo"
          )
        
        cat("\n\n4. GRUPOS COMPUTACIONAIS POR SESSAO\n")
        
        cat(
          "\nNumero de sessoes:",
          dplyr::n_distinct(
            grupos_sessao_auditoria$unidade
          ),
          "\n"
        )
        
        cat("\nSessoes contendo cada grupo:\n")
        
        print(
          grupos_sessao_auditoria %>%
            dplyr::summarise(
              n_sessoes =
                dplyr::n_distinct(unidade),
              .by = grupo_computacional
            )
        )
        
        cat("\nSessoes contendo cada grupo por sala:\n")
        
        print(
          grupos_sessao_auditoria %>%
            dplyr::summarise(
              n_sessoes =
                dplyr::n_distinct(unidade),
              .by = c(
                grupo_computacional,
                sala
              )
            )
        )
        
        
        # ---------------------------------------------------------------------------
        # 5. PROPORCAO DOS GRUPOS DENTRO DE CADA SESSAO
        # ---------------------------------------------------------------------------
        
        composicao_grupos_sessao_auditoria <- dados %>%
          dplyr::count(
            unidade,
            sala,
            ambiente,
            data,
            grupo_computacional,
            name = "n_grupo"
          ) %>%
          dplyr::group_by(
            unidade,
            sala,
            ambiente,
            data
          ) %>%
          dplyr::mutate(
            total_sessao = sum(n_grupo),
            proporcao = n_grupo / total_sessao
          ) %>%
          dplyr::ungroup()
        
        cat("\n\n5. COMPOSICAO DOS GRUPOS POR SESSAO\n")
        
        cat(
          "\nSessoes representadas:",
          dplyr::n_distinct(
            composicao_grupos_sessao_auditoria$unidade
          ),
          "\n"
        )
        
        cat("\nSoma das proporcoes dentro das sessoes:\n")
        
        print(
          composicao_grupos_sessao_auditoria %>%
            dplyr::summarise(
              soma = sum(proporcao),
              .by = unidade
            )
        )
        
        
        cat("\n============================================================\n")
        cat("FIM DA AUDITORIA\n")
        cat("============================================================\n")
        