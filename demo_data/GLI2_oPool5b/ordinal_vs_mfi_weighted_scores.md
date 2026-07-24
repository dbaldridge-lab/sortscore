# Sortscore MFI VS Ordinal Weights


- [Setup](#setup)
- [Score variants using ordinal values entered into the experimental
  setup
  file](#score-variants-using-ordinal-values-entered-into-the-experimental-setup-file)
- [Figures](#figures)

This notebook compares outputs scores produced using sortscore after
changing bin weights from median MFI (AFU) to ordinal values.

## Setup

Demo scores used for comparison were previously produced using
demo_data/single_experiment_demo.ipynb.

``` r
aa_scores <- readr::read_csv(
  "output/scores/demo_aa_scores.csv",
  show_col_types = FALSE
)
dna_scores <- readr::read_csv(
  "output/scores/demo_dna_scores.csv",
  show_col_types = FALSE
)
aa_scores_ordinal <- readr::read_csv(
  "output_ordinal/scores/opool5b_ordinal_aa_scores.csv",
  show_col_types = FALSE
)
dna_scores_ordinal <- readr::read_csv(
  "output_ordinal/scores/opool5b_ordinal_dna_scores.csv",
  show_col_types = FALSE
)

sortscore_python <- "venv/bin/python" # path to an interpreter with sortscore installed
experiment_name <- "opool5b_ordinal"
experiment_setup_path <- "demo_data/GLI2_oPool5b/experiment_setup_ordinal.csv"
config_path <- "demo_data/GLI2_oPool5b/config_ordinal.json"
```

## Score variants using ordinal values entered into the experimental setup file

``` r
library(here)
```

    here() starts at /Users/c.chitwood/code/sortscore

``` r
here::i_am("demo_data/GLI2_oPool5b/ordinal_vs_mfi_weighted_scores.qmd")
```

    here() starts at /Users/c.chitwood/code/sortscore

``` r
args <- c(
  "-m", "sortscore", "score",
  "-n", experiment_name,
  "-e", here::here(experiment_setup_path),
  "-c", here::here(config_path)
)

sortscore_python_path <- if (grepl("^(/|~)", sortscore_python)) {
  path.expand(sortscore_python)
} else {
  here::here(sortscore_python)
}

system2(sortscore_python_path, args)
```

## Figures

``` r
library(ggplot2)

complete_rows <- complete.cases(
  aa_scores$score,
  aa_scores_ordinal$score
)

aa_corr <- cor.test(
  x = aa_scores$score[complete_rows],
  y = aa_scores_ordinal$score[complete_rows],
  method = "spearman",
  exact = FALSE
)

label <- sprintf(
  "Spearman's \u03c1 = %.2f\nn = %d",
  unname(aa_corr$estimate),
  sum(complete_rows)
)

aa_scatterplot <- ggplot() +
  aes(
    x = aa_scores$score,
    y = aa_scores_ordinal$score
  ) +
  geom_smooth(method = lm, se = FALSE, color = "red", linewidth = 1) +
  geom_point(size = 1.5) +
  annotate(
      "text",
      x = Inf,
      y = -Inf,
      label = label,
      hjust = 1.1,
      vjust = -0.5
      ) +
  labs(
    title = "AA scores using MFI and ordinal bin weights",
    x = "MFI Weighted",
    y = "Ordinal Weighted"
    )

print(aa_scatterplot)
```

![](ordinal_vs_mfi_weighted_scores_files/figure-commonmark/AA%20scores%20-%20ordinal%20vs%20MFI%20weighted%20comparison-1.png)

``` r
ggsave(
  filename = "aa_scatterplot.png",
  plot = aa_scatterplot,
  width = 6, height = 5,
  units = "in",
  dpi = 300
)
```

``` r
dna_complete_rows <- complete.cases(
  dna_scores$score,
  dna_scores_ordinal$score
)

dna_corr <- cor.test(
  x = dna_scores$score[dna_complete_rows],
  y = dna_scores_ordinal$score[dna_complete_rows],
  method = "spearman",
  exact = FALSE
)

dna_label <- sprintf(
  "Spearman's \u03c1 = %.2f\nn = %d",
  unname(dna_corr$estimate),
  sum(dna_complete_rows)
)

dna_scatterplot <- ggplot() +
  aes(
    x = dna_scores$score,
    y = dna_scores_ordinal$score
  ) +
  geom_smooth(method = lm, se = FALSE, color = "red", linewidth = 1) +
  geom_point(size = 1.5) +
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = dna_label,
    hjust = 1.1,
    vjust = -0.5
  ) +
  labs(
    title = "DNA scores using MFI and ordinal bin weights",
    x = "MFI Weighted",
    y = "Ordinal Bin Weighted"
  ) 

print(dna_scatterplot)
```

![](ordinal_vs_mfi_weighted_scores_files/figure-commonmark/DNA%20scores%20-%20ordinal%20vs%20MFI%20weighted%20comparison-1.png)

``` r
ggsave(
  filename = "dna_scatterplot.png",
  plot = dna_scatterplot,
  width = 6, height = 5,
  units = "in",
  dpi = 300
)
```
