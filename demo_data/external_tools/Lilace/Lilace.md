# Lilace Scoring Demo


- [Use sortscore Output as Lilace
  Input](#use-sortscore-output-as-lilace-input)
  - [Parameters](#parameters)
  - [Installation](#installation)
  - [Setup](#setup)
  - [Generate Input File](#generate-input-file)
    - [Use sortscore output to get Lilace input
      counts](#use-sortscore-output-to-get-lilace-input-counts)
    - [Load input counts](#load-input-counts)
    - [Prepare counts for Lilace](#prepare-counts-for-lilace)
  - [Run Lilace](#run-lilace)
  - [Compare sortscore and Lilace
    Scores](#compare-sortscore-and-lilace-scores)

# Use sortscore Output as Lilace Input

The code in this notebook is adapted from the Lilace vignette to
demonstrate how to utilize sortscore output as an input for this
workflow.

Cite: Freudenberg, J., Rao, J., Howard, M. K., Macdonald, C., Greenwald,
N. F., Coyote-Maestas, W., & Pimentel, H. (2026). Accurate variant
effect estimation in FACS-based deep mutational scanning data with
Lilace. Genome Biology, 27(1), 48.
https://doi.org/10.1186/s13059-026-03934-1

Reference vignette: [Introduction to
Lilace](https://pimentellab.com/lilace/articles/intro.html).

To recreate this analysis, first run `demo_data/tiled_demo.ipynb` to
generate the required input files.

0.  The `sortscore integrate` command is used to return a Lilace-ready
    counts table.
1.  Import scores and metadata from sortscore into Lilace’s input
    format.
2.  Normalize according to cell sorting percentages, if needed.
3.  Run Lilace.
4.  Review the output scores.

## Parameters

Edit these values before running the analysis for a new dataset.

``` r
# Paths below are repo-relative unless you provide an absolute path.
counts_path <- "demo_data/external_tools/Lilace_input_counts.csv"
working_dir <- "demo_data/external_tools/Lilace"
input_scores_path <- "demo_data/_test_outputs/multitile_output/normalized/zscore_2pole/scores/batch_aa_scores.csv"

# Optional filter for tiled experiments. Set to NULL to include all batches.
batch_name <- "tile1"

# Export type for `sortscore integrate lilace`.
mutagenesis_type <- "aa" # set to "codon" for dna score tables

sortscore_python <- "venv/bin/python" # path to an interpreter with sortscore installed

# Set to TRUE if you want Lilace counts normalized to sorting proportions.
use_sort_normalization <- FALSE

# Sorting proportions aligned to the integrated Lilace count columns (`c_0`, `c_1`, ...).
sort_prop_values <- c(0.25, 0.25, 0.25, 0.25)

control_label <- "synonymous"
n_parallel_chains <- 4
seed <- 1
```

## Installation

The Lilace vignette recommends installing `cmdstanr` first, then
`CmdStan`, and finally Lilace itself.

## Setup

Once the remaining required packages are installed, load them in the
setup chunk below.

## Generate Input File

### Use sortscore output to get Lilace input counts

Generates the Lilace input counts file from a sortscore scores table.

``` r
integrate_args <- c(
  "-m", "sortscore", "integrate", "lilace",
  "--input", here::here(input_scores_path),
  "--output", here::here(counts_path),
  "--mutagenesis-type", mutagenesis_type
)

if (!is.null(batch_name) && nzchar(batch_name)) {
  integrate_args <- c(integrate_args, "--batch", batch_name)
}

sortscore_python_path <- if (grepl("^(/|~)", sortscore_python)) {
  path.expand(sortscore_python)
} else {
  here::here(sortscore_python)
}

system2(sortscore_python_path, integrate_args)
```

### Load input counts

The `sortscore integrate lilace` preprocessing step produces one row per
variant per replicate, including `variant_id`, `mutation_type`,
`position`, and Lilace count columns corresponding to the FACS bins.

``` r
counts_file <- here::here(counts_path)
working_dir_path <- here::here(working_dir)

raw_counts <- readr::read_csv(counts_file, show_col_types = FALSE)

dplyr::glimpse(raw_counts)
```

    Rows: 2,082
    Columns: 10
    $ variant_id             <chr> "A.13.*", "A.13.*", "A.13.*", "A.13.=", "A.13.=…
    $ mutation_type          <chr> "nonsense", "nonsense", "nonsense", "synonymous…
    $ position               <dbl> 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,…
    $ replicate              <dbl> 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1,…
    $ c_0                    <dbl> 47161, 40633, 45787, 99927, 88360, 82385, 23329…
    $ c_1                    <dbl> 36737, 46858, 38992, 45673, 66177, 67338, 23001…
    $ c_2                    <dbl> 16111, 24719, 25885, 123360, 149455, 127665, 63…
    $ c_3                    <dbl> 2979, 5568, 2625, 173779, 198722, 205987, 11855…
    $ sortscore_score        <dbl> -3.8703863, -3.8703863, -3.8703863, -0.5436406,…
    $ sortscore_score_column <chr> "score", "score", "score", "score", "score", "s…

``` r
head(raw_counts)
```

    # A tibble: 6 × 10
      variant_id mutation_type position replicate   c_0   c_1    c_2    c_3
      <chr>      <chr>            <dbl>     <dbl> <dbl> <dbl>  <dbl>  <dbl>
    1 A.13.*     nonsense            13         1 47161 36737  16111   2979
    2 A.13.*     nonsense            13         2 40633 46858  24719   5568
    3 A.13.*     nonsense            13         3 45787 38992  25885   2625
    4 A.13.=     synonymous          13         1 99927 45673 123360 173779
    5 A.13.=     synonymous          13         2 88360 66177 149455 198722
    6 A.13.=     synonymous          13         3 82385 67338 127665 205987
    # ℹ 2 more variables: sortscore_score <dbl>, sortscore_score_column <chr>

### Prepare counts for Lilace

The integrated input file already includes Lilace’s `variant_id`,
`mutation_type`, `position`, `replicate`, and `c_0`, `c_1`, … count
columns. The step below detects those count columns and separates any
extra metadata.

``` r
stopifnot(all(c("variant_id", "mutation_type", "position", "replicate") %in% names(raw_counts)))

count_col_names <- names(raw_counts)[grepl("^c_[0-9]+$", names(raw_counts))]
stopifnot(length(count_col_names) > 0)

lilace_counts <- raw_counts |>
  dplyr::select(
    variant_id,
    mutation_type,
    position,
    replicate,
    dplyr::all_of(count_col_names)
  )

lilace_counts[count_col_names] <- lapply(
  lilace_counts[count_col_names],
  function(col) as.integer(round(col))
)

metadata_cols <- raw_counts |>
  dplyr::select(-replicate, -dplyr::all_of(count_col_names))

lilace_counts |>
  dplyr::arrange(position, variant_id, replicate) |>
  head()
```

    # A tibble: 6 × 8
      variant_id mutation_type position replicate    c_0    c_1   c_2   c_3
      <chr>      <chr>            <dbl>     <dbl>  <int>  <int> <int> <int>
    1 C.1.*      nonsense             1         1 172748 141422 40350  8546
    2 C.1.*      nonsense             1         2 157075 155724 83279 14289
    3 C.1.*      nonsense             1         3 171978 142359 75873  5369
    4 C.1.=      synonymous           1         1   1683   1051  2910 10219
    5 C.1.=      synonymous           1         2   1404   1545  5604  7472
    6 C.1.=      synonymous           1         3   1520   1085  5148 12823

``` r
metadata_cols |>
  dplyr::arrange(position, variant_id) |>
  head()
```

    # A tibble: 6 × 5
      variant_id mutation_type position sortscore_score sortscore_score_column
      <chr>      <chr>            <dbl>           <dbl> <chr>                 
    1 C.1.*      nonsense             1           -3.86 score                 
    2 C.1.*      nonsense             1           -3.86 score                 
    3 C.1.*      nonsense             1           -3.86 score                 
    4 C.1.=      synonymous           1            1.53 score                 
    5 C.1.=      synonymous           1            1.53 score                 
    6 C.1.=      synonymous           1            1.53 score                 

Convert the counts table into a Lilace object.

``` r
lilace_obj <- lilace::lilace_from_counts(
  variant_id = lilace_counts$variant_id,
  mutation_type = lilace_counts$mutation_type,
  position = lilace_counts$position,
  replicate = lilace_counts$replicate,
  counts = lilace_counts |> dplyr::select(dplyr::all_of(count_col_names)),
  metadata = metadata_cols
)

head(lilace_obj$data)
```

    # A tibble: 6 × 15
      variant variant_id mutation_type position sortscore_score
      <chr>   <chr>      <chr>            <dbl>           <dbl>
    1 C.1.*   C.1.*      nonsense             1           -3.86
    2 C.1.*   C.1.*      nonsense             1           -3.86
    3 C.1.*   C.1.*      nonsense             1           -3.86
    4 C.1.=   C.1.=      synonymous           1            1.53
    5 C.1.=   C.1.=      synonymous           1            1.53
    6 C.1.=   C.1.=      synonymous           1            1.53
    # ℹ 10 more variables: sortscore_score_column <chr>, type <chr>,
    #   position.1 <dbl>, rep <dbl>, c_0 <int>, c_1 <int>, c_2 <int>, c_3 <int>,
    #   n_counts <dbl>, total_counts <dbl>

## Run Lilace

Optionally normalize to the FACS sorting proportions.

``` r
if (use_sort_normalization) {
  stopifnot(length(sort_prop_values) == length(count_col_names))
  sort_props <- stats::setNames(sort_prop_values, count_col_names)

  lilace_obj <- lilace::lilace_sorting_normalize(
    lilace_obj,
    sort_props,
    rep_specific = FALSE
  )
}

if (use_sort_normalization) {
  head(lilace_obj$normalized_data)
}
```

Fit the Lilace model. Per the Lilace vignette, synonymous variants are
used as the negative control group by default.

``` r
lilace_obj <- lilace::lilace_fit_model(
  lilace_obj,
  output_dir = working_dir_path,
  control_label = control_label,
  control_correction = TRUE,
  use_positions = TRUE,
  pseudocount = TRUE,
  n_parallel_chains = n_parallel_chains,
  seed = seed
)
```

    Running Lilace on input counts

    Init values were only set for a subset of parameters. 
    Missing init values for the following parameters:
     - chain 1: q, sigma_syn, theta, z, theta_syn, a, b
     - chain 2: q, sigma_syn, theta, z, theta_syn, a, b
     - chain 3: q, sigma_syn, theta, z, theta_syn, a, b
     - chain 4: q, sigma_syn, theta, z, theta_syn, a, b

    To disable this message use options(cmdstanr_warn_inits = FALSE).

    Running MCMC with 4 parallel chains...

    Chain 1 Iteration:    1 / 2000 [  0%]  (Warmup) 

    Chain 2 Rejecting initial value:

    Chain 2   Log probability evaluates to log(0), i.e. negative infinity.

    Chain 2   Stan can't start sampling from this initial value.

    Chain 2 Iteration:    1 / 2000 [  0%]  (Warmup) 
    Chain 3 Iteration:    1 / 2000 [  0%]  (Warmup) 
    Chain 4 Iteration:    1 / 2000 [  0%]  (Warmup) 
    Chain 3 Iteration:  250 / 2000 [ 12%]  (Warmup) 
    Chain 2 Iteration:  250 / 2000 [ 12%]  (Warmup) 
    Chain 1 Iteration:  250 / 2000 [ 12%]  (Warmup) 
    Chain 4 Iteration:  250 / 2000 [ 12%]  (Warmup) 
    Chain 4 Iteration:  500 / 2000 [ 25%]  (Warmup) 
    Chain 2 Iteration:  500 / 2000 [ 25%]  (Warmup) 
    Chain 3 Iteration:  500 / 2000 [ 25%]  (Warmup) 
    Chain 1 Iteration:  500 / 2000 [ 25%]  (Warmup) 
    Chain 4 Iteration:  750 / 2000 [ 37%]  (Warmup) 
    Chain 2 Iteration:  750 / 2000 [ 37%]  (Warmup) 
    Chain 1 Iteration:  750 / 2000 [ 37%]  (Warmup) 
    Chain 3 Iteration:  750 / 2000 [ 37%]  (Warmup) 
    Chain 2 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
    Chain 2 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
    Chain 4 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
    Chain 4 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
    Chain 3 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
    Chain 3 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
    Chain 1 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
    Chain 1 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
    Chain 2 Iteration: 1250 / 2000 [ 62%]  (Sampling) 
    Chain 2 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
    Chain 3 Iteration: 1250 / 2000 [ 62%]  (Sampling) 
    Chain 2 Iteration: 1750 / 2000 [ 87%]  (Sampling) 
    Chain 2 Iteration: 2000 / 2000 [100%]  (Sampling) 
    Chain 2 finished in 211.6 seconds.
    Chain 3 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
    Chain 1 Iteration: 1250 / 2000 [ 62%]  (Sampling) 
    Chain 4 Iteration: 1250 / 2000 [ 62%]  (Sampling) 
    Chain 3 Iteration: 1750 / 2000 [ 87%]  (Sampling) 
    Chain 3 Iteration: 2000 / 2000 [100%]  (Sampling) 
    Chain 3 finished in 262.7 seconds.
    Chain 1 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
    Chain 4 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
    Chain 1 Iteration: 1750 / 2000 [ 87%]  (Sampling) 
    Chain 4 Iteration: 1750 / 2000 [ 87%]  (Sampling) 
    Chain 1 Iteration: 2000 / 2000 [100%]  (Sampling) 
    Chain 1 finished in 369.5 seconds.
    Chain 4 Iteration: 2000 / 2000 [100%]  (Sampling) 
    Chain 4 finished in 421.7 seconds.

    All 4 chains finished successfully.
    Mean chain execution time: 316.4 seconds.
    Total execution time: 421.9 seconds.

    Warning: 7 of 4000 (0.0%) transitions hit the maximum treedepth limit of 10.
    See https://mc-stan.org/misc/warnings for details.

    Warning: 7 of 4000 (0.0%) transitions hit the maximum treedepth limit of 10.
    See https://mc-stan.org/misc/warnings for details.

    Joining with `by = join_by(p_map)`
    Joining with `by = join_by(map)`
    Joining with `by = join_by(variant, rep)`

Read the exported scores table.

The variant annotated as `multiple_aa` is a pathogenic control from the
literature (identified by multiple synonymous edits in tiled data where
the missense variant is not directly sequenced) and is excluded from
this analysis.

``` r
scores <- readr::read_tsv(
  file.path(working_dir_path, "lilace_output", "variant_scores.tsv"),
  show_col_types = FALSE
) |>
  dplyr::filter(type != "multiple_aa")

head(scores)
```

    # A tibble: 6 × 13
      variant type        variant_id mutation_type position sortscore_score
      <chr>   <chr>       <chr>      <chr>            <dbl>           <dbl>
    1 C.1.A   missense_aa C.1.A      missense_aa          1           -4.12
    2 C.1.D   missense_aa C.1.D      missense_aa          1           -3.96
    3 C.1.E   missense_aa C.1.E      missense_aa          1           -4.13
    4 C.1.F   missense_aa C.1.F      missense_aa          1           -4.12
    5 C.1.G   missense_aa C.1.G      missense_aa          1           -3.60
    6 C.1.H   missense_aa C.1.H      missense_aa          1           -4.21
    # ℹ 7 more variables: sortscore_score_column <chr>, effect <dbl>,
    #   effect_se <dbl>, lfsr <dbl>, pos_mean <dbl>, pos_sd <dbl>,
    #   discovery05 <dbl>

Basic summary views.

``` r
scores |>
  dplyr::count(type, sort = TRUE)
```

    # A tibble: 3 × 2
      type            n
      <chr>       <int>
    1 missense_aa   627
    2 nonsense       33
    3 synonymous     33

``` r
scores |>
  dplyr::arrange(lfsr, dplyr::desc(abs(effect))) %>%
  dplyr::select(variant, type, position, effect, effect_se, lfsr, discovery05) %>%
  head(20)
```

    # A tibble: 20 × 7
       variant type        position effect effect_se    lfsr discovery05
       <chr>   <chr>          <dbl>  <dbl>     <dbl>   <dbl>       <dbl>
     1 Y.30.*  nonsense          30 -1.66      0.408 0.00025          -1
     2 T.23.H  missense_aa       23  1.00      0.409 0.0005            1
     3 K.18.E  missense_aa       18 -1.63      0.417 0.001            -1
     4 C.1.F   missense_aa        1 -1.55      0.420 0.001            -1
     5 Y.30.W  missense_aa       30  0.950     0.420 0.001             1
     6 Q.20.M  missense_aa       20  0.892     0.408 0.00125           1
     7 E.4.V   missense_aa        4  0.821     0.406 0.002             1
     8 R.22.P  missense_aa       22 -1.49      0.421 0.0035           -1
     9 N.21.E  missense_aa       21  0.755     0.426 0.00725           1
    10 K.18.M  missense_aa       18 -1.41      0.413 0.00775          -1
    11 T.23.M  missense_aa       23  0.723     0.415 0.0085            1
    12 K.18.Q  missense_aa       18 -1.42      0.417 0.00875          -1
    13 F.10.K  missense_aa       10 -1.39      0.414 0.00875          -1
    14 N.21.M  missense_aa       21  0.720     0.412 0.009             1
    15 D.15.K  missense_aa       15 -1.38      0.398 0.01             -1
    16 F.10.P  missense_aa       10 -1.35      0.408 0.0112           -1
    17 S.25.F  missense_aa       25  0.701     0.412 0.0112            1
    18 R.22.T  missense_aa       22 -1.37      0.415 0.0125           -1
    19 C.1.Q   missense_aa        1 -1.32      0.405 0.0125           -1
    20 K.33.T  missense_aa       33 -1.35      0.412 0.0135           -1

Optional diagnostic plots from Lilace.

``` r
lilace::lilace_score_density(
  scores,
  savedir = working_dir_path,
  score.col = "effect",
  name = "score_histogram",
  hist = TRUE,
  scale.free = TRUE,
  savepdf = FALSE
)
```

![](Lilace_files/figure-commonmark/effect%20size%20heatmap-1.png)

``` r
heatmap_scores <- scores |>
  dplyr::mutate(
    wildtype = stringr::str_extract(variant, "^[^.]+"),
    mutation = stringr::str_extract(variant, "[^.]+$")
  )

lilace::lilace_score_heatmap(
  heatmap_scores,
  savedir = working_dir_path,
  score.col = "effect",
  name = "score_heatmap",
  x.text = 9,
  y.text = 10,
  seq.text = 3,
  savepdf = FALSE
)
```

![](Lilace_files/figure-commonmark/effect%20size%20heatmap-2.png)

``` r
lilace::lilace_score_heatmap(
  heatmap_scores,
  savedir = working_dir_path,
  score.col = "discovery05",
  name = "discovery05_heatmap",
  x.text = 9,
  y.text = 10,
  seq.text = 3,
  savepdf = FALSE
)
```

![](Lilace_files/figure-commonmark/discovery%20heatmap-1.png)

## Compare sortscore and Lilace Scores

``` r
library(ggplot2)

scores_filtered <- scores[scores$type != "multiple_aa", ]

corr <- cor.test(x=scores_filtered$sortscore_score, y=scores_filtered$effect, method='spearman', exact=FALSE)
# Add n and correlation to the plot
complete_scores <- scores_filtered[complete.cases(scores_filtered$sortscore_score, scores_filtered$effect),]

label <- sprintf(
  "Spearman's \u03c1 = %.2f\nn = %d",
  unname(corr$estimate),
  nrow(complete_scores)
)

scatterplot <- ggplot(scores) +
  aes(x = sortscore_score, y = effect, color = type) +
  geom_smooth(method=lm, se = FALSE, color = 'black', weight = 3) +
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
    title = "sortscore and Lilace scoring shows strong concordance (GLI2 DMS positions 518-550)",
    x = "Sortscore score",
    y = "Lilace Effect",
    color = "Type"
    )
  print(scatterplot)
```

![](Lilace_files/figure-commonmark/sortscore%20comparison-1.png)
