
# NMAforest

<!-- badges: start -->
<!-- badges: end -->

**NMAforest** is an R package for generating detailed forest plots in network meta-analysis (NMA). It visualizes direct, indirect, and network meta-analysis treatment effects, along with study- and path-level contribution proportions. The visualization is based on the evidence flow decomposition method by Papakonstantinou et al. (2018).


## Acknowledgments

This package relies on key infrastructure from the `netmeta`, `igraph`, and `ggplot2` R packages.  
It also adapts methods and code presented by **Papakonstantinou et al. (2018)** for evidence flow decomposition in network meta-analysis, including the `comparisonStreams()` function from the [flow_contribution GitHub repository](https://github.com/esm-ispm-unibe-ch/flow_contribution).

## Installation

The stable release of **NMAforest** can be installed from CRAN:

```r
install.packages("NMAPropForest")
```


## Required Data Format

The following columns are required (either with these names or specified via function arguments):

| Column     | Required For       | Description                                                  | Type                  |
|------------|--------------------|--------------------------------------------------------------|-----------------------|
|  `treat`    | All analyses        | Treatment label for each arm                                 | character or factor   |
|  `event`    | Binary outcomes     | Number of events in the arm                                  | numeric               |
|  `n`        | All analyses        | Sample size in each arm                                      | numeric               |
|  `mean`, `sd` | Continuous outcomes | Mean and standard deviation (used instead of `event`)        | numeric               |
|  `study`    | All analyses        | Study label or grouping variable                             | character or numeric  |
|  `study_id`       | Optional            | Unique numeric study identifier (auto-generated if missing)  | integer               |

##
**Note:** We recommend that users include an explicit `study_id` column where each value uniquely corresponds to a study label in the `study` column.  
If the `study_id` column is not present in the dataset, the function will automatically generate one and return the updated data frame with this column added.
