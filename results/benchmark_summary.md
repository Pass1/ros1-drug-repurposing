# ROS1 Benchmark Summary

This summary is generated from the benchmark registry in `data/benchmarks/`.

## `ros1_known_actives_background`

- Description: Known ROS1 inhibitors recovered from the screened library against the remainder of the scored library as presumed background.
- Mode: `active_vs_background`
- Target: ROS1
- Mutation/state: G2032R
- Declared labeled rows: 15
- Limitation note: Provisional benchmark only: negatives are presumed background rather than validated inactives or designed decoys.

### Metrics

| score_column    |   n_entities |   n_actives |   n_negatives |   n_explicit_negatives |   n_background |   roc_auc |   average_precision |    ef1 |   nef1 |    ef5 |   nef5 |   ef10 |   nef10 |   median_active_rank |   median_active_percentile |   hits_top_10 |   hits_top_25 |   hits_top_50 |   hits_top_100 |
|:----------------|-------------:|------------:|--------------:|-----------------------:|---------------:|----------:|--------------------:|-------:|-------:|-------:|-------:|-------:|--------:|---------------------:|---------------------------:|--------------:|--------------:|--------------:|---------------:|
| g2032r_score    |         2794 |           9 |          2785 |                      0 |           2785 |     0.927 |               0.056 | 11.087 |  0.111 | 11.087 |  0.556 |  6.652 |   0.667 |                  114 |                      0.041 |             1 |             1 |             2 |              3 |
| composite_score |         2794 |           9 |          2785 |                      0 |           2785 |     0.942 |               0.077 | 22.175 |  0.222 | 11.087 |  0.556 |  7.761 |   0.778 |                   72 |                      0.026 |             1 |             2 |             2 |              5 |

### Active Ranks: `composite_score`

| representative_name       |   score |   rank | screen_fraction   |
|:--------------------------|--------:|-------:|:------------------|
| entrectinib               | -11.162 |      4 | 0.14%             |
| Lorlatinib?(PF-6463922)   | -10.735 |     12 | 0.43%             |
| taletrectinib             |  -9.52  |     63 | 2.25%             |
| Ceritinib (LDK378)        |  -9.447 |     67 | 2.40%             |
| zidesamtinib              |  -9.4   |     72 | 2.58%             |
| repotrectinib             |  -8.891 |    194 | 6.94%             |
| Crizotinib hydrochloride  |  -8.772 |    250 | 8.95%             |
| Cabozantinib (BMS-907351) |  -8.62  |    321 | 11.49%            |
| Brigatinib (AP26113)      |  -8.33  |    518 | 18.54%            |

### Active Ranks: `g2032r_score`

| representative_name       |   score |   rank | screen_fraction   |
|:--------------------------|--------:|-------:|:------------------|
| entrectinib               | -10.662 |      4 | 0.14%             |
| Lorlatinib?(PF-6463922)   |  -9.735 |     40 | 1.43%             |
| zidesamtinib              |  -9.1   |     84 | 3.01%             |
| taletrectinib             |  -9.02  |    104 | 3.72%             |
| Ceritinib (LDK378)        |  -8.947 |    114 | 4.08%             |
| Cabozantinib (BMS-907351) |  -8.62  |    218 | 7.80%             |
| Crizotinib hydrochloride  |  -8.472 |    286 | 10.24%            |
| Brigatinib (AP26113)      |  -8.33  |    381 | 13.64%            |
| repotrectinib             |  -7.891 |    664 | 23.77%            |

