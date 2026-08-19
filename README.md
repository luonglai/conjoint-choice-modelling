# Choice-Based Conjoint Analysis: Skincare Attribute Preferences, Willingness-to-Pay, and Segmentation

A discrete choice experiment estimating what drives consumer choice of a premium
skincare serum, first with a pooled (aggregate) model, then with Hierarchical
Bayes (HB) to recover individual-level preferences and latent segments.

## Business question

Which product attributes most influence a consumer's choice of skincare serum,
what is each attribute worth in price terms, and is the market really one 
audience, or several with different preferences?

Specifically:

- Does a **dermatologist-recommended** positioning justify a price premium over
  a mass-prestige or luxury-heritage story?
- Which **active ingredient claim** (retinol, vitamin C, peptide complex) drives
  the most preference?
- Does **packaging format** (standard bottle vs. airless pump) matter enough to
  justify the cost of switching?
- Do different consumers actually want different things, and if so, does a
  single "average" model hide that, or mislead a launch decision?

## Data

**This is simulated data, not a real consumer panel.** Real conjoint studies
use proprietary respondent data that isn't publicly available for a portfolio
project. To build and demonstrate the full pipeline honestly, I instead:

1. Defined **three latent preference segments** with genuinely different,
   known part-worth utilities (see table below). A single shared set of 
   "true" preferences gives an aggregate model nothing to get wrong and
   HB nothing real to recover. Three segments let the project show what each
   method actually can and cannot see.
2. Assigned 300 simulated respondents to segments (35% / 30% / 35%), each with
   individual-level noise around their segment's mean utilities.
3. Generated 10 choice tasks per respondent x 3 product alternatives + a
   "None of these" option (12,000 rows total), with attribute levels assigned
   at random to each alternative.
4. Simulated each respondent's choice using a random-utility-maximisation
   process (their own true utility + Gumbel-distributed noise). This reflects 
   the same data-generating assumption underlying the multinomial logit model itself.

Because the true parameters are known at both the population and individual
level, the analysis can do something a real study cannot, which is to **validate 
what each model recovers correctly, and what it gets wrong.**

| Attribute | Levels |
|---|---|
| Positioning | Mass-Prestige (ref.), Dermatologist-Recommended, Luxury Heritage |
| Active ingredient | Retinol (ref.), Vitamin C, Peptide Complex |
| Format | Standard Bottle (ref.), Airless Pump |
| Price | €45 / €65 / €85 / €110 (continuous) |

### The three simulated segments

| Segment | Share | Defining preference |
|---|---:|---|
| Clinical/Dermatologist-Focused | 35% | Strong pull to dermatologist-recommended positioning and peptide-complex claims |
| Luxury Heritage Seekers | 30% | Strong pull to luxury-heritage positioning and premium format; low price sensitivity |
| Price-Sensitive Pragmatists | 35% | Weak preferences on positioning/ingredient; much more price-sensitive |

## Method

- **Design:** choice-based conjoint (CBC), 3 alternatives + "None" per task,
  attribute levels randomised per alternative (a production study would use a
  D-efficient/orthogonal design; see *Limitations*).
- **Step 1 — pooled estimation:** multinomial logit via `survival::clogit()`
  (conditional logistic regression on matched strata), giving one set of
  "average" part-worths for the whole market.
- **Step 2 — Hierarchical Bayes:** `bayesm::rhierMnlRwMixture()`, the standard
  R implementation of the HB-MNL method from Rossi, Allenby & McCulloch's
  *Bayesian Statistics and Marketing* — the reference text for HB conjoint
  estimation in market research. This estimates a **separate set of
  part-worths for every respondent**, shrunk toward a population distribution
  via a Bayesian hierarchical prior (single-component / random-coefficients
  specification. Segments are recovered afterward by clustering, not assumed
  by the model).
- **Step 3 — segmentation:** k-means (k=3) on the individual HB part-worths,
  then validated against the true simulated segment labels.

Reproduce with (R 4.3, packages: `survival`, `bayesm`, `dplyr`, `tidyr`,
`ggplot2`, `scales`):

```bash
Rscript scripts/01_simulate_data.R
Rscript scripts/02_estimate_model.R
Rscript scripts/03_business_outputs.R
Rscript scripts/04_hb_estimation.R      # ~20 seconds for 8,000 MCMC draws
Rscript scripts/05_hb_visualization.R
```

## Results

### Finding 1 — the pooled model is biased under heterogeneity, not just "unable to segment"

This is the headline finding, and it's a real, well-documented phenomenon: 
when a population has genuinely different preferences, a pooled
multinomial logit doesn't just fail to see the segments, its *aggregate*
estimates can be biased, because average choice probability across
heterogeneous individuals is not the same as choice probability evaluated at
the average individual (a nonlinearity in the logit link, sometimes called
aggregation bias). This shows up clearly here:

| Parameter | True population avg. | Pooled clogit | HB population mean |
|---|---:|---:|---:|
| Dermatologist-Recommended | 0.790 | 0.880 | 1.113 |
| Luxury Heritage | 0.460 | 0.615 | 0.699 |
| Vitamin C | 0.290 | 0.272 | 0.357 |
| Peptide Complex | 0.475 | 0.305 | 0.412 |
| Airless Pump | 0.320 | 0.376 | 0.459 |
| Price (per €) | -0.026 | **-0.017** | -0.042 |
| None (baseline) | -0.920 | **-0.047** | -1.905 |

The price coefficient and the "None" baseline are hit hardest: the pooled
model underestimates price sensitivity by roughly a third, and estimates a
"None" baseline nearly 20x too weak. Since willingness-to-pay is calculated
as *attribute utility ÷ price coefficient*, this bias directly inflates the
pooled model's WTP estimates (see Finding 2). This is an example of why a
launch decision based only on a pooled model can be misleading. Note that HB
doesn't perfectly recover the "None" parameter either (see *Limitations*).
Both approaches struggle with it, just in different directions.

### Finding 2 — aggregate view: attribute importance and WTP (pooled model)

![Attribute importance](outputs/attribute_importance.png)

Positioning (56.4%) dominates the pooled model's importance ranking, ahead of
format (24.1%) and active ingredient (19.6%).

![Willingness to pay](outputs/willingness_to_pay.png)

The pooled model estimates a dermatologist-recommended positioning is worth
**+€51** over a mass-prestige story. Given Finding 1, this is likely inflated
by the underestimated price coefficient. A real study would treat this
pooled WTP figure with real caution, which is also the argument for
running HB before finalizing a pricing recommendation.

### Finding 3 — HB recovers real, separable heterogeneity

![HB heterogeneity by segment](outputs/hb_heterogeneity_by_segment.png)

Individual-level HB part-worths, grouped by each respondent's *true* (unknown
to the model) segment, separate cleanly, especially on price sensitivity,
where the three segments barely overlap. This is the core value HB adds over
a pooled model: it recovers structure that the pooled model cannot see at all.

**Individual-level recovery quality** (correlation between each respondent's
HB-estimated and true part-worth, across all 300 respondents):

| Parameter | Correlation | Mean abs. error |
|---|---:|---:|
| Price (per €) | 0.808 | 0.025 |
| Luxury Heritage | 0.631 | 0.426 |
| Dermatologist-Recommended | 0.615 | 0.502 |
| Peptide Complex | 0.442 | 0.324 |
| Airless Pump | 0.252 | 0.301 |
| Vitamin C | 0.240 | 0.304 |
| None (baseline) | 0.305 | 1.028 |

Price is recovered best (it's the only continuous variable, observed on every
row, giving the model the most information per respondent). The categorical
attribute dummies are recovered more modestly, reflecting a realistic and expected
result given only 10 choice tasks per respondent (see *Limitations*).

### Finding 4 — segment recovery from HB + k-means

![Segment profiles](outputs/hb_segment_profiles.png)

Clustering the individual HB part-worths (k-means, k=3) and comparing against
the true simulated segment labels:

| True segment | Recovered as Cluster 1 | Cluster 2 | Cluster 3 |
|---|---:|---:|---:|
| Clinical/Dermatologist-Focused | **80** | 8 | 20 |
| Luxury Heritage Seekers | 10 | **71** | 13 |
| Price-Sensitive Pragmatists | 25 | 8 | **65** |

**Overall recovery accuracy: 72%** (best cluster-to-segment matching). The
recovered segment profiles are directionally correct and business-usable.
Cluster 1 clearly over-indexes on dermatologist positioning and peptide
complex, Cluster 2 on luxury heritage, and Cluster 3 is the flattest/most
price-driven profile. However, recovery isn't perfect, which is an
expected result given the sample design (see *Limitations*).

## Limitations

- **Simulated, not real, respondents.** Results describe a simulated market
  built to validate the pipeline, not actual consumer preferences.
- **Random rather than efficient design.** Choice sets were generated by
  random assignment rather than a D-efficient/orthogonal design, which would
  reduce standard errors for the same sample size in a real study.
- **Only 10 choice tasks per respondent.** This is on the lean side for CBC
  (production studies often use 12–20). It's a real driver of the modest
  individual-level correlations in Finding 3 and the imperfect segment
  recovery in Finding 4. HB needs enough within-respondent variation to
  separate a person's true preferences from noise, and 10 tasks is a
  meaningful constraint on that.
- **The "None" baseline is hard for both models to pin down.** Pooled clogit
  underestimates its magnitude roughly 20-fold; HB overestimates it relative
  to the true population average. Alternative-specific constants like this are
  known to be harder to identify precisely than attribute part-worths, and
  this project doesn't fully resolve why the two methods err in opposite
  directions — worth flagging honestly rather than picking whichever number
  looks better.
- **HB run with a single mixture component (ncomp = 1),** i.e. continuous
  individual heterogeneity via one multivariate normal, with segments
  recovered afterward via k-means. This mirrors common commercial practice 
  (HB first, cluster after) but means the model itself doesn't "know" there are 
  3 segments, that structure is recovered, not assumed.
- **No formal MCMC convergence diagnostics** (trace plots, Gelman-Rubin,
  effective sample size) were run. Burn-in was fixed at 50% of draws as a
  simple, conservative default rather than something empirically checked.
- **Linear price utility and no interaction effects,** same as the earlier
  pooled-only version of this project.

## Repo structure

```
├── data/
│   ├── simulated_choices.csv          # respondent-level choice data (long format)
│   ├── respondent_true_utilities.csv  # ground-truth individual-level utilities + segment
│   ├── true_utilities.rds             # population-weighted-average ground truth
│   ├── true_segment_utilities.rds     # ground-truth utilities per segment
│   └── true_segment_shares.csv        # true segment sizes
├── scripts/
│   ├── 01_simulate_data.R             # 3-segment experiment design + choice simulation
│   ├── 02_estimate_model.R            # pooled clogit estimation + validation vs. truth
│   ├── 03_business_outputs.R          # pooled-model importance, WTP, preference share
│   ├── 04_hb_estimation.R             # Hierarchical Bayes MNL (bayesm) + validation + k-means
│   └── 05_hb_visualization.R          # heterogeneity and segment-profile charts
├── outputs/                           # all CSVs, PNGs, and saved model objects
└── README.md
```
