# Breast cancer SPP1 fibrovascular context

This repository provides the reproducible analysis code, frozen specifications,
software requirements, provenance records, derived statistical tables, and
publication figures supporting the study:

> **Cross-cohort single-cell and spatial analyses reveal a fibrovascular,
> immune-infiltrated context for SPP1-associated macrophages in breast cancer**

Authors: Yin Xiaohui, Liu Ruihong, Zhang Ganlin, Zhu He, Guo Yinuo, and Zhang Qing.

## Scope and principal finding

The study integrates four independently processed single-cell cohorts, three
spatial transcriptomic series, longitudinal neoadjuvant immunotherapy data, and
a clinically annotated 94-patient TNBC cohort. A locked 13-gene macrophage
program recurred across independent datasets and occupied reproducible
fibrovascular, CD8/cytotoxic-infiltrated tissue neighborhoods. The available
data do not establish universal tumor induction, physical immune exclusion,
treatment resistance, or adverse prognosis.

## Repository contents

- `config/`: frozen genes, covariates, thresholds, null models, and clinical specifications.
- `scripts/R/` and `scripts/python/`: analysis and validation stages.
- `environment/`: Python and R dependency specifications.
- `results/tables/`: machine-readable statistical outputs used for the figures and claims.
- `results/figures/`: final main and supplementary figures in PNG format.
- `provenance/`: public-input checksums and software-version records.
- `data_manifest.tsv`: accession numbers, data roles, access links, and audited sample counts.
- `analysis_preregistration.md`: frozen analysis gates and dated amendments.
- `REPRODUCIBILITY.md`: complete execution order and verification instructions.

## Data availability

All inputs are publicly available. Raw datasets are not redistributed here.
Accession numbers and source URLs are listed in `data_manifest.tsv` and
`DATA_AVAILABILITY.md`. Download the public inputs and place them at the paths
documented in `REPRODUCIBILITY.md` before running the full pipeline.

## Reproduction levels

Install the Python dependencies:

```bash
python3 -m venv .venv
.venv/bin/pip install -r environment/requirements.txt
```

R 4.5.3 and the packages listed in `environment/R_requirements.txt` were used.
The complete stage order, computational requirements, and expected outputs are
documented in `REPRODUCIBILITY.md`.

The machine-readable tables permit verification and figure reconstruction
without rerunning the largest single-cell and spatial preprocessing stages.

## Reproducibility safeguards

- Outcome labels were withheld during discovery and gene locking.
- Patients or tissue sections, rather than cells or spots, are the inferential units.
- Spatial null models preserve pathology or malignant-expression composition.
- Frozen YAML files record primary signatures, thresholds, covariates, and prohibited inferences.
- Negative clinical boundary tests and leave-one-gene-out analyses are retained.

## Citation

Please cite the associated article after publication. A versioned Zenodo archive
and DOI will be added when the repository release is deposited.

## License

No reuse license has yet been selected. Copyright remains with the authors until
a license file is added.
