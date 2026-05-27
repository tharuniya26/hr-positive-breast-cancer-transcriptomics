# Data provenance and redistribution

Retrieved on **2026-10-16** from public resources.

| Artifact | Source | SHA-256 | Use |
|---|---|---|---|
| GSE45827 series matrix | [NCBI GEO](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE45nnn/GSE45827/matrix/GSE45827_series_matrix.txt.gz) | `741184cd3b3354ed3ea45fc6a428ed36401cadc20505e5cc5fb381d8aed0ad35` | 155 samples, 29,873 probe sets |
| GPL570 annotation, 2016-08-09 | [NCBI GEO](https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz) | `d7cd44352127b1e34f3a720ebea86093ef255a38f1612a85a2962b71bde8f394` | Probe-to-gene mapping |
| Reactome GMT, release 97 | [Official download](https://www.reactome.org/download/current/ReactomePathways.gmt.zip) | `89983d5c1f0af11c52edfeee7323eb425580ac6281d387a528562ab1787ce56b` (extracted GMT) | ORA pathway definitions |

The official Reactome site identifies [release 97](https://reactome.org/about/news/295-v97-released) as the current June 2026 release. The extracted official GMT checksum matches the independently recorded V97 checksum in the [THEMA provenance manifest](https://github.com/AviyahP/thema). The download URL is a moving `current` endpoint, so the script validates the exact SHA-256 and fails if the file changes. The source server presented an expired certificate in this environment; the download script falls back to an unverified TLS connection only for that certificate error, then requires the pinned checksum. The downloaded ZIP has SHA-256 `8c1dbc8578431da5d2d5118262718c60b553a9be3398e93658daa069e4a9afd4`. Raw GMT and ZIP are Git-ignored.

NCBI GEO is publicly accessible; dataset-level redistribution rights for original raw and processed files were not separately established. Thus raw matrix, annotation and GMT are Git-ignored. Derived sample IDs, statistical summaries, code and figures can be reviewed in the public repository. Users download upstream files with `python -m scripts.download_data` and the GMT link above. Publication of the original academic documents has not been authorized; both remain in `private_sources/`.

No identifying patient data were used. Source subtype labels are original GEO annotations. The original article is Gruosso T et al. (2016), *EMBO Molecular Medicine* 8:527–549, DOI [10.15252/emmm.201505891](https://doi.org/10.15252/emmm.201505891).
