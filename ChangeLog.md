# ChangeLog

## 3.0.0 — release preparation, 2026-10-22

Addresses [issue #4](https://github.com/BackofenLab/CopraRNA/issues/4), based on and
targeting `issue-1-python3-bugfixes`. This release completes the transition from
the original 2.1.3 scripts to the Python-only CopraRNA2 runtime.

- **Correct genomic accessibility boundaries ([issue #6](https://github.com/BackofenLab/CopraRNA/issues/6)).**
  Extend native target inputs on both sides by the target accessibility base-pair
  span (`tAccL`, 100 nt by default), honoring shared `accL` and the native zero/window
  semantics. Wrap circular genomes and stop at physical linear ends; preserve
  strand, spliced CDS and short-truncated-UTR selection conventions. Restrict
  interactions with `--tRegion` and native indexing so all reported positions,
  seeds and structures refer to the requested sequence. Trim target CSV sequences
  and normalize `E_norm` before filtering/statistics; reuse context for web exports.
  Keep context sidecars in prediction signatures and retained job artifacts, and
  invalidate old predictions. This intentionally changes accessibility energies and
  potentially rankings; independently captured fixtures for unchanged extraction,
  statistics and reports remain untouched. Test explicit boundary contracts and
  real IntaRNAsTar interactions, including the effect of genomic folding context.
- **Simplify missing report fields and update DAVID access.** Following review of
  [PR #7](https://github.com/BackofenLab/CopraRNA/pull/7) and the convention in
  [upstream issue #22](https://github.com/PatrickRWright/CopraRNA/issues/22), write
  unavailable values as empty cells in published result and web tables instead of
  the internal `NA` marker. This makes CSV imports unambiguous without changing
  statistical missing-value handling. Authenticate DAVID requests with the shared
  `rna@informatik.uni-freiburg.de` account instead of the former personal address.
- **Remove static organism catalogs.** `CopraRNA_available_organisms.txt` supplied
  display names and an input allow-list; `kegg2refseqnew.csv` supplied assembly
  membership and arbitrary organism codes. Neither is necessary: resolve current
  complete RefSeq assemblies and every chromosome/plasmid through NCBI
  E-utilities and Datasets v2, and obtain names from metadata/GenBank. No global
  table generation, KEGG identifiers or release-bound availability updates remain.
- **Versioned, reusable resolution.** Accept `NC_`, `NZ_` and `GCF_` IDs with optional
  versions. Cache validated metadata under the shared genome database using locks
  and atomic publication. Reject ambiguous/missing/draft assemblies, incomplete
  sequence reports and multiple inputs representing the same assembly. Follow
  all report pages. Use short, deterministic job-local species codes to respect
  DomClust's native name limits; snapshot those codes and exact genome identities
  in `genome_metadata.json`. Resolve assembly/plasmid/versioned FASTA headers to
  primary-chromosome IDs in a separate working copy for alignment and reporting.
  Scientific table labels stay unversioned.
- **Inspect, prefetch and refresh genomes.** Add `coprarna genomes`, `--download`
  and `--refresh`. Cache exact-version gzip GenBank files, fetch accessions directly
  without arbitrary first-hit searches, and retain old files on failed refresh.
  Pace NCBI requests with timeouts and bounded transient retries. Reuse shared
  compressed paths directly; prefetch metadata plus genomes for offline jobs.
- **Remove the original code base.** Delete the Perl/R/Python 2/shell/PhantomJS
  programs, old HTML wrappers, per-script mapping documents, migration inventory,
  catalog maintenance code/resources, optional interpreter environments and
  legacy source-capture tooling. Remove defect-emulation branches, script-name
  aliases and `coprarna status`. Keep independently captured scientific data and
  provenance with checksum protection; tests no longer execute original sources.
- **Refresh all 17 reference genomes.** Download current NCBI GenBank annotations
  and sequence versions, replacing the old bzip2 snapshots with versioned gzip
  files. The archive manifest records source URLs, retrieval time, record dates,
  lengths, CDS counts and SHA-256 hashes. A Python refresh tool stages and checks
  all records before publishing them. Reference tests check identities and all
  three target-region modes without invoking legacy scripts.
- **Publish one complete guide.** Replace the old README with the expanded Python
  guide covering installation, input choice, examples, supported assemblies,
  database updates, native settings, every prediction option, outputs, batches,
  migration, testing and biological citations. Remove README.2.md. Update the
  architecture ledger, wheel/sdist contents, package/CLI version and job-option
  version consistently to 3.0.0; include the current HTML report resource.

The validated numerical dependency/kernel pins are unchanged. Captured GEV,
correlation, aggregation, consensus and Jalview expectations are not regenerated.
Native integration verifies the current scientific interfaces, complete predictions,
file/inline equivalence, shared caches, preserved inputs and vector reports.
Live NCBI access was exercised for the reference refresh; offline request tests do
not certify future NCBI or DAVID availability. Cached metadata remains fixed until
explicit refresh; NCBI annotation updates can change predictions even when the
sequence accession version is unchanged. Use fresh job directories for new data.

## 3.0 development history — issues #1 and #3

The `3.0.0.dev0` development series addressed
[issue #1](https://github.com/BackofenLab/CopraRNA/issues/1).
The rewrite and subsequent corrections are reviewed separately in
[PR #2](https://github.com/BackofenLab/CopraRNA/pull/2) and
[PR #3](https://github.com/BackofenLab/CopraRNA/pull/3).

### Rewrite foundation — PR #2

The source workflow identifies itself as **CopraRNA 2.1.3**. Its controller and
auxiliary Perl, R, Python 2, shell and JavaScript programs were translated into a
Python 3 package. During that rewrite, original tracked files and bundled biological data were kept
byte unchanged; the Python runtime never invoked those scripts. The initial
rewrite provided counterparts for all 41 scripts, with adjacent mapping documents
and a machine-readable migration inventory. At the end of PR #3, the runtime retained 39
active counterparts and marked two obsolete scripts retired.

The architecture separates the command-line interface and isolated batch runner
(`cli.py`, `pipeline.py`) from the scientific stages:

- Genome retrieval, region/protein extraction and annotation: `ncbi.py`,
  `genome.py`, `annotation.py`; `homology.py` coordinates native protein clustering.
- Interaction prediction and statistical inference: `intarna.py`, `energy.py`,
  `combination.py`, `statistics.py`, `numerics.py`; native IntaRNA, BLAST+, DomClust
  and phylogeny tools remain scientific executables.
- Alignment and consensus mapping: `alignments.py`, `conservation.py`, with
  persisted R-compatible matrices and Jalview annotation archives.
- Reports and service adapters: `plots.py`, `phylogeny.py`, `enrichment.py`,
  `david.py`, `webserver.py`; Python plotting replaces R graphics, and the optional
  browser export helper replaces PhantomJS.
- Shared contracts: `helpers.py` dispatches active legacy-name aliases,
  `rcompat.py` preserves retained table/object conventions, and the reviewed
  `io.py`/`parameters.py` centralize compression and native parameter handling.

Independent jobs run in separate directories with preflight validation, bounded
concurrency, ordered result records and per-job failures. Packaging includes the
organism tables and report resources, with build/install checks outside the source
checkout. Scientific fixtures were captured from original programs rather than
regenerated from translated output. Retained comparisons cover text bytes,
statistical fits, matrix contents, plot geometry and Jalview exports. Whole native
workflows additionally exercise real external-tool integration.

### Correctness fixes — PR #3

- **Cluster p-value headers:** name the p-value column explicitly. R-style implicit
  first-row/header interpretation could drop or mislabel an observation.
- **Multi-record GenBank regions:** extract a feature from its own record sequence
  and topology. The old path reused the first record's sequence for later records.
- **Ambiguous reverse-strand bases:** use the complete IUPAC reverse complement,
  including lowercase symbols; a four-base transliteration left ambiguous bases
  incorrect in extracted flanks.
- **Organism-code lookup:** derive an accession from the GenBank basename, including
  compressed suffixes. Absolute directory prefixes previously broke mapping.
- **Suboptimal interactions:** retain the second-best sorted interaction rather
  than overwriting it with later, worse interactions for the same target.
- **Enrichment JSON:** escape quotes and backslashes using JSON serialization so
  biological labels cannot produce invalid report data.
- **Historical correlation scoring:** corrected histogram-density scoring in the
  earlier follow-up. The affected obsolete mode, its implementation and its tests
  were subsequently removed under the current review scope.
- **DomClust argument order:** put options before positional files; the executable
  could otherwise ignore them and fail to emit the expected output.
- **DomClust output interface:** replace removed table mode 5 with supported mode 9.
  Reconstruct the cluster table in Python, preserving species order, multiple
  homologs, missing species and split-domain suffixes.
- **Lowercase CDS plots:** accept the documented `cds` setting without relying on
  an uninitialized plot-center variable.
- **Short prediction lists:** bound plotting to available rows and include partial
  final pages, avoiding padding failures and out-of-range indexing.
- **Numerical reproducibility:** retain R-compatible reductions, simplex fitting
  and one-column QR rounding. Select the validated OpenBLAS Haswell kernel on
  supported Linux x86-64/AVX2 systems before NumPy import, respecting explicit
  caller overrides. CPU/BLAS differences outside that setup remain a parity limit.

PR #3 supported explicit thread-local reproduction of retained defects for
historical helper comparisons. That comparison mode was removed for 3.0.0.

### Review-driven behavior and integration changes — PR #3

Following [Martin Raden's review](https://github.com/BackofenLab/CopraRNA/pull/3#issuecomment-6324702688):

- Simplify prediction to the current comparative algorithm. Remove the obsolete
  mode's dispatch, sampling/correlation implementation, options, fixtures and
  user documentation; mark its two original scripts retired in the inventory.
- Require **IntaRNA >=3** and default to **IntaRNAsTar**. Parse native CSV rather
  than requesting the old output emulator. The frontend's duplex text is formatted
  from native CSV in Python. Treat printed native errors as failures even when
  IntaRNA returns exit status zero.
- Accept `--intarna-params FILE` or `--intarna-params-text CONTENTS`. Forward the
  same file to primary and web predictions, preserve a settings snapshot, and
  remove temporary inline files after success or failure. Remove duplicated
  native window/distance settings from the CopraRNA CLI.
- Check current dependency releases and interfaces. Adopt BLAST+ 2.17.0
  (`makeblastdb`, `blastp` format 6), Biopython 1.88, NumPy 2.5.4, SciPy 1.18.1,
  pandas 3.0.6 and Matplotlib 3.11.2 in the Python 3.12 environment. Keep compatible
  NumPy/SciPy pins for Python 3.11. Verify IntaRNA 3.4.1, ViennaRNA 2.7.2 and the
  supported DomClust, alignment and phylogeny interfaces. Install Matplotlib's
  wheel to avoid incompatible conda ICU rendering requirements.
- Write primary and web IntaRNA CSV directly as gzip and process compressed
  optimal/suboptimal tables without plain intermediate copies. Keep raw prediction
  archives and filtering audits after cleanup.
- Add a compressed shared genome database via `--genome-db`. Request compressed
  HTTP responses or compress while streaming; parse gzip directly, publish cache
  files atomically, and lock each accession to share concurrent downloads. Reuse
  cached paths directly rather than copying genomes into every job directory.
- Preserve the original input FASTA in place and byte unchanged. Normalize into a
  separate working file, protect the original during cleanup, and reject names
  that collide with required generated files before work starts. Cleanup also
  protects supplied parameter files and their reproducibility snapshot, including
  when enrichment cleanup uses broad intermediate-file patterns.
- Remove stacking-only helices spanning at least 90% of either RNA's effective
  sequence/accessibility-window length **before** fitting energy distributions.
  Apply the same filter to native web exports before selecting a retained
  interaction, so discarded best hits cannot reappear in reports.
  Record removed IDs/coordinates and retained counts. Reject insufficient variable
  data for GEV fitting rather than continuing with an unusable distribution.
- Validate native prediction reuse with a hash of tool version, personality,
  input contents, parameter contents and filtering policy. Prevent stale files
  from silently bypassing new scientific settings.
- Render `Phylogeny/compatible.treefile.svg` in Python. Export every diagram as
  PDF and SVG, use vector heatmap cells, and include all extended-list pages in
  both formats. Produce SVG thumbnails and Python functional-enrichment graphics;
  production reports no longer require a browser or format converter.
- Preserve `NA` comparative p-values in the web table. The more selective sTar
  predictions can leave a valid per-organism interaction without enough homologs
  for a comparative score; report generation must not crash on that case.
- Add function and inline documentation covering file lifetimes, scientific matrix
  and coordinate conventions, numerical choices and native-tool constraints.
  Record maintenance requirements in `AGENTS.md`; update the archive README,
  user guide, environment files, mappings and CI.

### Validation and release limits

Retained scientific/helper expectations still come from original sources.
New behavior is checked with focused contracts and complete native runs, including
file/inline parameter equivalence, shared gzip genomes, input preservation, real
BLAST+ clustering, vector reports and installed-package execution. The reviewed
sTar defaults and artifact filter intentionally change predictions; historical
whole-workflow equality is no longer the acceptance criterion.

Live NCBI/DAVID service availability is not certified by offline request tests.
The conda environment targets Linux; it is not a published Bioconda recipe or a
portable lockfile. PR #3 retained original sources and separate optional Perl/R oracle environments
for auditing. Version 3.0.0 removes them from the distribution; captured results
and provenance remain in the test suite, and the sources remain in Git history.


### Follow-up completeness audit — 2026-10-22

Rechecked the original issue, Martin's full review and all review/inline-comment
endpoints. The scientific changes above remain covered. The audit additionally:

- Extends input preservation to generated-output basename collisions and symbolic
  aliases, and prevents cleanup destinations from overwriting nested inputs.
  Invalid input/output combinations fail before native prediction writes files.
- Enforces IntaRNA >=3 independently at the homology, primary prediction and web
  export entry points. Standalone helper use now has the same release requirement
  as the complete workflow. Prediction caching uses the detected release even
  when an old option file records another version.
- Adds the remaining public helper/API documentation, including table schemas,
  numerical precision, maintenance RNG scope and service-adapter side effects.
  AGENTS.md records a requirement-by-requirement implementation/test map.
