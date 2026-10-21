# ChangeLog

## Unreleased — Python 3 CopraRNA2

Development release `3.0.0.dev0`, addressing
[issue #1](https://github.com/BackofenLab/CopraRNA/issues/1).
The rewrite and subsequent corrections are reviewed separately in
[PR #2](https://github.com/BackofenLab/CopraRNA/pull/2) and
[PR #3](https://github.com/BackofenLab/CopraRNA/pull/3).

### Rewrite foundation — PR #2

The source workflow identifies itself as **CopraRNA 2.1.3**. Its controller and
auxiliary Perl, R, Python 2, shell and JavaScript programs were translated into a
Python 3 package. Original tracked files and bundled biological data remain
byte unchanged; the Python runtime never invokes those scripts. The initial
rewrite provided counterparts for all 41 scripts, with adjacent mapping documents
and a machine-readable migration inventory. The reviewed runtime now retains 39
active counterparts and marks two obsolete scripts retired.

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

Historical helper comparisons can explicitly reproduce retained defects in a
thread-local context. Complete production predictions always use corrected code.

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
  that collide with required generated files before work starts.
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
portable lockfile. Original sources remain for auditing, with optional interpreter
requirements separated in `environment-legacy.yml` and `environment-legacy-r.yml`.
Keep these optional interpreters in separate environments to avoid incompatible
XML/ICU libraries. The retained RobustRankAggreg 1.1 oracle installs from its
archived pure-R source, avoiding obsolete conda ABI pins.
