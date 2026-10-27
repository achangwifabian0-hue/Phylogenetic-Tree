# Zero to Tree: Bacterial Phylogenetics for Absolute Beginners

*A beginner-friendly teaching workflow for comparative genomics of* Klebsiella pneumoniae *using Ubuntu/WSL, NCBI Datasets, Snippy, IQ-TREE 2, and iTOL.*

> **Scope:** This repository is a learning exercise, not a clinical or outbreak-investigation pipeline. A tree alone cannot prove direct transmission or establish transmission direction.

## Learning objectives

By the end, you should be able to:

- Navigate files and directories in Ubuntu or Windows Subsystem for Linux (WSL).
- Download and inspect public bacterial genome assemblies.
- Run reference-based variant calling with Snippy.
- Build a core-genome alignment with `snippy-core`.
- Infer a preliminary maximum-likelihood tree with IQ-TREE 2.
- Visualize a Newick tree in iTOL.
- Identify quality checks needed before interpreting results.

## Workflow overview

1. Download ten genome assemblies: one KPPR1 reference and nine comparison strains.
2. Inspect the assemblies and verify their identities.
3. Run Snippy on each of the nine comparison strains against KPPR1.
4. Combine the Snippy results into a core-genome alignment.
5. Inspect and clean the alignment.
6. Infer a preliminary tree with IQ-TREE 2.
7. Inspect labels and visualize the tree in iTOL.

---

## Section 0. Set up the environment

Run the commands below in an Ubuntu terminal or WSL terminal, not in the Windows Command Prompt.

### 0.1 Activate your Conda environment

If you already created the environment:

```bash
conda activate snippy-env
```

If you have not created it, install the required tools with Conda:

```bash
conda create -n snippy-env -c conda-forge -c bioconda snippy iqtree seqkit ncbi-datasets-cli unzip
conda activate snippy-env
```

Check that the tools are available:

```bash
snippy --version
iqtree2 --version
seqkit version
datasets --version
unzip -v | head -n 2
```

If a command returns `command not found`, check that the environment is active and that the relevant package installed successfully. The exact output and version numbers will vary.

### 0.2 Create the project folders

```bash
mkdir -p ~/phylogenetic_tree_project/results/genomes
mkdir -p ~/phylogenetic_tree_project/results/genome_downloads
mkdir -p ~/phylogenetic_tree_project/results/genome_downloads/extracted
cd ~/phylogenetic_tree_project
pwd
```

The `pwd` command should show a path ending in `phylogenetic_tree_project`.

---

## Section 1. Download the genomes

This tutorial uses ten *K. pneumoniae* assemblies: KPPR1 is the reference, and the other nine are comparison genomes. The reference is included in the final alignment as a reference sequence; it is **not** a tenth clinical sample.

### 1.1 Assembly accessions

| Role | Strain label used in this tutorial | Assembly accession |
|---|---|---|
| Reference | KPPR1 | `GCF_000742755.1` |
| Sample 1 | MS6671 | `GCF_001455995.1` |
| Sample 2 | NTUH-K2044 | `GCF_000009885.1` |
| Sample 3 | HS11286 | `GCF_000240185.1` |
| Sample 4 | MGH78578 (MGH 78578) | `GCF_000016305.1` |
| Sample 5 | KPNIH10 | `GCA_000281435.2` |
| Sample 6 | KPNIH1 | `GCA_000281535.2` |
| Sample 7 | KPNIH24 | `GCA_000714675.1` |
| Sample 8 | KPNIH27 | `GCA_000695935.1` |
| Sample 9 | KPNIH29 | `GCA_000784945.1` |

**Accession check:** The accession-to-strain mappings in this table were checked against NCBI assembly/BioProject records during this review. For example, NCBI identifies `GCF_000742755.1` as KPPR1, `GCF_001455995.1` as MS6671, `GCA_000281535.2` as KPNIH1, and the listed `GCA_` records as KPNIH10, KPNIH24, KPNIH27, and KPNIH29. Recheck the NCBI record if you intentionally update an accession. The `GCF_` prefix denotes a RefSeq assembly record and `GCA_` denotes a GenBank assembly record; keep the version suffix because it identifies a specific assembly version.

### 1.2 Download the ten NCBI Datasets packages

Run these commands from the project directory. Each command saves a ZIP package with a distinct filename, so downloads do not overwrite one another. but first do the checks below 


### 1.2.a Check whether Conda is available
In your Ubuntu terminal, run:

```Bash
conda --version
```
- If you see a version number, Conda is available. Continue to Step 2.
- If you see conda: command not found, skip to Step 3.

### Step 2  Install NCBI Datasets
If Step 1 showed that Conda is available, run:

```Bash
conda activate snippy-env
```
Then install NCBI Datasets into that environment:

```Bash
conda install -c conda-forge ncbi-datasets-cli
```

If Conda asks you to confirm the installation, type y and press Enter.

### Step 3 — Confirm the installation
Once installation finishes, run:

```Bash
datasets --version
```

You should see a version number rather than command not found.
If that works, verify that your project directory is correct:

```Bash
cd ~/phylogenetic_tree_project
mkdir -p results/genome_downloads
```

### Now test just one genome download first:

```Bash
datasets download genome accession GCF_000742755.1 \
  --include genome \
  --filename results/genome_downloads/KPPR1.zip
```
### Check whether the ZIP file exists:
``` Bash
ls -lh results/genome_downloads/KPPR1.zip
```
### If the file is listed with a nonzero size, the first download appears to have succeeded. We can then download the remaining nine genomes.

### if all works, then run 

```bash
cd ~/phylogenetic_tree_project
mkdir -p results/genome_downloads
datasets download genome accession GCF_000742755.1 --include genome --filename results/genome_downloads/KPPR1.zip
datasets download genome accession GCF_001455995.1 --include genome --filename results/genome_downloads/MS6671.zip
datasets download genome accession GCF_000009885.1 --include genome --filename results/genome_downloads/NTUH-K2044.zip
datasets download genome accession GCF_000240185.1 --include genome --filename results/genome_downloads/HS11286.zip
datasets download genome accession GCF_000016305.1 --include genome --filename results/genome_downloads/MGH78578.zip
datasets download genome accession GCA_000281435.2 --include genome --filename results/genome_downloads/KPNIH10.zip
datasets download genome accession GCA_000281535.2 --include genome --filename results/genome_downloads/KPNIH1.zip
datasets download genome accession GCA_000714675.1 --include genome --filename results/genome_downloads/KPNIH24.zip
datasets download genome accession GCA_000695935.1 --include genome --filename results/genome_downloads/KPNIH27.zip
datasets download genome accession GCA_000784945.1 --include genome --filename results/genome_downloads/KPNIH29.zip
```

If any command fails, stop and resolve that download before continuing. Do not assume that all ten archives exist just because the commands were entered.

Check that there are ten ZIP files:

```bash
find results/genome_downloads -maxdepth 1 -type f -name '*.zip' | sort
find results/genome_downloads -maxdepth 1 -type f -name '*.zip' | wc -l
```

Expected count: `10`.

### 1.3 Extract the genome FASTA files

NCBI Datasets packages normally contain a genomic FASTA file under `ncbi_dataset/data/`. The loop below extracts each archive into a separate folder, locates its genomic FASTA, and copies it to `results/genomes/` with a consistent name.

```bash
cd ~/phylogenetic_tree_project
mkdir -p results/genomes results/genome_downloads/extracted

for strain in KPPR1 MS6671 NTUH-K2044 HS11286 MGH78578 \
              KPNIH10 KPNIH1 KPNIH24 KPNIH27 KPNIH29
do
    archive="results/genome_downloads/${strain}.zip"
    extract_dir="results/genome_downloads/extracted/${strain}"

    if [ ! -f "$archive" ]; then
        echo "ERROR: Missing archive: $archive"
        continue
    fi

    mkdir -p "$extract_dir"
    unzip -oq "$archive" -d "$extract_dir"

    fasta=$(find "$extract_dir/ncbi_dataset/data" \
        -type f -name '*_genomic.fna' -print -quit)

    if [ -n "$fasta" ]; then
        cp "$fasta" "results/genomes/${strain}.fna"
        echo "Saved results/genomes/${strain}.fna"
    else
        echo "ERROR: No *_genomic.fna found for ${strain}"
    fi
done
```

Check the resulting files:

```bash
ls -lh results/genomes/
find results/genomes -maxdepth 1 -type f -name '*.fna' | wc -l
```

Expected count: `10`. If the count is not ten, do not proceed until the missing file(s) are resolved.

The expected filenames are:

```text
KPPR1.fna
MS6671.fna
NTUH-K2044.fna
HS11286.fna
MGH78578.fna
KPNIH10.fna
KPNIH1.fna
KPNIH24.fna
KPNIH27.fna
KPNIH29.fna
```

### 1.4 Manual download alternative

If you cannot use the command-line tool, find each accession in the [NCBI Assembly database](https://www.ncbi.nlm.nih.gov/assembly/) and download its genomic FASTA. Save all ten files in `~/phylogenetic_tree_project/results/genomes/` using exactly the filenames listed above.

Do not mix assemblies downloaded under different accessions without recording which accession produced each file.

---

## Section 2. Inspect the assemblies

### 2.1 Check sequence counts and lengths

```bash
cd ~/phylogenetic_tree_project
seqkit stats results/genomes/*.fna
```

**Checkpoint:** Confirm that ten files are listed. Genome assemblies can contain multiple sequences, including plasmids, so the number of sequences and total length can vary among assemblies.

### 2.2 Inspect FASTA headers

```bash
for file in results/genomes/*.fna; do
    echo "=== $file ==="
    grep '^>' "$file" | head -n 5
done
```

FASTA headers start with `>`. Use them to check the assembly identity and see whether the file contains chromosome and/or plasmid sequences.

### 2.3 Check that the reference is the expected assembly

The KPPR1 reference accession listed above is a single-chromosome complete assembly. Check its length and headers rather than assuming the file is correct:

```bash
seqkit stats results/genomes/KPPR1.fna
grep '^>' results/genomes/KPPR1.fna
```

The reference sequence should be consistent with the selected assembly record. If it contains unexpected extra records or appears to be a different strain, resolve the discrepancy before continuing.

**Quality note:** These checks confirm basic file presence and structure; they do not establish species purity or rule out contamination.

---

## Section 3. Run reference-based variant calling with Snippy

### 3.1 The basic idea

KPPR1 is the coordinate reference. Snippy compares each of the nine other assemblies against KPPR1 and records variants relative to that reference. The reference itself is not run as a sample.

Run the following from the project directory. Each sample has a separate output folder.

```bash
cd ~/phylogenetic_tree_project

snippy --ref results/genomes/KPPR1.fna --ctgs results/genomes/MS6671.fna --outdir results/snippy_MS6671 --cpus 2
snippy --ref results/genomes/KPPR1.fna --ctgs results/genomes/NTUH-K2044.fna --outdir results/snippy_NTUH-K2044 --cpus 2
snippy --ref results/genomes/KPPR1.fna --ctgs results/genomes/HS11286.fna --outdir results/snippy_HS11286 --cpus 2
snippy --ref results/genomes/KPPR1.fna --ctgs results/genomes/MGH78578.fna --outdir results/snippy_MGH78578 --cpus 2
snippy --ref results/genomes/KPPR1.fna --ctgs results/genomes/KPNIH10.fna --outdir results/snippy_KPNIH10 --cpus 2
snippy --ref results/genomes/KPPR1.fna --ctgs results/genomes/KPNIH1.fna --outdir results/snippy_KPNIH1 --cpus 2
snippy --ref results/genomes/KPPR1.fna --ctgs results/genomes/KPNIH24.fna --outdir results/snippy_KPNIH24 --cpus 2
snippy --ref results/genomes/KPPR1.fna --ctgs results/genomes/KPNIH27.fna --outdir results/snippy_KPNIH27 --cpus 2
snippy --ref results/genomes/KPPR1.fna --ctgs results/genomes/KPNIH29.fna --outdir results/snippy_KPNIH29 --cpus 2
```

These commands may take time. Review each run's terminal output and log files. Do not continue if Snippy reports an error.

### 3.2 Verify all nine Snippy output folders

```bash
cd ~/phylogenetic_tree_project

for strain in MS6671 NTUH-K2044 HS11286 MGH78578 \
              KPNIH10 KPNIH1 KPNIH24 KPNIH27 KPNIH29
do
    if [ -d "results/snippy_${strain}" ]; then
        echo "FOUND: results/snippy_${strain}"
        test -f "results/snippy_${strain}/snps.vcf" \
          && echo "  snps.vcf found" \
          || echo "  WARNING: snps.vcf missing"
    else
        echo "ERROR: Output directory missing for ${strain}"
    fi
done
```

**Checkpoint:** All nine output directories should exist, and each should contain expected Snippy output files. A directory existing by itself does not guarantee that the run completed successfully.

---

## Section 4. Build the core-genome alignment

### 4.1 Run `snippy-core`

Combine the nine Snippy result folders. Use the same KPPR1 reference used for all Snippy runs.

```bash
cd ~/phylogenetic_tree_project

snippy-core \
  --ref results/genomes/KPPR1.fna \
  --prefix results/core10 \
  results/snippy_MS6671 \
  results/snippy_NTUH-K2044 \
  results/snippy_HS11286 \
  results/snippy_MGH78578 \
  results/snippy_KPNIH10 \
  results/snippy_KPNIH1 \
  results/snippy_KPNIH24 \
  results/snippy_KPNIH27 \
  results/snippy_KPNIH29
```

The prefix `core10` is used because the resulting alignment is expected to include ten sequences: KPPR1 reference plus nine comparison strains.

### 4.2 Check the alignment files

```bash
ls -lh results/core10.full.aln results/core10.aln
grep -c '^>' results/core10.full.aln
seqkit stats results/core10.full.aln results/core10.aln
```

The full alignment should contain ten sequence records. The exact alignment length and SNP count depend on the assemblies and analysis results; do not copy an expected SNP count from another run as if it were guaranteed.

If a file is missing or the sequence count is not ten, stop and check the Snippy runs and `snippy-core` output.

---

## Section 5. Clean the alignment and infer a preliminary tree

### 5.1 Clean the full core alignment

Use Snippy's alignment-cleaning utility to remove problematic alignment positions according to the tool's documented behavior:

```bash
cd ~/phylogenetic_tree_project
snippy-clean_full_aln < results/core10.full.aln > results/core10.clean.full.aln
```

Confirm that the cleaned file exists and inspect its statistics:

```bash
ls -lh results/core10.clean.full.aln
seqkit stats results/core10.clean.full.aln
```

If `snippy-clean_full_aln` is not found, check that Snippy is installed and the correct Conda environment is active.

### 5.2 Infer the tree from the cleaned full alignment

```bash
iqtree2 \
  -s results/core10.clean.full.aln \
  -m MFP \
  -B 1000 \
  -T 2 \
  --prefix results/klebsiella_tree
```

This asks IQ-TREE 2 to use ModelFinder (`-m MFP`) to choose a substitution model and to perform 1,000 ultrafast bootstrap replicates (`-B 1000`). `-T 2` requests two threads. The run may take time.

Expected outputs include:

- `results/klebsiella_tree.treefile` — inferred maximum-likelihood tree in Newick format.
- `results/klebsiella_tree.iqtree` — analysis report.
- `results/klebsiella_tree.contree` — consensus tree with branch support.

Confirm that these files exist and inspect the IQ-TREE log for successful completion, the selected model, and warnings:

```bash
ls -lh results/klebsiella_tree.treefile results/klebsiella_tree.iqtree results/klebsiella_tree.contree
cat results/klebsiella_tree.treefile
```

**Do not mix alignment types without changing the model appropriately.** The main workflow above uses the cleaned full core alignment. If you choose a SNP-only alignment instead, ascertainment-bias correction may be needed; consult the IQ-TREE documentation for the correct model syntax for your installed version. The SNP-only route is not part of the main commands in this tutorial.

### 5.3 Optional advanced extension: recombination filtering

*K. pneumoniae* may undergo homologous recombination. Recombination filtering is an advanced optional extension, not a required step in the basic workflow. If you choose to use Gubbins, install it separately, check its version-specific documentation, and inspect its output files before using any filtered alignment for tree inference. Do not run a command against a filename that has not actually been produced.

---

## Section 6. Visualize the tree in iTOL

### 6.1 Inspect tip labels first

```bash
cd ~/phylogenetic_tree_project
cat results/klebsiella_tree.treefile
```

Inspect the actual tip labels before renaming them. Do not blindly replace labels based on an assumption about how your Snippy version names them.

If the tree contains `snippy_` prefixes and the reference tip is exactly named `Reference`, you can make a presentation copy with simplified labels:

```bash
sed 's/snippy_//g; s/Reference/KPPR1/g' \
  results/klebsiella_tree.treefile \
  > results/klebsiella_tree_iTOL.tree
```

This command changes labels only; it does not change the inferred topology or branch lengths. If your labels differ, edit the substitutions to match the actual tree labels. Check the output:

```bash
cat results/klebsiella_tree_iTOL.tree
```

### 6.2 Upload to iTOL

1. Open [iTOL](https://itol.embl.de/).
2. Upload `results/klebsiella_tree_iTOL.tree` (or the original `.treefile` if no label changes are needed).
3. Check that the expected ten tips are present: the KPPR1 reference plus nine comparison strains.
4. Inspect branch lengths and support values.
5. Customize the display and export an image if needed.

On WSL, you can open the results folder in Windows Explorer with:

```bash
explorer.exe "$(wslpath -w "$PWD/results")"
```

**Interpretation reminder:** Genetic similarity alone does not prove direct transmission between patients or establish transmission direction.

---

## Section 7. Scientific interpretation: common traps

| Tree feature | What it means | Common trap |
|---|---|---|
| Topology | The branching pattern inferred from the alignment and model | Closely grouped tips do not, by themselves, prove direct transmission. |
| Branch length | Estimated evolutionary change under the fitted model | Branch lengths are not calendar time. |
| Bootstrap support | How consistently a clade is recovered under the resampling procedure | High support does not guarantee error-free data or prove a biological claim. |
| Reference genome | Provides the coordinate system for variant calling and is represented as a sequence in the alignment | KPPR1 is a reference, not an additional clinical sample. |

---

## Section 8. Quality-control checklist

Before presenting or publishing results, check that:

- [ ] Each assembly accession has been verified against the strain name in NCBI.
- [ ] All ten expected genome FASTA files exist and have been inspected.
- [ ] Assembly completeness and potential contamination have been assessed where appropriate.
- [ ] The KPPR1 reference is suitable for the intended comparison.
- [ ] All nine Snippy runs completed successfully and logs were reviewed.
- [ ] The core alignment contains ten records and has been inspected for missing or ambiguous data.
- [ ] The selected substitution model and tree-inference settings were recorded.
- [ ] The reference is clearly identified as a reference, not a clinical sample.
- [ ] Tree relationships are interpreted alongside biological and epidemiological context.
- [ ] Software versions and command parameters have been recorded.

---

## Section 9. Record software versions

Run this from the project directory to save the versions available in the active environment:

```bash
cd ~/phylogenetic_tree_project
mkdir -p results

{
  echo "Execution date: $(date -Is)"
  echo "---------------------------"
  echo "Snippy version:"
  snippy --version
  echo "IQ-TREE 2 version:"
  iqtree2 --version
  echo "SeqKit version:"
  seqkit version
  echo "NCBI Datasets version:"
  datasets --version
} > results/software_versions.txt

cat results/software_versions.txt
```

Keep this file with the results. Software output formats and command options can vary between versions.

---

## Section 10. Troubleshooting

<details>
<summary><strong>1. <code>snippy: command not found</code></strong></summary>

Activate the environment and check the executable:

```bash
conda activate snippy-env
which snippy
snippy --version
```
</details>

<details>
<summary><strong>2. <code>iqtree2: command not found</code></strong></summary>

Activate the environment and check IQ-TREE 2:

```bash
conda activate snippy-env
which iqtree2
iqtree2 --version
```

If it is missing, install it in the active environment:

```bash
conda install -c conda-forge -c bioconda iqtree
```
</details>

<details>
<summary><strong>3. <code>No such file or directory</code></strong></summary>

Check that you are in the project directory and that the expected files exist:

```bash
cd ~/phylogenetic_tree_project
pwd
ls -lh results/genomes/
```

File names are case-sensitive. Confirm that the file name in the command matches the actual file.
</details>

<details>
<summary><strong>4. NCBI download or ZIP extraction failed</strong></summary>

Check that `datasets` and `unzip` are installed, that the internet connection is working, and that the archive exists and is non-empty:

```bash
datasets --version
unzip -t results/genome_downloads/KPPR1.zip
```

Repeat the check with the archive that failed. Do not continue with an incomplete set of genomes.
</details>

<details>
<summary><strong>5. Snippy or snippy-core failed</strong></summary>

Review the relevant Snippy output directory and log files. Confirm that every sample FASTA is valid and that all nine Snippy runs used the same KPPR1 reference file. Re-run only after correcting the reported issue.
</details>

---

## Official documentation

- [NCBI Datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/)
- [Snippy](https://github.com/tseemann/snippy)
- [IQ-TREE 2](https://www.iqtree.org/)
- [SeqKit](https://bioinf.shenwei.me/seqkit/)
- [iTOL](https://itol.embl.de/)
- [AliView](https://ormbunkar.se/aliview/)

## Completion checklist

You have completed the basic exercise when you can:

1. Download and inspect all ten assemblies.
2. Run Snippy on the nine comparison strains.
3. Generate and inspect the core-genome alignment.
4. Infer a preliminary tree with IQ-TREE 2.
5. Display the tree in iTOL.
6. Explain topology, branch lengths, and support values.
7. Describe the quality checks needed before drawing scientific conclusions.

This course is intended for learning and reproducible bioinformatics practice. Do not use a preliminary tree as the sole basis for clinical decisions, outbreak declarations, or transmission claims.

## Acknowledgements

This course was developed from learning and practical experience in phylogenetic tree construction and analysis during the genomics course of the African Stars Fellowship at Institut Pasteur de Dakar and CARE.

## Producer

- Fabian Achangwi

## Course facilitators

- Ndeye Awa Ndiaye
- Abdoul Wahab Soumare
