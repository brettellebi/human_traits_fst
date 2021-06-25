# For snakemake

module load singularity-3.7.0-gcc-9.3.0-dp5ffrp
bsub -Is bash
cd /hps/software/users/birney/ian/repos/human_traits_fst
conda activate snakemake_6.4.1
snakemake \
  --jobs 5000 \
  --latency-wait 100 \
  --cluster-config config/cluster.yaml \
  --cluster 'bsub -g /snakemake_bgenie -J {cluster.name} -n {cluster.n} -M {cluster.memory} -outdir {cluster.outdir} -o {cluster.outfile} -e {cluster.error}' \
  --keep-going \
  --rerun-incomplete \
  --use-conda \
  --use-singularity \
  -s workflow/Snakefile \
  -p

# For R

