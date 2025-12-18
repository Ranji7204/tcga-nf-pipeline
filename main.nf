nextflow.enable.dsl=2

params.input  = "${projectDir}/input"
params.outdir = "${projectDir}/results"

Channel
  .fromPath("${params.input}/*", type: 'dir')
  .map { dir -> tuple(dir.baseName, dir) }
  .set { REGION_INPUTS }


process BATCH_CHECK {

  tag "$region"

  publishDir "${params.outdir}/${region}", mode: 'copy'

  input:
  tuple val(region), path(region_dir)

  output:
  tuple val(region), path("batch_median.txt"),
                     path("logCPM_normalized.txt"),
                     path("metadata_with_batch_info.txt"),
                     path("PCA_before_batch_correction.png"),
                     path("heatmap_before_batch_correction.png")

  script:
  """
  cp ${region_dir}/* .
  Rscript ${projectDir}/scripts/tcga_barcode_batch_pipeline_auto.R
  """
}


process BATCH_CORRECTION {

  tag "$region"

  publishDir "${params.outdir}/${region}", mode: 'copy'

  input:
  tuple val(region),
      path("batch_median.txt"),
      path("logCPM_normalized.txt"),
      path("metadata_with_batch_info.txt"),
      path("PCA_before_batch_correction.png"),
      path("heatmap_before_batch_correction.png")

  output:
  tuple val(region),
        path("logCPM_batch_corrected.txt"),
        path("metadata_batch_corrected.txt"),
        path("PCA_after_batch_correction.png"),
        path("heatmap_after_batch_correction.png")

  script:
  """
  Rscript ${projectDir}/scripts/batch_correction.R
  """
}


process MERGE_ALL {

  publishDir "${params.outdir}/merged", mode: 'copy'

  input:
  val regions

  output:
  path "logCPM_merged_all_regions.txt"
  path "metadata_merged_all_regions.txt"
  path "PCA_merged_regions.png"
  path "heatmap_merged_regions.png"

  script:
  """
  Rscript ${projectDir}/scripts/merge_all_files.R ${params.outdir}
  """
}


workflow {

  checked   = BATCH_CHECK(REGION_INPUTS)
  corrected = BATCH_CORRECTION(checked)

  MERGE_ALL(corrected.map{ it[0] }.collect())
}
