import yaml as yml
import logging as logg
import os
import subprocess
import glob
import pandas as pd

def load_config(config_file):
    with open(config_file) as file:
        config = yml.safe_load(file)
    return config

def setup_logging(log_file):
    logg.basicConfig(filename=log_file, level=logg.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    logg.info("Pipeline started")

def run_command(tool):
    logg.info(f"command is starting: {tool}")
    result = subprocess.run(tool, capture_output=True, text=True)
    returncode = result.returncode
    if returncode == 0:
        logg.info(f"command is successful: {returncode}")
    else:
        logg.error(f"command failed: {returncode}")
    return result

def run_fastqc(fastq_file, output_dir):
    run_command(["fastqc", fastq_file, "-o", output_dir])
    
def run_trimmomatic(r1, r2, trimmed_dir, adapter_file):
    r1_trimmed = os.path.basename(r1).replace(".fastq.gz", "_paired.fastq.gz")
    r2_trimmed = os.path.basename(r2).replace(".fastq.gz", "_paired.fastq.gz")
    r1_trimmed_unpaired = os.path.basename(r1).replace(".fastq.gz", "_unpaired.fastq.gz")
    r2_trimmed_unpaired = os.path.basename(r2).replace(".fastq.gz", "_unpaired.fastq.gz")

    run_command([config["trimmomatic"], "PE", r1, r2, os.path.join(trimmed_dir, r1_trimmed), os.path.join(trimmed_dir, r1_trimmed_unpaired), os.path.join(trimmed_dir, r2_trimmed), os.path.join(trimmed_dir, r2_trimmed_unpaired), f"ILLUMINACLIP:{adapter_file}:2:30:10"])

def download_sample(url, input_dir):
    run_command(["wget", url, "-P", input_dir])

def run_star_index(fasta, gtf, index_dir, threads):
    run_command([config["star"], "--runMode", "genomeGenerate", "--genomeDir", index_dir, "--genomeFastaFiles", fasta, "--sjdbGTFfile", gtf, "--runThreadN", threads, "--genomeSAindexNbases", "11"])

def run_star_align(r1, r2, index_dir, output_dir, threads):
    run_command([config["star"], "--runMode", "alignReads", "--genomeDir", index_dir, "--readFilesIn", r1, r2, "--outFileNamePrefix", output_dir, "--runThreadN", threads, "--readFilesCommand", "gunzip -c", "--outSAMtype", "BAM", "SortedByCoordinate"])

def run_samtools_index(bam_file):
    run_command([config["samtools"], "index", bam_file])

def run_featurecounts(bam_file, gtf, output_dir):
    output_file = os.path.join(output_dir, os.path.basename(bam_file).replace(".bam", ".txt"))
    run_command([config["featurecounts"], "-a", gtf, "-o", output_file, bam_file])

def run_deseq2(counts_dir, output_dir, script_path, metadata):
    run_command(["Rscript", script_path, counts_dir, output_dir, metadata])
    
if __name__== "__main__":
    config = load_config("./config/config.yaml")
    setup_logging(os.path.join(config["logs_dir"], "pipeline.log"))
    
    # for d in ["fastqc", "trimmed_fastqc", "bam", "featurecounts"]:
    #     os.makedirs(os.path.join(config["output_dir"], d), exist_ok=True)
    # os.makedirs(config["trimmed_dir"], exist_ok=True)

    # metadata = pd.read_csv(config["metadata"])
    # url_combined = metadata["r1_url"].tolist() + metadata["r2_url"].tolist()
    
    # for i in url_combined:
    #     download_sample(i, config["input_dir"])

    # for fastq_file in glob.glob(os.path.join(config["input_dir"], "*.fastq.gz")):
    #     run_fastqc(fastq_file, os.path.join(config["output_dir"], "fastqc"))
    
    #logg.info("=" * 50)

    # for r1 in glob.glob(os.path.join(config["input_dir"], "*_R1.fastq.gz")):
    #     r2 = r1.replace("_R1", "_R2")
    #     run_trimmomatic(r1, r2, config["trimmed_dir"], config["adapter_file"])

    #logg.info("=" * 50)

    # for fastq_file in glob.glob(os.path.join(config["trimmed_dir"], "*_paired.fastq.gz")):
    #     run_fastqc(fastq_file, os.path.join(config["output_dir"], "trimmed_fastqc"))
    
    #run_star_index(config["fasta_file"], config["annotation_file"], config["STAR_index"], str(config["threads"]))

    logg.info("=" * 50)

    for r1 in glob.glob(os.path.join(config["trimmed_dir"], "*_R1_paired.fastq.gz")):
        r2 = r1.replace("_R1", "_R2")
        run_star_align(r1, r2, config["STAR_index"], os.path.join(config["output_dir"], "bam", os.path.basename(r1).replace("_R1_paired.fastq.gz", "_")), str(config["threads"]))

    # for sample_bam in glob.glob(os.path.join(config["output_dir"], "bam/*.bam")):
    #     run_samtools_index(sample_bam)

    # for sample_bam in glob.glob(os.path.join(config["output_dir"], "bam/*.bam")):
    #     run_featurecounts(sample_bam, 
    #                       config["annotation_file"], 
    #                       os.path.join(config["output_dir"], "featurecounts"))

    # run_deseq2(
    #     os.path.join(config["output_dir"], "featurecounts"), 
    #     os.path.join(config["output_dir"], "deseq2"), 
    #     config["deseq2_script"], config["metadata"]
    #     )