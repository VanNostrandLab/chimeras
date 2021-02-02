#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
A pipeline designed to identify miRNA–target chimeras.
"""

import os
import sys
import logging
import subprocess
import argparse

from datetime import datetime

import ruffus


logger = logging.getLogger("CHIMERAS")
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler(sys.stderr)
handler.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] %(message)s',
                                       datefmt='%Y-%m-%d %H:%M:%S'))
handler.setLevel(logging.DEBUG)
logger.addHandler(handler)

FASTQS = []
CORES = 4
GENOME = '/storage/vannostrand/genomes/hg19/from_yeolab/star_sjdb_encode'
RNA = '/storage/vannostrand/software/chimeras/data/mature.U.to.T.fa'


def run_cmd(cmd, output_mode='wt', **kwargs):
    """
    Run cmd or throw exception if run fails.
    """
    
    def cmding(cmd):
        cmd = [str(c) for c in cmd]
        if '<' in cmd:
            raise ValueError('Invalid cmd, standard input via "<" not supported yet.')
        message = ' '.join(cmd).replace(' -', ' \\\n  -').replace(' >', ' \\\n  >')
        if '>' in cmd:
            output = cmd[cmd.index('>') + 1]
            cmd = cmd[:cmd.index('>')]
        else:
            output = ''
        logger.debug(f'\n{message}')
        return cmd, output
    
    cmd, output = cmding(cmd)
    if output:
        text_mode = True if 't' in output_mode else None
        with open(output, output_mode) as out:
            process = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=text_mode, **kwargs)
        if process.returncode:
            message = process.stderr or open(output).read()
            raise Exception(f'Failed to run {cmd[0]} (exit code {process.returncode}):\n{message}')
    else:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, **kwargs)
        stdout, stderr = process.communicate()
        if process.returncode:
            raise Exception(f'Failed to run {cmd[0]} (exit code {process.returncode}):{stderr or stdout}')


@ruffus.jobs_limit(1)
@ruffus.follows(ruffus.mkdir('logs', 'data', 'results', 'scripts'))
@ruffus.transform(FASTQS,
                  ruffus.formatter(r'.+/(?P<BASENAME>.*).f[ast]*q.gz$'),
                  '{BASENAME[0]}.fastq.gz')
def soft_link_fastq(fastq, link):
    link = os.path.abspath(link)
    if fastq == link:
        logger.warning("No symbolic was link made. "
                       "You are directly working on the original data files.")
    else:
        logger.debug(f'Soft link {os.path.basename(fastq)}:\n  ln -s {fastq} {link}')
        os.symlink(fastq, link)


@ruffus.transform(soft_link_fastq,
                  ruffus.formatter(r'.+/(?P<BASENAME>.*).fastq.gz$'),
                  '{BASENAME[0]}.umi.fastq.gz')
def extract_umi(fastq, umi_extracted_fastq):
    cmd = ['umi_tools',
           'extract',
           '--random-seed', 1,
           '--stdin', fastq,
           '--bc-pattern', 'NNNNNNNNNN',
           '--log', umi_extracted_fastq.replace('.fastq.gz', '.extract.metrics'),
           '--stdout', umi_extracted_fastq]
    logger.debug(f'Extracting UMIs for {os.path.basename(fastq)} ...')
    run_cmd(cmd)
    logger.debug(f'Extracting UMIs for {os.path.basename(fastq)} completed.')


@ruffus.jobs_limit(1)
@ruffus.transform(extract_umi, ruffus.suffix('.fastq.gz'), '.trim.fastq')
def cut_adapter(fastq, trimmed_fastq):
    trim_tmp = trimmed_fastq.replace('.trim.fastq', '.trim.tmp.fastq')
    trim_metrics = trimmed_fastq.replace('.trim.fastq', '.trim.metrics')
    trim_tmp_metrics = trimmed_fastq.replace('.trim.fastq', '.trim.tmp.metrics')
    trim_umi_metrics = trimmed_fastq.replace('.trim.fastq', '.trim.umi.metrics')
    cmd = ['cutadapt',
           '-j', CORES,
           '-O', 1,
           '--match-read-wildcards',
           '--times', 1,
           '-e', 0.1,
           '--quality-cutoff', 6,
           '-m', 18,
           '-o', trim_tmp,
           # '-a', 'NNNNNNNNNNAGATC',
           # '-a', 'NNNNNNNNNAGATCG',
           # '-a', 'NNNNNNNNAGATCGG',
           # '-a', 'NNNNNNNAGATCGGA',
           # '-a', 'NNNNNNAGATCGGAA',
           # '-a', 'NNNNNAGATCGGAAG',
           # '-a', 'NNNNAGATCGGAAGA',
           # '-a', 'NNNAGATCGGAAGAG',
           # '-a', 'NNAGATCGGAAGAGC',
           # '-a', 'NAGATCGGAAGAGCA',
           '-a', 'AGATCGGAAGAGCAC',
           '-a', 'GATCGGAAGAGCACA',
           '-a', 'ATCGGAAGAGCACAC',
           '-a', 'TCGGAAGAGCACACG',
           '-a', 'CGGAAGAGCACACGT',
           '-a', 'GGAAGAGCACACGTC',
           '-a', 'GAAGAGCACACGTCT',
           '-a', 'AAGAGCACACGTCTG',
           '-a', 'AGAGCACACGTCTGA',
           '-a', 'GAGCACACGTCTGAA',
           '-a', 'AGCACACGTCTGAAC',
           '-a', 'GCACACGTCTGAACT',
           '-a', 'CACACGTCTGAACTC',
           '-a', 'ACACGTCTGAACTCC',
           '-a', 'CACGTCTGAACTCCA',
           '-a', 'ACGTCTGAACTCCAG',
           '-a', 'CGTCTGAACTCCAGT',
           '-a', 'GTCTGAACTCCAGTC',
           '-a', 'TCTGAACTCCAGTCA',
           '-a', 'CTGAACTCCAGTCAC',
           fastq, '>', trim_tmp_metrics]
    logger.debug(f'Trimming adapters for {fastq} ...')
    run_cmd(cmd)
    logger.debug(f'Trimming adapters {fastq} completed.')

    cmd = ['cutadapt',
           '-j', CORES,
           '-u', -9,
           '-o', trimmed_fastq,
           fastq, '>', trim_umi_metrics]
    logger.debug(f'Trimming UMIs for {trim_tmp} ...')
    run_cmd(cmd)
    logger.debug(f'Trimming UMIs {trim_tmp} completed.')
    with open(trim_metrics, 'w') as o:
        with open(trim_tmp_metrics) as f:
            o.write(f.read())
        with open(trim_umi_metrics) as f:
            o.write(f.read())
    os.unlink(trim_tmp)
    os.unlink(trim_tmp_metrics)
    os.unlink(trim_umi_metrics)


@ruffus.transform(cut_adapter, ruffus.suffix('.trim.fastq'), '.trim.sorted.fastq')
def sort_fastq(fastq, sorted_fastq):
    cmd = ['fastq-sort', '--id', fastq, '>', sorted_fastq]
    logger.debug(f'Sorting adapters trimmed fastq for {fastq} ...')
    run_cmd(cmd)
    logger.debug(f'Sorting adapters trimmed fastq for {fastq} completed.')


@ruffus.transform(sort_fastq, ruffus.suffix('.trim.sort.fastq'), '.trim.sort.fasta')
def fastq_to_fasta(fastq, fasta):
    cmd = ['seqtk', 'seq', '-A', fastq, '>', fasta]
    logger.debug(f'Converting fastq to fasta for {fastq} ...')
    run_cmd(cmd)
    logger.debug(f'Converting fastq to fasta for {fastq} completed.')


@ruffus.transform(fastq_to_fasta, ruffus.suffix('.trim.sort.fasta'), '.trim.sort.collapse.fasta')
def collapse_fasta(fasta, collapsed_fasta):
    cmd = ['fasta2collapse.pl', fasta, collapsed_fasta]
    logger.debug(f'Collapsing identical fasta sequences for fasta file {fasta} ...')
    env = os.environ.copy()
    path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    env['PATH'] = f'{os.path.join(path, "venv", "ctk-1.1.4")}:{env["PATH"]}'
    env['PERL5LIB'] = os.path.join(path, "venv", "czplib-1.0.8")
    run_cmd(cmd, env=env)
    logger.debug(f'Collapsing identical fasta sequences for fasta file {fasta} completed.')


@ruffus.transform(collapse_fasta, ruffus.suffix('.umi.trim.sort.collapse.fasta'), '.bowtie.index.log')
def bowtie_index(collapsed_fasta, log):
    cmd = ['bowtie-build',
           '--offrate', 2,
           collapsed_fasta,
           log.replace('.bowtie.index.log', ''),
           '>', log]
    logger.debug(f'Building bowtie index using fasta file {collapsed_fasta} ...')
    run_cmd(cmd)
    logger.debug(f'Building bowtie index using fasta file {collapsed_fasta} completed.')


@ruffus.jobs_limit(1)
@ruffus.transform(bowtie_index, ruffus.suffix('.bowtie.index.log'), '.reverse.mir.alignments.tsv')
def reverse_map_mirna(bowtie_index_log, reversed_alignment):
    cmd = ['bowtie',
           '-a', '-e', 35, '-f', '-l', 8, '-n', 1, '-p', CORES,
           bowtie_index_log.replace('.bowtie.index.log', ''),
           RNA,
           reversed_alignment,
           '>',  reversed_alignment.replace('.alignments.tsv', '.alignments.log')]
    collapsed_fasta = bowtie_index_log.replace('.bowtie.index.log', '.umi.trim.sort.collapsed.fasta')
    logger.debug(f'Reverse mapping miRNA to collapsed fasta file {collapsed_fasta} ...')
    run_cmd(cmd)
    logger.debug(f'Reverse mapping miRNA to collapsed fasta file {collapsed_fasta} completed.')


@ruffus.transform(reverse_map_mirna, ruffus.suffix('.alignments.tsv'), '.alignments.filtered.tsv')
def filter_alignment(alignment, filtered_alignment):
    cmd = ['python', 'collapse_bowtie_results.py',
           '--bowtie_align', alignment,
           '--out_file', filtered_alignment]
    logger.debug(f'Filtering alignment {alignment} ...')
    run_cmd(cmd)
    logger.debug(f'Filtering alignment {alignment} completed.')


@ruffus.transform(filter_alignment, ruffus.suffix('.alignments.filtered.tsv'), '.chimeric.candidates.fasta')
def candidate_chimeras(filtered_alignment, candidates_fasta):
    cmd = ['python', 'find_candidate_chimeric_seqs_from_mir_alignments.py',
           '--bowtie_align', filtered_alignment,
           '--fa_file', filtered_alignment.replace('.alignments.filtered.tsv', '.umi.trim.sort.fasta'),
           '--metrics_file', filtered_alignment.replace('.alignments.filtered.tsv', '.chimeric.candidates.metrics'),
           '--out_file', candidates_fasta]
    logger.debug(f'Find candidate chimeric sequences for alignment {filtered_alignment} ...')
    run_cmd(cmd)
    logger.debug(f'Find candidate chimeric sequences for alignment {filtered_alignment} completed.')


@ruffus.jobs_limit(1)
@ruffus.transform(candidate_chimeras,
                  ruffus.suffix('.chimeric.candidates.fasta'),
                  '.Aligned.sortedByCoord.out.bam')
def genomic_map(candidates_fasta, bam):
    cmd = ['STAR',
           '--runMode', 'alignReads',
           '--runThreadN', CORES,
           '--alignEndsType', 'EndToEnd',
           '--genomeDir', GENOME,
           '--genomeLoad', 'NoSharedMemory',
           '--outBAMcompression', 10,
           '--outFileNamePrefix', bam.replace('Aligned.sortedByCoord.out.bam', ''),
           '--outFilterMatchNminOverLread', 0.66,
           '--outFilterMultimapNmax', 1,
           '--outFilterMultimapScoreRange', 1,
           '--outFilterScoreMin', 10,
           '--outFilterScoreMinOverLread', 0.66,
           '--outFilterType', 'BySJout',
           '--outReadsUnmapped', 'Fastx',
           '--outSAMattrRGline', 'ID:foo',
           '--outSAMattributes', 'All',
           '--outSAMmode', 'Full',
           '--outSAMtype', 'BAM', 'SortedByCoordinate',
           '--outSAMunmapped', 'Within',
           '--outStd', 'Log',
           '--readFilesIn', candidates_fasta]
    logger.debug(f'Mapping candidate chimeras to genome for {candidates_fasta} ...')
    run_cmd(cmd)
    logger.debug(f'Mapping candidate chimeras to genome for {candidates_fasta} completed.')


@ruffus.transform(genomic_map,
                  ruffus.suffix('.Aligned.sortedByCoord.out.bam'),
                  '.Aligned.sortedByCoord.out.bam.bai')
def index_bam(bam, bai):
    cmd = ['samtools', 'index', bam, bai]
    logger.debug(f'Indexing {os.path.basename(bam)} ...')
    run_cmd(cmd)
    logger.debug(f'Indexing {os.path.basename(bam)} completed.')


@ruffus.transform(genomic_map, ruffus.suffix('.Aligned.sortedByCoord.out.bam'), '.dedup.bam')
def dedup_bam(bam, deduped_bam):
    cmd = ['umi_tools', 'dedup',
           '--random-seed', 1,
           '-I', bam,
           '--method', 'unique',
           '--output-stats', deduped_bam.replace('.deduped.bam', '.dedup'),
           '--log', deduped_bam.replace('.deduped.bam', '.dedup.metrics'),
           '-S', deduped_bam]
    logger.debug(f'Indexing {os.path.basename(bam)} ...')
    run_cmd(cmd)
    logger.debug(f'Indexing {os.path.basename(bam)} completed.')


@ruffus.transform(dedup_bam, ruffus.suffix('.dedup.bam'), '.dedup.bg')
def chimeras(deduped_bam, erich_regions):
    cmd = ['genomeCoverageBed', '-ibam', deduped_bam, '-bg', '>', erich_regions]
    logger.debug(f'Finding chimeric enriched regions for {deduped_bam} ...')
    run_cmd(cmd)
    logger.debug(f'Finding chimeric enriched regions for {deduped_bam} completed.')


def schedule(scheduler, fastqs, outdir, rna, cores, memory, hours, jobname, email):
    sbatch = """#!/usr/bin/env bash

#SBATCH -n {cores}                  # Number of cores (-n)
#SBATCH -N 1                        # Ensure that all cores are on one Node (-N)
#SBATCH -t {runtime}                # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH --mem={memory}G             # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --job-name={jobname}        # Short name for the job
"""
    sbatch_email = """
#SBATCH --mail-user={email}
#SBATCH --mail-type=ALL
"""
    pbs = """ #!/usr/bin/env bash

#PBS -l nodes=1:ppn={cores}
#PBS -l walltime={runtime}
#PBS -l vmem={memory}gb
#PBS -j oe
#PBS -N {jobname}
"""
    pbs_email = """
#PBS -M {email}
#PBS -m abe
"""
    code = r"""
export TMPDIR={tmpdir}
export TEMP={tmpdir}
export TMP={tmpdir}
source CHIMERAS_ENVIRONMENT

echo [$(date +"%m-%d-%Y %H:%M:%S")] "Chimeras start."
source ENVIRONMENT

chimeras \\
    -o {outdir} \\
    -r {rna} \\
    {fastqs}
echo [$(date +"%m-%d-%Y %H:%M:%S")] "Chimeras complete."
"""
    if scheduler.upper() in ('PBS', 'QSUB'):
        runtime, directive, exe, mail = f'{hours}:00:00', pbs, 'qsub', pbs_email
    elif scheduler.upper() in ('SLURM', 'SBATCH'):
        days, hours = divmod(hours, 24)
        runtime, directive, exe, mail = f'{days}-{hours:02}:00', sbatch, 'sbatch', sbatch_email
    else:
        raise ValueError(f'Unsupported scheduler: {scheduler}, see help for supported schedulers.')
    tmpdir = os.path.join(os.getcwd(), 'tmp')
    if not os.path.isdir(tmpdir):
        os.mkdir(tmpdir)
    data = {'cores': cores, 'runtime': runtime, 'memory': memory, 'jobname': jobname, 'fastqs': fastqs,
            'outdir': outdir, 'tmpdir': tmpdir, 'rna': rna, 'email': email}
    text = [directive, mail, code] if email else [directive, code]
    text = ''.join(text).format(**data)
    
    submitter = 'submit.sh'
    with open(submitter, 'w') as o:
        o.write(text)
    
    print(f'Job submit script was saved to {submitter}')
    subprocess.run([exe, submitter])
    print(f'Job {jobname} was successfully submitted with the following settings:')
    data = {'Job name:': jobname, 'Output directory:': os.getcwd(),
            'Number of cores:': cores, 'Job memory:': memory,
            'Job runtime:': f'{runtime} (D-HH:MM)'}
    for k, v in data.items():
        print(f'{k:>20} {v}')


def main():
    parser = argparse.ArgumentParser(description=__doc__, prog='clip')
    parser.add_argument('FASTQ', type=str, nargs='+',
                        help='Path to one or multiple FASTQ files (separated by space).')
    parser.add_argument('-o', type=str, dest='outdir',
                        help="Path to the output directory (if not exist, will try to create it first).",
                        default=os.getcwd())
    parser.add_argument('-g', type=str, dest='genome',
                        help="Path to STAR genome index directory.")
    parser.add_argument('-r', type=str, dest='rna',
                        help="Path to mature RNA FASTA file (generated using miRBase and "
                             "U needs to be replaced by T).")
    parser.add_argument('-j', type=str, dest='jobname',
                        help="Name of your job, default: chimeras",
                        default='eCLIP')
    parser.add_argument('-e', type=str, dest='email',
                        help='Email address for notifying you the start, end, and abort of you job.')
    parser.add_argument('-s', type=str, dest='scheduler',
                        help='Name of the scheduler on your cluster, '
                             'e.g., PBS (or QSUB) or SBATCH (or SLURM), case insensitive.')
    parser.add_argument('-t', type=int, dest='time',
                        help='Time (in integer hours) for running your job, default: 4.',
                        default=4)
    parser.add_argument('-m', type=int, dest='memory',
                        help='Amount of memory (in GB) for all cores needed for your job, default: 32.',
                        default=32)
    parser.add_argument('-n', type=int, dest='cores',
                        help='Number of cores can be used for your job, default: 8.',
                        default=8)
    parser.add_argument('--dry-run', action='store_true', dest='dry',
                        help='Print out steps and files involved in each step without actually '
                             'running the pipeline.')
    
    args = parser.parse_args()
    global FASTQS
    FASTQS = [os.path.abspath(fastq) for fastq in args.FASTQ]
    if args.dry:
        ruffus.pipeline_printout()
    else:
        outdir = args.outdir
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        os.chdir(outdir)
    
        global logger
        log = f'clip_{datetime.today().strftime("%Y-%m-%d_%H:%M:%S")}.log'
        handler = logging.FileHandler(filename=log, mode='w')
        handler.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] %(message)s',
                                               datefmt='%Y-%m-%d %H:%M:%S'))
        handler.setLevel(logging.DEBUG)
        logger.addHandler(handler)
    
        global CORES
        CORES = args.cores if args.cores else CORES
        
        global GENOME
        GENOME = args.genome if args.genome else GENOME
        
        global RNA
        GENOME = args.genome if args.genome else GENOME
        
        scheduler = args.scheduler
        if scheduler:
            fastqs = ' \\\n    '.join(args.fastqs)
            schedule(scheduler, fastqs, outdir, args.rna,
                     args.cores, args.memory, args.time, args.jobname, args.email)
        else:
            ruffus.pipeline_run(multiprocess=CORES, verbose=1)


if __name__ == '__main__':
    main()
