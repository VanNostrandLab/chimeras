# chimeras

## Usage

Inside chimeras directory, issue the following command to check the usage:

```shell
$ ./chimeras -h

usage: chimeras [-h] [-o OUTDIR] [-g GENOME] [-r RNA] [-j JOBNAME] 
                [-e EMAIL] [-s SCHEDULER] [-t TIME] [-m MEMORY] 
                [-n CORES] [--dry-run] FASTQ [FASTQ ...]

A pipeline designed to identify miRNA-target chimeras.

positional arguments:
  FASTQ         Path to one or multiple FASTQ files (separated by space).

optional arguments:
  -h, --help    show this help message and exit
  -o OUTDIR     Path to the output directory (if not exist, will try to 
                create it first).
  -g GENOME     Path to STAR genome index directory.
  -r RNA        Path to mature RNA FASTA file (generated using miRBase and 
                U needs to be replaced by T).
  -j JOBNAME    Name of your job, default: chimeras
  -e EMAIL      Email address for notifying you the start, end, and abort 
                of you job.
  -s SCHEDULER  Name of the scheduler on your cluster, e.g., PBS (or QSUB) 
                or SBATCH (or SLURM), case insensitive.
  -t TIME       Time (in integer hours) for running your job, default: 4.
  -m MEMORY     Amount of memory (in GB) for all cores needed for your job, 
                default: 32.
  -n CORES      Number of cores can be used for your job, default: 8.
  --dry-run     Print out steps and files involved in each step without 
                actually running the pipeline.

```