## 01: Install parallel-fastq-dump, fastp, psamtools, hisat2
## 02: conda install -c bioconda -y parallel-fastq-dump
## 03: conda install -c bioconda fastp
## 04: do not install samtools from bioconda since the following error sometimes cannot be resolved:
# dyld: Library not loaded: @rpath/libcrypto.1.0.0.dylib ##
#  04: brew install samtools
# see also https://github.com/conda-forge/conda-forge.github.io/issues/701 ##
## 05: conda install -c bioconda star

import os

is_sra = 'yes'
is_paired = 'yes'
is_qualitycheck = 'yes'
is_delete_FASTQ = 'no'
thread_num = 12 # thread number
species = 'hg38'
annotation_file = species+'.gtf'
annotation = os.path.join('/Users/petadimensionlab/ws/ref/annotation',annotation_file)
STARidx_dir = '/Users/petadimensionlab/ws/ref/STARidx'

root_dir = '/Users/petadimensionlab/Downloads/test01/'
input_dir = os.path.join(root_dir,'input')
output_dir = os.path.join(root_dir,'output')

fr = open(os.path.join(root_dir,'fastq_files.txt'),'r').readlines()
for line in fr:
    line = line.replace('\n','')
    lst = line.split(',')
    filename = lst[0]
    samplename = filename
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    msg = '%s is now processing...' % (samplename)
    print( msg )
    cmd = 'python exec_sra2readcount_STAR.py %s %s %s %s %s %s %s %s %s %s %s' % (is_sra,is_paired,is_qualitycheck,is_delete_FASTQ,thread_num,samplename,species,annotation,input_dir,output_dir,STARidx_dir)
    print( cmd )
    os.system(cmd)

