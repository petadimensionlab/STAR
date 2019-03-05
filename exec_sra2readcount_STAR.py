import os
import sys
import re
import shutil
REPanno = '\.g[t|f]f'

argvs = sys.argv
argc = len(argvs)
if argc!=12:
    msg = 'Usage: # python %s <is_sra: yes or no> <is_paired: yes or no> <is_qualitycheck: yes or no> <is_delete_FASTQ: yes or no> <thread-number> <sample name> <species> <annotation> <input-directory> <output-directory> <STAR index directory)' % (argvs[0])
    print( msg )
    quit()

## args ##
is_sra = argvs[1] # check the type of input-data (sra or fastq)
is_paired = argvs[2]
is_qualitycheck = argvs[3]
is_delete_FASTQ = argvs[4]
thread_num = argvs[5] # thread number
samplename = argvs[6]
species = argvs[7]
annotation = argvs[8]
input_dir = argvs[9]
output_dir = argvs[10]
STARidx_dir = argvs[11]

## tools ##
fastq_dump = 'parallel-fastq-dump'
fastp = 'fastp'
samtools = 'samtools'
ws_dir = os.getcwd()
bam2rc = os.path.join(ws_dir,'exec_bam2readcount_STAR.R')

## directories ##
sample_dir = os.path.join(output_dir,samplename)
genome_dir =  os.path.join(STARidx_dir,species)

## STAR options ##
runThreadN = str(thread_num)
readFilesCommand = 'gunzip -c'
outFilterMultimapScoreRange = str(1)
outFilterMultimapNmax = str(10)
outFilterMismatchNmax = str(5)
alignIntronMax = str(500000)
alignMatesGapMax = str(1000000)
sjdbScore = str(2)
alignSJDBoverhangMin = str(1)
genomeLoad = 'NoSharedMemory'
limitBAMsortRAM = str(0)
outFilterMatchNminOverLread = str(0.33) 
outFilterScoreMinOverLread = str(0.33) 
sjdbOverhang = str(100)
outSAMstrandField = 'intronMotif'
outSAMattributes = 'NH HI NM MD AS XS'
outSAMunmapped = 'Within'
outSAMtype = 'BAM Unsorted'
outSAMheaderHD = '@HD VN:1.4'

## create a directory for given sample ID ##
if os.path.exists(sample_dir):
    msg = 'Directory %s already exists.' % (samplename)
    print( msg )
else:
    os.mkdir(sample_dir)

#### Execution ####
if is_sra == 'yes':
    ## 1st step : SRA to FASTQ ##
    msg = '%s : SRA to FASTQ...' % (samplename)
    print( msg )
    if is_paired=='yes':
        #cmd = '%s --sra-id %s/%s.sra' % (fastq_dump,input_dir,samplename)
        #print( cmd )
        cmd = '%s --sra-id %s/%s.sra --threads %s --split-files --gzip' % (fastq_dump,input_dir,samplename,thread_num)
        print( cmd )
        os.system(cmd)
    else:
        cmd = '%s --sra-id %s/%s.sra --threads %s --gzip' % (fastq_dump,input_dir,samplename,thread_num)
        print( cmd )
        os.system(cmd)
else:
    msg = '%s : Skip conversion of SRA to FASTQ.' % (samplename)
    print( msg )
    msg = '%s : Copy FASTQ files into an output directory.' % (samplename)
    print( msg )
    if is_paired=='yes':
        file1 = '%s_1.fastq.gz' % (samplename)
        file2 = '%s_2.fastq.gz' % (samplename)
        if os.path.exists(os.path.join(input_dir,file1)):
            shutil.copy(os.path.join(input_dir,file1),os.path.join(sample_dir,file1))
        if os.path.exists(os.path.join(input_dir,file2)):
            shutil.copy(os.path.join(input_dir,file2),os.path.join(sample_dir,file2))
    else:
        file = '%s.fastq.gz' % (samplename)
        if os.path.exists(os.path.join(input_dir,file)):
            shutil.copy(os.path.join(input_dir,file),os.path.join(sample_dir,file))

## 2nd step : fastp ##
if is_qualitycheck == 'yes':
    msg = '%s : fastp...' % (samplename)
    print( msg )
    report_file = os.path.join(sample_dir,samplename+'.html')
    if is_paired =='yes':
        fwd_fastq = os.path.join(ws_dir,samplename+'_1.fastq.gz')
        rev_fastq = os.path.join(ws_dir,samplename+'_2.fastq.gz')
        trim_fwd_fastq = os.path.join(sample_dir,samplename+'_trim_paired_1.fastq.gz')
        trim_rev_fastq = os.path.join(sample_dir,samplename+'_trim_paired_2.fastq.gz')
        cmd = '%s -w %s -h %s -i %s -I %s -o %s -O %s' % (fastp,int(thread_num),report_file,fwd_fastq,rev_fastq,trim_fwd_fastq,trim_rev_fastq)
        os.system(cmd)
    else:
        fastq = os.path.join(sample_dir,samplename+'.fastq.gz')
        trim_fastq = os.path.join(ws_dir,samplename+'_trim.fastq.gz')
        cmd = '%s -w %s -h %s -i %s -o %s' % (fastp,int(thread_num),report_file,fastq,trim_fastq)
        print( cmd )
        os.system(cmd)

## 3rd step : STAR ##
msg = '%s : STAR...' % (samplename)
print( msg )
outFileNamePrefix = os.path.join(output_dir,samplename,samplename+'.')
if is_paired=='yes':
        fwd_fastq = os.path.join(output_dir,samplename,samplename+'_trim_paired_1.fastq.gz')
        rev_fastq = os.path.join(output_dir,samplename,samplename+'_trim_paired_2.fastq.gz')
        readFilesIn = fwd_fastq+' '+rev_fastq
        #print(readFilesIn)
else:
        fastq_file = os.path.join(output_dir,samplename+'trim_paired.fastq.gz')
        readFilesIn = fastq_file
        
cmd = 'STAR --runThreadN %s --readFilesCommand %s --genomeDir %s --readFilesIn %s --outFileNamePrefix %s --outFilterMultimapScoreRange %s --outFilterMultimapNmax %s --outFilterMismatchNmax %s --alignIntronMax %s --alignMatesGapMax %s --sjdbScore %s --alignSJDBoverhangMin %s --genomeLoad %s --limitBAMsortRAM %s --outFilterMatchNminOverLread %s --outFilterScoreMinOverLread %s --sjdbOverhang %s --outSAMstrandField %s --outSAMattributes %s --outSAMunmapped %s --outSAMtype %s --outSAMheaderHD %s' % (runThreadN,readFilesCommand,genome_dir,readFilesIn,outFileNamePrefix,outFilterMultimapScoreRange,outFilterMultimapNmax,outFilterMismatchNmax,alignIntronMax,alignMatesGapMax,sjdbScore,alignSJDBoverhangMin,genomeLoad,limitBAMsortRAM,outFilterMatchNminOverLread,outFilterScoreMinOverLread,sjdbOverhang,outSAMstrandField,outSAMattributes,outSAMunmapped,outSAMtype,outSAMheaderHD)
print(cmd)
os.system(cmd)

## 4th step : Extract read count from bam file ##
os.chdir(sample_dir)
msg = '%s : Bam to read count...' % (samplename)
print( msg )
cmd = 'R --vanilla --slave --args %s %s %s %s < %s' % (annotation,samplename,thread_num,sample_dir,bam2rc)
print( cmd )
os.system(cmd)

## 5th step : Delete copied FASTQ file ##
if is_delete_FASTQ == 'yes':
    os.chdir(sample_dir) # confirm that we are in the output directory
    msg = '%s : Delete copied/dumped FASTQ files...' % (samplename)
    print( msg )
    if is_paired =='yes':
        fastq1 = os.path.join(ws_dir,samplename+'_1.fastq.gz')
        fastq2 = os.path.join(ws_dir,samplename+'_2.fastq.gz')
        cmd = 'rm %s %s' % (fastq1,fastq2)
    else:
        fastq = '%s.fastq.gz' % (samplename)
        cmd = 'rm %s' % (fastq)
    os.system(cmd)
    ## delete intermediate sam file ##
    cmd = 'rm %s' % (samfile)
    os.system(cmd)

## finish ##
msg = '%s : finish!' % (samplename)
print( msg )
