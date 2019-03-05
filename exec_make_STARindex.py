#### conda config --add channels defaults ####
#### conda config --add channels bioconda ####
#### conda config --add channels conda-forge ####
#### conda install -c bioconda star ####
#### conda install -c bioconda samtools ####
import os
import glob

## configuration ##
species_list = ['hg38','mm10','sarsca1.0']
species = species_list[0]

root_dir = os.getcwd()
STARidx_dir =  os.path.join('/Users/petadimensionlab/ws/ref/STARidx/'+species)
genome_dir = '/Users/petadimensionlab/ws/ref/genome/'

if not os.path.exists(STARidx_dir):
	cmd = 'mkdir %s' % (STARidx_dir)
	os.system(cmd)

## options ##
thread_num = 12
genomeFastaFile = os.path.join(genome_dir,species+'.fa')
runMode = 'genomeGenerate'
sjdbGTFfile = os.path.join('/Users/petadimensionlab/ws/ref/annotation',species+'.gtf') 
sjdbOverhang = str(100)


## build index ##
cmd = 'STAR --runThreadN %s --runMode %s --genomeDir %s --genomeFastaFiles %s --sjdbGTFfile %s --sjdbOverhang %s' % (int(thread_num),runMode,STARidx_dir,genomeFastaFile,sjdbGTFfile,sjdbOverhang)
#print(cmd)
os.system(cmd)

## Change the access permission ##
items = glob.glob('*.ht2')
for item in items:
    cmd = 'chmod +x %s' % (item)
    os.system(cmd)
    srcf = os.path.join(root_dir,item)
    destf = os.path.join(STARidx_dir,item)
    cmd = 'mv %s %s' % (srcf,destf)
    os.system(cmd)
