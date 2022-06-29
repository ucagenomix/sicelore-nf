#! /usr/bin/env nextflow

params.nbcpu 	= "20"
params.project 	= "sicelore"
params.baminput = "/export/data2/000-10x-visium/000-MOB/GEUS10xAttributes.umifound_onepercent.bam"
params.samtools = "/export/apps/bin/samtools"
params.java 	= "/bin/java"
params.sicelore = "/export/data2/lebrigand/sicelore/Jar/Sicelore-2.0.jar"
params.classpath= "-cp </export/data2/lebrigand/sicelore/Jar/lib/*>"
params.minimap2 = "/export/apps/bin/minimap2"
workDir 		= '/export/data/analysis/NEXTFLOW'
csv				= ''
meta			= ''
params.juncbed	= '/export/data/references/mapper_indexes/GRCm39/gencode/gencode.vM27.annotation.bed'
params.mmi		= '/export/data/references/mapper_indexes/GRCm39/minimap/GRCm39.mmi'
params.refflat	= '/export/data/references/mapper_indexes/GRCm39/gencode/gencode.vM27.annotation.refflat.txt'
params.cellscsv = '/export/data2/000-10x-visium/mob164315.barcodes.tsv'

workflow {
	
	CHROMNAMES(params.baminput) | splitText | map{it -> it.trim()} | SPLITBAM | CONSENSUS_CALL | CONCAT | collectFile | CONSENSUS_DEDUP | CONSENSUS_MAPPING | CONSENSUS_ADDTAGS | CONSENSUS_ADDGENE | ISOMATRIX 	
	
	//CONSENSUS_METRICS(CONSENSUS_MAPPING.out.molbam)

}

process CHROMNAMES {
    input:
	path(bam)
    
    output:
    path("chromo.csv")
    
    """
	$params.samtools view -H $params.baminput | grep SQ | awk '{ print \$2 }' | sed 's/SN://' | grep -v 'ERCC\\|SIRV\\|phiX174' > chromo.csv
	"""
}

process SPLITBAM {
    
    input:
 	val(chromo)
 	
 	output:
 	path("*.bam")	, emit: bam
 	
    """
    echo processing chromosome $chromo
    
    $params.samtools view -Sb $params.baminput $chromo -o chromosome.bam
	$params.samtools index -@ 10 chromosome.bam
    """
}

process CONSENSUS_CALL {
    
    input:
 	path(bam)
 	
 	output:
 	path("*.fq")	, emit: fq
 	
    """
    echo CONSENSUS CALLING
    $params.java -jar -Xmx44g $params.sicelore ComputeConsensus T=$params.nbcpu I=$bam O=perchrom.fq VALIDATION_STRINGENCY=SILENT
    """
}


process CONCAT {

  publishDir '/export/data/analysis/NEXTFLOW/'

  input:
  path x
  
  output:
  path 'all.fq'	, emit: cons
  
  script:
  """
  < $x cat > all.fq
  """
}

process CONSENSUS_DEDUP {
    
    input:
 	path(fq)
 	
 	output:
 	path("*.fq")	, emit: dedup
 	
    """
    echo DeduplicateMolecule
 	$params.java -jar -Xmx44g $params.sicelore DeduplicateMolecule I=$fq O=dedup.fq SELECT=true  VALIDATION_STRINGENCY=SILENT
    """
}

process CONSENSUS_MAPPING {
    
    input:
 	path(fq)
 	
 	output:
 	path("*.bam")	, emit: molbam
 	
    """
    echo MAPPING
	$params.minimap2 -ax splice -uf --sam-hit-only -t $params.nbcpu --junc-bed $params.juncbed $params.mmi $fq > molecules.sam
	$params.samtools view -Sb -@ 20 molecules.sam -o deleted.bam
	$params.samtools sort -@ 20 deleted.bam -o molecules.bam
	$params.samtools index -@ 20 molecules.bam
	rm -f deleted.bam molecules.sam
    """
}

process CONSENSUS_ADDTAGS {
    
    input:
 	path(bam)
 	
 	output:
 	path("*.bam")	, emit: tagbam
 	
    """
    echo AddBamMoleculeTags
	$params.java -jar -Xmx44g $params.sicelore AddBamMoleculeTags I=$bam O=tags.bam
	$params.samtools index -@ 10 tags.bam
    """
}

process CONSENSUS_METRICS {
    
    input:
 	path(bam)
 	
 	output:
 	path("*.txt")	, emit: metrics
 	
    """
    echo GetMoleculeMetrics
    $params.java -jar -Xmx44g $params.sicelore GetMoleculeMetrics I=$bam O=metrics.txt
    """
}

process CONSENSUS_ADDGENE {
    
    input:
 	path(bam)
 	
 	output:
 	path("*.bam")	, emit: genebam
 	
    """
    echo AddBamMoleculeTags
	$params.java -jar -Xmx44g $params.sicelore AddGeneNameTag I=$bam O=molecules.tags.GE.bam REFFLAT=$params.refflat GENETAG=GE ALLOW_MULTI_GENE_READS=true USE_STRAND_INFO=true VALIDATION_STRINGENCY=SILENT
	$params.samtools index -@ 10 molecules.tags.GE.bam
    """
}

process ISOMATRIX {
    
    input:
 	path(bam)
 	
 	output:
 	path("*.txt")	, emit: isotxt
 	
    """
    echo ISOMATRIX
	$params.java -jar -Xmx44g $params.sicelore IsoformMatrix DELTA=2 METHOD=STRICT ISOBAM=true GENETAG=GE I=$bam REFFLAT=$params.refflat CSV=$params.cellscsv OUTDIR=./ PREFIX=sicmol VALIDATION_STRINGENCY=SILENT
    """
}














