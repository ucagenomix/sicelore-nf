#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

workflow {
	CHROMNAMES(params.baminput) | splitText | map{it -> it.trim()} | SPLITBAM | CONSENSUS_CALL | CONCAT | collectFile | CONSENSUS_DEDUP | CONSENSUS_MAPPING | CONSENSUS_ADDTAGS | CONSENSUS_ADDGENE | ISOMATRIX 	
	CONSENSUS_METRICS(CONSENSUS_MAPPING.out.molbam)
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
    $params.samtools view -Sb $params.baminput $chromo -o chromosome.bam
	$params.samtools index -@ 10 chromosome.bam
    """
}

process CONSENSUS_CALL {
    input:
 	path(bam)
 	
 	output:
 	path 'chr.fq'	, emit: fq
 	
    """
    $params.java -jar -Xmx44g $params.sicelore ComputeConsensus T=$params.max_cpus I=$bam O=chr.fq VALIDATION_STRINGENCY=SILENT
    """
}

process CONCAT {
  input:
  path x
  
  output:
  path 'consensus_all.fq'	, emit: cons
  
  script:
  """
  < $x cat > consensus_all.fq
  """
}

process CONSENSUS_DEDUP {
    publishDir "${params.outdir}/${params.consdir}", mode: 'copy'
    
    input:
 	path(fq)
 	
 	output:
 	path 'consensus_dedup.fq'	, emit: dedup
 	
    """
 	$params.java -jar -Xmx44g $params.sicelore DeduplicateMolecule I=$fq O=consensus_dedup.fq SELECT=true  VALIDATION_STRINGENCY=SILENT
    """
}

process CONSENSUS_MAPPING {
    publishDir "${params.outdir}/${params.bamdir}", mode: 'copy'
    
    input:
 	path(fq)
 	
 	output:
 	path 'molecules.bam'		, emit: molbam
 	path 'molecules.bam.bai'	, emit: molbai
 	
    """
	$params.minimap2 -ax splice -uf --sam-hit-only -t $params.max_cpus --junc-bed $params.juncbed $params.mmi $fq > molecules.sam
	$params.samtools view -Sb -@ 20 molecules.sam -o deleted.bam
	$params.samtools sort -@ 20 deleted.bam -o molecules.bam
	$params.samtools index -@ 20 molecules.bam
    """
}

process CONSENSUS_ADDTAGS {
    publishDir "${params.outdir}/${params.bamdir}", mode: 'copy'
   
    input:
 	path(bam)
 	path(bai)
 	
 	output:
 	path 'molecules.tags.bam'		, emit: tagbam
 	path 'molecules.tags.bam.bai'	, emit: tagbai
 	
    """
	$params.java -jar -Xmx44g $params.sicelore AddBamMoleculeTags I=$bam O=molecules.tags.bam
	$params.samtools index -@ 10 molecules.tags.bam
    """
}

process CONSENSUS_METRICS {
    input:
 	path(bam)
 	
 	output:
 	path("*.txt")	, emit: metrics
 	
    """
    $params.java -jar -Xmx44g $params.sicelore GetMoleculeMetrics I=$bam O=metrics.txt
    """
}

process CONSENSUS_ADDGENE {
    publishDir "${params.outdir}/${params.bamdir}", mode: 'copy'
    
    input:
 	path(bam)
 	path(bai)
 	
 	output:
 	path 'molecules.tags.GE.bam'		, emit: genebam
 	path 'molecules.tags.GE.bam.bai'	, emit: genebai
 	
    """
	$params.java -jar -Xmx44g $params.sicelore AddGeneNameTag I=$bam O=molecules.tags.GE.bam REFFLAT=$params.refflat GENETAG=GE ALLOW_MULTI_GENE_READS=true USE_STRAND_INFO=true VALIDATION_STRINGENCY=SILENT
	$params.samtools index -@ 10 molecules.tags.GE.bam
    """
}

process ISOMATRIX {
    publishDir "${params.outdir}/${params.matrixdir}", mode: 'copy'
    
    input:
 	path(bam)
 	path(bai)
 	
 	output:
 	path("*.txt")	, emit: isotxt
 	
    """
	$params.java -jar -Xmx44g $params.sicelore IsoformMatrix DELTA=2 METHOD=STRICT ISOBAM=true GENETAG=GE I=$bam REFFLAT=$params.refflat CSV=$params.cellscsv OUTDIR=./ PREFIX=sicmol VALIDATION_STRINGENCY=SILENT
    """
}














