#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

workflow {
	CELLBARCODES()
	CHROMNAMES(params.readTagged) | splitText | map{it -> it.trim()} | SPLITBAM | CONSENSUS_CALL | CONCAT | collectFile | CONSENSUS_DEDUP | CONSENSUS_MAPPING | CONSENSUS_ADDTAGS | CONSENSUS_ADDGENE
	ISOMATRIX(CONSENSUS_ADDGENE.out.genebam, CONSENSUS_ADDGENE.out.genebai, CELLBARCODES.out.csv)
}

process CELLBARCODES {
    cpus 2
    
    output:
    path 'cellbarcodes.csv'	, emit: csv
    
    // alternative solution to test
    //java -jar picard.jar SplitSamByNumberOfReads \
    // I=paired_unmapped_input.bam \
    // OUTPUT=out_dir \ 
    // TOTAL_READS_IN_INPUT=800000000 \ 
    // SPLIT_TO_N_READS=48000000
    
    
    """
	$params.java -jar -Xmx4g $params.sicelore SelectValidCellBarcode I=$params.barcodeassigned O=cellbarcodes.csv MINUMI=$params.MINUMI ED0ED1RATIO=$params.ED0ED1RATIO
	"""
}

process CHROMNAMES {
    input:
	path(bam)
    
    output:
    path("chromo.csv")
    
    // alternative solution to test
    //java -jar picard.jar SplitSamByNumberOfReads \
    // I=paired_unmapped_input.bam \
    // OUTPUT=out_dir \ 
    // TOTAL_READS_IN_INPUT=800000000 \ 
    // SPLIT_TO_N_READS=48000000
    
    
    """
	$params.samtools view -H $params.readTagged | grep SQ | awk '{ print \$2 }' | sed 's/SN://' | grep -v 'ERCC\\|SIRV\\|phiX174' > chromo.csv
	"""
}

process SPLITBAM {
    cpus 2
    
    input:
 	val(chromo)
 	
 	output:
 	path 'chromosome.bam'	, emit: bam
 	
    """
    $params.samtools view -Sb $params.readTagged $chromo -o chromosome.bam
	$params.samtools index -@ 10 chromosome.bam
    """
}

process CONSENSUS_CALL {  
    input:
 	path(bam)
 	
 	output:
 	path 'chr.fq'	, emit: fq
 	
    """
    $params.java -jar -Xmx44g $params.sicelore ComputeConsensus T=$params.max_cpus I=$bam O=chr.fq CELLTAG=$params.CELLTAG UMITAG=$params.UMITAG GENETAG=$params.GENETAG TSOENDTAG=$params.TSOENDTAG POLYASTARTTAG=$params.POLYASTARTTAG CDNATAG=$params.CDNATAG USTAG=$params.USTAG RNTAG=$params.RNTAG MAPQV0=$params.MAPQV0 TMPDIR=$params.TMPDIR VALIDATION_STRINGENCY=SILENT MAXREADS=$params.MAXREADS MINPS=$params.MINPS MAXPS=$params.MAXPS DEBUG=$params.DEBUG
    """
}

process CONCAT {
    cpus 2
    
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
 	$params.java -jar -Xmx44g $params.sicelore DeduplicateMolecule I=$fq O=consensus_dedup.fq SELECT=true VALIDATION_STRINGENCY=SILENT
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
	$params.java -jar -Xmx44g $params.sicelore AddBamMoleculeTags I=$bam O=molecules.tags.bam CELLTAG=$params.CELLTAG UMITAG=$params.UMITAG RNTAG=$params.RNTAG
	$params.samtools index -@ 10 molecules.tags.bam
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
	$params.java -jar -Xmx44g $params.sicelore AddGeneNameTag I=$bam O=molecules.tags.GE.bam REFFLAT=$params.REFFLAT GENETAG=$params.GENETAG ALLOW_MULTI_GENE_READS=$params.ALLOW_MULTI_GENE_READS USE_STRAND_INFO=$params.USE_STRAND_INFO VALIDATION_STRINGENCY=SILENT
	$params.samtools index -@ 10 molecules.tags.GE.bam
    """
}

process ISOMATRIX {
    publishDir "${params.outdir}/${params.matrixdir}", mode: 'copy'
    
    input:
 	path(bam)
 	path(bai)
	path(csv)
 	
 	output:
 	path("*.txt")			, emit: isotxt
 	path("*_isobam.bam")	, emit: isobam
 	
    """
	$params.java -jar -Xmx44g $params.sicelore IsoformMatrix I=$bam REFFLAT=$params.REFFLAT CSV=$csv OUTDIR=./ PREFIX=$params.PREFIX CELLTAG=$params.CELLTAG UMITAG=$params.UMITAG GENETAG=$params.GENETAG TSOENDTAG=$params.TSOENDTAG POLYASTARTTAG=$params.POLYASTARTTAG CDNATAG=$params.CDNATAG USTAG=$params.USTAG RNTAG=$params.RNTAG MAPQV0=$params.MAPQV0 DELTA=$params.DELTA METHOD=$params.METHOD ISOBAM=$params.ISOBAM AMBIGUOUS_ASSIGN=$params.AMBIGUOUS_ASSIGN VALIDATION_STRINGENCY=SILENT
    """
}














