rule index:
	message: "Indexing Ct references"
	input:
		reference = CTREF
	output:
		status = OUTDIR / "status" / "ctReferences" / "ctReference.status"
	params:
		prefix = "resources/ctReferences/all"
	threads: config["threads"]["bowtieindex"]
	log: OUTDIR / "log" / "ctReference.index.log"
	conda: "../envs/bowtie.yaml"
	shell:"""
	bowtie2-build --threads {threads} {input.reference} {params.prefix} &> {log}

	touch {output.status}
	"""

rule bowtie:
	input:
		r1 = rules.scrub.output.r1,
		r2 = rules.scrub.output.r2,
		status = rules.index.output.status,
	output:
		bam = OUTDIR / "{sample}" / "aligned" / "{sample}.bam",
		coverage = OUTDIR / "{sample}" / "coverage.{sample}.tsv",
		status = OUTDIR / "status" / "aligned.{sample}.txt",
	params:
		prefix = rules.index.params.prefix
	threads: config["threads"]["bowtie"]
	log: OUTDIR / "{sample}" / "log" / "bowtie2.{sample}.log"
	conda: "../envs/bowtie.yaml"
	shell:"""
	bowtie2 -x \
	{params.prefix} \
	-1 {input.r1} \
	-2 {input.r2}\
	--threads {threads} | \
	samtools view -bSh - | \
	samtools sort -@{threads} \
	-o {output.bam} 2> {log} 1> {log}

	samtools index {output.bam}
	samtools coverage {output.bam} > {output.coverage}

	touch {output.status}
	"""

rule bamtofastq:
	input:
		bam = rules.bowtie.output.bam,
	output:
		r1 = OUTDIR / "{sample}" / "aligned" / "{sample}_refg_r1.fastq.gz",
		r2 = OUTDIR / "{sample}" / "aligned" / "{sample}_refg_r2.fastq.gz",
		status = OUTDIR / "status" / "aligned.{sample}.txt",
	threads: config["threads"]["bedtools"]
	log: OUTDIR / "{sample}" / "log" / "bowtie2.{sample}.log"
	conda: "../envs/bowtie.yaml"
	shell:"""
	bamToFastq \
	-i {input.bam} \
	-fq {output.r1} \
	-fq2 {output.r2}

	touch {output.status}
	"""

use rule shovill as mapping_shovill with:
	input:
		r1 = rules.bamtofastq.output.r1,
		r2 = rules.bamtofastq.output.r2,
	output:
		status = OUTDIR / "status" / "mapping.shovill.{sample}.txt",
		outdir = directory(OUTDIR / "{sample}" / "mapping.shovill"),
		contig = OUTDIR / "{sample}" / "mapping.shovill" / "contigs.fa",
