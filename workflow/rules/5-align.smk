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

rule mapping_shovill:
	message: "Assembling Sample {wildcards.sample}"
	input:
		r1 = rules.bamtofastq.output.r1,
		r2 = rules.bamtofastq.output.r2,
	output:
		status = OUTDIR / "status" / "mapping.shovill.{sample}.txt",
		outdir = directory(OUTDIR / "{sample}" / "mapping.shovill"),
		contig = OUTDIR / "{sample}" / "mapping.shovill" / "contigs.fa",
	params:
		gsize = "1.04M",
		depth = 200,
	conda: "../envs/shovill.yaml"
	log: OUTDIR / "{sample}" / "log" / "shovill.{sample}.log"
	benchmark: OUTDIR / "{sample}" / "benchmark" / "shovill.{sample}.txt"
	threads: config["threads"]["shovill"]
	shell:"""
	shovill \
	--R1 {input.r1} \
	--R2 {input.r2} \
	--outdir {output.outdir} \
	--gsize {params.gsize} \
	--depth {params.depth} \
	--force \
	--cpus {threads} 2> {log}

	touch {output.status}
	"""

rule mapping_scaffold:
	input: 
		contigs = rules.mapping_shovill.output.contig
	output:
		outdir = directory(OUTDIR / "{sample}" / "mapping_scaffold"),
		scaffold = OUTDIR / "{sample}" / "mapping_scaffold" / "ragtag.scaffold.fasta",
		status = OUTDIR / "status" / "mapping_scaffold.{sample}.txt",
	params:
		ref = SCAFFOLDREF,
		min_len = 500
	threads: config['threads']['ragtag']
	conda: "../envs/scaffold.yaml"
	log: OUTDIR / "{sample}" / "log" / "mapping_scaffold.{sample}.log"
	benchmark: OUTDIR / "{sample}" / "benchmark" / "mapping_scaffold.{sample}.txt"
	shell:"""
	ragtag.py scaffold \
	{params.ref} \
	{input.contigs} \
	-o {output.outdir} \
	-f {params.min_len} \
	-r \
	-u \
	-t {threads} 2> {log}

    touch {output.status}
	"""

rule gapfiller:
	input:
		scaffold = rules.mapping_scaffold.output.scaffold,
		r1 = rules.scrub.output.r1,
		r2 = rules.scrub.output.r2,
	output:
		filled = OUTDIR / "{sample}" / "mapping_gap2seq" / "filled.fasta",
		status = OUTDIR / "status" / "mapping_gap2seq.{sample}.txt",
	threads: config['threads']['gap2seq']
	conda: "../envs/scaffold.yaml"
	log: OUTDIR / "{sample}" / "log" / "mapping_gap2seq.{sample}.log"
	benchmark: OUTDIR / "{sample}" / "benchmark" / "mapping_gap2seq.{sample}.txt"
	shell:"""
	EXITCODE=$(Gap2Seq.sh \
	--scaffolds {input.scaffold} \
	--filled {output.filled} \
	--reads {input.r1},{input.r2} \
	--nb-core {threads} 2>&1)

	# gap2seq will fail if filling fails.
	# capture error, copy scaffolded as filled.
	if [ $EXITCODE -eq 1 ]
    then
        cp {input.scaffold} {output.filled}
    fi

	touch {output.status}
	"""
