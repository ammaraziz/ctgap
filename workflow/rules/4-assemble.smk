rule shovill:
	input:
		r1 = rules.scrub.output.r1,
		r2 = rules.scrub.output.r2,
	output:
		status = OUTDIR / "status" / "shovill.{sample}.txt",
		outdir = directory(OUTDIR / "{sample}" / "shovill"),
		contig = OUTDIR / "{sample}" / "shovill" / "contigs.fa",
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

# are the reads here interleaved? bbsplit outputs r1 and r2.
rule spades:
	input:
		r1 = rules.scrub.output.r1,
		r2 = rules.scrub.output.r2,
	output:
		status = OUTDIR / "status" / "spades.{sample}.txt",
		assembly = directory(OUTDIR / "{sample}" / "spades"),
	threads: config["threads"]["spades"]
	conda: "../envs/shovill.yaml"
	log: OUTDIR / "{sample}" / "log" / "spades.{sample}.logs"
	benchmark: OUTDIR / "{sample}" / "benchmark" / "spades.{sample}.txt"
	shell:"""
	spades.py \
	--careful \
	--threads {threads} \
	--pe1-1 {params.clean_r1} \
	--pe1-2 {params.clean_r2} \
	-o {output.assembly} 2> {log}
	"""