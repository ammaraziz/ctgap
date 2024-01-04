# currently not used
# replaced by shovill/seqkt
rule normalise:
	input:
		r1 = rules.scrub.output.r1,
		r2 = rules.scrub.output.r2,
	output:
		r1 = OUTDIR / "{sample}" / "normalised" / "{sample}_trim_scrub_norm_R1.fastq.gz",
		r2 = OUTDIR / "{sample}" / "normalised" / "{sample}_trim_scrub_norm_R2.fastq.gz",
		status = OUTDIR / "status" / "normalise.{sample}.txt"
	params:
		target = 100,
		minimum = 5,
	conda: "../envs/normalise.yaml"
	log: OUTDIR / "{sample}" / "log" / "normalise.{sample}.log"
	threads: 10
	shell:"""
	bbnorm.sh \
	in1={input.r1} \
	in2={input.r2} \
	out1={output.r1} \
	out2={output.r2} \
	target={params.target} \
	threads={threads} \
	min={params.minimum} 2> {log}

	touch {output.status}
	"""