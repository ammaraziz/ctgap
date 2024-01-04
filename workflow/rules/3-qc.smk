rule fastp_qc:
	input:
		r1 = rules.scrub.output.r1,
		r2 = rules.scrub.output.r2,
	output:
		html = OUTDIR / "{sample}" / "scrub" / "{sample}_scrub.html",
		json = OUTDIR / "{sample}" / "scrub" / "{sample}_scrub.json",
		status = OUTDIR / "status" / "trim.scrub.{sample}.txt",
	threads: config["threads"]["fastp"]
	conda: "../envs/trim.yaml"
	log: OUTDIR / "{sample}" / "log" / "trim.scrub.{sample}.log"
	shell:"""
	fastp \
	-i {input.r1} \
	-I {input.r2} \
	-h {output.html} \
	-j {output.json} \
	--thread {threads} 2> {log}

	touch {output.status}
	"""

rule multiqc:
	input:
		status_spades = expand(OUTDIR / "status" / "spades.{sample}.txt", sample = SAMPLES),
		status_shovill = expand(OUTDIR / "status" / "shovill.{sample}.txt", sample = SAMPLES)
	output:
		report = OUTDIR / "multiqc_report.html"
	params:
		fastp_dir = OUTDIR / "qc"
	log: OUTDIR / "log" / "multiqc.log",
	conda: "../envs/misc.yaml",
	shell:"""
	multiqc {params.fastp_dir} -n {output.report} 2> {log}
	"""

