rule collate_coverage:
	input:
		status = expand(OUTDIR / "status" / "aligned.{sample}.txt", sample = SAMPLES),
		coverages = expand(OUTDIR / "{sample}" / "coverage.{sample}.tsv", sample = SAMPLES),
	output:
		coverages = OUTDIR / "all.coverage.tsv",
		status = OUTDIR / "status" / "collage.coverage.txt",
	threads: 1
	conda: "../envs/misc.yaml"
	shell:"""
	csvtk concat -C $ {input.coverages} >  {output.coverages}

	touch {output.status}
	"""
	
rule collate_blast:
	input:
		status_shovill = expand(OUTDIR / "status" / "shovill.{sample}.txt", sample = SAMPLES),
	output:
		tsv = OUTDIR / "all.blast.tsv",
		status = OUTDIR / "status" / "collate.blast.txt",
	params:
		outdir = OUTDIR,
		pattern = "**/blast/*.tab",
	threads: 1
	shell:"""
	echo -e "query\tsubject\tpident\tlength\tmismatch\tgapopen\tquery_start\tquery_end\tsubject_start\tsubject_end\tevalue\tbitscore" > {output.tsv}
	grep "" {params.outdir}/{params.pattern} >> {output.tsv}

	touch {output.status}
	"""