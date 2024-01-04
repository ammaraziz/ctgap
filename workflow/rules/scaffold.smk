rule scaffold:
	input: 
		scaffold = rules.shovill.output.contig
	output:
		outdir = directory(OUTDIR / "{sample}" / "scaffold"),
		scaffold = OUTDIR / "{sample}" / "scaffold" / "ragtag.scaffold.fasta",
		status = OUTDIR / "status" / "scaffold.{sample}.txt",
	params:
		ref = SCAFFOLDREF,
		min_len = 500
	threads: config['threads']['ragtag']
	conda: "../envs/scaffold.yaml"
	log: OUTDIR / "{sample}" / "log" / "scaffold.{sample}.log"
	benchmark: OUTDIR / "{sample}" / "benchmark" / "scaffold.{sample}.txt"
	shell:"""
	ragtag.py scaffold \
	{params.ref} \
	{input.scaffold} \
	-o {output.outdir} \
	-f {params.min_len} \
	-r \
	-u \
	-t {threads} 2> {log}

    touch {output.status}
	"""

