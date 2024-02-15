rule scaffold:
	input: 
		contigs = rules.shovill.output.contig
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
		scaffold = rules.scaffold.output.scaffold,
		r1 = rules.scrub.output.r1,
		r2 = rules.scrub.output.r2,
	output:
		filled = OUTDIR / "{sample}" / "gap2seq" / "filled.fasta",
		status = OUTDIR / "status" / "gap2seq.{sample}.txt",
	threads: config['threads']['gap2seq']
	conda: "../envs/scaffold.yaml"
	log: OUTDIR / "{sample}" / "log" / "gap2seq.{sample}.log"
	benchmark: OUTDIR / "{sample}" / "benchmark" / "gap2seq.{sample}.txt"
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
