rule blast_ompa:
	input:
		contig = rules.shovill.output.contig,
	output:
		tab = OUTDIR / "{sample}" / "blast" / "blast.ompa.tab",
		status = OUTDIR / "status" / "blastn.{sample}.txt"
	log: OUTDIR / "{sample}" / "log" / "blastompa.{sample}.log"
	params:
		outfmt = 6,
		db = OMPABLASTDB,
		targets = 1
	threads: config["threads"]["blast"]
	conda: "../envs/misc.yaml"
	shell:"""
	blastn \
	-query {input.contig} \
	-db {params.db} \
	-max_target_seqs {params.targets} \
	-html  \
	-outfmt {params.outfmt} \
	-out {output.tab}

	touch {output.status}
	"""