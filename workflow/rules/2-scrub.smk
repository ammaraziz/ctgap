rule scrub:
	input:
		r1 = rules.fastp.output.r1,
		r2 = rules.fastp.output.r2,
	output:
		r1 = OUTDIR / "{sample}" / "scrub" / "{sample}_trim_scrub_R1.fastq.gz",
		r2 = OUTDIR / "{sample}" / "scrub" / "{sample}_trim_scrub_R2.fastq.gz",
		r1tmp = OUTDIR / "{sample}" / "scrub" / "{sample}_scrub_first_R1.fastq.gz",
		r2tmp = OUTDIR / "{sample}" / "scrub" / "{sample}_scrub_first_R2.fastq.gz",
		status = OUTDIR / "{sample}" / "status" / "scrub.{sample}.txt"
	params:
		db = KRAKENDB,
		dbname = KRAKENDB.stem,
		human = HUMANREF,
		minlen = 50,
		kraken_taxa_extract = 51291, #Chlamydiales
		workdir = directory(OUTDIR / "{sample}" / "scrub" / "scrubby_temp"),
	threads: config["threads"]["scrubby"]
	conda: "../envs/scrub.yaml"
	log: OUTDIR / "{sample}" / "log" / "scrub.{sample}.log"
	benchmark: OUTDIR / "{sample}" / "benchmark" / "scrub.{sample}.txt"
	shell:"""
	scrubby scrub-reads \
	-i {input.r1} {input.r2} \
	-o {output.r1tmp} {output.r2tmp} \
	--kraken-db {params.db} \
	--kraken-taxa "Archaea Eukaryota Holozoa Nucletmycea" \
	--min-len {params.minlen} \
	--minimap2-index {params.human} \
	--kraken-threads {threads} \
	--keep \
	--workdir {params.workdir:q} 2> {log}	

	echo -e "\nScrubby Kraken Extract \n" >> {log}

	scrubby scrub-kraken \
	-i {output.r1tmp} {output.r2tmp} \
	-o {output.r1} {output.r2} \
	--extract \
	--kraken-taxa {params.kraken_taxa_extract} \
	--kraken-reads {params.workdir}/0-{params.dbname}.kraken \
	--kraken-report {params.workdir}/0-{params.dbname}.report 2>> {log}

	touch {output.status}
	"""
