rule shovill:
	message: "Assembling Sample {wildcards.sample}"
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
	--force --cpus {threads} 2>> {log} 1>> {log}

	touch {output.status}
	"""

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

rule blast_ompa:
	input:
		contig = rules.gapfiller.output.filled,
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
	-html \
	-outfmt {params.outfmt} \
	-out {output.tab} 2> {log}

	touch {output.status}
	"""

rule mlst:
	input:
		rules.gapfiller.output.filled
	output:
		generic = OUTDIR / "{sample}" / "mlst" / "{sample}.genome.chlamydiales.mlst.txt",
		ct = OUTDIR / "{sample}" / "mlst" / "{sample}.genome.ctrachomatis.mlst.txt",
		plasmid =  OUTDIR / "{sample}" / "mlst" / "{sample}.genome.plasmid.mlst.txt",
		status = OUTDIR / "status" / "mlst.{sample}.txt",
	log: OUTDIR / "{sample}" / "log" / "mlst.{sample}.log"
	benchmark: OUTDIR / "{sample}" / "benchmark" / "mlst.{sample}.txt"
	conda: "../envs/mlst.yaml"
	params:
		dbgeneric = MLSTDB / "chlamydiales",
		dbct = MLSTDB / "c.trachomatis",
		dbplasmid = MLSTDB / "plasmid",
	threads: config["threads"]["mlst"]
	shell:"""
	echo -e "chlamydiales\n" >> {log}
	claMLST search \
	{params.dbgeneric} \
	{input} > {output.generic} 2>> {log}

	echo -e "\nctrachomatis\n" >> {log}
	claMLST search \
	{params.dbct} \
	{input} > {output.ct} 2>> {log}

	echo -e "\nplasmid\n" >> {log}
	claMLST search \
	{params.dbplasmid} \
	{input} > {output.plasmid} 2>> {log}

	touch {output.status}
	"""

rule denovo_collate_blast:
	input:
		status_blast = expand(OUTDIR / "status" / "blastn.{sample}.txt", sample = SAMPLES),
	output:
		tsv = OUTDIR / "denovo.blast.tsv",
		status = OUTDIR / "status" / "denovo.collate.blast.txt",
	params:
		outdir = OUTDIR,
		pattern = "**/blast/*.tab",
	threads: 1
	shell:"""
	echo -e "query\tsubject\tpident\tlength\tmismatch\tgapopen\tquery_start\tquery_end\tsubject_start\tsubject_end\tevalue\tbitscore" > {output.tsv}
	grep "" {params.outdir}/{params.pattern} >> {output.tsv}

	touch {output.status}
	"""

rule denovo_collate_mlst:
	input:
		generic = expand(OUTDIR / "{sample}" / "mlst" / "{sample}.genome.chlamydiales.mlst.txt", sample = SAMPLES),
		ct = expand(OUTDIR / "{sample}" / "mlst" / "{sample}.genome.ctrachomatis.mlst.txt", sample = SAMPLES),
		plasmid =  expand(OUTDIR / "{sample}" / "mlst" / "{sample}.genome.plasmid.mlst.txt", sample = SAMPLES),
	output:
		generic = OUTDIR / "mlst.generic.results.tsv",
		cd = OUTDIR / "mlst.ct.results.tsv",
		plasmid = OUTDIR / "mlst.plasmid.results.tsv",
		status = OUTDIR / "status" / "mlst.collate.txt"
	conda: "../envs/misc.yaml"
	threads: 1
	shell:"""
	csvtk concat {input.generic} -o {output.generic}
	csvtk concat {input.ct} -o {output.ct}
	csvtk concat {input.plasmid} -o {output.plasmid}

	touch {output.status}
	"""

