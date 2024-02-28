rule shovill:
	message: "Assembling Sample {wildcards.sample}"
	input:
		r1 = rules.scrub.output.r1,
		r2 = rules.scrub.output.r2,
	output:
		status = OUTDIR / "status" / "shovill.{sample}.txt",
		outdir = directory(OUTDIR / "{sample}" / "denovo" / "shovill"),
		contig = OUTDIR / "{sample}" / "denovo" / "shovill" / "contigs.fa",
	params:
		gsize = config['shovill']['gsize'],
		depth = config['shovill']['downsample'],
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
	--cpus {threads} 2>> {log} 1>> {log}

	touch {output.status}
	"""

rule scaffold:
	input: 
		contigs = rules.shovill.output.contig
	output:
		outdir = directory(OUTDIR / "{sample}" / "denovo" / "scaffold"),
		scaffold = OUTDIR / "{sample}" / "denovo" / "scaffold" / "ragtag.scaffold.fasta",
		status = OUTDIR / "status" / "denovo.scaffold.{sample}.txt",
	params:
		ref = SCAFFOLDREF,
		min_len = config['ragtag']['min_len']
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
		filled = OUTDIR / "{sample}" / "denovo" / "gap2seq" / "filled.fasta",
		status = OUTDIR / "status" / "denovo.gap2seq.{sample}.txt",
	shadow: "minimal"
	threads: config['threads']['gap2seq']
	conda: "../envs/scaffold.yaml"
	log: OUTDIR / "{sample}" / "log" / "denovo.gap2seq.{sample}.log"
	benchmark: OUTDIR / "{sample}" / "benchmark" / "denovo.gap2seq.{sample}.txt"
	shell:"""
	set +e
	Gap2Seq \
	--scaffolds {input.scaffold} \
	--filled {output.filled} \
	--reads {input.r1},{input.r2} \
	--threads {threads}	> {log} 2>&1
	exitcode=$?

	# gap2seq will fail if filling fails.
	# capture error, copy scaffolded as filled.
	if [ $exitcode -ne 0 ]
	then
		cp {input.scaffold} {output.filled}
	fi

	touch {output.status}
	"""

rule rename:
	input:
		filled = rules.gapfiller.output.filled
	output:
		renamed = OUTDIR / "{sample}" / "denovo" / "{sample}.final.fasta",
	params:
		sample_name = lambda w: w.sample
	conda: "../envs/misc.yaml"
	shell:"""
	seqkit replace -p "(.*)" -r '{params.sample_name}_contig{{nr}}' {input.filled} > {output.renamed} 
	"""

rule blast_ompa:
	input:
		contig = rules.rename.output.renamed,
	output:
		tab = OUTDIR / "{sample}" / "denovo" / "blast" / "blast.ompa.tab",
		status = OUTDIR / "status" / "denovo.blastn.{sample}.txt"
	log: OUTDIR / "{sample}" / "log" / "denovo.blastompa.{sample}.log"
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
	-out {output.tab} > {log} 2>&1

	touch {output.status}
	"""

rule mlst:
	input:
		contig = rules.rename.output.renamed,
	output:
		generic = OUTDIR / "{sample}" / "denovo" / "mlst" / "{sample}.genome.chlamydiales.mlst.txt",
		ct = OUTDIR / "{sample}" / "denovo" / "mlst" / "{sample}.genome.ctrachomatis.mlst.txt",
		plasmid =  OUTDIR / "{sample}" / "denovo" / "mlst" / "{sample}.genome.plasmid.mlst.txt",
		status = OUTDIR / "status" / "denovo.mlst.{sample}.txt",
	log: OUTDIR / "{sample}" / "log" / "denovo.mlst.{sample}.log"
	benchmark: OUTDIR / "{sample}" / "benchmark" / "denovo.mlst.{sample}.txt"
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
		status_blast = expand(OUTDIR / "status" / "denovo.blastn.{sample}.txt", sample = SAMPLES),
	output:
		tsv = OUTDIR / "denovo.blast.tsv",
		status = OUTDIR / "status" / "denovo.collate.blast.txt",
	params:
		outdir = OUTDIR,
		pattern = "**/denovo/blast/*.tab",
	threads: 1
	shell:"""
	echo -e "query\tsubject\tpident\tlength\tmismatch\tgapopen\tquery_start\tquery_end\tsubject_start\tsubject_end\tevalue\tbitscore" > {output.tsv}
	grep "" {params.outdir}/{params.pattern} >> {output.tsv}

	touch {output.status}
	"""

rule denovo_collate_mlst:
	input:
		generic = expand(OUTDIR / "{sample}" / "denovo" / "mlst" / "{sample}.genome.chlamydiales.mlst.txt", sample = SAMPLES),
		ct = expand(OUTDIR / "{sample}" / "denovo" / "mlst" / "{sample}.genome.ctrachomatis.mlst.txt", sample = SAMPLES),
		plasmid =  expand(OUTDIR / "{sample}" / "denovo" / "mlst" / "{sample}.genome.plasmid.mlst.txt", sample = SAMPLES),
	output:
		generic = OUTDIR / "denovo.mlst.generic.results.tsv",
		ct = OUTDIR / "denovo.mlst.ct.results.tsv",
		plasmid = OUTDIR / "denovo.mlst.plasmid.results.tsv",
		status = OUTDIR / "status" / "denovo.mlst.collate.txt"
	conda: "../envs/misc.yaml"
	threads: 1
	shell:"""
	csvtk concat {input.generic} -o {output.generic}
	csvtk concat {input.ct} -o {output.ct}
	csvtk concat {input.plasmid} -o {output.plasmid}

	touch {output.status}
	"""

rule ska_input_prep:
	input:
		filled = expand(OUTDIR / "{sample}" / "denovo" / "gap2seq" / "filled.fasta", sample = SAMPLES)
	output:
		glist = OUTDIR / "denovo" / "tree" / "input.list",
	params:
		samples = SAMPLES
	threads: 10
	run:
		print(input.filled)
		print(params.samples)
		for f,s in zip(input.filled, params.samples):
			with open(output.glist, "w") as handle:
				handle.write(f"{s}\t{f}\n")

rule ska_alignment:
	input:
		glist = rules.ska_input_prep.output.glist
	output:
		alignment = OUTDIR / "denovo" / "tree" / "alignment.fasta",
		tree = OUTDIR / "denovo" / "tree" / f"{ORG}.final_tree.tre"
	params:
		reference = SCAFFOLDREF,
		organism = ORG,
	conda: "../envs/tree.yaml"
	log: OUTDIR / "denovo" / "tree" / "tree.log"
	benchmark: OUTDIR / "denovo" / "tree" / "benchmark.tree.txt"
	threads: config["threads"]["tree"]
	shell:"""
	scripts/generate_ska_alignment.py \
	--threads {threads} \
	--reference {params.reference} \
	--input {input.glist} \
	--out {output.alignment}
	"""