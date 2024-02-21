rule index:
	message: "Indexing references"
	input:
		reference = REF,
		refset24 = REFSET,
	output:
		status = OUTDIR /"status" / "reference.status"
	params:
		prefix_ref = REF.stem,
		prefix_ref24 =  REFSET.stem,
		loc_ref = REF.parent,
		loc_ref24 = REFSET.parent,
	threads: config["threads"]["bowtieindex"]
	conda: "../envs/bowtie.yaml"
	shell:"""
	bowtie2-build --threads {threads} {input.reference} {params.loc_ref}/{params.prefix_ref} &> /dev/null
	bowtie2-build --threads {threads} {input.refset24} {params.loc_ref24}/{params.prefix_ref24} &> /dev/null

	touch {output.status}
	"""

rule bowtie:
	message: "Aligning to specified reference"
	input:
		r1 = rules.scrub.output.r1,
		r2 = rules.scrub.output.r2,
		status = rules.index.output.status,
	output:
		bam = OUTDIR / "{sample}" / "ref-denovo" / "{sample}.bam",
		status = OUTDIR / "status" / "bowtie2.{sample}.txt",
	params:
		prefix = rules.index.params.prefix_ref,
		loc = rules.index.params.loc_ref,
	threads: config["threads"]["bowtie"]
	log: OUTDIR / "{sample}" / "log" / "bowtie2.{sample}.log"
	conda: "../envs/bowtie.yaml"
	shell:"""
	bowtie2 -x \
	{params.loc}/{params.prefix} \
	-1 {input.r1} \
	-2 {input.r2} \
	--threads {threads} | \
	samtools view -bSh - | \
	samtools sort -@{threads} \
	-o {output.bam} 2> {log}

	samtools index {output.bam}

	touch {output.status}
	"""

rule bowtie_ref24:
	message: "Align against reference set"
	input:
		r1 = rules.scrub.output.r1,
		r2 = rules.scrub.output.r2,
		status = rules.index.output.status,
	output:
		bam = OUTDIR / "{sample}" / "ref-denovo" / "{sample}.ref24.bam",
		coverage = OUTDIR / "{sample}" / "ref-denovo" / "coverage.ref24.{sample}.tsv",
		status = OUTDIR / "status" / "bowtie2.ref24.{sample}.txt",
	params:
		prefix = rules.index.params.prefix_ref24,
		loc = rules.index.params.loc_ref24,
		sample_name = lambda w: w.sample,
		tmp = lambda w: OUTDIR / f"{w.sample}" / "ref-denovo" / "tmp.cov.tsv",
	threads: config["threads"]["bowtie"]
	log: OUTDIR / "{sample}" / "log" / "bowtie2.ref24.{sample}.log"
	conda: "../envs/bowtie.yaml"
	shell:"""
	bowtie2 -x \
	{params.loc}/{params.prefix} \
	-1 {input.r1} \
	-2 {input.r2} \
	--threads {threads} | \
	samtools view -bSh - | \
	samtools sort -@{threads} \
	-o {output.bam} 2> {log}

	samtools index {output.bam}
	samtools coverage {output.bam} > {params.tmp}

	csvtk mutate2 -t -C $ {params.tmp} -n "sample_name" --at 1 -e " '{params.sample_name}' " > {output.coverage}

	touch {output.status}
	"""

rule bamtofastq:
	input:
		bam = rules.bowtie.output.bam,
	output:
		r1 = OUTDIR / "{sample}" / "ref-denovo" / "{sample}_refg_r1.fastq.gz",
		r2 = OUTDIR / "{sample}" / "ref-denovo" / "{sample}_refg_r2.fastq.gz",
		status = OUTDIR / "status" / "bamtofastq.{sample}.txt",
	threads: config["threads"]["bedtools"]
	log: OUTDIR / "{sample}" / "log" / "bowtie2.{sample}.log"
	conda: "../envs/misc.yaml"
	shell:"""
	bedtools bamtofastq \
	-i {input.bam} \
	-fq {output.r1} \
	-fq2 {output.r2} 2> {log}

	touch {output.status}
	"""

rule ref_shovill:
	message: "Assembling Sample {wildcards.sample}"
	input:
		r1 = rules.bamtofastq.output.r1,
		r2 = rules.bamtofastq.output.r2,
	output:
		status = OUTDIR / "status" / "ref-denovo.shovill.{sample}.txt",
		outdir = directory(OUTDIR / "{sample}" / "ref-denovo" / "shovill"),
		contig = OUTDIR / "{sample}" / "ref-denovo" / "shovill" / "contigs.fa",
	params:
		gsize = "1.04M",
		depth = 200,
	conda: "../envs/shovill.yaml"
	log: OUTDIR / "{sample}" / "log" / "ref-denovo.shovill.{sample}.log"
	benchmark: OUTDIR / "{sample}" / "benchmark" / "ref-denovo.shovill.{sample}.txt"
	threads: config["threads"]["shovill"]
	shell:"""
	shovill \
	--R1 {input.r1} \
	--R2 {input.r2} \
	--outdir {output.outdir} \
	--gsize {params.gsize} \
	--depth {params.depth} \
	--force \
	--cpus {threads} > {log} 2>&1

	touch {output.status}
	"""

rule ref_scaffold:
	input: 
		contigs = rules.ref_shovill.output.contig
	output:
		outdir = directory(OUTDIR / "{sample}" / "ref-denovo" / "scaffold"),
		scaffold = OUTDIR / "{sample}" / "ref-denovo" / "scaffold" / "ragtag.scaffold.fasta",
		status = OUTDIR / "status" / "ref-denovo.scaffold.{sample}.txt",
	params:
		ref = SCAFFOLDREF,
		min_len = 500
	threads: config['threads']['ragtag']
	conda: "../envs/scaffold.yaml"
	log: OUTDIR / "{sample}" / "log" / "ref-denovo_scaffold.{sample}.log"
	benchmark: OUTDIR / "{sample}" / "benchmark" / "ref-denovo.scaffold.{sample}.txt"
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

rule ref_gapfiller:
	input:
		scaffold = rules.ref_scaffold.output.scaffold,
		r1 = rules.scrub.output.r1,
		r2 = rules.scrub.output.r2,
	output:
		filled = OUTDIR / "{sample}" / "ref-denovo" / "gap2seq" / "{sample}.fasta",
		status = OUTDIR / "status" / "ref-denovo.gap2seq.{sample}.txt",
	params:
		sample_name = lambda w: w.sample,
	shadow: "copy-minimal"
	threads: config['threads']['gap2seq']
	conda: "../envs/scaffold.yaml"
	log: OUTDIR / "{sample}" / "log" / "ref-denovo.gap2seq.{sample}.log"
	benchmark: OUTDIR / "{sample}" / "benchmark" / "ref-denovo.gap2seq.{sample}.txt"
	shell:"""
	EXITCODE=$(Gap2Seq \
	--scaffolds {input.scaffold} \
	--filled {output.filled} \
	--reads {input.r1},{input.r2} \
	--threads {threads} > {log} 2>&1)

	# gap2seq will fail if filling fails.
	# capture error, copy scaffolded as filled.
	if [ "$EXITCODE" = "1" ]
	then
		cp {input.scaffold} {output.filled}
	fi

	touch {output.status}
	"""

rule ref_rename:
	input:
		filled = rules.ref_gapfiller.output.filled
	output:
		renamed = OUTDIR / "{sample}" / "ref-denovo" / "{sample}.final.fasta",
	conda: "../envs/misc.yaml"
	shell:"""
	seqkit replace -p "(.*)" -r "F3_contig{nr}" {input.filled} > {output.renamed} 
	"""


rule ref_blast_ompa:
	input:
		contig = rules.ref_rename.output.rename,
	output:
		tab = OUTDIR / "{sample}" / "ref-denovo" / "blast" / "blast.ompa.tab",
		status = OUTDIR / "status" / "ref-denovo.blastn.{sample}.txt"
	log: OUTDIR / "{sample}" / "log" / "ref-denovo.blastompa.{sample}.log"
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

rule ref_mlst:
	input:
		rules.ref_rename.output.rename,
	output:
		generic = OUTDIR / "{sample}" / "ref-denovo" / "mlst" / "{sample}.genome.chlamydiales.mlst.txt",
		ct = OUTDIR / "{sample}" / "ref-denovo" / "mlst" / "{sample}.genome.ctrachomatis.mlst.txt",
		plasmid =  OUTDIR / "{sample}" / "ref-denovo" / "mlst" / "{sample}.genome.plasmid.mlst.txt",
		status = OUTDIR / "{sample}" / "status" / "ref-denovo.mlst.{sample}.txt",
	log: OUTDIR / "{sample}" / "log" / "ref-denovo.mlst.{sample}.log"
	benchmark: OUTDIR / "{sample}" / "benchmark" / "ref-denovo.mlst.{sample}.txt"
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

rule ref_collate_coverage:
	input:
		coverages = expand(OUTDIR / "{sample}" / "ref-denovo" / "coverage.ref24.{sample}.tsv", sample = SAMPLES),
	output:
		coverages = OUTDIR / "ref-denovo.coverage.tsv",
		status = OUTDIR / "status" / "ref-denovo.collage.coverage.txt",
	threads: 1
	conda: "../envs/misc.yaml"
	shell:"""
	csvtk concat -C $ {input.coverages} >  {output.coverages}

	touch {output.status}
	"""

rule ref_collate_blast:
	input:
		blast_status = expand(OUTDIR / "status" / "ref-denovo.blastn.{sample}.txt", sample = SAMPLES),
	output:
		tsv = OUTDIR / "ref-denovo.blast.tsv",
		status = OUTDIR / "status" / "ref-denovo.collate.blast.txt",
	params:
		outdir = OUTDIR,
		pattern = "**/ref-denovo/blast/*.tab",
	threads: 1
	shell:"""
	echo -e "query\tsubject\tpident\tlength\tmismatch\tgapopen\tquery_start\tquery_end\tsubject_start\tsubject_end\tevalue\tbitscore" > {output.tsv}
	grep "" {params.outdir}/{params.pattern} >> {output.tsv}

	touch {output.status}
	"""

rule ref_collate_mlst:
	input:
		generic = expand(OUTDIR / "{sample}" / "ref-denovo" / "mlst" / "{sample}.genome.chlamydiales.mlst.txt", sample = SAMPLES),
		ct = expand(OUTDIR / "{sample}" / "ref-denovo" / "mlst" / "{sample}.genome.ctrachomatis.mlst.txt", sample = SAMPLES),
		plasmid =  expand(OUTDIR / "{sample}" / "ref-denovo" / "mlst" / "{sample}.genome.plasmid.mlst.txt", sample = SAMPLES),
	output:
		generic = OUTDIR / "ref-denovo.mlst.generic.results.tsv",
		ct = OUTDIR / "ref-denovo.mlst.ct.results.tsv",
		plasmid = OUTDIR / "ref-denovo.mlst.plasmid.results.tsv",
		status = OUTDIR / "status" / "mlst.collate.txt"
	conda: "../envs/misc.yaml"
	threads: 1
	shell:"""
	csvtk concat {input.generic} -o {output.generic}
	csvtk concat {input.ct} -o {output.ct}
	csvtk concat {input.plasmid} -o {output.plasmid}

	touch {output.status}
	"""
