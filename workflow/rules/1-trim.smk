rule fastp:
	input:
		r1 = INDIR / "{sample}_R1.fastq.gz",
		r2 = INDIR / "{sample}_R2.fastq.gz",
	output:
		r1 = OUTDIR / "{sample}" / "trim" / "{sample}_trim_R1.fastq.gz",
		r2 = OUTDIR / "{sample}" / "trim" / "{sample}_trim_R2.fastq.gz",
		html = OUTDIR / "{sample}" / "trim" /"{sample}.html",
		json = OUTDIR / "{sample}" / "trim" /"{sample}.json",
		status = OUTDIR / "status" / "trim.{sample}.txt",
	params:
		adapters = ADAPTERS,
		barcodes = BARCODES,
		window_size = 1,
		mean_quality = 3,
		right_window_size = 4,
		right_mean_quality = 15,
	threads: config["threads"]["fastp"]
	conda: "../envs/trim.yaml"
	log: OUTDIR / "{sample}" / "log" / "trim.{sample}.log"
	shell:"""
	fastp \
	-i {input.r1} \
	-I {input.r2} \
	-o {output.r1} \
	-O {output.r2} \
	-h {output.html} \
	-j {output.json} \
	--detect_adapter_for_pe \
	--cut_front \
	--cut_front_window_size={params.window_size} \
	--cut_front_mean_quality={params.mean_quality} \
	--cut_tail \
	--cut_tail_window_size={params.window_size} \
	--cut_tail_mean_quality={params.mean_quality} \
	--cut_right \
	--cut_right_window_size={params.right_window_size} \
	--cut_right_mean_quality={params.right_mean_quality} \
	--thread {threads} 2> {log}

	touch {output.status}
	"""

