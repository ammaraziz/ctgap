rule mlst_genome:
    input:
        rules.shovill.output.contig
    output:
        generic = OUTDIR / "{sample}" / "mlst" / "{sample}.genome.chlamydiales.mlst.txt",
        ct = OUTDIR / "{sample}" / "mlst" / "{sample}.genome.ctrachomatis.mlst.txt",
        plasmid =  OUTDIR / "{sample}" / "mlst" / "{sample}.genome.plasmid.mlst.txt",
    log: OUTDIR / "{sample}" / "log" / "mlst.{sample}.log"
    benchmark: OUTDIR / "{sample}" / "benchmark" / "mlst.{sample}.txt"
    params:
        dbgeneric = MLSTDBLOC / "chlamydiales",
        dbct = MLSTDBLOC / "c.trachomatis",
        plasmid = MLSTDBLOC / "plasmid",
    threads: config["threads"]["mlst"]
    shell:"""
    echo -e "chlamydiales\n"
    claMLST search \
    {params.dbgeneric} \ 
    {input} > {output.generic} 2> {log}

    echo -e "\nctrachomatis\n"
    claMLST search \
    {params.ct} \ 
    {input} > {output.ct} 2>> {log}

    echo -e "\nplasmid\n"
    claMLST search \
    {params.plasmid} \ 
    {input} > {output.plasmid} 2>> {log}
    """
