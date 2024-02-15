rule mlst:
    input:
        rules.scaffold.output.scaffold
    output:
        generic = OUTDIR / "{sample}" / "mlst" / "{sample}.genome.chlamydiales.mlst.txt",
        ct = OUTDIR / "{sample}" / "mlst" / "{sample}.genome.ctrachomatis.mlst.txt",
        plasmid =  OUTDIR / "{sample}" / "mlst" / "{sample}.genome.plasmid.mlst.txt",
        status = OUTDIR / "status" / "mlst.{sample}.txt",
    log: OUTDIR / "{sample}" / "log" / "mlst.{sample}.log"
    benchmark: OUTDIR / "{sample}" / "benchmark" / "mlst.{sample}.txt"
    conda: "../envs/mlst.yaml"
    params:
        dbgeneric = MLSTDBLOC / "chlamydiales",
        dbct = MLSTDBLOC / "c.trachomatis",
        dbplasmid = MLSTDBLOC / "plasmid",
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