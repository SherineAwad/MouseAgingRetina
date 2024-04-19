with open(config['SAMPLES']) as fp:
    SAMPLES= fp.read().splitlines()
print(SAMPLES)



rule all: 
    input:
        expand("{sample}_gex.pileup", sample = SAMPLES),
        expand("{sample}_atac.pileup", sample = SAMPLES), 
        "12wk1_gex.pileup", 
        "12wk1_atac.pileup",
        expand("{sample}_gex.snp.vcf", sample = SAMPLES),
        expand("{sample}_gex.indel.vcf", sample = SAMPLES), 
        expand("{sample}_atac.snp.vcf", sample = SAMPLES),
        expand("{sample}_atac.indel.vcf", sample = SAMPLES),
        expand("{sample}_gex.filtered.vcf",sample =SAMPLES),
        expand("{sample}_atac.filtered.vcf", sample =SAMPLES),
        expand("{sample}_gex.filtered.Germline.hc.vcf", sample =SAMPLES),
        expand("{sample}_gex.filtered.LOH.hc.vcf", sample =SAMPLES),
        expand("{sample}_gex.filtered.Somatic.hc.vcf", sample =SAMPLES),
        expand("{sample}_atac.filtered.Germline.hc.vcf", sample =SAMPLES),
        expand("{sample}_atac.filtered.LOH.hc.vcf", sample =SAMPLES),
        expand("{sample}_atac.filtered.Somatic.hc.vcf", sample =SAMPLES),
        expand("{sample}_gex_hc_germline.mm10_multianno.vcf", sample =SAMPLES),
        expand("{sample}_gex_hc_LOH.mm10_multianno.vcf", sample =SAMPLES),
        expand("{sample}_gex_hc_somatic.mm10_multianno.vcf",sample =SAMPLES),
        expand("{sample}_atac_hc_germline.mm10_multianno.vcf",sample =SAMPLES),
        expand("{sample}_atac_hc_LOH.mm10_multianno.vcf",sample =SAMPLES),
        expand("{sample}_atac_hc_somatic.mm10_multianno.vcf",sample =SAMPLES)


rule pileup:
    input:
       "{sample}_gex_possorted.bam",
       "{sample}_atac_possorted.bam"
    output:
       "{sample}_gex.pileup",
       "{sample}_atac.pileup"
    params:
       genome  = config['GENOME'],
    shell:
        """
        samtools mpileup -q 20 -Q 25 -B -d 1000 -f {params.genome} {input[0]} > {output[0]}   
        samtools mpileup -q 20 -Q 25 -B -d 1000 -f {params.genome} {input[1]} > {output[1]} 
        """
rule varscan:
        input: 
             "{sample}_gex.pileup",
             "{sample}_atac.pileup",
             "12wk1_gex.pileup",
             "12wk1_atac.pileup", 
        output: 
            "{sample}_gex.snp.vcf",
            "{sample}_gex.indel.vcf", 
            "{sample}_atac.snp.vcf", 
            "{sample}_atac.indel.vcf" 
        params: 
            "{sample}_gex",
            "{sample}_atac",
           
        shell:
            """
            varscan somatic {input[2]} {input[0]} {params[0]} --strand-filter 1 --output-vcf 1  
            varscan somatic {input[3]} {input[1]} {params[1]} --strand-filter 1 --output-vcf 1         
            """            
rule filter: 
        input: 
            "{sample}_gex.snp.vcf",
            "{sample}_gex.indel.vcf", 
            "{sample}_atac.snp.vcf",
            "{sample}_atac.indel.vcf"
        output: 
            "{sample}_gex.filtered.vcf", 
            "{sample}_atac.filtered.vcf" 
	shell:  
            """
            varscan somaticFilter {input[0]} --indel-file {input[1]} --output-file {output[0]} --p-value 0.01 --min-strands2 2 --output-vcf 1
            varscan somaticFilter {input[2]} --indel-file {input[3]} --output-file {output[1]} --p-value 0.01 --min-strands2 2 --output-vcf 1
            """


rule process_somatic:
    input:
        "{sample}_gex.filtered.vcf",
        "{sample}_atac.filtered.vcf"
    output:
        "{sample}_gex.filtered.Germline.hc.vcf",
        "{sample}_gex.filtered.LOH.hc.vcf",
        "{sample}_gex.filtered.Somatic.hc.vcf",
        "{sample}_atac.filtered.Germline.hc.vcf",
        "{sample}_atac.filtered.LOH.hc.vcf",
        "{sample}_atac.filtered.Somatic.hc.vcf"
    shell:
       """
       varscan processSomatic {input[0]} --p-value 0.05
       varscan processSomatic {input[1]} --p-value 0.05
       """


rule annotate: 
      input: 
          "{sample}_gex.filtered.Germline.hc.vcf",
          "{sample}_gex.filtered.LOH.hc.vcf",
          "{sample}_gex.filtered.Somatic.hc.vcf",
          "{sample}_atac.filtered.Germline.hc.vcf",
          "{sample}_atac.filtered.LOH.hc.vcf",
          "{sample}_atac.filtered.Somatic.hc.vcf"
      params: 
          "{sample}_gex_hc_germline",
          "{sample}_gex_hc_LOH", 
          "{sample}_gex_hc_somatic", 
          "{sample}_atac_hc_germline",
          "{sample}_atac_hc_LOH",
          "{sample}_atac_hc_somatic",
      output: 
          "{sample}_gex_hc_germline.mm10_multianno.vcf",
          "{sample}_gex_hc_LOH.mm10_multianno.vcf", 
          "{sample}_gex_hc_somatic.mm10_multianno.vcf", 
          "{sample}_atac_hc_germline.mm10_multianno.vcf",
          "{sample}_atac_hc_LOH.mm10_multianno.vcf",
          "{sample}_atac_hc_somatic.mm10_multianno.vcf",
      shell:
          """
           {config[ANNOVAR]}/table_annovar.pl {input[0]}  {config[ANNOVARDB]} -buildver mm10 -out {params[0]} -remove -protocol refGene -operation g   -nastring . -vcfinput
           {config[ANNOVAR]}/table_annovar.pl {input[1]}  {config[ANNOVARDB]} -buildver mm10 -out {params[1]} -remove -protocol refGene -operation g   -nastring . -vcfinput
           {config[ANNOVAR]}/table_annovar.pl {input[2]}  {config[ANNOVARDB]} -buildver mm10 -out {params[2]} -remove -protocol refGene -operation g   -nastring . -vcfinput
           
           {config[ANNOVAR]}/table_annovar.pl {input[3]}  {config[ANNOVARDB]} -buildver mm10 -out {params[3]} -remove -protocol refGene -operation g   -nastring . -vcfinput
           {config[ANNOVAR]}/table_annovar.pl {input[4]}  {config[ANNOVARDB]} -buildver mm10 -out {params[4]} -remove -protocol refGene -operation g   -nastring . -vcfinput
           {config[ANNOVAR]}/table_annovar.pl {input[5]}  {config[ANNOVARDB]} -buildver mm10 -out {params[5]} -remove -protocol refGene -operation g   -nastring . -vcfinput
          """

