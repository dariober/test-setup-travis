from collections import OrderedDict
import re
import os
import pandas 
import sys

# Dir where slurm logs will go. It must match the sbatch option `--output=<dir>/...`
os.makedirs('slurm', exist_ok=True)

REF= config['ref']
CHROMS= config['chroms']

manifest= pandas.read_csv(config['manifest'], sep= '\t')

# TODO for devel: Add some sanity check to manifest dataframe.

manifest= manifest.assign(fastq_dir= [os.path.abspath(os.path.dirname(x)) for x in manifest['fastq']])
manifest= manifest.assign(fastq_base= [os.path.basename(x).replace('.fastq.gz', '') for x in manifest['fastq']])

pon_library= list(OrderedDict.fromkeys(manifest[(manifest['sample_type'] == 'PON')]['library']))

tumour_lib= list(OrderedDict.fromkeys(manifest[(manifest['sample_type'] == 'Tumour')]['library']))

# Do not interpret wildcards as regexes. This is to prevent {library}.bam to
# match also {library}.tmp.bam See also
# https://groups.google.com/forum/#!topic/snakemake/wVlJW9X-9EU
wildcard_constraints:
    fastq_base= '|'.join([re.escape(x) for x in manifest.fastq_base]),
    library= '|'.join([re.escape(x) for x in manifest.library]),
    tumour_lib= '|'.join([re.escape(x) for x in tumour_lib]),
    pon_library= '|'.join([re.escape(x) for x in pon_library]),
    chrom= '|'.join([re.escape(x) for x in CHROMS])

os.chdir(config['output_dir'])

# Dir where slurm logs will go. It must match the sbatch option `--output=<dir>/...`
os.makedirs('slurm', exist_ok=True)

def normalForTumourLib(tumour_lib):
    """Find the matched normal library for the tumour library in input
    """
    patient= manifest[(manifest['library'] == tumour_lib)]['patient']
    patient= list(set(patient))
    if len(patient) != 1:
        raise Exception('Expected exactly 1 patient assigned to library %s. Got: %s' %(tumour_lib, ', '.join(patient)))
    lib= manifest[(manifest['patient'] == patient[0]) & (manifest['sample_type'] == 'BloodNormal')]['library']
    lib= list(set(lib))
    if len(lib) != 1:
        raise Exception('Expected exactly 1 matched normal. Got: %s' % ', '.join(lib))
    return lib[0]

rule all:
    input:
        expand('pindel/{tumour_lib}.vep.vcf.gz.tbi', tumour_lib= tumour_lib),
        expand('facets/{tumour_lib}.cn.bed', tumour_lib= tumour_lib),
        expand('manta/{tumour_lib}.sv.vcf.gz.tbi', tumour_lib= tumour_lib),
        expand('gatk4/{tumour_lib}.vep.vcf.gz.tbi', tumour_lib= tumour_lib),
        expand('samstats/{library}.stats', library= manifest.library),
        expand('tdf/{library}.tdf', library= manifest.library),
        expand('bam/{library}.bam', library= manifest.library),
        expand('bam/{library}.bam.bai', library= manifest.library),
        'multiqc/multiqc_report.html',
        expand('fastqc/{fastq_base}_fastqc.zip', fastq_base= manifest.fastq_base),

rule filter_manta:
    input:
        vcf= 'manta/{tumour_lib}/results/variants/somaticSV.vcf.gz',
        tbi= 'manta/{tumour_lib}/results/variants/somaticSV.vcf.gz.tbi',
        exclude= 'manta/{tumour_lib}.ids_to_exclude.txt',
    output:
        vcf= 'manta/{tumour_lib}.sv.vcf.gz',    
        tbi= 'manta/{tumour_lib}.sv.vcf.gz.tbi',
    params:
        mem= 4000,
        jobname= '{tumour_lib}',
        cpus_per_task= 1,
        bin_dir= config['bin_dir'],
    shell:
        """
        {params.bin_dir}/filter_manta.py --tumour_sample {wildcards.tumour_lib} --exclude_ids {input.exclude} --min_support 7 -i {input.vcf} \\
        | bcftools view -f PASS --output-type z > {output.vcf}
        tabix -f {output.vcf}
        """

rule manta_ids_in_pon:
    # Collect the manta IDs (and MATEIDs for BNDs) whose start or end region
    # overlaps with the PON.
    input:
        vcf= 'manta/{tumour_lib}/results/variants/somaticSV.vcf.gz',
        tbi= 'manta/{tumour_lib}/results/variants/somaticSV.vcf.gz.tbi',
        pon= 'manta/pon/panelOfNormals.bed.gz',
        poi= 'manta/pon/panelOfNormals.bed.gz.tbi',
    output:
        temp('manta/{tumour_lib}.ids_to_exclude.txt'),
    params:
        mem= 4000,
        jobname= '{tumour_lib}',
        cpus_per_task= 1,
    shell:
        """
        bcftools view {input.vcf} \\
        | vcfToBedpe \\
        | grep -v '^#' \\
        | awk -v FS='\\t' -v OFS='\\t' '{{
        if($13 ~ "^MATEID=" || $13 ~ ";MATEID=") {{
            # INFO field (13th in bedpe) contains tag 'MATEID': get its value.
            sub(".*MATEID=", "", $13);
            sub(";.*", "", $13);
            mate=$13
        }} else {{
            mate=$7
        }};
        print $1, $2, $3, $7
        print $4, $5, $6, mate}}' \\
        | awk -v OFS='\\t' '{{
            # Because of issue https://github.com/Illumina/manta/issues/127
            if($2 < 0) $2=0        
            print $0
        }}' \\
        | intersectBed -a - -b {input.pon} \\
        | cut -f 4 \\
        | sort \\
        | uniq > {output}
        """

mem= 16000
memGb= int(mem/1000)
rule manta:
    input:
        tumour= 'bam/{tumour_lib}.bam',
        tbai= 'bam/{tumour_lib}.bam.bai',
        normal= lambda wildcards: 'bam/%s.bam' % normalForTumourLib(wildcards.tumour_lib),
        nbai= lambda wildcards: 'bam/%s.bam.bai' % normalForTumourLib(wildcards.tumour_lib),
    output:
        'manta/{tumour_lib}/results/variants/somaticSV.vcf.gz',
        'manta/{tumour_lib}/results/variants/somaticSV.vcf.gz.tbi',
    params:
        ref= REF,
        mem= mem,
        baseDir= 'manta/{tumour_lib}',
        memGb= memGb,
        jobname= '{tumour_lib}',
        cpus_per_task= 8,
    shell:
        """
        rm -rf {params.baseDir}
        
        ~/applications/manta/manta-1.3.2.centos6_x86_64/bin/configManta.py \\
            --tumorBam {input.tumour} \\
            --normalBam {input.normal}  \\
            --referenceFasta {params.ref} \\
            --runDir {params.baseDir}

        {params.baseDir}/runWorkflow.py --jobs {params.cpus_per_task} --mode local --memGb {params.memGb}
        rm -r {params.baseDir}/workspace 
        """

rule manta_excludable:
    input:
        vcf= 'manta/pon/panelOfNormals.vcf.gz',
        tbi= 'manta/pon/panelOfNormals.vcf.gz.tbi',
        ref= REF + '.fai',
    output:
        bed= 'manta/pon/panelOfNormals.bed.gz',
        tbi= 'manta/pon/panelOfNormals.bed.gz.tbi',
    params:
        mem= 1000,
        jobname= 'panelOfNormals',
        cpus_per_task= 1,
    shell:
        """
        bcftools view {input.vcf} | vcfToBedpe | grep -v '^#' | cut -f 1-3 > {output.bed}.tmp
        bcftools view {input.vcf} | vcfToBedpe | grep -v '^#' | cut -f 4-6 >> {output.bed}.tmp
        sort -k1,1 -k2,2n {output.bed}.tmp \\
        | awk -v OFS='\\t' '{{
            # Because of issue https://github.com/Illumina/manta/issues/127
            if($2 < 0) $2=0        
            print $0
        }}' \\
        | slopBed -b 50 -i - -g {input.ref} \\
        | mergeBed \\
        | bgzip > {output.bed}
        tabix -f {output.bed}
        rm {output.bed}.tmp
        """

rule merge_manta_pon:
    input:
        vcf= lambda wildcards: expand('manta/pon/{pon_library}.diploidSV.vcf.gz', pon_library= pon_library),
        tbi= lambda wildcards: expand('manta/pon/{pon_library}.diploidSV.vcf.gz.tbi', pon_library= pon_library),
    output:
        vcf= 'manta/pon/panelOfNormals.vcf.gz',
        tbi= 'manta/pon/panelOfNormals.vcf.gz.tbi',
    params:
        mem= 1000,
        jobname= 'panelOfNormals',
        cpus_per_task= 1,
    shell:
        """
        bcftools concat --allow-overlaps {input.vcf} \\
        | bcftools view --output-type z --include 'FILTER == "PASS"' > {output.vcf}
        tabix -f {output.vcf}
        """

rule manta_pon_rename:
    input:
        vcf= 'manta/pon/{pon_library}/results/variants/diploidSV.vcf.gz',
        tbi= 'manta/pon/{pon_library}/results/variants/diploidSV.vcf.gz.tbi',
    output:
        vcf= temp('manta/pon/{pon_library}.diploidSV.vcf.gz'),
        tbi= temp('manta/pon/{pon_library}.diploidSV.vcf.gz.tbi'),
        sname= temp('manta/pon/{pon_library}.tmp'),
    params:
        mem= 1000,
        jobname= '{pon_library}',
        cpus_per_task= 1,
    shell:
        """
        echo "normal" > {output.sname}
        bcftools reheader --samples {output.sname} {input.vcf} > {output.vcf}
        tabix -f {output.vcf}
        """

mem= 16000
memGb= int(mem/1000)
rule manta_pon:
    input:
        bam= 'bam/{pon_library}.bam',
        bai= 'bam/{pon_library}.bam.bai',
    output:
        vcf= 'manta/pon/{pon_library}/results/variants/diploidSV.vcf.gz',
        tbi= 'manta/pon/{pon_library}/results/variants/diploidSV.vcf.gz.tbi'
    params:
        ref= REF,
        mem= mem,
        baseDir= 'manta/pon/{pon_library}',
        memGb= memGb,
        jobname= '{pon_library}',
        cpus_per_task= 8,
    shell:
        """
        rm -rf {params.baseDir}
        
        ~/applications/manta/manta-1.3.2.centos6_x86_64/bin/configManta.py \\
            --normalBam {input.bam}  \\
            --referenceFasta {params.ref} \\
            --runDir {params.baseDir}

        {params.baseDir}/runWorkflow.py --jobs {params.cpus_per_task} --mode local --memGb {params.memGb}
        rm -r {params.baseDir}/workspace
        """


rule vepPindel:
    input:
        vcf= 'pindel/{tumour_lib}.vcf.gz',
        tbi= 'pindel/{tumour_lib}.vcf.gz.tbi',
    output:
        vcf= 'pindel/{tumour_lib}.vep.vcf.gz',
        tbi= 'pindel/{tumour_lib}.vep.vcf.gz.tbi',
        stats= 'pindel/{tumour_lib}.vep.summary.html'
    params:
        jobname= '{tumour_lib}',
        cpus_per_task= 4,
        mem= 4000,
        bin_dir= config['bin_dir'],
        dir_cache= config['vep_dir_cache'],
    shell:
        """
        ~/applications/vep/ensembl-vep-release-91.2/vep --flag_pick --offline --cache --dir_cache {params.dir_cache} --fork 4 --vcf --force_overwrite -i {input.vcf} -o STDOUT --stats_file {output.stats} \\
        | {params.bin_dir}/vepReorder.py \\
        | bgzip > {output.vcf}
        sleep 120
        tabix {output.vcf}
        """

rule pindelConcat:
    input:
        vcf= lambda wildcards: expand('pindel/{tumour_lib}.{chrom}.intx.vcf.gz', chrom= CHROMS, tumour_lib= wildcards.tumour_lib),
        tbi= lambda wildcards: expand('pindel/{tumour_lib}.{chrom}.intx.vcf.gz.tbi', chrom= CHROMS, tumour_lib= wildcards.tumour_lib),
    output:
        vcf= temp('pindel/{tumour_lib}.vcf.gz'),
        tbi= temp('pindel/{tumour_lib}.vcf.gz.tbi'),
    params:
        jobname= '{tumour_lib}',
        cpus_per_task= 1,
        mem= 1000,
    shell:
        """
        # Exclude VCFs with no records as they don't have sample names.
        vcf=''
        for x in {input.vcf}
        do
            n=`bcftools view -H $x | wc -l`
            if [[ $n > 0 ]]
            then
                vcf="${{vcf}} $x"
            fi
        done

        bcftools concat --output-type z ${{vcf}} > {output.vcf}
        tabix -f {output.vcf} 
        """

rule pindelChromFilter:
    input:
        vcf= 'pindel/{tumour_lib}.{chrom}.vcf.gz',
        tbi= 'pindel/{tumour_lib}.{chrom}.vcf.gz.tbi',
        pon= 'pindel/pon/pon.bed.gz',
        poi= 'pindel/pon/pon.bed.gz.tbi',
        repeats= 'pindel/pon/simpleRepeat.bed.gz',
        rei= 'pindel/pon/simpleRepeat.bed.gz.tbi',
    output:
        vcf= temp('pindel/{tumour_lib}.{chrom}.intx.vcf.gz'),
        tbi= temp('pindel/{tumour_lib}.{chrom}.intx.vcf.gz.tbi'),
    params:
        jobname= '{tumour_lib}.{chrom}',
        cpus_per_task= 1,
        mem= 2000,
    shell:
        """
        ## Prepare header
        bcftools view -h {input.vcf} | grep -v '#CHROM' > {input.vcf}.tmp.vcf

        echo '##FILTER=<ID=panel_of_normals,Description="Found in panel of normals">' >> {input.vcf}.tmp.vcf
        echo '##FILTER=<ID=simple_repeat,Description="Overlaps a simple repeat">' >> {input.vcf}.tmp.vcf
        bcftools view -h {input.vcf} | grep '#CHROM' >> {input.vcf}.tmp.vcf

        intersectBed -a {input.vcf} -b {input.pon} -wa -loj \\
        | awk -v OFS='\\t' '{{
                            if($NF == "panelOfNormals") {{if($7 == "PASS") {{$7="panel_of_normals"}} else {{$7=$7";panel_of_normals"}}}}
                            else if($NF == "."){{}}
                            else {{exit 1}}
                            print $0
                            }}' \\
        | cut -f 1-11 >> {input.vcf}.tmp.vcf

        intersectBed -a {input.vcf}.tmp.vcf -b {input.repeats} -wa -loj -header \\
        | awk -v OFS='\t' '{{
                            if($0 ~ "^#") {{}}
                            else if($12 ~ "chr") {{if($7 == "PASS") {{$7="simple_repeat"}} else {{$7=$7";simple_repeat"}}}}
                            else if($12 == "."){{}}
                            else {{exit 1}}
                            print $0
                            }}' \\
        | cut -f 1-11 \\
        | bgzip > {output.vcf}
        tabix -f {output.vcf}

        rm {input.vcf}.tmp.vcf
        """

pindel_prefix= 'pindel/{tumour_lib}.{chrom}'
rule pindelChromCandidate:
    input:
        'pindel/{tumour_lib}.{chrom}.done',
    output:
        vcf= temp('pindel/{tumour_lib}.{chrom}.vcf.gz'),
        tbi= temp('pindel/{tumour_lib}.{chrom}.vcf.gz.tbi'),
    params:
        bin_dir= config['bin_dir'],
        prefix= pindel_prefix,
        ref= REF,
        jobname= '{tumour_lib}.{chrom}',
        cpus_per_task= 1,
        mem= 1000,
    shell:
        """
        {params.bin_dir}/filterSomaticPindel.py --input-prefix {params.prefix} \\
            --type D SI \\
            --output-prefix {params.prefix}_filter \\
            --normal normal \\
            --tumour tumour \\
            --tumour-support 4 \\
            --tumour-freq 0.03

        pindel2vcf --pindel_output_root {params.prefix}_filter \\
                    --vcf {params.prefix}_filter.vcf \\
                    --reference {params.ref} \\
                    --reference_name GRCh38 \\
                    --reference_date 19000101

        sed 's/##FORMAT=<ID=PL,Number=3,/##FORMAT=<ID=PL,Number=G,/' {params.prefix}_filter.vcf \\
        | bcftools norm -f {params.ref} --threads 8 - \\
        | bcftools view -i 'SVLEN < 3000 && SVLEN > -3000' \\
        | bgzip > {output.vcf}
        tabix -f {output.vcf}

        rm {params.prefix}_filter.vcf
        rm {params.prefix}_filter_*
        """

rule concatPindelPon:
    input:
        vcf= expand('pindel/pon/{pon_library}.vcf.gz', pon_library= pon_library),
        tbi= expand('pindel/pon/{pon_library}.vcf.gz.tbi', pon_library= pon_library),
    output:
        vcf= 'pindel/pon/pon.bed.gz',
        tbi= 'pindel/pon/pon.bed.gz.tbi',
    params:
        jobname= 'pindelPon',
        cpus_per_task= 1,
        mem= 4000,
    shell:
        """
        # Exclude VCFs with no records as they don't have sample names.
        vcf=''
        for x in {input.vcf}
        do
            n=`bcftools view -H $x | wc -l`
            if [[ $n > 0 ]]
            then
                vcf="${{vcf}} $x"
            fi
        done

        bcftools concat --allow-overlaps ${{vcf}} \\
        | mergeBed \\
        | awk -v OFS='\\t' '{{print $0, "panelOfNormals"}}' \\
        | bgzip > {output.vcf}

        tabix -f {output.vcf}
        """

rule concatPindelNormal:
    input:
        vcf= lambda wildcards: expand("pindel/pon/{pon_library}.{chrom}.vcf.gz", pon_library= wildcards.pon_library, chrom= CHROMS),
        tbi= lambda wildcards: expand("pindel/pon/{pon_library}.{chrom}.vcf.gz.tbi", pon_library= wildcards.pon_library, chrom= CHROMS)
    output:
        vcf= 'pindel/pon/{pon_library}.vcf.gz',
        tbi= 'pindel/pon/{pon_library}.vcf.gz.tbi',
    params:
        jobname= '{pon_library}',
        cpus_per_task= 1,
        mem= 1000,
    shell:
        """
        # Exclude VCFs with no records as they don't have sample names.
        vcf=''
        for x in {input.vcf}
        do
            n=`bcftools view -H $x | wc -l`
            if [[ $n > 0 ]]
            then
                vcf="${{vcf}} $x"
            fi
        done

        bcftools concat --output-type z ${{vcf}} > {output.vcf}
        tabix -f {output.vcf}
        """

pindel_prefix_pon= 'pindel/{pon_library}.{chrom}'
rule pindelChromNormalToVcf:
    input:
        'pindel/{pon_library}.{chrom}.done',
    output:
        vcf= temp('pindel/pon/{pon_library}.{chrom}.vcf.gz'),
        tbi= temp('pindel/pon/{pon_library}.{chrom}.vcf.gz.tbi'),
    params:
        bin_dir= config['bin_dir'],
        prefix= pindel_prefix_pon,
        ref= REF,
        jobname= '{pon_library}.{chrom}',
        cpus_per_task= 1,
        mem= 1000,
    shell:
        """
        {params.bin_dir}/filterSomaticPindel.py --input-prefix {params.prefix} \\
            --type D SI \\
            --output-prefix {params.prefix}_nfilter \\
            --tumour normal \\
            --tumour-support 2 \\
            --tumour-freq 0.02

        pindel2vcf --pindel_output_root {params.prefix}_nfilter \\
                    --vcf {params.prefix}_nfilter.vcf \\
                    --reference {params.ref} \\
                    --reference_name GRCh38 \\
                    --reference_date 19000101

        sed 's/##FORMAT=<ID=PL,Number=3,/##FORMAT=<ID=PL,Number=G,/' {params.prefix}_nfilter.vcf \\
        | bcftools norm -f {params.ref} --threads 8 - \\
        | bcftools view -i 'SVLEN < 3000 && SVLEN > -3000' \\
        | bgzip > {output.vcf}
        tabix -f {output.vcf}

        rm {params.prefix}_nfilter.vcf
        rm {params.prefix}_nfilter_*
        """

rule pindelChrom:
    input:
        tumour= 'bam/{tumour_lib}.bam',
        normal= lambda wildcards: 'bam/%s.bam' % normalForTumourLib(wildcards.tumour_lib),
        config= 'pindel/{tumour_lib}.config',
        done= 'pindel/{tumour_lib}.{chrom}.exclude.done',
        exclude= 'pindel/{tumour_lib}.{chrom}.exclude.bed',
    output:
        done= touch('pindel/{tumour_lib}.{chrom}.done')
    params:
        prefix= pindel_prefix,
        ref= REF,
        chrom= '{chrom}',
        jobname= '{tumour_lib}.{chrom}',
        cpus_per_task= 8,
        mem= lambda wildcards: 24000 if not wildcards.chrom in ['chr3', 'chr16'] else 48000,
    shell:
        """
        pindel \\
            --exclude {input.exclude} \\
            --report_inversions 0 \\
            --report_duplications 0 \\
            --anchor_quality 10 \\
            --minimum_support_for_event 3 \\
            --number_of_threads 8 \\
            --fasta {params.ref} \\
            --config-file {input.config} \\
            --chromosome {params.chrom} \\
            --output-prefix {params.prefix}
        """

rule pindelChromExcludable:
    input:
        normal= lambda wildcards: 'bam/%s.bam' % normalForTumourLib(wildcards.tumour_lib),
    output:
        exclude= 'pindel/{tumour_lib}.{chrom}.exclude.bed',
        done= touch('pindel/{tumour_lib}.{chrom}.exclude.done'),
    params:
        prefix= pindel_prefix,
        chrom= '{chrom}',
        jobname= '{tumour_lib}.{chrom}',
        cpus_per_task= 1,
        mem= 1000
    shell:
        """
        samtools depth -r {params.chrom} {input.normal} \\
        | awk -v OFS='\\t' '$3 > 1000 {{print $1, $2-1, $2}}' \\
        | mergeBed > {output.exclude}
        """

localrules: pindelConfigTumNorm
rule pindelConfigTumNorm:
    input:
        t_stat= 'samstats/{tumour_lib}.stats',
        t_bam= 'bam/{tumour_lib}.bam',
        n_stat= lambda wildcards: 'samstats/%s.stats' % normalForTumourLib(wildcards.tumour_lib),
        n_bam= lambda wildcards: 'bam/%s.bam' % normalForTumourLib(wildcards.tumour_lib),
    output:
        'pindel/{tumour_lib}.config'
    shell:
        # Insert size must be an INT so we need to round the float in input.
        # If you insert a float, pindel interprets the decimal part as the sample name (?!)
        """
        n_ins=`grep -P '^SN\\tinsert size average' {input.n_stat} | cut -f 3 | xargs printf "%.*f\\n" 0`
        echo -e "{input.n_bam}\\t${{n_ins}}\\tnormal" > {output}

        t_ins=`grep -P '^SN\\tinsert size average' {input.t_stat} | cut -f 3 | xargs printf "%.*f\\n" 0`
        echo -e "{input.t_bam}\\t${{t_ins}}\\ttumour" >> {output}
        """

rule pindelChromPon:
    input:
        normal= 'bam/{pon_library}.bam',
        config= 'pindel/{pon_library}.config',
        done= 'pindel/{pon_library}.{chrom}.exclude.done',
        exclude= 'pindel/{pon_library}.{chrom}.exclude.bed',
    output:
        done= touch('pindel/{pon_library}.{chrom}.done')
    params:
        prefix= pindel_prefix_pon,
        ref= REF,
        chrom= '{chrom}',
        jobname= '{pon_library}.{chrom}',
        cpus_per_task= 8,
        mem= 24000
    shell:
        """
        pindel \\
            --exclude {input.exclude} \\
            --report_inversions 0 \\
            --report_duplications 0 \\
            --anchor_quality 10 \\
            --minimum_support_for_event 3 \\
            --number_of_threads 8 \\
            --fasta {params.ref} \\
            --config-file {input.config} \\
            --chromosome {params.chrom} \\
            --output-prefix {params.prefix}
        """

rule pindelChromExcludablePon:
    input:
        normal= 'bam/{pon_library}.bam',
    output:
        exclude= 'pindel/{pon_library}.{chrom}.exclude.bed',
        done= touch('pindel/{pon_library}.{chrom}.exclude.done'),
    params:
        prefix= pindel_prefix_pon,
        chrom= '{chrom}',
        jobname= '{pon_library}.{chrom}',
        cpus_per_task= 1,
        mem= 1000
    shell:
        """
        samtools depth -r {params.chrom} {input.normal} \\
        | awk -v OFS='\\t' '$3 > 1000 {{print $1, $2-1, $2}}' \\
        | mergeBed > {output.exclude}
        """

localrules: pindelConfigPon
rule pindelConfigPon:
    input:
        n_stat= 'samstats/{pon_library}.stats',
        n_bam= 'bam/{pon_library}.bam',
    output:
        'pindel/{pon_library}.config'
    shell:
        # Insert size must be an INT so we need to round the float in input.
        # If you insert a float, pindel interprets the decimal part as the sample name (?!)
        """
        n_ins=`grep -P '^SN\\tinsert size average' {input.n_stat} | cut -f 3 | xargs printf "%.*f\\n" 0`
        echo -e "{input.n_bam}\\t${{n_ins}}\\tnormal" > {output}
        """

rule facets_run:
    input:
        'facets/{tumour_lib}.csv.gz'
    output:
        bed= 'facets/{tumour_lib}.cn.bed',
        pdf= 'facets/{tumour_lib}.cn.pdf',
        spider= 'facets/{tumour_lib}.spider.pdf'
    params:
        jobname= '{tumour_lib}',
        cpus_per_task= 16,
        mem= 48000,
        bin_dir= config['bin_dir'],
    shell:
        """
        {params.bin_dir}/run_facets.R \\
            --pileup {input} \\
            --out {output.bed} \\
            --pdf {output.pdf} \\
            --pdf_spider {output.spider} \\
            --chrom \\
            --snp_nbhd 500 \\
            --cval 400 \\
            --rnd_seed 2222
        """

rule facets_pileup:
    input:
        tumour= 'bam/{tumour_lib}.bam',
        tbai= 'bam/{tumour_lib}.bam.bai',
        normal= lambda wildcards: 'bam/%s.bam' % normalForTumourLib(wildcards.tumour_lib),
        nbai= lambda wildcards: 'bam/%s.bam.bai' % normalForTumourLib(wildcards.tumour_lib),
    output:
        'facets/{tumour_lib}.csv.gz' # This can be made temp()
    params:
        jobname= '{tumour_lib}',
        cpus_per_task= 1,
        mem= 2000,
        dbsnp= config['dbsnp'],
    shell:
        """
        ~/local/lib64/R/library/facets/extcode/snp-pileup \\
            --gzip \\
            --min-map-quality 15 \\
            --min-base-quality 20 \\
            --pseudo-snps 100 \\
            --min-read-counts 20,0 \\
           {params.dbsnp} {output} {input.normal} {input.tumour}
        """


rule vepMutect:
    input:
        vcf= 'gatk4/{tumour_lib}.vcf.gz',
        tbi= 'gatk4/{tumour_lib}.vcf.gz.tbi'
    output:
        vcf= 'gatk4/{tumour_lib}.vep.vcf.gz',
        tbi= 'gatk4/{tumour_lib}.vep.vcf.gz.tbi',
        stats= 'gatk4/{tumour_lib}.vep.summary.html'
    params:
        jobname= '{tumour_lib}',
        cpus_per_task= 4,
        mem= 4000,
        bin_dir= config['bin_dir'],
        dir_cache= config['vep_dir_cache'],
    shell:
        """
        ~/applications/vep/ensembl-vep-release-91.2/vep --flag_pick --offline --cache --dir_cache {params.dir_cache} --fork 4 --vcf --force_overwrite -i {input.vcf} -o STDOUT --stats_file {output.stats} \\
        | {params.bin_dir}/vepReorder.py \\
        | bgzip > {output.vcf}
        sleep 120
        tabix {output.vcf}
        """

rule filterMutectCalls:
    input:
        vcf= 'gatk4/{tumour_lib}.tmp.vcf.gz',
        tbi= 'gatk4/{tumour_lib}.tmp.vcf.gz.tbi',
    output:
        vcf= temp('gatk4/{tumour_lib}.vcf.gz'),
        tbi= temp('gatk4/{tumour_lib}.vcf.gz.tbi')
    params:
        jobname= '{tumour_lib}',
        cpus_per_task= 2,
        mem= 4000,
    shell:
        """
        gatk FilterMutectCalls --variant {input.vcf} --output {output.vcf}
        """

rule mergeMutect:
    input:
        vcf= lambda wildcards: expand("gatk4/{tumour_lib}.{chrom}.tmp.vcf.gz", tumour_lib= wildcards.tumour_lib, chrom= CHROMS),
        tbi= lambda wildcards: expand("gatk4/{tumour_lib}.{chrom}.tmp.vcf.gz.tbi", tumour_lib= wildcards.tumour_lib, chrom= CHROMS),
    output:
        vcf= temp('gatk4/{tumour_lib}.tmp.vcf.gz'),
        tbi= temp('gatk4/{tumour_lib}.tmp.vcf.gz.tbi')
    params:
        jobname= '{tumour_lib}',
        cpus_per_task= 1,
        mem= 4000,
    shell:
        """
        bcftools concat {input.vcf} | bgzip > {output.vcf}
        tabix -f {output.vcf}
        """

rule mutect:
    input:
        tumour= 'bam/{tumour_lib}.bam',
        tbai= 'bam/{tumour_lib}.bam.bai',
        normal= lambda wildcards: 'bam/%s.bam' % normalForTumourLib(wildcards.tumour_lib),
        nbai= lambda wildcards: 'bam/%s.bam.bai' % normalForTumourLib(wildcards.tumour_lib),
        pon= 'gatk4/pon/panelOfNormals.vcf.gz',
        idx= 'gatk4/pon/panelOfNormals.vcf.gz.tbi'
    output:
        vcf= temp('gatk4/{tumour_lib}.{chrom}.tmp.vcf.gz'),
        tbi= temp('gatk4/{tumour_lib}.{chrom}.tmp.vcf.gz.tbi'),
    params:
        ref= REF,
        chrom= '{chrom}',
        tumour_lib= '{tumour_lib}',
        normal_lib= lambda wildcards: '%s' % normalForTumourLib(wildcards.tumour_lib),
        jobname= '{tumour_lib}.{chrom}',
        cpus_per_task= 1,
        mem= 4000,
    shell:
        """
        gatk --java-options '-Xmx3500m' Mutect2 \\
            -R {params.ref} \\
            -I {input.tumour} \\
            -tumor {params.tumour_lib} \\
            -I {input.normal} \\
            -normal {params.normal_lib} \\
            --panel-of-normals {input.pon} \\
            -L {params.chrom} \\
            -O {output.vcf}
        """

normal_variants= expand('gatk4/pon/{pon_library}.vcf.gz', pon_library= pon_library)
rule combinePanelOfNormals:
    input:
        normal_variants
    output:
        vcf= 'gatk4/pon/panelOfNormals.vcf.gz',
        idx= 'gatk4/pon/panelOfNormals.vcf.gz.tbi'
    params:
        ref= REF,
        vcfs= ' '.join(['-vcfs ' + re.sub('.tbi$', '', x) for x in normal_variants]),
        jobname= 'panelOfNormals',
        cpus_per_task= 2,
        mem= 8000,
    shell:
        """
        gatk CreateSomaticPanelOfNormals \\
            {params.vcfs} \\
            --output {output.vcf} \\
            --TMP_DIR gatk4
        """

rule mergeMutectNormals:
    input:
        # See https://groups.google.com/forum/#!topic/Snakemake/eYHSn3EoA6Y for explanation of this lambda
        vcf= lambda wildcards: expand("gatk4/pon/{pon_library}.{chrom}.vcf.gz", pon_library= wildcards.pon_library, chrom= CHROMS),
        tbi= lambda wildcards: expand("gatk4/pon/{pon_library}.{chrom}.vcf.gz.tbi", pon_library= wildcards.pon_library, chrom= CHROMS)
    output:
        vcf= 'gatk4/pon/{pon_library}.vcf.gz',
        tbi= 'gatk4/pon/{pon_library}.vcf.gz.tbi'
    params:
        jobname= '{pon_library}',
        cpus_per_task= 1,
        mem= 2000,
    shell:
        """
        bcftools concat {input.vcf} | bgzip > {output.vcf}
        tabix -f {output.vcf}
        """

rule mutectNormal:
    input:
        bam= 'bam/{pon_library}.bam',
        bai= 'bam/{pon_library}.bam.bai',
    output:
        vcf= temp('gatk4/pon/{pon_library}.{chrom}.vcf.gz'),
        tbi= temp('gatk4/pon/{pon_library}.{chrom}.vcf.gz.tbi')
    params:
        chrom= '{chrom}',
        ref= REF,
        pon_library= '{pon_library}',
        jobname= '{pon_library}.{chrom}',
        cpus_per_task= 1,
        mem= 3000,
    shell:
        """
        gatk --java-options '-Xmx2500m' Mutect2 \\
            -R {params.ref} \\
            -I {input.bam} \\
            -tumor {params.pon_library} \\
            -L {params.chrom} \\
            -O {output.vcf} \\
            --TMP_DIR gatk4
        """

chrom_sizes= re.sub('\.fa$', '.chrom.sizes', REF)
rule bamToTdf:
    input:
        bam= 'bam/{library}.bam',
        genome= chrom_sizes 
    output:
        'tdf/{library}.tdf'
    params:
        jobname= '{library}',
        cpus_per_task= 1,
        mem= 4000,
    shell:
        """
        igvtools count --minMapQuality 5 {input.bam} {output} {input.genome}
        """

rule samtoolsStats:
    input:
        'bam/{library}.bam',
    output:
        'samstats/{library}.stats'
    params:
        jobname= '{library}',
        cpus_per_task= 2,
        mem= 2000,
    shell:
        """
        samtools stats -@ {params.cpus_per_task} {input} > {output}
        """

def get_fq1(wildcards):
    fq= list(manifest[(manifest['library'] == wildcards.library) & (manifest['rmate'] == 'R1')]['fastq'])
    fq.sort()
    return fq

def get_fq2(wildcards):
    fq= list(manifest[(manifest['library'] == wildcards.library) & (manifest['rmate'] == 'R2')]['fastq'])
    fq.sort()
    return fq

rule bamsort:
    input:
        'bam/{library}.tmp.bam'    
    output:
        bam= 'bam/{library}.bam',
        bai= 'bam/{library}.bam.bai'
    params:
        dupmetrics= 'bam/{library}.dupmetrics.txt',
        jobname= '{library}',
        cpus_per_task= 16,
        mem= 24000,
    shell:
        # bamsort: We mark duplicates but we don't remove them. If removed, you
        # get "Mate not found for paired read" from picard ValidateSamFile.
        """
        bamsort inputformat=bam markduplicates=1 rmdup=0 fixmates=1 inputthreads=12 outputthreads=12 \\
          M={params.dupmetrics} I={input} O={output.bam} index=1 indexfilename={output.bai}
        sleep 60
        touch -c {output.bai}
        """

rule align:
    input:
        REF + '.bwt',
        fq1= get_fq1,
        fq2= get_fq2,
    output:
        bam= temp('bam/{library}.tmp.bam'),
    params:
        ref= REF,
        trimlog= 'bam/{library}.trim.log',
        RG= r'@RG\tID:{id}\tSM:{sm}\tLB:{lb}\tPL:ILLUMINA\tPU:NA'.format(id= '{library}', sm= '{library}', lb= '{library}'),
        jobname= '{library}',
        cpus_per_task= 8,
        mem= 24000
    shell:
        """
        catInterleaveFastq.sh -1 {input.fq1} \\
                              -2 {input.fq2} \\
        | bbduk.sh in=stdin.fq out=stdout.fq interleaved=true qtrim=rl minlength=15 trimq=12 2> {params.trimlog} \\
        | bwa mem -R '{params.RG}' -t 10 -p {params.ref} - \\
        | samtools view -b - > {output}
        """

rule multiqc:
    input:
        expand('fastqc/{fastq_base}_fastqc.zip', fastq_base= manifest.fastq_base)
    output:
        'multiqc/multiqc_report.html'
    params:
        jobname= 'multiqc',
        cpus_per_task= 2,
        mem= 2000,
    shell:
        """
        multiqc --force --zip-data-dir --filename {output} fastqc/
        """

rule fastqc:
    input:
        # {fastq_dir} is not used in output so we use lambda expression to have
        # only {fastq_base} used as wildcard.
        # Note the zip function to have fastq_dir paired to fastq_base in a
        # ono-to-one way rather than all vs all.
        lambda wildcards: expand('{fastq_dir}/{fastq_base}.fastq.gz', zip, fastq_dir= manifest.fastq_dir, fastq_base= wildcards.fastq_base)
    output:
        'fastqc/{fastq_base}_fastqc.zip'
    params:
        jobname= '{fastq_base}',
        cpus_per_task= 2,
        mem= 2000,
    shell:
        """
        fastqc --noextract --outdir fastqc/ {input}
        """

localrules: igvChromSizes
rule igvChromSizes:
    input:
        REF + '.fai'
    output:
        temp(chrom_sizes)
    shell:
        """
        cp {input} {output}
        """

localrules: faidx
rule faidx:
    input:
        REF
    output:
        REF + '.fai'
    shell:
        """
        samtools faidx {input}
        """

rule bwaIndex:
    input:
        ancient(REF)
    output:
        REF + '.bwt'
    params:
        jobname= 'bwaIndex',
        cpus_per_task= 16, # You don't need many threads but put this high enough to prevent slurm overbooking MEM
        mem= 32000,
    shell:
        """
        bwa index {input}
        """

localrules: simpleRepeat
rule simpleRepeat:
    output:
        bed= 'pindel/pon/simpleRepeat.bed.gz',
        tbi= 'pindel/pon/simpleRepeat.bed.gz.tbi',
    shell:
        """
        curl -s -S http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz \\
        | zcat \\
        | awk '$6 <= 6' \\
        | cut -f 2- \\
        | bgzip > {output.bed}
        tabix -f {output.bed}
    """
