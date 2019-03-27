#!/usr/bin/env nextflow

// Nextflow pipeline for running Mixcr / V'DJer / RSEM

params.star_ref = '/datastore/nextgenout4/share/labs/bioinformatics/seqware/ref/hg38/star_hg38_d1'
params.mixcr_species = 'hsa'
params.transcript_fasta = '/datastore/nextgenout4/seqware-analysis/lmose/ref/transcriptome/hg19_transcripts.hg19.20091027.fa'
params.vdjer_igh = '/datastore/nextgenout4/seqware-analysis/lmose/ref/vdjer/igh'
params.vdjer_extra_params = ''

//params.star_ref = '/datastore/nextgenout4/seqware-analysis/lmose/ref/mouse/mm10/star'
//params.mixcr_species = 'mmu'
//params.transcript_fasta = '/datastore/nextgenout4/seqware-analysis/lmose/ref/mouse/transcript/gencode.vM13.transcripts.fa'
//params.vdjer_igh = '/datastore/nextgenout4/seqware-analysis/lmose/vdjer/ref/mouse/igh'

// Mouse + sensitive mode
//params.vdjer_extra_params = '--vr chr12:113572929-116009954 --cr chr12:113256204-116009954 --mq 60 --mf 2 --rs 25 --ms 2 --mcs -5.5'

params.output_dir = '/datastore/nextgenout4/seqware-analysis/lmose/nextflow/mixcr_vdjer/anti_pd1_test'
params.input = '/datastore/nextgenout4/seqware-analysis/lmose/nextflow/mixcr_vdjer/anti_pd1_test/inventory.txt'

Channel
    .fromPath(params.input)
    .splitCsv(header:true, sep:'\t')
    .map{ row-> tuple(row.sample, file(row.fq1), file(row.fq2)) }
    .set { input_channel }

input_channel.into {
    input_channel1
    input_channel2
    input_channel3
    input_channel4
}


//star_read_command = "cat"
//if (params.fq1.endsWith(".gz")) {
    star_read_command = "zcat"
//}

//
// MiXCR pipeline (Runs MiXCR, FastQC, VDJTools)
//
process mixcr {
    cpus 8
    memory '32 GB'

    input:
    set sample, file(fq1), file(fq2) from input_channel1

    output:
    file 'mixcr_output/*' into mixcr_output
    file 'mixcr.log' into mixcr_log

    """
    hostname;
    mkdir mixcr_output;
    WD=`pwd`;
    docker run --rm=true  -v /datastore:/datastore:shared dockerreg.bioinf.unc.edu:5000/mixcr_2.1.9:4 bash -c \
        "cd \$WD;
        source /import/run_mixcr.sh \
        --import_dir=/import \
        --chains=ALL \
        --rna_seq=true \
        --use_existing_vdjca=true \
        --species=${params.mixcr_species} \
        --threads=8 \
        --r1_path=../${fq1} \
        --r2_path=../${fq2} \
        --output_dir=mixcr_output \
        --sample_name=${sample} \
        --separate_by_c=false
        " > mixcr.log 2>&1
    """
}

//
// 1-pass align with STAR (include unmapped reads)
//
process star_map {
    cpus 8
    memory '40 GB'

    input:
    set sample, file(fq1), file(fq2) from input_channel2
    
    output:
    file 'STAR_Aligned.out.bam' into star_genome

    """
    echo "Processing ${sample} ${fq1} ${fq2}"
    WD=`pwd`;
    docker run --rm=true  -v /datastore:/datastore:shared dockerreg.bioinf.unc.edu:5000/mozack/vdjer_test:0.3 sh -c \
        "cd \$WD;
        STAR \
        --genomeDir ${params.star_ref} \
        --readFilesIn ${fq1} ${fq2} \
        --readFilesCommand ${star_read_command} \
        --outSAMunmapped Within \
        --outSAMtype BAM Unsorted \
        --runThreadN 16 \
        --outFileNamePrefix STAR_
        "
    """
}

//
// Sort/index STAR output
//
process sort_and_index {
    cpus 8
    memory '40 GB'
    
    input:
    file input_bam from star_genome

    output:
    file 'star.sort.bam' into sorted_star_genome_bam, sorted_star_genome_bam2
    file 'star.sort.bam.bai' into sorted_star_genome_index, sorted_star_genome_index2

    """
    WD=`pwd`;
    docker run --rm -v /datastore:/datastore:shared biocontainers/samtools:v1.7.0_cv3 sh -c \
        "cd \$WD;
        samtools sort -@ 8 -m 3G ${input_bam} -o star.sort.bam;
        samtools index -@ 8 star.sort.bam"
    """
}

//
// Use BWA to estimate insert length
//
process estimate_insert_len {
    cpus 8
    memory '16 GB'

    input:
    set sample, file(fq1), file(fq2) from input_channel4

    output:
    file 'insert_size.txt' into insert_size

    """
    WD=`pwd`;
    docker run --rm=true  -v /datastore:/datastore:shared dockerreg.bioinf.unc.edu:5000/mozack/vdjer_test:0.3 sh -c \
        "cd \$WD;
        bwa mem -t 8 -S /datastore/nextgenout4/seqware-analysis/lmose/ref/transcriptome/hg19_transcripts.hg19.20091027.fa \
             '<zcat ${fq1} | head -4000000' '<zcat ${fq2} | head -4000000' 2> bwa.head.log | samtools view -1 -bS -F 0xC -f 0x02 - > bwa.head.bam
        samtools view bwa.head.bam | python /vdjer-0.12//insert_len2.py > insert_length.txt
        tail -n 1 insert_length.txt | cut -f 7 > insert_size.txt
        "
    """
}

//
// Run V'DJer on IGH chain
//
process vdjer_igh {
    cpus 8
    memory '64 GB'

    input:
    file input_bam from sorted_star_genome_bam
    file input_bai from sorted_star_genome_index
    file input_insert_size from insert_size

    output:
    file 'vdj_contigs.fa' into igh_vdjer_contigs, igh_vdjer_contigs2
    file 'vdjer_IGH.sam' into igh_vdjer_sam
    file 'vdjer_IGH.log' into igh_vdjer_log

    """
    WD=`pwd`;
    docker run --rm=true  -v /datastore:/datastore:shared dockerreg.bioinf.unc.edu:5000/mozack/vdjer_test:0.3 sh -c \
        "cd \$WD;
        vdjer --in ${input_bam} --ins `cat insert_size.txt` --t 8 --chain IGH --ref-dir /datastore/nextgenout4/seqware-analysis/lmose/ref/vdjer/igh ${params.vdjer_extra_params} 2> vdjer_IGH.log > vdjer_IGH.sam
        "
    """
}

//
// Quantify V'DJer output using RSEM
//
process rsem_igh {
    cpus 8
    memory '32 GB'

    input:
    file igh_contigs_input from igh_vdjer_contigs
    file igh_sam_input from igh_vdjer_sam

    output:
    file 'rsem_results.isoforms.results' into rsem_igh_quant

    """
    if [ -s ${igh_contigs_input} ]; then
      WD=`pwd`;
      docker run --rm=true  -v /datastore:/datastore:shared dockerreg.bioinf.unc.edu:5000/mozack/vdjer_test:0.3 sh -c \
        "cd \$WD;
	rsem.bash ${igh_sam_input}
	"
    else
      touch rsem_results.isoforms.results
    fi
    """
}

//
// Copy files to final output destination
//
process collect_output {
    cpus 1
    memory '2 GB'

    input:
    set sample, file(fq1), file(fq2) from input_channel3
    file rsem_igh_quant
    file igh_vdjer_contigs2
    file igh_vdjer_log
    file sorted_star_genome_bam2
    file sorted_star_genome_index2
    file mixcr_output
    file mixcr_log

    """
    mkdir -p ${params.output_dir}/${sample}
    cp ${rsem_igh_quant} ${params.output_dir}/${sample}
    cp ${igh_vdjer_contigs2} ${params.output_dir}/${sample}
    cp ${igh_vdjer_log} ${params.output_dir}/${sample}
    cp ${sorted_star_genome_bam2} ${params.output_dir}/${sample}
    cp ${sorted_star_genome_index2} ${params.output_dir}/${sample}
    mkdir -p ${params.output_dir}/${sample}/mixcr_output
    cp ${mixcr_output} ${params.output_dir}/${sample}/mixcr_output
    cp ${mixcr_log} ${params.output_dir}/${sample}
    """
}
