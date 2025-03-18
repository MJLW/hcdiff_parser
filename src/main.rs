use std::io::Write;
use std::{env, fs::File};

use polars::prelude::*;
use polars_rows_iter::*;
use rust_htslib::{bcf::Read, faidx::Reader};

struct Region<'a> {
    contig: &'a str,
    start: usize,
    end: usize, // Open ended (+1)
}

#[derive(Debug, FromDataFrameRow, PartialEq)] // for assert_eq
struct MyRow<'a> {
    #[column("Chromosome")]
    contig: &'a str,
    #[column("Start position")]
    position: i64,
    #[column("Reference")]
    reference: &'a str,
    #[column("Variant")]
    alt: &'a str,
    #[column("De novo assessment")]
    denovo: &'a str,
}

fn main() {
    let args: Vec<String> = env::args().collect();

    /*
     * Read in reference file
     */
    // let reference_path = "data/hs_ref_GRCh37.p5_all_contigs.fa";
    let reference_path = &args[1];
    let reader = Reader::from_path(reference_path)
        .expect(format!("Could not read reference file at {}.", reference_path).as_str());

    /*
     * Read in HCDiff file
     */
    let hcdiff_columns: Vec<Expr> = vec![
        "Chromosome",
        "Start position",
        "Reference",
        "Variant",
        "De novo assessment",
    ]
    .iter()
    .map(|name| col(*name))
    .collect();

    // let hcdiff_path = "data/DNA-003114B_hcdiffs.denovo.txt";
    let hcdiff_path = &args[2];
    let hcdiff = LazyCsvReader::new(hcdiff_path)
        .with_has_header(true)
        .with_infer_schema_length(Some(10000))
        .with_separator(b'\t')
        .finish()
        .expect(format!("Could not read HCDiff file at {}", hcdiff_path).as_str())
        .select(hcdiff_columns)
        .with_columns([
            when(col("Variant").is_null())
                .then(col("Start position") - lit(1))
                .otherwise(col("Start position"))
                .alias("Start position"),
            when(col("Reference").is_null())
                .then(lit(""))
                .otherwise(col("Reference"))
                .alias("Reference"),
            when(col("Variant").is_null())
                .then(lit(""))
                .otherwise(col("Variant"))
                .alias("Variant"),
            when(col("De novo assessment").is_null())
                .then(lit(""))
                .otherwise(col("De novo assessment"))
                .alias("De novo assessment"),
        ])
        .collect()
        .expect("Could not find target columns in HCDiff, this is likely not a trio.");

    let hcdiff_clone = hcdiff.clone();
    let hcdiff_rows: Vec<MyRow> = hcdiff_clone
        .rows_iter::<MyRow>()
        .unwrap()
        .map(|row| row.unwrap())
        .collect();

    let reference_base: Vec<String> = hcdiff_rows
        .iter()
        .map(|row| {
            let region = Region {
                contig: row.contig,
                start: (row.position - 1) as usize,
                end: (row.position - 1) as usize,
            };

            reader
                .fetch_seq_string(region.contig, region.start, region.end)
                .expect(
                    format!(
                        "Could not read region {}:{}-{} from reference",
                        region.contig, region.start, region.end
                    )
                    .as_str(),
                )
        })
        .collect();

    let preceding_reference = Series::new(PlSmallStr::from("Reference base"), reference_base);

    let hcdiff_with_preceding = hcdiff
        .clone()
        .with_column(preceding_reference)
        .unwrap()
        .clone()
        .lazy()
        .with_columns([
            when(
                col("Reference")
                    .str()
                    .len_chars()
                    .eq(lit(0))
                    .or(col("Variant").str().len_chars().eq(lit(0))),
            )
            .then(col("Reference base") + col("Reference"))
            .otherwise(col("Reference"))
            .alias("Reference"),
            when(
                col("Reference")
                    .str()
                    .len_chars()
                    .eq(lit(0))
                    .or(col("Variant").str().len_chars().eq(lit(0))),
            )
            .then(col("Reference base") + col("Variant"))
            .otherwise(col("Variant"))
            .alias("Variant"),
        ])
        .collect()
        .unwrap();

    let mut rows: Vec<MyRow> = hcdiff_with_preceding
        .rows_iter::<MyRow>()
        .unwrap()
        .map(|row| row.unwrap())
        .collect();

    println!("{}", hcdiff_with_preceding);

    // let vcf_path = "data/DNA-003114B_haplotypecaller.vcf.gz";
    let vcf_path = &args[3];
    let mut in_vcf = rust_htslib::bcf::Reader::from_path(vcf_path).expect("Could not VCF");
    let header_view = in_vcf.header().clone();

    rows.sort_by(|a, b| {
        header_view
            .name2rid(a.contig.as_bytes())
            .unwrap()
            .cmp(&header_view.name2rid(b.contig.as_bytes()).unwrap())
            .then(a.position.cmp(&b.position))
    });

    let header = r#"##fileformat=VCFv4.1
##contig=<ID=chr1,length=249250621>
##contig=<ID=chr2,length=243199373>
##contig=<ID=chr3,length=198022430>
##contig=<ID=chr4,length=191154276>
##contig=<ID=chr5,length=180915260>
##contig=<ID=chr6,length=171115067>
##contig=<ID=chr7,length=159138663>
##contig=<ID=chr8,length=146364022>
##contig=<ID=chr9,length=141213431>
##contig=<ID=chr10,length=135534747>
##contig=<ID=chr11,length=135006516>
##contig=<ID=chr12,length=133851895>
##contig=<ID=chr13,length=115169878>
##contig=<ID=chr14,length=107349540>
##contig=<ID=chr15,length=102531392>
##contig=<ID=chr16,length=90354753>
##contig=<ID=chr17,length=81195210>
##contig=<ID=chr18,length=78077248>
##contig=<ID=chr19,length=59128983>
##contig=<ID=chr20,length=63025520>
##contig=<ID=chr21,length=48129895>
##contig=<ID=chr22,length=51304566>
##contig=<ID=chrX,length=155270560>
##contig=<ID=chrY,length=59373566>
##contig=<ID=chrMT,length=16569>
##reference=file:///ifs/data/lib/genomes/human/GRCh37.p5/hs_ref_GRCh37.p5_all_contigs.fa
##INFO=<ID=DENOVO,Number=1,Type=String,Description="De novo field">
"#;

    let output_path = "output.vcf";
    let mut file =
        File::create(output_path).expect(format!("Could not open file: {}", output_path).as_str());

    let _ = file.write(header.as_bytes());
    let _ = file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n".as_bytes());

    let mut hcdiff_row_index: usize = 0;
    for record_result in in_vcf.records() {
        let record = record_result.expect("Failed to read record.");

        if record.allele_count() > 2 {
            let alt_allele_count = record.allele_count() - 1;
            hcdiff_row_index += alt_allele_count as usize;
            continue;
        }

        let chr =
            std::str::from_utf8(header_view.rid2name(record.rid().unwrap()).unwrap()).unwrap();
        let position = record.pos();
        let reference = std::str::from_utf8(record.alleles()[0]).unwrap();

        let mut annotated_record = record.clone();
        annotated_record.unpack();

        for j in 1..record.allele_count() {
            let alt = std::str::from_utf8(record.alleles()[j as usize]).unwrap();

            let row = rows.get(hcdiff_row_index).unwrap();

            if (chr != row.contig)
                | (position != row.position - 1)
                | (reference != row.reference)
                | (alt != row.alt)
            {
                println!("MISMATCH~!");
                println!(
                    "HCDIFF: {}\t{}\t{}\t{}",
                    row.contig,
                    row.position - 1,
                    row.reference,
                    row.alt
                );
                println!("VCF: {}\t{}\t{}\t{}", chr, position, reference, alt);
                continue;
            }

            let denovo = row.denovo.replace(" ", "_");

            let line = format!(
                "{}\t{}\t.\t{}\t{}\t.\t.\tDENOVO={}\n",
                chr,
                position + 1,
                reference,
                alt,
                denovo
            );

            let _ = file.write(line.as_bytes());

            // let tag = "DENOVO".as_bytes();
            // let values = &(vec![row.denovo.as_bytes()] as Vec<&[u8]>);
            // let mut out_record = out_vcf.empty_record();
            // out_record.set_rid(annotated_record.rid());
            // out_record.set_pos(annotated_record.pos());
            // let _ = out_record
            //     .set_alleles(&annotated_record.alleles())
            //     .expect("Failed to set alleles");
            // out_record.set_qual(annotated_record.qual());
            //
            // out_record
            //     .push_info_string(tag, values)
            //     .expect("Could not add to tag INFO/DENOVO");

            // in_vcf.header().header_records().iter().map(f)
            // in_vcf.header()

            // out_record
            //     .push_format_string(tag, values)
            //     .expect("Could not add to tag DENOVO");

            // out_vcf.write(&out_record).unwrap();

            hcdiff_row_index += 1;
        }
    }
}
