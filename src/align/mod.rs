pub mod cli;
pub mod engine;
use engine::{align_with_config, AlignPlusConfig};

pub fn run(matches: &clap::ArgMatches) -> anyhow::Result<()> {
    let reference = matches.get_one::<String>("reference").unwrap();
    let queries: Vec<String> = matches
        .get_many::<String>("query")
        .unwrap()
        .cloned()
        .collect();
    let output = matches.get_one::<String>("output").unwrap();
    let is_hpc = matches.get_flag("hpc");
    let kmer = matches.get_one::<i16>("kmer").copied();
    let window = matches.get_one::<i16>("window").copied();
    let mid_occ_frac = matches
        .get_one::<(f32, Option<f32>)>("mid_occ_frac")
        .copied();
    let mask_level = matches.get_one::<f32>("mask_level").copied();
    let bounds_of_occurrence = matches
        .get_one::<(i32, Option<i32>)>("bounds_of_occurrence")
        .copied();
    let max_gap = matches.get_one::<i32>("max_gap").copied();
    let max_gap_ref = matches.get_one::<i32>("max_gap_ref").copied();
    let max_frag_len = matches.get_one::<i32>("max_frag_len").copied();
    let bw = matches.get_one::<(i32, Option<i32>)>("bw").copied();
    let min_cnt = matches.get_one::<i32>("min_cnt").copied();
    let min_chain_score = matches.get_one::<i32>("min_chain_score").copied();

    let matching_score = matches.get_one::<i32>("matching_score").copied();
    let mismatch_penalty = matches.get_one::<i32>("mismatch_penalty").copied();
    let gap_open = matches.get_one::<(i32, Option<i32>)>("gap_open").copied();
    let gap_extension = matches
        .get_one::<(i32, Option<i32>)>("gap_extension")
        .copied();
    let z_drop = matches.get_one::<(i32, Option<i32>)>("z_drop").copied();
    let min_dp_max = matches.get_one::<i32>("min_dp_max").copied();

    let max_qlen = matches.get_one::<i32>("max_qlen").copied();
    let batch_size = matches.get_one::<u64>("batch_size").unwrap();
    let soft_clip = matches.get_flag("soft_clip");
    let output_bam_format = matches.get_flag("output_bam");
    let output_cigar = matches.get_flag("output_cigar");
    let eqx = matches.get_flag("eqx");
    let cs = matches.get_one::<String>("cs").cloned();
    let secondary = matches.get_one::<String>("secondary").unwrap();
    let best_n = matches.get_one::<i32>("best_n").copied();
    let pri_ratio = matches.get_one::<f32>("pri_ratio").copied();
    let mini_batch_size = matches.get_one::<i64>("mini_batch_size").unwrap();
    let seed = matches.get_one::<i32>("seed").copied();
    let preset = matches.get_one::<String>("preset").unwrap();
    let porec_gap_rescue = matches.get_flag("porec_gap_rescue");
    let rescue_k = matches.get_one::<i16>("rescue_k").copied();
    let rescue_w = matches.get_one::<i16>("rescue_w").copied();
    let min_gap_len = matches.get_one::<usize>("min_gap_len").copied().unwrap_or(100);

    let threads = matches.get_one::<usize>("threads").unwrap();
    let realign = !matches.get_flag("no_realign");
    let rescue_mode = matches.get_one::<String>("rescue_mode").cloned();
    let mapq_rescue = matches.get_one::<u8>("mapq_rescue").copied();
    let porec_mapq_calibrate = matches.get_flag("porec_mapq_calibrate");
    let candidate_probability = matches.get_flag("candidate_probability");
    let graph_assignment = matches.get_flag("graph_assignment");
    let re_site = matches.get_one::<String>("re_site").cloned();
    let secondary = match secondary.as_str() {
        "yes" => true,
        "no" => false,
        _ => false,
    };



    let config = AlignPlusConfig {
        reference: reference.clone(),
        queries,
        output: output.clone(),
        output_bam_format,
        output_cigar,
        eqx,
        cs,
        kmer,
        window,
        is_hpc,
        batch_size: *batch_size,
        soft_clip,
        secondary,
        mid_occ_frac,
        bounds_of_occurrence,
        max_gap,
        max_gap_ref,
        max_frag_len,
        mask_level,
        bw,
        min_cnt,
        min_chain_score,
        matching_score,
        mismatch_penalty,
        gap_open,
        gap_extension,
        z_drop,
        min_dp_max,
        best_n,
        pri_ratio,
        max_qlen,
        mini_batch_size: *mini_batch_size,
        seed,
        porec_gap_rescue,
        rescue_k,
        rescue_w,
        min_gap_len,
        preset: preset.clone(),
        threads: *threads,
        realign,
        rescue_mode,
        mapq_rescue,
        porec_mapq_calibrate,
        candidate_probability,
        graph_assignment,
        re_site,
        methylation: matches.get_one::<String>("meth_bed").map(|bed| crate::methalign::MethConfig {
            bed: bed.clone(),
            match_score: *matches.get_one::<i32>("meth_match").unwrap(),
            ref_penalty: *matches.get_one::<i32>("meth_ref_penalty").unwrap(),
            read_penalty: *matches.get_one::<i32>("meth_read_penalty").unwrap(),
            ref_prob_cutoff: *matches.get_one::<f64>("meth_ref_cutoff").unwrap(),
            prob_cutoff: *matches.get_one::<u8>("meth_cutoff").unwrap(),
            designate_mapq: *matches.get_one::<u8>("meth_mapq").unwrap(),
            cpg: matches.get_flag("meth_cpg"),
        }),
    };

    align_with_config(config)
}
