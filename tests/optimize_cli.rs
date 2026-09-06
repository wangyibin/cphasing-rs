use cphasing::cli;
use cphasing::clm::{ClmbWriter, encode_endpoint};
use std::fs;
use std::process::Command;

fn orientation_cli_fixture(n: usize) -> tempfile::TempDir {
    let directory = tempfile::tempdir().unwrap();
    let names: Vec<_> = (0..n).map(|id| format!("ctg{id}")).collect();
    let counts = directory.path().join("group.txt");
    fs::write(&counts, names.iter().map(|name| format!("{name}\t10\t100000\n")).collect::<String>()).unwrap();
    let clmb = directory.path().join("group.clmb");
    let mut writer = ClmbWriter::create_synchronous(&clmb, &names, 4096, None, None).unwrap();
    for u in 0..n {
        for v in u + 1..n {
            let distances = if v == u + 1 { [20_000, 100_000, 100_000, 180_000] } else { [100_000; 4] };
            for (index, distance) in distances.into_iter().enumerate() {
                writer.write_record(
                    encode_endpoint(u as u32, (index / 2) as u8).unwrap(),
                    encode_endpoint(v as u32, (index % 2) as u8).unwrap(),
                    &[distance; 64],
                ).unwrap();
            }
        }
    }
    writer.finish().unwrap();
    directory
}

fn run_orientation_cli(directory: &std::path::Path, extra: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .current_dir(directory)
        .env("RUST_LOG", "info")
        .args(["optimize", "group.txt", "group.clmb", "--skipGA", "--threads", "1"])
        .args(extra).output().unwrap()
}

fn run_resumed_orientation_fixture(n: usize, input: &str, extra: &[&str]) -> (Vec<String>, String) {
    let directory = orientation_cli_fixture(n);
    fs::write(directory.path().join("group.tour"), format!("{input}\n")).unwrap();
    let mut args = vec!["--resume"];
    args.extend_from_slice(extra);
    let output = run_orientation_cli(directory.path(), &args);
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(output.status.success(), "{stderr}");
    let tour = fs::read_to_string(directory.path().join("group.tour")).unwrap();
    (tour.split_whitespace().map(str::to_string).collect(), stderr)
}

#[test]
fn resume_input_errors_preserve_the_tour_and_existing_backups() {
    for (input, message) in [
        (None, "cannot resume from"),
        (Some("\n \n"), "tour is empty"),
        (Some("ctg0+ unknown-\n"), "unknown contig"),
        (Some("ctg0+ ctg0-\n"), "duplicate contig"),
        (Some("ctg0+\n"), "expected all 2 contigs"),
    ] {
        let directory = orientation_cli_fixture(2);
        let tour = directory.path().join("group.tour");
        if let Some(input) = input { fs::write(&tour, input).unwrap(); }
        let backup = directory.path().join("group.tour.sav");
        fs::write(&backup, "older backup\n").unwrap();
        let output = run_orientation_cli(directory.path(), &["--resume"]);
        assert_eq!(output.status.code(), Some(2));
        let stderr = String::from_utf8(output.stderr).unwrap();
        assert!(stderr.contains(message), "{stderr}");
        assert!(!stderr.contains("panicked"), "{stderr}");
        assert_eq!(fs::read_to_string(&tour).ok().as_deref(), input);
        assert_eq!(fs::read_to_string(&backup).unwrap(), "older backup\n");
        assert!(!directory.path().join("group.tour.sav.1").exists());
    }
}

#[test]
fn invalid_clmb_does_not_move_a_resumed_tour_or_replace_its_backup() {
    let directory = orientation_cli_fixture(2);
    let tour = directory.path().join("group.tour");
    let backup = directory.path().join("group.tour.sav");
    fs::write(&tour, "ctg0+ ctg1-\n").unwrap();
    fs::write(&backup, "older backup\n").unwrap();
    fs::write(directory.path().join("group.clmb"), b"invalid-clmb").unwrap();
    let output = run_orientation_cli(directory.path(), &["--resume"]);
    assert_eq!(output.status.code(), Some(2));
    assert!(String::from_utf8(output.stderr).unwrap().contains("invalid CLMB optimize input"));
    assert_eq!(fs::read_to_string(tour).unwrap(), "ctg0+ ctg1-\n");
    assert_eq!(fs::read_to_string(backup).unwrap(), "older backup\n");
    assert_eq!(fs::read_dir(directory.path()).unwrap().count(), 4);
}

#[test]
fn failure_after_loading_clmb_preserves_resumed_files() {
    let directory = orientation_cli_fixture(2);
    let tour = directory.path().join("group.tour");
    let backup = directory.path().join("group.tour.sav");
    fs::write(&tour, "ctg0+ ctg1-\n").unwrap();
    fs::write(&backup, "older backup\n").unwrap();
    let output = run_orientation_cli(directory.path(), &[
        "--resume", "--endpoint-multiscale", "--split-contacts", "missing.contacts",
    ]);
    assert!(!output.status.success());
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stderr.contains("missing.contacts"), "{stderr}");
    assert_eq!(fs::read_to_string(tour).unwrap(), "ctg0+ ctg1-\n");
    assert_eq!(fs::read_to_string(backup).unwrap(), "older backup\n");
    assert_eq!(fs::read_dir(directory.path()).unwrap().count(), 4);
}

#[test]
fn successful_resume_keeps_numbered_backups_and_end_greedy_input_signs() {
    let directory = orientation_cli_fixture(2);
    let original = ">old history\nctg0- ctg1-\n\n \n";
    let tour = directory.path().join("group.tour");
    fs::write(&tour, original).unwrap();
    fs::write(directory.path().join("group.tour.sav"), "older backup\n").unwrap();
    fs::write(directory.path().join("group.tour.sav.1"), "another backup\n").unwrap();
    let args = ["--resume", "--initializer", "end-greedy", "--split-contacts", "unused",
        "--orientation-method", "banded", "--orientation-min-confidence", "1"];
    let first = run_orientation_cli(directory.path(), &args);
    assert!(first.status.success(), "{}", String::from_utf8_lossy(&first.stderr));
    assert_eq!(fs::read_to_string(&tour).unwrap().trim(), "ctg0- ctg1-");
    assert_eq!(fs::read_to_string(directory.path().join("group.tour.sav.2")).unwrap(), original);
    let first_result = fs::read(&tour).unwrap();
    let second = run_orientation_cli(directory.path(), &args);
    assert!(second.status.success(), "{}", String::from_utf8_lossy(&second.stderr));
    assert_eq!(fs::read(directory.path().join("group.tour.sav.3")).unwrap(), first_result);
    assert_eq!(fs::read_to_string(directory.path().join("group.tour.sav")).unwrap(), "older backup\n");
    assert_eq!(fs::read_to_string(directory.path().join("group.tour.sav.1")).unwrap(), "another backup\n");
    assert_eq!(fs::read_dir(directory.path()).unwrap().count(), 7);
}

#[test]
fn banded_default_gate_accepts_decisive_evidence_with_a_nonzero_prior() {
    let directory = orientation_cli_fixture(21);
    let names: Vec<_> = (0..21).map(|id| format!("ctg{id}")).collect();
    let mut writer = ClmbWriter::create_synchronous(
        directory.path().join("group.clmb"), &names, 4096, None, None,
    ).unwrap();
    // Coordinates x=100000, y=50000 on two 100000 bp contigs.
    // Only ctg0 has a directional bias; ctg1 is centered.
    for (index, distance) in [50_000, 50_000, 150_000, 150_000].into_iter().enumerate() {
        writer.write_record(
            encode_endpoint(0, (index / 2) as u8).unwrap(),
            encode_endpoint(1, (index % 2) as u8).unwrap(), &[distance; 64],
        ).unwrap();
    }
    writer.finish().unwrap();
    let input = names.iter().map(|name| format!("{name}-")).collect::<Vec<_>>().join(" ");
    for method in ["banded", "banded-contact"] {
        fs::write(directory.path().join("group.tour"), &input).unwrap();
        let output = run_orientation_cli(directory.path(), &["--resume", "--orientation-method", method]);
        let stderr = String::from_utf8(output.stderr).unwrap();
        assert!(output.status.success(), "{stderr}");
        let result = fs::read_to_string(directory.path().join("group.tour")).unwrap();
        let tokens: Vec<_> = result.split_whitespace().collect();
        assert_eq!(tokens[0], "ctg0+");
        assert!(tokens[1..].iter().all(|token| token.ends_with('-')));
        assert!(stderr.contains("rejected low-confidence changes=0"), "{stderr}");
    }
}

#[test]
fn default_solver_repairs_signs_even_when_the_entire_scaffold_must_flip() {
    let (tour, _) = run_resumed_orientation_fixture(2, "ctg0- ctg1-", &[]);
    assert_eq!(tour, ["ctg0+", "ctg1+"]);
    let (guarded, _) = run_resumed_orientation_fixture(2, "ctg0- ctg1-", &["--orientation-method", "banded"]);
    assert_eq!(guarded, ["ctg0-", "ctg1-"]);
}

#[cfg(unix)]
#[test]
fn atomic_tour_publication_preserves_normal_and_existing_permissions() {
    use std::os::unix::fs::PermissionsExt;
    let directory = orientation_cli_fixture(2);
    let tour = directory.path().join("group.tour");
    let mode = |path: &std::path::Path| fs::metadata(path).unwrap().permissions().mode() & 0o777;
    let output = run_orientation_cli(directory.path(), &[]);
    assert!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
    assert_eq!(mode(&tour), mode(&directory.path().join("group.txt")));
    fs::set_permissions(&tour, fs::Permissions::from_mode(0o640)).unwrap();
    let output = run_orientation_cli(directory.path(), &["--resume"]);
    assert!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
    assert_eq!(mode(&tour), 0o640);
    assert_eq!(mode(&directory.path().join("group.tour.sav")), 0o640);
}

#[test]
fn default_solver_repairs_a_signed_block_and_block_span_zero_disables_it() {
    let input = "ctg0+ ctg1+ ctg5- ctg4- ctg3- ctg2- ctg6+ ctg7+";
    let (tour, stderr) = run_resumed_orientation_fixture(8, input, &[]);
    assert_eq!(tour, (0..8).map(|id| format!("ctg{id}+")).collect::<Vec<_>>());
    assert!(stderr.contains("Using historical banded signed-block refinement"));
    let (without_blocks, _) = run_resumed_orientation_fixture(8, input, &["--orientation-block-span", "0"]);
    let names = |tokens: Vec<String>| tokens.into_iter().map(|token| token[..token.len()-1].to_string()).collect::<Vec<_>>();
    assert_eq!(names(without_blocks), names(input.split_whitespace().map(str::to_string).collect()));
}

#[test]
fn block_context_is_opt_in_and_rejects_incompatible_modes() {
    let matches = cli::cli().try_get_matches_from(["cphasing", "optimize", "group.txt", "group.clmb"]).unwrap();
    assert_eq!(matches.subcommand_matches("optimize").unwrap().get_one::<f64>("ORIENTATION_BLOCK_CONTEXT_WEIGHT"), Some(&0.0));
    for extra in [vec!["--orientation-method", "banded"], vec!["--orientation-block-span", "0"]] {
        let output = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
            .args(["optimize", "missing.txt", "missing.clmb", "--orientation-block-context-weight", "0.01"])
            .args(extra).output().unwrap();
        assert_eq!(output.status.code(), Some(2));
        assert!(String::from_utf8(output.stderr).unwrap().contains("--orientation-block-context-weight requires"));
    }
    let input = "ctg0+ ctg1+ ctg5- ctg4- ctg3- ctg2- ctg6+ ctg7+";
    let (_, stderr) = run_resumed_orientation_fixture(8, input, &["--orientation-block-context-weight", "0.01"]);
    assert!(stderr.contains("Signed block context:"), "{stderr}");
}

#[test]
fn joint_block_comparison_is_bounded_opt_in_and_checks_modes_before_inputs() {
    let matches = cli::cli().try_get_matches_from(["cphasing", "optimize", "group.txt", "group.clmb"]).unwrap();
    assert_eq!(matches.subcommand_matches("optimize").unwrap().get_one::<usize>("ORIENTATION_BLOCK_CANDIDATES"), Some(&0));
    for value in ["-1", "17", "1.5", "NaN"] {
        assert!(cli::cli().try_get_matches_from([
            "cphasing", "optimize", "group.txt", "group.clmb", "--orientation-block-candidates", value,
        ]).is_err());
    }
    for extra in [vec!["--orientation-method", "banded"], vec!["--orientation-method", "robust"], vec!["--orientation-block-span", "0"]] {
        let output = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
            .args(["optimize", "missing.txt", "missing.clmb", "--orientation-block-candidates", "4"])
            .args(extra).output().unwrap();
        assert_eq!(output.status.code(), Some(2));
        assert!(String::from_utf8(output.stderr).unwrap().contains("--orientation-block-candidates requires"));
    }
    let input = "ctg0+ ctg1+ ctg5- ctg4- ctg3- ctg2- ctg6+ ctg7+";
    let (expected, _) = run_resumed_orientation_fixture(8, input, &[]);
    for k in ["0", "1", "4", "16"] {
        let (tour, stderr) = run_resumed_orientation_fixture(8, input, &["--orientation-block-candidates", k]);
        assert_eq!(tour, expected);
        assert_eq!(stderr.contains("Signed block joint comparison:"), k != "0");
    }
    let (_, stderr) = run_resumed_orientation_fixture(8, input, &[
        "--orientation-block-candidates", "4", "--orientation-block-context-weight", "0.01",
    ]);
    assert!(stderr.contains("Signed block joint comparison:"));
    assert!(stderr.contains("Signed block context:"));
}

#[test]
fn historical_mode_rejects_unsupported_gates_before_opening_inputs() {
    for (option, value) in [
        ("--orientation-min-confidence", "0.2"),
        ("--orientation-max-flip-bp-fraction", "0.05"),
        ("--orientation-block-max-bp-fraction", "0.05"),
    ] {
        let output = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
            .args(["optimize", "missing-input.txt", "missing-input.clmb", option, value])
            .output().unwrap();
        assert_eq!(output.status.code(), Some(2));
        let stderr = String::from_utf8(output.stderr).unwrap();
        assert!(stderr.contains("use --orientation-method banded for confidence or bp limits"), "{stderr}");
    }
}

#[test]
fn experimental_endpoint_mode_has_a_separate_effect_default() {
    let matches = cli::cli().try_get_matches_from([
        "cphasing", "optimize", "group.txt", "group.clmb", "--orientation-method", "robust",
        "--orientation-audit", "audit.tsv",
    ]).unwrap();
    let args = matches.subcommand_matches("optimize").unwrap();
    assert_eq!(args.get_one::<f64>("ORIENTATION_MIN_CONFIDENCE"), Some(&0.2));
    assert_eq!(args.get_one::<usize>("ORIENTATION_BLOCK_SPAN"), Some(&0));
    assert!(!args.get_flag("ORIENTATION_TRUST_INPUT"));
    assert!(cli::cli().try_get_matches_from([
        "cphasing", "optimize", "group.txt", "group.clmb", "--orientation-method", "robust",
        "--orientation-trust-input",
    ]).is_err());
    let matches = cli::cli().try_get_matches_from([
        "cphasing", "optimize", "group.txt", "group.clmb", "--orientation-method", "robust",
        "--resume", "--orientation-trust-input", "--orientation-min-confidence", "0.4",
    ]).unwrap();
    let args = matches.subcommand_matches("optimize").unwrap();
    assert!(args.get_flag("ORIENTATION_TRUST_INPUT"));
    assert_eq!(args.get_one::<f64>("ORIENTATION_MIN_CONFIDENCE"), Some(&0.4));
}

#[test]
fn backbone_mode_is_enabled_by_default() {
    let matches = cli::cli()
        .try_get_matches_from(["cphasing", "optimize", "group.txt", "group.clmb"])
        .unwrap();
    let optimize = matches.subcommand_matches("optimize").unwrap();

    assert!(!optimize.get_flag("NO_BACKBONE"));
    assert_eq!(
        optimize
            .get_one::<String>("INITIALIZER")
            .map(String::as_str),
        Some("random")
    );
    assert_eq!(
        optimize
            .get_one::<String>("ORIENTATION_METHOD")
            .map(String::as_str),
        Some("banded-legacy")
    );
    assert_eq!(optimize.get_one::<usize>("ORIENTATION_WINDOW"), Some(&3));
    assert_eq!(
        optimize
            .get_one::<String>("ORIENTATION_PAIR_WEIGHT")
            .map(String::as_str),
        Some("sqrt-links")
    );
    assert_eq!(optimize.get_one::<usize>("ORIENTATION_MIN_LINKS"), Some(&3));
    assert_eq!(optimize.get_one::<f64>("ORIENTATION_PRIOR"), Some(&0.05));
    assert_eq!(
        optimize.get_one::<f64>("ORIENTATION_MIN_CONFIDENCE"),
        Some(&0.0)
    );
    assert_eq!(
        optimize.get_one::<f64>("ORIENTATION_MAX_FLIP_BP_FRACTION"),
        Some(&1.0)
    );
    assert_eq!(
        optimize.get_one::<usize>("ORIENTATION_BLOCK_SPAN"),
        Some(&32)
    );
    assert_eq!(
        optimize.get_one::<f64>("ORIENTATION_BLOCK_MAX_BP_FRACTION"),
        Some(&1.0)
    );
    assert_eq!(
        optimize.get_one::<f64>("ORIENTATION_BLOCK_MIN_GAIN"),
        Some(&0.0001)
    );
    assert_eq!(optimize.get_one::<usize>("ORIENTATION_BLOCK_PASSES"), Some(&4));
}

#[test]
fn conservative_banded_mode_retains_its_previous_defaults() {
    for method in ["banded", "banded-contact"] {
        let matches = cli::cli().try_get_matches_from([
            "cphasing", "optimize", "group.txt", "group.clmb", "--orientation-method", method,
        ]).unwrap();
        let args = matches.subcommand_matches("optimize").unwrap();
        assert_eq!(args.get_one::<f64>("ORIENTATION_MIN_CONFIDENCE"), Some(&0.95));
        assert_eq!(args.get_one::<f64>("ORIENTATION_MAX_FLIP_BP_FRACTION"), Some(&0.05));
        assert_eq!(args.get_one::<usize>("ORIENTATION_BLOCK_SPAN"), Some(&0));
        assert_eq!(args.get_one::<f64>("ORIENTATION_BLOCK_MAX_BP_FRACTION"), Some(&0.05));
        assert_eq!(args.get_one::<f64>("ORIENTATION_BLOCK_MIN_GAIN"), Some(&0.05));
    }
}

#[test]
fn robust_orientation_options_are_selectable_and_validated() {
    let matches = cli::cli()
        .try_get_matches_from([
            "cphasing",
            "optimize",
            "group.txt",
            "group.clmb",
            "--orientation-method",
            "banded",
            "--orientation-window",
            "4",
            "--orientation-pair-weight",
            "equal-pair",
            "--orientation-min-links",
            "7",
            "--orientation-prior",
            "0.1",
            "--orientation-min-confidence",
            "0.3",
            "--orientation-max-flip-bp-fraction",
            "0.02",
            "--orientation-block-span",
            "0",
            "--orientation-block-max-bp-fraction",
            "0.2",
            "--orientation-block-passes",
            "2",
            "--orientation-block-min-gain",
            "0.001",
        ])
        .unwrap();
    let optimize = matches.subcommand_matches("optimize").unwrap();

    assert_eq!(optimize.get_one::<usize>("ORIENTATION_WINDOW"), Some(&4));
    assert_eq!(
        optimize
            .get_one::<String>("ORIENTATION_PAIR_WEIGHT")
            .map(String::as_str),
        Some("equal-pair")
    );
    assert_eq!(optimize.get_one::<usize>("ORIENTATION_MIN_LINKS"), Some(&7));
    assert_eq!(optimize.get_one::<f64>("ORIENTATION_PRIOR"), Some(&0.1));
    assert_eq!(
        optimize.get_one::<f64>("ORIENTATION_MIN_CONFIDENCE"),
        Some(&0.3)
    );
    assert_eq!(
        optimize.get_one::<f64>("ORIENTATION_MAX_FLIP_BP_FRACTION"),
        Some(&0.02)
    );
    assert_eq!(
        optimize.get_one::<usize>("ORIENTATION_BLOCK_SPAN"),
        Some(&0)
    );
    assert_eq!(
        optimize.get_one::<f64>("ORIENTATION_BLOCK_MAX_BP_FRACTION"),
        Some(&0.2)
    );

    let invalid_window = cli::cli().try_get_matches_from([
        "cphasing",
        "optimize",
        "group.txt",
        "group.clmb",
        "--orientation-window",
        "17",
    ]);
    assert!(invalid_window.is_err());

    let invalid_block_span = cli::cli().try_get_matches_from([
        "cphasing",
        "optimize",
        "group.txt",
        "group.clmb",
        "--orientation-block-span",
        "1",
    ]);
    assert!(invalid_block_span.is_err());

    let invalid_block_bp_fraction = cli::cli().try_get_matches_from([
        "cphasing",
        "optimize",
        "group.txt",
        "group.clmb",
        "--orientation-block-max-bp-fraction",
        "1.1",
    ]);
    assert!(invalid_block_bp_fraction.is_err());

    let invalid_confidence = cli::cli().try_get_matches_from([
        "cphasing",
        "optimize",
        "group.txt",
        "group.clmb",
        "--orientation-min-confidence",
        "1.1",
    ]);
    assert!(invalid_confidence.is_err());

    let invalid_flip_bp_fraction = cli::cli().try_get_matches_from([
        "cphasing",
        "optimize",
        "group.txt",
        "group.clmb",
        "--orientation-max-flip-bp-fraction",
        "0",
    ]);
    assert!(invalid_flip_bp_fraction.is_err());
}

#[test]
fn seriation_initializer_is_selectable() {
    let matches = cli::cli()
        .try_get_matches_from([
            "cphasing",
            "optimize",
            "group.txt",
            "group.clmb",
            "--initializer",
            "seriation",
        ])
        .unwrap();
    let optimize = matches.subcommand_matches("optimize").unwrap();

    assert_eq!(
        optimize
            .get_one::<String>("INITIALIZER")
            .map(String::as_str),
        Some("seriation")
    );
}

#[test]
fn split_contacts_default_to_hierarchical_end_initialization() {
    let matches = cli::cli()
        .try_get_matches_from([
            "cphasing",
            "optimize",
            "group.txt",
            "group.clmb",
            "--split-contacts",
            "group.split.contacts.gz",
        ])
        .unwrap();
    let optimize = matches.subcommand_matches("optimize").unwrap();

    assert_eq!(
        optimize
            .get_one::<String>("INITIALIZER")
            .map(String::as_str),
        Some("end-hierarchical")
    );
}

#[test]
fn explicit_initializer_overrides_the_split_contacts_default() {
    let matches = cli::cli()
        .try_get_matches_from([
            "cphasing",
            "optimize",
            "group.txt",
            "group.clmb",
            "--split-contacts",
            "group.split.contacts.gz",
            "--initializer",
            "random",
        ])
        .unwrap();
    let optimize = matches.subcommand_matches("optimize").unwrap();

    assert_eq!(
        optimize
            .get_one::<String>("INITIALIZER")
            .map(String::as_str),
        Some("random")
    );
}

#[test]
fn end_tsp_initializer_accepts_split_contacts() {
    let matches = cli::cli()
        .try_get_matches_from([
            "cphasing",
            "optimize",
            "group.txt",
            "group.clmb",
            "--initializer",
            "end-tsp",
            "--split-contacts",
            "group.split.contacts.gz",
        ])
        .unwrap();
    let optimize = matches.subcommand_matches("optimize").unwrap();

    assert_eq!(
        optimize
            .get_one::<String>("INITIALIZER")
            .map(String::as_str),
        Some("end-tsp")
    );
    assert_eq!(
        optimize
            .get_one::<String>("SPLIT_CONTACTS")
            .map(String::as_str),
        Some("group.split.contacts.gz")
    );
}

#[test]
fn end_tsp_initializer_requires_split_contacts() {
    let result = cli::cli().try_get_matches_from([
        "cphasing",
        "optimize",
        "group.txt",
        "group.clmb",
        "--initializer",
        "end-tsp",
    ]);

    assert!(result.is_err());
}

#[test]
fn end_hierarchical_initializer_accepts_split_contacts() {
    let matches = cli::cli()
        .try_get_matches_from([
            "cphasing",
            "optimize",
            "group.txt",
            "group.clmb",
            "--initializer",
            "end-hierarchical",
            "--split-contacts",
            "group.split.contacts.gz",
            "--skipGA",
        ])
        .unwrap();
    let optimize = matches.subcommand_matches("optimize").unwrap();

    assert_eq!(
        optimize
            .get_one::<String>("INITIALIZER")
            .map(String::as_str),
        Some("end-hierarchical")
    );
    assert!(optimize.get_flag("SKIPGA"));
}

#[test]
fn end_beam_initializer_accepts_split_contacts() {
    let matches = cli::cli()
        .try_get_matches_from([
            "cphasing",
            "optimize",
            "group.txt",
            "group.clmb",
            "--initializer",
            "end-beam",
            "--split-contacts",
            "group.split.contacts.gz",
            "--skipGA",
        ])
        .unwrap();
    let optimize = matches.subcommand_matches("optimize").unwrap();

    assert_eq!(
        optimize
            .get_one::<String>("INITIALIZER")
            .map(String::as_str),
        Some("end-beam")
    );
    assert!(optimize.get_flag("SKIPGA"));
}

#[test]
fn end_beam_initializer_requires_split_contacts() {
    let result = cli::cli().try_get_matches_from([
        "cphasing",
        "optimize",
        "group.txt",
        "group.clmb",
        "--initializer",
        "end-beam",
    ]);

    assert!(result.is_err());
}

#[test]
fn no_backbone_disables_backbone_mode() {
    let matches = cli::cli()
        .try_get_matches_from([
            "cphasing",
            "optimize",
            "group.txt",
            "group.clmb",
            "--no-backbone",
        ])
        .unwrap();
    let optimize = matches.subcommand_matches("optimize").unwrap();

    assert!(optimize.get_flag("NO_BACKBONE"));
}

#[test]
fn length_tiered_objective_is_selectable_and_excludes_log_distance() {
    let matches = cli::cli()
        .try_get_matches_from([
            "cphasing",
            "optimize",
            "group.txt",
            "group.clmb",
            "--length-tiered",
        ])
        .unwrap();
    let optimize = matches.subcommand_matches("optimize").unwrap();
    assert!(optimize.get_flag("LENGTH_TIERED"));

    let conflict = cli::cli().try_get_matches_from([
        "cphasing",
        "optimize",
        "group.txt",
        "group.clmb",
        "--length-tiered",
        "--logDist",
    ]);
    assert!(conflict.is_err());
}

#[test]
fn endpoint_multiscale_requires_split_contacts() {
    let matches = cli::cli()
        .try_get_matches_from([
            "cphasing",
            "optimize",
            "group.txt",
            "group.clmb",
            "--split-contacts",
            "group.split.contacts.gz",
            "--endpoint-multiscale",
        ])
        .unwrap();
    let optimize = matches.subcommand_matches("optimize").unwrap();
    assert!(optimize.get_flag("ENDPOINT_MULTISCALE"));

    let missing = cli::cli().try_get_matches_from([
        "cphasing",
        "optimize",
        "group.txt",
        "group.clmb",
        "--endpoint-multiscale",
    ]);
    assert!(missing.is_err());
}
