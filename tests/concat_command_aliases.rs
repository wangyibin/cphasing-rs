use std::process::Command;

#[test]
fn concat_subcommands_keep_porec_aliases() {
    let commands: [(&[&str], &[&str], &str); 14] = [
        (
            &["paf2concat"],
            &["paf2porec"],
            "Usage: cphasing-rs paf2concat",
        ),
        (
            &["concat2depth"],
            &["porec2depth"],
            "Usage: cphasing-rs concat2depth",
        ),
        (
            &["concat2pqs"],
            &["porec2pqs"],
            "Usage: cphasing-rs concat2pqs",
        ),
        (
            &["concat-split"],
            &["porec-split"],
            "Usage: cphasing-rs concat-split",
        ),
        (
            &["concat2pairs"],
            &["porec2pairs"],
            "Usage: cphasing-rs concat2pairs",
        ),
        (
            &["concat-break"],
            &["porec-break"],
            "Usage: cphasing-rs concat-break",
        ),
        (
            &["concat-dup"],
            &["porec-dup"],
            "Usage: cphasing-rs concat-dup",
        ),
        (
            &["concat-merge"],
            &["porec-merge"],
            "Usage: cphasing-rs concat-merge",
        ),
        (
            &["concat2reads"],
            &["porec2reads"],
            "Usage: cphasing-rs concat2reads",
        ),
        (
            &["concat-intersect"],
            &["porec-intersect"],
            "Usage: cphasing-rs concat-intersect",
        ),
        (
            &["concat-downsample"],
            &["porec-downsample"],
            "Usage: cphasing-rs concat-downsample",
        ),
        (
            &["concat-chr2ctg"],
            &["porec-chr2ctg"],
            "Usage: cphasing-rs concat-chr2ctg",
        ),
        (
            &["concatbamstat"],
            &["porecbamstat"],
            "Usage: cphasing-rs concatbamstat",
        ),
        (
            &["simulator", "concat"],
            &["simulator", "porec"],
            "Usage: cphasing-rs simulator concat",
        ),
    ];

    for (canonical, legacy_alias, canonical_usage) in commands {
        for command in [canonical, legacy_alias] {
            let result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
                .args(command)
                .arg("--help")
                .output()
                .unwrap();
            assert!(
                result.status.success(),
                "{} --help failed: {}",
                command.join(" "),
                String::from_utf8_lossy(&result.stderr)
            );
            assert!(
                String::from_utf8_lossy(&result.stdout).contains(canonical_usage),
                "{} --help did not use canonical command name",
                command.join(" ")
            );
        }
    }
}

#[test]
fn root_help_groups_visible_commands() {
    let result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("--help")
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let help = String::from_utf8_lossy(&result.stdout);
    assert!(
        !help.contains("\nCommands:\n"),
        "the default ungrouped command list is still present"
    );
    assert!(
        !help.contains("OTHER COMMANDS:"),
        "every visible command must belong to a named group"
    );
    assert!(
        help.contains("\nOPTIONS:\n  -h, --help"),
        "the native help options are retained after the grouped command list"
    );

    let groups = [
        (
            "PAF:",
            &[
                "methalign",
                "paf2concat",
                "paf2depth",
                "paf2pairs",
                "paf-downsample",
            ][..],
        ),
        (
            "CONCATEMER:",
            &[
                "concat2pqs",
                "concat2pairs",
                "concat2depth",
                "concat2reads",
                "concat-split",
                "concat-merge",
                "concat-intersect",
                "concat-downsample",
                "concat-break",
                "concat-dup",
                "concat-chr2ctg",
            ][..],
        ),
        (
            "PAIRS:",
            &[
                "pairs-downsample",
                "pairs-merge",
                "pairs-filter",
                "pairs-break",
                "pairs-dup",
                "pairs-split",
                "pairs2contacts",
                "pairs2clm",
                "pairs2depth",
                "pairs2mnd",
                "pairs2bam",
                "pairs2porec",
                "pairs-chr2ctg",
                "pairs-prune",
                "pairs-intersect",
                "splitcontacts",
                "cool2mcool",
            ][..],
        ),
        (
            "BAM:",
            &[
                "bam2paf",
                "bam2pairs",
                "bam-chr2ctg",
                "bam-prune",
                "bam2fastq",
                "bam2fasta",
                "modbam2fq",
                "bamstat",
                "splitbam",
            ][..],
        ),
        (
            "ASSEMBLY & GRAPH:",
            &[
                "kprune",
                "clm",
                "mergeclm",
                "splitclm",
                "gfa",
                "optimize",
                "phase-reads",
            ][..],
        ),
        (
            "SEQUENCE, RESTRICTION ENZYME & SIMULATION:",
            &[
                "splitfastq",
                "extract-fasta",
                "slidefastq",
                "slidefasta",
                "slide2raw",
                "digest",
                "count_re",
                "cutsite",
                "chromsizes",
                "modfa",
                "simulator",
            ][..],
        ),
    ];
    let positions = groups
        .iter()
        .map(|(heading, _)| help.find(heading).expect("group heading is listed"))
        .collect::<Vec<_>>();
    assert!(positions.windows(2).all(|pair| pair[0] < pair[1]));

    for ((heading, commands), start) in groups.iter().zip(positions.iter()) {
        let end = positions
            .iter()
            .copied()
            .find(|position| position > start)
            .unwrap_or(help.len());
        let section = &help[*start..end];
        for command in *commands {
            assert!(
                section.contains(command),
                "{command} is not listed under {heading}"
            );
        }
    }
}
