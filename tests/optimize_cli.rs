use cphasing::cli;

#[test]
fn backbone_mode_is_enabled_by_default() {
    let matches = cli::cli()
        .try_get_matches_from(["cphasing", "optimize", "group.txt", "group.clmb"])
        .unwrap();
    let optimize = matches.subcommand_matches("optimize").unwrap();

    assert!(!optimize.get_flag("NO_BACKBONE"));
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
