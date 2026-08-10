use cphasing::cli;

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
