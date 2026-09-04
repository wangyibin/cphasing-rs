use std::process::Command;

fn run_cool2mcool(arguments: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("cool2mcool")
        .args(arguments)
        .output()
        .expect("run cphasing-rs cool2mcool")
}

#[test]
fn cool2mcool_help_exposes_native_generation_options() {
    let output = run_cool2mcool(&["--help"]);

    assert!(output.status.success());
    let help = String::from_utf8_lossy(&output.stdout);
    assert!(help.contains("fixed 1 kb COOL input"));
    assert!(help.contains("--level-parallelism <N>"));
    assert!(help.contains("--aggregation-mode <MODE>"));
    assert!(help.contains("--compression-level <LEVEL>"));
    assert!(help.contains("--resolutions <BP[,BP...]>"));
    assert!(help.contains("--kr-min-resolution <BP>"));
}

#[test]
fn cool2mcool_rejects_invalid_options_before_reading_input() {
    let threads = run_cool2mcool(&["--threads", "0", "input.cool", "output.mcool"]);
    assert!(!threads.status.success());
    assert!(String::from_utf8_lossy(&threads.stderr).contains("value must be greater than 0"));

    let duplicate_resolutions =
        run_cool2mcool(&["--resolutions", "1000,1000", "input.cool", "output.mcool"]);
    assert!(!duplicate_resolutions.status.success());
    assert!(
        String::from_utf8_lossy(&duplicate_resolutions.stderr)
            .contains("--resolutions cannot contain duplicate values")
    );
}
