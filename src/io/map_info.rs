//! map_info.json writer — summary statistics for a mapping run.

use std::path::Path;

use anyhow::Result;

use crate::index::reference_index::RefSigInfo;

/// All parameters needed to write a map_info.json file.
pub struct MapInfoParams<'a> {
    pub path: &'a Path,
    /// Mapping mode: "bulk", "sc-rna", or "sc-atac".
    pub mode: &'a str,
    pub num_reads: u64,
    pub num_mapped: u64,
    pub num_poisoned: u64,
    pub elapsed_secs: f64,
    /// Mapping pipeline only, excluding index load and final file backpatching.
    pub mapping_elapsed_secs: f64,
    pub sig_info: Option<&'a RefSigInfo>,
    pub piscem_rs_version: &'a str,
    pub num_threads: usize,
    /// The effective mapping/decode allocation after resource capping and
    /// decoder-handle reconciliation.
    pub execution_plan: Option<&'a crate::io::fastx::ExecutionPlan>,
    pub broker_report: Option<&'a thread_broker::BrokerReport>,
    pub broker_failure: Option<&'a crate::io::broker::BrokerFailure>,
    pub producer_measurement: Option<&'a thread_broker::ProducerMeasurementStats>,
    pub consumer_measurement: Option<&'a crate::io::threads::ConsumerMeasurement>,
    pub pipeline_tuning: Option<&'a crate::io::fastx::PipelineTuning>,
    pub index_path: &'a Path,
    pub k: usize,
    pub m: usize,
    pub num_refs: usize,
    /// Skipping strategy: "permissive", "strict", or "every-kmer".
    pub skipping_strategy: &'a str,
}

/// Write a map_info.json file with mapping statistics and run metadata.
pub fn write_map_info(params: &MapInfoParams) -> Result<()> {
    let percent_mapped = if params.num_reads > 0 {
        (params.num_mapped as f64 / params.num_reads as f64) * 100.0
    } else {
        0.0
    };

    let cmdline: Vec<String> = std::env::args().collect();

    let mut info = serde_json::json!({
        "mode": params.mode,
        "piscem_rs_version": params.piscem_rs_version,
        "index_path": params.index_path.display().to_string(),
        "k": params.k,
        "m": params.m,
        "num_refs": params.num_refs,
        "num_threads": params.num_threads,
        "skipping_strategy": params.skipping_strategy,
        "num_reads": params.num_reads,
        "num_mapped": params.num_mapped,
        "num_poisoned": params.num_poisoned,
        "percent_mapped": format!("{:.2}", percent_mapped),
        "runtime_seconds": format!("{:.2}", params.elapsed_secs),
        "mapping_seconds": params.mapping_elapsed_secs,
        "cmdline": cmdline.join(" "),
    });

    if let Some(sig) = params.sig_info {
        info["signatures"] = serde_json::json!({
            "sha256_names": sig.sha256_names,
            "sha256_seqs":  sig.sha256_seqs,
            "sha512_names": sig.sha512_names,
            "sha512_seqs":  sig.sha512_seqs,
        });
    }

    if let Some(plan) = params.execution_plan {
        info["execution_plan"] = serde_json::to_value(plan)?;
    }
    if let Some(report) = params.broker_report {
        info["thread_broker"] = serde_json::to_value(report)?;
    }
    if let Some(failure) = params.broker_failure {
        info["thread_broker_failure"] = serde_json::to_value(failure)?;
    }
    if let Some(measurement) = params.producer_measurement {
        info["producer_measurement"] = serde_json::to_value(measurement)?;
    }
    if let Some(measurement) = params.consumer_measurement {
        info["consumer_measurement"] = serde_json::to_value(measurement)?;
    }
    if let Some(tuning) = params.pipeline_tuning {
        info["pipeline_tuning"] = serde_json::to_value(tuning)?;
    }

    let file = std::fs::File::create(params.path)?;
    let writer = std::io::BufWriter::new(file);
    serde_json::to_writer_pretty(writer, &info)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_map_info() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("map_info.json");
        let index_path = std::path::PathBuf::from("/tmp/test_index");

        write_map_info(&MapInfoParams {
            path: &path,
            mode: "bulk",
            num_reads: 1000,
            num_mapped: 800,
            num_poisoned: 10,
            elapsed_secs: 42.5,
            mapping_elapsed_secs: 40.25,
            sig_info: None,
            piscem_rs_version: "0.1.0",
            num_threads: 8,
            execution_plan: None,
            broker_report: None,
            broker_failure: None,
            producer_measurement: None,
            consumer_measurement: None,
            pipeline_tuning: None,
            index_path: &index_path,
            k: 31,
            m: 19,
            num_refs: 500,
            skipping_strategy: "permissive",
        })
        .unwrap();

        let content = std::fs::read_to_string(&path).unwrap();
        let val: serde_json::Value = serde_json::from_str(&content).unwrap();
        assert_eq!(val["num_reads"], 1000);
        assert_eq!(val["num_mapped"], 800);
        assert_eq!(val["num_poisoned"], 10);
        assert_eq!(val["percent_mapped"], "80.00");
        assert_eq!(val["mode"], "bulk");
        assert_eq!(val["k"], 31);
        assert_eq!(val["m"], 19);
        assert_eq!(val["num_refs"], 500);
        assert_eq!(val["num_threads"], 8);
        assert_eq!(val["skipping_strategy"], "permissive");
        assert_eq!(val["piscem_rs_version"], "0.1.0");
        assert_eq!(val["mapping_seconds"], 40.25);
    }

    #[test]
    fn writes_the_effective_execution_plan() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("map_info.json");
        let index_path = std::path::PathBuf::from("/tmp/test_index");
        let plan = crate::io::fastx::ExecutionPlan {
            requested_budget: 32,
            effective_budget: 8,
            map_threads: 6,
            decode_slots: 2,
            allocation: crate::io::fastx::DecodeAllocation::Adaptive,
            requested_decode_slots: None,
        };

        write_map_info(&MapInfoParams {
            path: &path,
            mode: "bulk",
            num_reads: 1,
            num_mapped: 1,
            num_poisoned: 0,
            elapsed_secs: 1.0,
            mapping_elapsed_secs: 0.8,
            sig_info: None,
            piscem_rs_version: "0.1.0",
            num_threads: plan.effective_budget,
            execution_plan: Some(&plan),
            broker_report: None,
            broker_failure: None,
            producer_measurement: None,
            consumer_measurement: None,
            pipeline_tuning: None,
            index_path: &index_path,
            k: 31,
            m: 19,
            num_refs: 1,
            skipping_strategy: "permissive",
        })
        .unwrap();

        let value: serde_json::Value =
            serde_json::from_slice(&std::fs::read(path).unwrap()).unwrap();
        assert_eq!(value["num_threads"], 8);
        assert_eq!(value["execution_plan"]["requested_budget"], 32);
        assert_eq!(value["execution_plan"]["effective_budget"], 8);
        assert_eq!(value["execution_plan"]["map_threads"], 6);
        assert_eq!(value["execution_plan"]["decode_slots"], 2);
        assert_eq!(value["execution_plan"]["allocation"]["kind"], "adaptive");
    }

    #[test]
    fn writes_broker_and_measurement_diagnostics() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("map_info.json");
        let index_path = std::path::PathBuf::from("/tmp/test_index");
        let measurement = crate::io::threads::ConsumerMeasurement {
            busy_nanos: 11,
            callback_setup_nanos: 3,
            output_flush_nanos: 2,
        };
        let report = thread_broker::BrokerReport {
            final_consumer_threads: 6,
            final_producer_limit: 2,
            final_phase: thread_broker::BrokerPhase::Steady,
            final_model: Some(thread_broker::Model {
                producer_cost_share: 0.25,
                producer_cost_share_uncertainty: 0.01,
                ideal_producer_slots: 2,
                useful_cap: 3,
                useful_cap_reason: thread_broker::ProducerCapReason::Slack,
                effective_deadband_threads: 1,
                effective_resurvey_distance: 3,
            }),
            ..thread_broker::BrokerReport::default()
        };
        let producer_measurement = thread_broker::ProducerMeasurementStats {
            busy_nanos: 17,
            accounted_busy_nanos: Some(19),
            completed_worker_cpu_nanos: Some(13),
            completed_auxiliary_cpu_nanos: Some(4),
            cpu_accounting_failures: Some(0),
            sampler_cpu_nanos: Some(3),
            sampler_cpu_accounting_failures: 0,
            calibration_samples: 5,
            monitoring_samples: 2,
            mode_changes: 1,
            final_mode: thread_broker::ProducerMeasurementMode::Monitoring,
            observation_nanos: 7,
            calibration_interval_micros: 2_000,
            monitoring_interval_micros: 25_000,
        };
        let broker_failure = crate::io::broker::BrokerFailure {
            stage: crate::io::broker::BrokerFailureStage::ControllerRuntime,
            message: "injected resize timeout".to_string(),
        };
        let tuning = crate::io::fastx::PipelineTuning {
            reader_batch_size: 1024,
            progress_flush_every: 64,
        };

        write_map_info(&MapInfoParams {
            path: &path,
            mode: "bulk",
            num_reads: 1,
            num_mapped: 1,
            num_poisoned: 0,
            elapsed_secs: 1.0,
            mapping_elapsed_secs: 0.8,
            sig_info: None,
            piscem_rs_version: "0.1.0",
            num_threads: 8,
            execution_plan: None,
            broker_report: Some(&report),
            broker_failure: Some(&broker_failure),
            producer_measurement: Some(&producer_measurement),
            consumer_measurement: Some(&measurement),
            pipeline_tuning: Some(&tuning),
            index_path: &index_path,
            k: 31,
            m: 19,
            num_refs: 1,
            skipping_strategy: "permissive",
        })
        .unwrap();

        let value: serde_json::Value =
            serde_json::from_slice(&std::fs::read(path).unwrap()).unwrap();
        assert_eq!(value["thread_broker"]["final_phase"], "steady");
        assert_eq!(
            value["thread_broker_failure"]["stage"],
            "controller_runtime"
        );
        assert_eq!(
            value["thread_broker_failure"]["message"],
            "injected resize timeout"
        );
        assert_eq!(
            value["thread_broker"]["final_model"]["useful_cap_reason"],
            "slack"
        );
        assert_eq!(value["consumer_measurement"]["busy_nanos"], 11);
        assert_eq!(value["producer_measurement"]["final_mode"], "monitoring");
        assert_eq!(value["producer_measurement"]["calibration_samples"], 5);
        assert_eq!(value["producer_measurement"]["accounted_busy_nanos"], 19);
        assert_eq!(value["consumer_measurement"]["callback_setup_nanos"], 3);
        assert_eq!(value["consumer_measurement"]["output_flush_nanos"], 2);
        assert_eq!(value["pipeline_tuning"]["reader_batch_size"], 1024);
        assert_eq!(value["pipeline_tuning"]["progress_flush_every"], 64);
    }
}
