"""
Automated acceptance tests for Bounty #3: DryLab4 & SiLA 2 Robotic Bridge.
Run: pytest tests/test_bounty3_sila2.py -v
"""
import json
import pathlib
import subprocess
import time
import statistics
import pytest

MANIFEST = pathlib.Path("results/sila2_manifest.json")

GRPC_LATENCY_P99_MS_MAX = 50.0
PTP_UTC_OFFSET_MAX_MS = 1.0  # IEEE 1588 PTP vs GPS-disciplined UTC reference
CLOCK_JITTER_MAX_MS = 1.0
DRYLAB4_RT_ERROR_MAX_FRACTION = 0.02  # 2%
ICH_Q14_REQUIRED_FIELDS = {"actor", "timestamp", "delta", "operation_id"}


@pytest.fixture(scope="module")
def manifest():
    if not MANIFEST.exists():
        pytest.skip("results/sila2_manifest.json not found — run stack first")
    with open(MANIFEST) as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# Test 1: SiLA2 XML schema validation
# ---------------------------------------------------------------------------
def test_sila2_schema_validation(manifest):
    feature_xml = manifest.get("sila2_feature_descriptor_path")
    assert feature_xml, "sila2_feature_descriptor_path not set"
    result = subprocess.run(
        ["xmllint", "--noout", "--schema", "schemas/sila2_core_v1.0.0.xsd", feature_xml],
        capture_output=True, text=True
    )
    assert result.returncode == 0, (
        f"SiLA2 feature descriptor failed schema validation:\n{result.stderr}"
    )


# ---------------------------------------------------------------------------
# Test 2: gRPC roundtrip latency p99 < 50 ms
# ---------------------------------------------------------------------------
def test_grpc_latency(manifest):
    latencies = manifest["grpc_benchmark"]["latencies_ms"]
    assert len(latencies) >= 1000, f"Need ≥ 1000 latency samples, got {len(latencies)}"
    sorted_lat = sorted(latencies)
    p99_idx = int(len(sorted_lat) * 0.99)
    p99 = sorted_lat[p99_idx]
    assert p99 < GRPC_LATENCY_P99_MS_MAX, (
        f"gRPC p99 latency {p99:.2f} ms ≥ {GRPC_LATENCY_P99_MS_MAX} ms"
    )


# ---------------------------------------------------------------------------
# Test 3: DryLab4 retention-time prediction error < 2%
# ---------------------------------------------------------------------------
def test_drylab4_rt_prediction(manifest):
    predictions = manifest["drylab4"]["retention_time_predictions"]
    for entry in predictions:
        predicted = entry["predicted_min"]
        reference = entry["reference_min"]
        error = abs(predicted - reference) / reference
        assert error <= DRYLAB4_RT_ERROR_MAX_FRACTION, (
            f"RT prediction error {error:.2%} > 2% for compound '{entry.get('compound', '?')}'"
        )


# ---------------------------------------------------------------------------
# Test 4: IEEE 1588 PTP clock offset < 1 ms vs. GPS-disciplined UTC reference
# ---------------------------------------------------------------------------
def test_master_clock_jitter(manifest):
    clock_data = manifest["master_clock"]
    utc_offset_ms = clock_data.get("utc_offset_ms", clock_data.get("jitter_ms_60s_window"))
    assert utc_offset_ms < PTP_UTC_OFFSET_MAX_MS, (
        f"PTP UTC offset {utc_offset_ms:.3f} ms ≥ {PTP_UTC_OFFSET_MAX_MS} ms (IEEE 1588 requirement)"
    )


# ---------------------------------------------------------------------------
# Test 5: ICH Q14 audit trail completeness
# ---------------------------------------------------------------------------
def test_ich_q14_audit_trail(manifest):
    audit_log_path = manifest.get("ich_q14_audit_log")
    assert audit_log_path, "ich_q14_audit_log not set in manifest"
    with open(audit_log_path) as f:
        entries = [json.loads(line) for line in f if line.strip()]
    assert entries, "Audit log is empty"
    for i, entry in enumerate(entries):
        missing = ICH_Q14_REQUIRED_FIELDS - set(entry.keys())
        assert not missing, (
            f"Audit entry #{i} missing required ICH Q14 fields: {missing}"
        )


# ---------------------------------------------------------------------------
# Test 6: Docker Compose stack starts cleanly
# ---------------------------------------------------------------------------
def test_docker_compose_start():
    result = subprocess.run(
        ["docker", "compose", "up", "--wait", "--timeout", "120"],
        capture_output=True, text=True, timeout=150
    )
    assert result.returncode == 0, (
        f"docker compose up --wait failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
    )
    # Cleanup
    subprocess.run(["docker", "compose", "down"], capture_output=True)


# ---------------------------------------------------------------------------
# Test 7: End-to-end integration (mock instruments)
# ---------------------------------------------------------------------------
def test_e2e_mock_integration(manifest):
    e2e_cmd = manifest.get("e2e_test_command")
    if not e2e_cmd:
        pytest.skip("e2e_test_command not specified in manifest")
    result = subprocess.run(e2e_cmd, capture_output=True, text=True, timeout=300)
    assert result.returncode == 0, (
        f"E2E integration test failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
    )
