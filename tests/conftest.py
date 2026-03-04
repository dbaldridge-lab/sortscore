
import os
import tempfile
import json
import pytest
import shutil
import copy

CONFIG_JSON_PATH = os.path.join(os.path.dirname(__file__), '../demo_data/GLI2_oPool5b/config.json')
BATCH_CONFIG_JSON_PATH = os.path.join(os.path.dirname(__file__), '../demo_data/batch_config.json')


@pytest.fixture(scope="function")
def config_dict(request):
    """
    Loads config from JSON file for tests. To override, use @pytest.mark.parametrize('config_dict', [your_dict], indirect=True)
    """
    with open(CONFIG_JSON_PATH, 'r') as f:
        config = json.load(f)
    if hasattr(request, "param"):
        config.update(request.param)
    return config



@pytest.fixture(scope="function")
def batch_config_dict(request):
    """
    Loads batch config from JSON file for tests. To override, use @pytest.mark.parametrize('batch_config_dict', [your_dict], indirect=True)
    """
    with open(BATCH_CONFIG_JSON_PATH, 'r') as f:
        config = json.load(f)
    if hasattr(request, "param"):
        config.update(request.param)
    return config

@pytest.fixture(scope="function")
def config_path(config_dict):
    # Write config dict to a temp file for CLI use
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as tf:
        json.dump(config_dict, tf)
        tf.flush()
        yield tf.name
    os.remove(tf.name)


@pytest.fixture(scope="function")
def batch_config_path(batch_config_dict):
    # Write batch config dict to a temp file for CLI/use in tests
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as tf:
        json.dump(batch_config_dict, tf)
        tf.flush()
        yield tf.name
    os.remove(tf.name)

@pytest.fixture(scope="function")
def cleanup_outputs(config_dict):
    data_dir = config_dict["output_dir"]
    yield
    import datetime
    # Cleanup function to remove only files and directories generated today
    if os.path.exists(data_dir):
        today_str = datetime.datetime.now().strftime('%Y%m%d')
        for filename in os.listdir(data_dir):
            if today_str in filename:
                file_path = os.path.join(data_dir, filename)
                try:
                    if os.path.isfile(file_path) or os.path.islink(file_path):
                        os.unlink(file_path)
                    elif os.path.isdir(file_path):
                        shutil.rmtree(file_path, ignore_errors=True)
                except Exception:
                    pass


@pytest.fixture(scope="function")
def isolated_runtime_env(tmp_path, monkeypatch):
    runtime_dir = tmp_path / "runtime"
    runtime_dir.mkdir()
    monkeypatch.setenv("MPLBACKEND", "Agg")
    monkeypatch.setenv("MPLCONFIGDIR", str(runtime_dir))
    monkeypatch.setenv("XDG_CACHE_HOME", str(runtime_dir))
    monkeypatch.setenv("HOME", str(runtime_dir))
    return os.environ.copy()


# Parameters to test individually for boundary cases
PARAMS_TO_TEST = [
    ("bins_required", 1),
    ("bins_required", 100),
    ("reps_required", 1),
    ("reps_required", 100),
    ("avg_method", "simple-avg"),
    ("avg_method", "rep-weighted"),
    ("minread_threshold", 0),
    ("minread_threshold", 10000),
    ("max_cv", 0.0),
    ("max_cv", 100.0),
    ("read_count", []),
    ("read_count", [10000] * 10),
    ("output_dir", "_test_outputs"),
    ("mutagenesis_variants", ["W", "F", "Y", "P", "M", "I", "L", "V", "A", "G", "C", "S", "T", "Q", "N", "D", "E", "H", "R", "K"])
]

@pytest.mark.parametrize("param,value", PARAMS_TO_TEST)
def test_config_single_param_boundary(config_dict, param, value):
    # Update one parameter at a time
    config = config_dict.copy()
    config[param] = value
    # Basic checks
    assert config[param] == value
    # TODO: #21 add more assertions for each param
