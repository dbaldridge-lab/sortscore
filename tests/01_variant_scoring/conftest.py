import os
import tempfile
import json
import pytest
import shutil

DEFAULT_CONFIG_DICT = {
    "bins_required": 3,
    "reps_required": 1,
    "avg_method": "rep-weighted",
    "minread_threshold": 0,
    "wt_seq": "ATCCCTGGCTGCACCAAGAGATACACCGACCCTAGCAGCCTGAGGAAGCACGTGAAGACCGTGCACGGCCCTGACGCCCACGTGACCAAGAAGCAGAGG",
    "min_pos": 551,
    "max_pos": 583,
    "output_dir": '/Users/c.chitwood/code/sortscore/tests/scratch/_test_output',
    "biophysical_prop": True
}

@pytest.fixture(scope="function")
def config_dict(request):
    """
    Default config dict fixture, parameterizable via pytest.mark.parametrize.
    To override, use @pytest.mark.parametrize('config_dict', [your_dict], indirect=True)
    """
    config = DEFAULT_CONFIG_DICT.copy()
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
    ("output_dir", "tests/scratch/_test_output"),
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
