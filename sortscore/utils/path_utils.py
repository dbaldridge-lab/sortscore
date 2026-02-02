import os

def resolve_path(path, config_file_dir):
    """Return absolute path, resolving relative to config file if needed."""
    if os.path.isabs(path):
        return path
    return os.path.abspath(os.path.join(config_file_dir, path))