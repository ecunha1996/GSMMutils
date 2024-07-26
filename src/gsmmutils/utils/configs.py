import configparser
import importlib
from os import getenv
from os.path import join, normpath, realpath
import importlib.resources

config = None


def get_config():
    """
    Function to get base configurations
    Returns
    -------
    config: configparser.ConfigParser
    """
    global config

    if config is not None:
        return config

    config = read_config()
    return config


def get_defaults():
    """
    Function to get default configurations
    Returns
    -------
    cfg: configparser.ConfigParser
    """
    default_package_root = join(str(importlib.resources.files("gsmmutils").parent), "../")
    normalized_path = normpath(default_package_root)
    resolved_path = realpath(normalized_path)
    GSMMUTILS_ROOT = getenv("GSMMUTILS_ROOT", resolved_path)
    cfg = configparser.ConfigParser()
    cfg.add_section("PATHS")
    cfg.set("PATHS", "GSMMUTILS_ROOT", GSMMUTILS_ROOT)
    cfg.set("PATHS", "DATA_PATH", join(GSMMUTILS_ROOT, 'data'))
    cfg.set("PATHS", "SRC_PATH", join(GSMMUTILS_ROOT, 'src'))
    cfg.set("PATHS", "CONFIG_PATH", join(GSMMUTILS_ROOT, "config"))
    return cfg


def read_config():
    """
    Function to read configurations from file "configurations.ini"
    Returns
    -------
    cfg: configparser.ConfigParser
    """
    cfg = get_defaults()
    cfg.read(join(cfg["PATHS"]["CONFIG_PATH"], 'configurations.ini'))
    #update_config(cfg)
    return cfg


def update_config(cfg):
    """
    Function to write the configurations in the "configurations.ini" file
    Parameters
    ----------
    cfg: configparser.ConfigParser

    Returns
    -------
        None
    """
    try:
        with open(join(cfg["PATHS"]["CONFIG_PATH"], 'configurations.ini'), "w") as f:
            cfg.write(f)
            f.flush()
            f.close()
    except Exception as e:
        pass
