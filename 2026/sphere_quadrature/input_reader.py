
import configparser
import os


class ConfigError(Exception):
    pass


class Config:
    pass


_VALID_CHECKS = (
    "polar_polynomial_exactness",
    "spherical_harmonic_orthonormality",
    "orthonormality_beyond_resolution",
    "gauss_chebyshev_trapezoidal_equivalence",
    "quadrature_nodes_3d",
)

_VALID_AZIMUTHAL = ("trapezoidal", "gauss_chebyshev")
_VALID_COLOR_BY = ("weight", "function")
_VALID_TEST_FUNCTION = ("constant", "mu_power", "spherical_harmonic")


def _get_int(cp, section, key):
    try:
        return cp.getint(section, key)
    except (configparser.NoOptionError, configparser.NoSectionError):
        raise ConfigError(f"Missing value '{key}' in section [{section}] of input.txt")
    except ValueError:
        raw = cp.get(section, key)
        raise ConfigError(
            f"Value '{key}' in section [{section}] must be an integer, found: '{raw}'")


def _get_str(cp, section, key):
    try:
        return cp.get(section, key).strip()
    except (configparser.NoOptionError, configparser.NoSectionError):
        raise ConfigError(f"Missing value '{key}' in section [{section}] of input.txt")


def _get_str_optional(cp, section, key):
    """Like _get_str, but returns '' instead of raising if the section/key is missing."""
    try:
        return cp.get(section, key).strip()
    except (configparser.NoOptionError, configparser.NoSectionError):
        return ""


def _get_choice(cp, section, key, choices):
    val = _get_str(cp, section, key).lower()
    if val not in choices:
        raise ConfigError(
            f"'{key}' in section [{section}] must be one of {choices}, found: '{val}'")
    return val


def _require_positive_int(name, value):
    if value <= 0:
        raise ConfigError(f"'{name}' must be a positive integer")


def read_config(path="input.txt"):
    if not os.path.isfile(path):
        raise ConfigError(f"Input file not found: '{path}'")

    cp = configparser.ConfigParser(inline_comment_prefixes=("#",))
    try:
        with open(path, "r", encoding="utf-8") as f:
            cp.read_file(f)
    except configparser.Error as e:
        raise ConfigError(f"File '{path}' is not correctly formatted: {e}")

    cfg = Config()

    
    cfg.which = _get_choice(cp, "ACTION", "which", _VALID_CHECKS)

    
    if cfg.which == "polar_polynomial_exactness":
        section = "POLAR_POLYNOMIAL_EXACTNESS"
        cfg.n_mu = _get_int(cp, section, "n_mu")
        _require_positive_int("n_mu", cfg.n_mu)
        raw_k_max = _get_str_optional(cp, section, "k_max")
        if raw_k_max:
            try:
                cfg.k_max = int(raw_k_max)
            except ValueError:
                raise ConfigError(f"'k_max' in section [{section}] must be an integer, "
                                   f"found: '{raw_k_max}'")
        else:
            cfg.k_max = None  #

    elif cfg.which == "spherical_harmonic_orthonormality":
        section = "SPHERICAL_HARMONIC_ORTHONORMALITY"
        cfg.n_mu = _get_int(cp, section, "n_mu")
        cfg.n_phi = _get_int(cp, section, "n_phi")
        cfg.n_max = _get_int(cp, section, "n_max")
        cfg.azimuthal = _get_choice(cp, section, "azimuthal", _VALID_AZIMUTHAL)
        _require_positive_int("n_mu", cfg.n_mu)
        _require_positive_int("n_phi", cfg.n_phi)
        if cfg.n_max < 0:
            raise ConfigError("'n_max' must be >= 0")

    elif cfg.which == "orthonormality_beyond_resolution":
        section = "ORTHONORMALITY_BEYOND_RESOLUTION"
        cfg.n_mu = _get_int(cp, section, "n_mu")
        cfg.n_phi = _get_int(cp, section, "n_phi")
        cfg.azimuthal = _get_choice(cp, section, "azimuthal", _VALID_AZIMUTHAL)
        _require_positive_int("n_mu", cfg.n_mu)
        _require_positive_int("n_phi", cfg.n_phi)

    elif cfg.which == "gauss_chebyshev_trapezoidal_equivalence":
        section = "GAUSS_CHEBYSHEV_TRAPEZOIDAL_EQUIVALENCE"
        cfg.n_phi = _get_int(cp, section, "n_phi")
        _require_positive_int("n_phi", cfg.n_phi)

    else:  
        section = "QUADRATURE_NODES_3D"
        cfg.n_mu = _get_int(cp, section, "n_mu")
        cfg.n_phi = _get_int(cp, section, "n_phi")
        cfg.azimuthal = _get_choice(cp, section, "azimuthal", _VALID_AZIMUTHAL)
        _require_positive_int("n_mu", cfg.n_mu)
        _require_positive_int("n_phi", cfg.n_phi)

        cfg.color_by = _get_choice(cp, section, "color_by", _VALID_COLOR_BY)

        cfg.test_function = None
        cfg.k = None
        cfg.n = None
        cfg.m = None
        if cfg.color_by == "function":
            cfg.test_function = _get_choice(cp, section, "test_function", _VALID_TEST_FUNCTION)
            if cfg.test_function == "mu_power":
                cfg.k = _get_int(cp, section, "k")
            elif cfg.test_function == "spherical_harmonic":
                cfg.n = _get_int(cp, section, "n")
                cfg.m = _get_int(cp, section, "m")
                if abs(cfg.m) > cfg.n:
                    raise ConfigError(f"'m' must satisfy |m| <= n in section [{section}]")

    # --- OUTPUT ---
    cfg.out_dir = _get_str(cp, "OUTPUT", "output_folder")

    return cfg
