import configparser
import os


class ConfigError(Exception):
    """Error while reading or interpreting the input file."""
    pass


class Config:
    """Simple container for all data read from input.txt."""
    pass


_MISSING = object()  # sentinel: "no default given" (required key)

_GETTERS = {
    float: ("getfloat", "a number"),
    int: ("getint", "an integer"),
    bool: ("getboolean", "yes/no"),
}


def _get(cp, section, key, kind=str, default=_MISSING):
    if kind is str:
        try:
            return cp.get(section, key).strip()
        except (configparser.NoOptionError, configparser.NoSectionError):
            if default is not _MISSING:
                return default
            raise ConfigError(f"Missing value '{key}' in section [{section}] of input.txt")

    method = getattr(cp, _GETTERS[kind][0])
    try:
        return method(section, key)
    except (configparser.NoOptionError, configparser.NoSectionError):
        if default is not _MISSING:
            return default
        raise ConfigError(f"Missing value '{key}' in section [{section}] of input.txt")
    except ValueError:
        raw = cp.get(section, key)
        raise ConfigError(
            f"Value '{key}' in section [{section}] must be {_GETTERS[kind][1]}, found: '{raw}'")


def _get_choice(cp, section, key, choices, default=_MISSING):
    val = _get(cp, section, key, kind=str, default=default).lower()
    if val not in choices:
        raise ConfigError(
            f"'{key}' in section [{section}] must be one of {choices}, found: '{val}'")
    return val


def _parse_int_list(raw, section, key):
    items = [x.strip() for x in raw.split(",") if x.strip()]
    try:
        return [int(x) for x in items]
    except ValueError:
        raise ConfigError(
            f"'{key}' in section [{section}] must be a comma-separated list of integers, "
            f"found: '{raw}'")


def _parse_bc_pairs(raw, section, key):
    valid = ("void", "reflect")
    pairs = []
    for chunk in raw.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        parts = chunk.split("-")
        if len(parts) != 2:
            raise ConfigError(
                f"'{key}' in section [{section}]: each combination must have the form "
                f"'left-right' (e.g. 'void-reflect'), found: '{chunk}'")
        left, right = parts[0].strip(), parts[1].strip()
        if left not in valid or right not in valid:
            raise ConfigError(
                f"'{key}' in section [{section}]: allowed values are {valid}, "
                f"found: '{chunk}'")
        pairs.append((left, right))
    if not pairs:
        raise ConfigError(f"'{key}' in section [{section}] is empty: at least one combination is required")
    return pairs


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

    # geometry
    cfg.d = _get(cp, "GEOMETRY", "thickness_cm", kind=float)
    if cfg.d <= 0:
        raise ConfigError("'thickness_cm' must be positive")
    cfg.I = _get(cp, "GEOMETRY", "number_of_cells", kind=int)
    if cfg.I <= 0:
        raise ConfigError("'number_of_cells' must be a positive integer")

    # material
    cfg.material_mode = _get_choice(cp, "MATERIAL", "mode", ("homogeneous", "heterogeneous"))
    cfg.sigma_t = _get(cp, "MATERIAL", "sigma_t_cm-1", kind=float)
    if cfg.sigma_t <= 0:
        raise ConfigError("'sigma_t_cm-1' must be positive")
    if cfg.material_mode == "homogeneous":
        cfg.c = _get(cp, "MATERIAL", "c", kind=float)
        if not (0.0 <= cfg.c < 1.0):
            raise ConfigError("'c' must be between 0 (included) and 1 (excluded)")
        cfg.c_left = cfg.c_right = cfg.c
    else:
        cfg.c_left = _get(cp, "MATERIAL", "c_left", kind=float)
        cfg.c_right = _get(cp, "MATERIAL", "c_right", kind=float)
        for name, val in [("c_left", cfg.c_left), ("c_right", cfg.c_right)]:
            if not (0.0 <= val < 1.0):
                raise ConfigError(f"'{name}' must be between 0 (included) and 1 (excluded)")
        cfg.c = None

    # linear-anisotropic (P1) extension
    cfg.g = _get(cp, "MATERIAL", "g", kind=float, default=0.0)
    if not (-1.0 / 3.0 <= cfg.g <= 1.0 / 3.0):
        raise ConfigError(
            "'g' in section [MATERIAL] must satisfy -1/3 <= g <= 1/3 "
            "(otherwise the scattering phase function p(mu0)=1/2(1+3 g mu0) "
            "would be negative for some mu0)")

    # source
    cfg.source_type = _get_choice(cp, "SOURCE", "type", ("uniform", "peaked"))
    cfg.sext = _get(cp, "SOURCE", "intensity", kind=float)

    # linear-anisotropic extension
    cfg.source_angular_shape = _get_choice(
        cp, "SOURCE", "angular_shape", ("isotropic", "forward_peaked"), default="isotropic")
    if cfg.source_type == "peaked":
        cfg.peak_width = _get(cp, "SOURCE", "peaked_width_cm", kind=float)
        if not (0 < cfg.peak_width <= cfg.d):
            raise ConfigError("'peaked_width_cm' must be between 0 and the slab thickness")
    else:
        cfg.peak_width = None

    # boundary conditions
    cfg.bc_left = _get_choice(cp, "BOUNDARY_CONDITIONS", "left", ("void", "reflect"))
    cfg.bc_right = _get_choice(cp, "BOUNDARY_CONDITIONS", "right", ("void", "reflect"))

    # angular quadrature
    cfg.N_quad = _get(cp, "ANGULAR_QUADRATURE", "order_N", kind=int)
    if cfg.N_quad <= 0 or cfg.N_quad % 2 != 0:
        raise ConfigError("'order_N' must be a positive even integer")

    # solver
    cfg.tol = _get(cp, "SOLVER", "tolerance", kind=float)
    cfg.max_iter = _get(cp, "SOLVER", "max_iterations", kind=int)

    # actions
    cfg.do_solve_and_plot = _get(cp, "ACTIONS", "solve_and_plot", kind=bool)
    cfg.do_compare_analytical = _get(cp, "ACTIONS", "compare_with_analytical", kind=bool)
    cfg.do_grid_convergence = _get(cp, "ACTIONS", "grid_convergence_study", kind=bool)
    cfg.do_iterations_vs_c = _get(cp, "ACTIONS", "iterations_vs_c_study", kind=bool)
    cfg.compare_bc = _get(cp, "ACTIONS", "compare_boundary_conditions", kind=bool)

    # linear-anisotropic extension
    cfg.do_plot_anisotropic = _get(cp, "ACTIONS", "plot_anisotropic_functions", kind=bool, default=False)

    if not any([cfg.do_solve_and_plot, cfg.do_compare_analytical,
                cfg.do_grid_convergence, cfg.do_iterations_vs_c,
                cfg.do_plot_anisotropic]):
        raise ConfigError(
            "All actions in [ACTIONS] are set to 'no': nothing to do. "
            "Set at least one of them to 'yes'.")

    # compare_boundary_conditions
    if cfg.compare_bc:
        if cfg.do_compare_analytical:
            raise ConfigError(
                "'compare_boundary_conditions = yes' cannot be combined with "
                "'compare_with_analytical = yes' (the analytical check assumes a single, "
                "fixed reflect/reflect configuration). Turn one of the two off.")
        if not (cfg.do_solve_and_plot or cfg.do_iterations_vs_c):
            raise ConfigError(
                "'compare_boundary_conditions = yes' has no effect unless 'solve_and_plot' "
                "and/or 'iterations_vs_c_study' is also set to 'yes'.")
        raw_bc_list = _get(cp, "BOUNDARY_CONDITIONS", "compare_list")
        cfg.bc_compare_list = _parse_bc_pairs(raw_bc_list, "BOUNDARY_CONDITIONS", "compare_list")

    # compare_with_analytical requires an infinite-homogeneous-medium setup
    if cfg.do_compare_analytical:
        if cfg.material_mode != "homogeneous":
            raise ConfigError(
                "'compare_with_analytical = yes' requires [MATERIAL] mode = homogeneous "
                "(the analytical formula assumes a single c everywhere)")
        if cfg.bc_left != "reflect" or cfg.bc_right != "reflect":
            raise ConfigError(
                "'compare_with_analytical = yes' requires both boundary conditions to be "
                "'reflect' (reflect/reflect approximates an infinite homogeneous medium)")
        if cfg.g != 0.0 or cfg.source_angular_shape != "isotropic":
            raise ConfigError(
                "'compare_with_analytical = yes' requires isotropic scattering and source "
                "(the analytical formula phi0 = s_ext/(sigma_t*(1-c)) assumes both); "
                "set [MATERIAL] g = 0.0 and [SOURCE] angular_shape = isotropic, "
                "or turn this action off")
        # optional convergence table across several mesh sizes
        raw_extra_I = _get(cp, "ANALYTICAL_COMPARISON", "extra_cell_counts", default="")
        if raw_extra_I:
            cfg.analytical_cell_counts = _parse_int_list(
                raw_extra_I, "ANALYTICAL_COMPARISON", "extra_cell_counts")
        else:
            cfg.analytical_cell_counts = None

    # grid convergence study
    if cfg.do_grid_convergence:
        raw_Ilist = _get(cp, "GRID_CONVERGENCE_STUDY", "cell_counts")
        cfg.gc_cell_counts = _parse_int_list(raw_Ilist, "GRID_CONVERGENCE_STUDY", "cell_counts")
        raw_Iplot = _get(cp, "GRID_CONVERGENCE_STUDY", "cell_counts_for_plot")
        cfg.gc_cell_counts_plot = _parse_int_list(
            raw_Iplot, "GRID_CONVERGENCE_STUDY", "cell_counts_for_plot")

    # iterations-vs-c study
    if cfg.do_iterations_vs_c:
        cfg.ic_c_min = _get(cp, "ITERATIONS_VS_C_STUDY", "c_min", kind=float)
        cfg.ic_c_threshold = _get(cp, "ITERATIONS_VS_C_STUDY", "c_threshold", kind=float)
        cfg.ic_c_max = _get(cp, "ITERATIONS_VS_C_STUDY", "c_max", kind=float)
        cfg.ic_step_coarse = _get(cp, "ITERATIONS_VS_C_STUDY", "coarse_step", kind=float)
        cfg.ic_step_fine = _get(cp, "ITERATIONS_VS_C_STUDY", "fine_step", kind=float)
        cfg.ic_max_iter = _get(cp, "ITERATIONS_VS_C_STUDY", "max_iterations", kind=int)

    # output
    cfg.out_dir = _get(cp, "OUTPUT", "output_folder")

    return cfg