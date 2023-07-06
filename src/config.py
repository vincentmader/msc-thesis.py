import os
import toml

from constants import M_sun, L_sun, AU

dirname = os.path.dirname(__file__)
path_to_config_toml = os.path.join(dirname, "..", "config.toml")

PATH_TO_DARKMODE = os.path.join(
    dirname, "visualization", "mpl-styles", "dark.mplstyle"
)
PATH_TO_FIGURES = os.path.join(dirname, "..", "figures")


class Config():

    def __init__(
        self,
        stellar_mass=None,
        stellar_luminosity=None,
        disk_mass_ratio=None,
        dust_to_gas_ratio=None,
        distance_to_star=None,
        flaring_angle=None,
        radial_min_value=None,
        radial_max_value=None,
        radial_resolution=None,
        radial_axis_scale=None,
        mass_min_value=None,
        mass_max_value=None,
        mass_resolution=None,
        mass_axis_scale=None,
        time_min_value=None,
        time_max_value=None,
        time_resolution=None,
        time_axis_scale=None,
        enable_coagulation=None,
        enable_fragmentation=None,
        enable_physical_gas_density=None,
        enable_physical_cross_sections=None,
        enable_physical_relative_velocities=None,
        enable_cancellation_handling=None,
        fragmentation_variants=None,
        collision_outcome_variant=None,
        solver_variant=None,
        mpl_dark_mode=None,
    ):
        cfg = toml.load(path_to_config_toml)

        cfg_i = cfg["disk"]
        if stellar_mass is None:
            stellar_mass = M_sun * cfg_i["stellar_mass"]
        if stellar_luminosity is None:
            stellar_luminosity = L_sun * cfg_i["stellar_luminosity"]
        if disk_mass_ratio is None:
            disk_mass_ratio = cfg_i["disk_mass_ratio"]
        disk_mass = stellar_mass * disk_mass_ratio
        if dust_to_gas_ratio is None:
            dust_to_gas_ratio = cfg_i["dust_to_gas_ratio"]
        if distance_to_star is None:
            distance_to_star = AU * cfg_i["distance_to_star"]
        if flaring_angle is None:
            flaring_angle = cfg_i["flaring_angle"]
        if enable_physical_gas_density is None:
            enable_physical_gas_density = cfg_i["enable_physical_gas_density"]

        cfg_i = cfg["radial_discretization"]
        if radial_min_value is None:
            radial_min_value = AU * cfg_i["radial_min_value"]
        if radial_max_value is None:
            radial_max_value = AU * cfg_i["radial_max_value"]
        if radial_resolution is None:
            radial_resolution = cfg_i["radial_resolution"]
        if radial_axis_scale is None:
            radial_axis_scale = cfg_i["radial_axis_scale"]

        cfg_i = cfg["mass_discretization"]
        if mass_min_value is None:
            mass_min_value = cfg_i["mass_min_value"]
        if mass_max_value is None:
            mass_max_value = cfg_i["mass_max_value"]
        if mass_resolution is None:
            mass_resolution = cfg_i["mass_resolution"]
        if mass_axis_scale is None:
            mass_axis_scale = cfg_i["mass_axis_scale"]

        cfg_i = cfg["time_discretization"]
        if time_min_value is None:
            time_min_value = cfg_i["time_min_value"]
        if time_max_value is None:
            time_max_value = cfg_i["time_max_value"]
        if time_resolution is None:
            time_resolution = cfg_i["time_resolution"]
        if time_axis_scale is None:
            time_axis_scale = cfg_i["time_axis_scale"]

        cfg_i = cfg["kernel"]
        if enable_coagulation is None:
            enable_coagulation = cfg_i["enable_coagulation"]
        if enable_fragmentation is None:
            enable_fragmentation = cfg_i["enable_fragmentation"]
        if enable_physical_cross_sections is None:
            enable_physical_cross_sections = cfg_i["enable_physical_cross_sections"]
        if enable_physical_relative_velocities is None:
            enable_physical_relative_velocities = cfg_i["enable_physical_relative_velocities"]
        if enable_cancellation_handling is None:
            enable_cancellation_handling = cfg_i["enable_cancellation_handling"]
        if fragmentation_variants is None:
            fragmentation_variants = cfg_i["fragmentation_variants"]
        if collision_outcome_variant is None:
            collision_outcome_variant = cfg_i["collision_outcome_variant"]

        cfg_i = cfg["solver"]
        if solver_variant is None:
            solver_variant = cfg_i["solver_variant"]

        cfg_i = cfg["plotting"]
        if mpl_dark_mode is None:
            mpl_dark_mode = cfg_i["mpl_dark_mode"]

        self.mass_min_value = mass_min_value
        self.mass_max_value = mass_max_value
        self.mass_resolution = mass_resolution
        self.mass_axis_scale = mass_axis_scale
        self.time_min_value = time_min_value
        self.time_max_value = time_max_value
        self.time_resolution = time_resolution
        self.time_axis_scale = time_axis_scale
        self.radial_min_value = radial_min_value
        self.radial_max_value = radial_max_value
        self.radial_resolution = radial_resolution
        self.radial_axis_scale = radial_axis_scale
        self.enable_coagulation = enable_coagulation
        self.enable_fragmentation = enable_fragmentation
        self.enable_physical_cross_sections = enable_physical_cross_sections
        self.enable_physical_relative_velocities = enable_physical_relative_velocities
        self.enable_cancellation_handling = enable_cancellation_handling
        self.enable_physical_gas_density = enable_physical_gas_density
        self.fragmentation_variants = fragmentation_variants
        self.collision_outcome_variant = collision_outcome_variant
        self.solver_variant = solver_variant
        self.stellar_mass = stellar_mass
        self.stellar_luminosity = stellar_luminosity
        self.dust_to_gas_ratio = dust_to_gas_ratio
        self.flaring_angle = flaring_angle
        self.mpl_dark_mode = mpl_dark_mode
        self.distance_to_star = distance_to_star
        self.disk_mass_ratio = disk_mass_ratio
        self.disk_mass = disk_mass
