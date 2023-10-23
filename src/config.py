import os
from typing import Optional

import toml

from constants import M_sun, L_sun, AU

dirname = os.path.dirname(__file__)
PATH_TO_CFG_TOML = os.path.join(dirname, "..", "config.toml")
PATH_TO_LIB      = os.path.join(dirname, "..", "lib")
PATH_TO_FIGURES  = os.path.join(dirname, "..", "out", "figures")
PATH_TO_DARKMODE = os.path.join(PATH_TO_LIB, "mpl-styles", "dark.mplstyle")
PATH_TO_COAG     = os.path.join(PATH_TO_LIB, "coag_py")


class Config():
    __slots__ = [
        "stellar_mass", "stellar_luminosity", "disk_mass_ratio", "dust_to_gas_ratio",
        "distance_to_star", "flaring_angle", 
        "radial_min_value", "radial_max_value", "radial_resolution", "radial_axis_scale",
        "mass_min_value", "mass_max_value", "mass_resolution", "mass_axis_scale",
        "time_min_value", "time_max_value", "time_resolution", "time_axis_scale",
        "enable_coagulation", "enable_fragmentation", "enable_physical_gas_density",
        "enable_physical_collisions", "enable_cancellation_handling", "enable_collision_sampling",
        "relative_velocity_components", "fragmentation_variant", "fragmentation_velocity",
        "collision_outcome_variant", "solver_variant", "mpl_dark_mode", 
        "dust_particle_density", "viscosity_alpha", "nr_of_samples", "disk_mass",
        "initialization_variant", "initial_mass_bin",
    ]

    def __init__(
        self,
        stellar_mass:                   Optional[float]     = None,
        stellar_luminosity:             Optional[float]     = None,
        disk_mass_ratio:                Optional[float]     = None,
        dust_to_gas_ratio:              Optional[float]     = None,
        distance_to_star:               Optional[float]     = None,
        flaring_angle:                  Optional[float]     = None,
        radial_min_value:               Optional[float]     = None,
        radial_max_value:               Optional[float]     = None,
        radial_resolution:              Optional[int]       = None,
        radial_axis_scale:              Optional[str]       = None,
        mass_min_value:                 Optional[float]     = None,
        mass_max_value:                 Optional[float]     = None,
        mass_resolution:                Optional[int]       = None,
        mass_axis_scale:                Optional[str]       = None,
        time_min_value:                 Optional[float]     = None,
        time_max_value:                 Optional[float]     = None,
        time_resolution:                Optional[int]       = None,
        time_axis_scale:                Optional[str]       = None,
        enable_coagulation:             Optional[bool]      = None,
        enable_fragmentation:           Optional[bool]      = None,
        enable_physical_gas_density:    Optional[bool]      = None,
        enable_physical_collisions:     Optional[bool]      = None,
        relative_velocity_components:   Optional[list[str]] = None,
        enable_cancellation_handling:   Optional[bool]      = None,
        fragmentation_variant:          Optional[str]       = None,
        fragmentation_velocity:         Optional[float]     = None,
        collision_outcome_variant:      Optional[str]       = None,
        solver_variant:                 Optional[str]       = None,
        mpl_dark_mode:                  Optional[bool]      = None,
        dust_particle_density:          Optional[float]     = None,
        viscosity_alpha:                Optional[float]     = None,
        enable_collision_sampling:      Optional[bool]      = None,
        nr_of_samples:                  Optional[int]       = None,
        initialization_variant:         Optional[str]       = None,
        initial_mass_bin:               Optional[int]       = None,
    ):
        cfg = toml.load(PATH_TO_CFG_TOML)

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
        if enable_physical_collisions is None:
            enable_physical_collisions = cfg_i["enable_physical_collisions"]
        if dust_particle_density is None:
            dust_particle_density = cfg_i["dust_particle_density"]
        if viscosity_alpha is None:
            viscosity_alpha = cfg_i["viscosity_alpha"]
        if relative_velocity_components is None:
            relative_velocity_components = cfg_i["relative_velocity_components"]
        if fragmentation_variant is None:
            fragmentation_variant = cfg_i["fragmentation_variant"]
        if fragmentation_velocity is None:
            fragmentation_velocity = cfg_i["fragmentation_velocity"]
        if collision_outcome_variant is None:
            collision_outcome_variant = cfg_i["collision_outcome_variant"]

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
        if enable_cancellation_handling is None:
            enable_cancellation_handling = cfg_i["enable_cancellation_handling"]
        if enable_collision_sampling is None:
            enable_collision_sampling = cfg_i["enable_collision_sampling"]
        if nr_of_samples is None:
            nr_of_samples = cfg_i["nr_of_samples"]

        cfg_i = cfg["solver"]
        if solver_variant is None:
            solver_variant = cfg_i["solver_variant"]

        cfg_i = cfg["plotting"]
        if mpl_dark_mode is None:
            mpl_dark_mode = cfg_i["mpl_dark_mode"]

        cfg_i = cfg["initialization"]
        if initialization_variant is None:
            initialization_variant = cfg_i["initialization_variant"]
        if initial_mass_bin is None:
            initial_mass_bin = cfg_i["initial_mass_bin"]

        self.collision_outcome_variant    = collision_outcome_variant
        self.disk_mass                    = disk_mass
        self.disk_mass_ratio              = disk_mass_ratio
        self.distance_to_star             = distance_to_star
        self.dust_particle_density        = dust_particle_density
        self.dust_to_gas_ratio            = dust_to_gas_ratio
        self.enable_cancellation_handling = enable_cancellation_handling
        self.enable_coagulation           = enable_coagulation
        self.enable_fragmentation         = enable_fragmentation
        self.enable_physical_collisions   = enable_physical_collisions
        self.enable_physical_gas_density  = enable_physical_gas_density
        self.flaring_angle                = flaring_angle
        self.fragmentation_variant        = fragmentation_variant
        self.fragmentation_velocity       = fragmentation_velocity
        self.mass_axis_scale              = mass_axis_scale
        self.mass_max_value               = mass_max_value
        self.mass_min_value               = mass_min_value
        self.mass_resolution              = mass_resolution
        self.mpl_dark_mode                = mpl_dark_mode
        self.radial_axis_scale            = radial_axis_scale
        self.radial_max_value             = radial_max_value
        self.radial_min_value             = radial_min_value
        self.radial_resolution            = radial_resolution
        self.relative_velocity_components = relative_velocity_components
        self.solver_variant               = solver_variant
        self.stellar_luminosity           = stellar_luminosity
        self.stellar_mass                 = stellar_mass
        self.time_axis_scale              = time_axis_scale
        self.time_max_value               = time_max_value
        self.time_min_value               = time_min_value
        self.time_resolution              = time_resolution
        self.viscosity_alpha              = viscosity_alpha
        self.enable_collision_sampling    = enable_collision_sampling
        self.nr_of_samples                = nr_of_samples
        self.initialization_variant       = initialization_variant
        self.initial_mass_bin             = initial_mass_bin
