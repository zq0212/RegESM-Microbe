# RegESM-Microbe

Source-code replacements for coupling CLM-Microbe to RegCM-Chem-YIBs, supporting the study *Divergent Responses of Methane Emissions to Nitrogen and Sulfur Deposition across Asian Wetlands*.

## Upstream models

- **RegCM-Chem-YIBs**: Xie, N., Wang, T., & Xie, X. (2023). *RegCM-Chem-YIBs: RegCM-Chem-YIBs_v1.0*. Zenodo. [doi:10.5281/zenodo.8186164](https://doi.org/10.5281/zenodo.8186164). [Source and installation instructions](https://github.com/honegege/RegCM-Chem-YIBs/tree/code).
- **CLM-Microbe**: Xu, X., He, L., & Wang, Y. (2022). *CLM-Microbe v1.0* (1.0). Zenodo. [doi:10.5281/zenodo.7439312](https://doi.org/10.5281/zenodo.7439312). [Source repository](https://github.com/email-clm/CLM-Microbe).

Please cite both upstream model releases when using this code. The online coupling is described in [Zhang et al. (2025), Earth's Future](https://doi.org/10.1029/2025EF005991).

## Installation and use

1. Prepare a working RegCM-Chem-YIBs installation using the upstream source, dependencies, and installation instructions above.
2. Back up the corresponding files in that installation, then copy the files in this repository to the destinations below, replacing existing files with matching names. Preserve the rest of the upstream source tree.

| Files in this repository | Destination in the RegCM-Chem-YIBs source tree |
| --- | --- |
| `mod_lm_interface.F90`, `mod_output.F90` | `Main/` |
| `chemlib/*.F90` | `Main/chemlib/` |
| `mpplib/*.F90` | `Main/mpplib/` |
| `radlib/*.F90` | `Main/radlib/` |
| `clm4.5/*.F90`, `clm4.5/Makefile.am`, `clm4.5/Makefile.in` | `Main/clmlib/clm4.5/` |

3. Reconfigure and rebuild the receiving model with CLM4.5 and the microbial module enabled. Ensure that the receiving build defines `MICROBE` and includes the microbial source modules. Compiler, MPI, and netCDF settings follow the upstream installation. Regenerate build files as required by that installation after applying the replacements.
4. Prepare the required forcing, initial conditions, surface data, and namelist for the intended experiment, then run the rebuilt model following the upstream workflow.

### Configuration of the supplied snapshot

The supplied `clm4.5/mod_clm_regcm.F90` currently sets `clm_a2l%forc_sdep(:)` to zero. The alternative deposition conversion is present as a comment on the same line. Check this assignment against the intended experiment before enabling S deposition. The supplied `mod_clm_microbeMod.F90` contains empirical S-dependent multipliers in methane-production calculations; inspect these locations when configuring or interpreting the S treatment.

This repository supplies source replacements. Simulation outputs and experiment input datasets are not included. Compilation and model execution must be checked in the receiving RegCM-Chem-YIBs environment; they have not been verified as part of this source publication. Compiled libraries, editor backups, and machine-generated Makefiles are excluded from version control.

## License and attribution

The combined distribution is provided under **GNU GPL version 3 or later**, consistent with the RegCM license notices in the included source files. See [LICENSE](LICENSE). Original file-level copyright and license notices are retained. Upstream MIT licenses and the CESM copyright notice are reproduced in [LICENSES/](LICENSES/); see [NOTICE.md](NOTICE.md) for attribution.

Maintainer: Qian Zhang, Laoshan Laboratory.
