# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres
to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.0.5] - 2021-05-11
### Added
- support for DACE opacities original data

## [1.0.4] - 2021-05-04
### Added
- added Rayleigh scattering converter. This class computes the Rayleigh scattering for the indicated atom and can be used as input data for the RAPOC models.
### Changed
- we removed "cubic" mode and added "loglinear" instead. This mode is the same as "linear" but uses the logarithm of pressure instead of simple pressure.
- we added ``self.exceptions`` variable to model. If it is not `None`, ``model.estimate`` and the validation procedures behave accordingly. The variable is given by the ``loader.read_content`` and is ``None`` by default.
- _ktable_ name replaced with _opacities_

## [1.0.0] - 2021-03-02
First RAPOC release


[Unreleased]: https://github.com/ExObsSim/Rapoc
[1.0.5]: https://github.com/ExObsSim/Rapoc-public/compare/v1.0.4...v1.0.5
[1.0.4]: https://github.com/ExObsSim/Rapoc-public/compare/v1.0.0...v1.0.4
[1.0.0]: https://github.com/ExObsSim/Rapoc-public/releases/tag/v1.0.0