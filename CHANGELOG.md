# dincalcilab/samurai: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.5 - [2024-12-12]

### Documentation

- Add usage docs (PR #31, @SaraPotente) (#14)

### Internal changes

- Rename `test_ascat` profile to `test`

## v1.0.4 - [2024-11-29]

### Bug Fixes

- Fix generating WisecondorX plots (@lbeltrame) (#28)

### Internal changes

- Initial work towards implementation of a test profile (@SaraPotente) (#6)

## v1.0.3 - [2024-10-16]

### Bug Fixes

- Fix downloading genomes from iGenomes (#19)
- Fix a pipeline hang when a panel of normals was built at the same time as an ichorCNA run (#20)

## v1.0.2 - [2024-09-13]

### Bug Fixes

- Fix CIN signature output file detection.

## v1.0.1 - [2024-09-05]

### Bug Fixes

- Bump container version for gistic-cli to fix an off-by-one error in upstream gistic-cli.

## v1.0.0 - [2024-08-08]

Initial release of dincalcilab/samurai, created with the [nf-core](https://nf-co.re/) template.
