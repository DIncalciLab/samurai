# dincalcilab/samurai: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.2.3 - "Nitta Yoshisada" - [2025-07-21]

This release fixes a critical bug in the logR correction for GISTIC in ichorCNA, which caused wrong logR values to be generated. In addition, there are other fixes in ichorCNA (in particular now PoNs in data without chromosome X can be generated) and a small new feature (absolute integer CN in ASCAT.sc). Nitta Yoshisada (新田義貞, 1301-1338) was a samurai lord during the Nanboku-cho period in Japan, also briefly involved in the wars following the so-called "Kenmu Restoration".

## New Features

- Add integer absolute copy numbers to ASCAT.sc (#38, @SaraPotente)

## Bug Fixes

- Fix GISTIC logR correction in ichorCNA (@lbeltrame, 002273536314d96f8446ee0e9cd84836af7ccf45)
- Make sure the ichorCNA PoN generation works even in absence of data on chrX (@lbeltrame, 5e292cef413676448f503f97660c395ee9ba0d99)

## v1.2.2 - "Shibata Katsuie" - [2025-07-10]

This release fixes a bug in the workflow caused by the new genomic plotting feature in the ichorCNA workflow.

Shibata Katsuie (柴田勝家; 1522-1583), also known with the name of Gonroku (権六) was a retainer of the Nobunaga clan, serving under Oda Nobuhide and Oda Nobunaga during the "Warring States" period in Japan.

### Bug Fixes

- Fix undefined output variable (9418166871a6504ea25470af5e8ef2393f216664, @SaraPotente)

## v1.2.1 - "Maeda Toshiie" - [2025-07-04]

This release fixes a bug in the new ichorCNA plotting introduced in version 1.2.0: plotting would not work due to a missing dependency in the container. A new container image has been made to correct this problem.

Maeda Toshiie (前田利家; 1538-1599), also known as "Yari no Mataza" (槍の又左) was a leading general under Oda Nobunaga during the "Warring States" period in Japan.

### Bug Fixes

- Add missing `svglite` dependency to the container (3a0f0760617bc1749614305b2883c4c19c114d00)

## v1.2.0 - "Ohori Tsuruhime" - [2025-06-30]

This release adds support for ichorCNA in the solid biopsy workflow (see PR #42 for caveats) and new ploidy-aware plots for ichorCNA (off by default), which can be enabled by specifying the option `--ichorcna_ploidy_aware_plot`.

This release is named after Ohori Tsuruhime (大祝鶴姫; 1526–1543), an "onna-musha" (female warrior) who lived and fought during the "Warring States" period in Japan. A suit of armor allegedly belonging to Tsuruhime is kept in the treasure hall of Oyamazumi Shrine, on the island of Omishima.

### New features

- New plots for ichorCNA output (PR #44, @SaraPotente)
- Allow ichorCNA also in the solid biopsy workflow (PR #42, @lbeltrame)

### Internal changes

- Move ichorCNA to a separate subworkflow

## v1.1.0 - "Takeda Shingen" - [2025-04-29]

### New features

- Integration of HRDCNA module for HRDCNA score detection in solid biopsy subworkflow (PR #35, @SaraPotente)

Takeda Shingen (1521-1573) was a warlord who ruled the province of Kai during the "Warring States" period.

## v1.0.6 - [2025-01-30]

### Bug Fixes

- Fix genomic plot concatenations when running under Docker (#32)
- Allow running WisecondorX with size selection (8824631a699a33fa66fbe9c6b0021e021a64272f)

### Internal changes

- Don't try to be smart with HMMcopy's ReadCounter (1bbae5755acda8574b508ab93481fb5d0268cfc5)

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
